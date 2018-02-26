#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <algorithm>

#include "arithmetic_encode.h"
#include "arithmetic_coder_ut.h"
#include "arithmetic_coder_cfg.h"

static const int      QH_SCALE  = 1 << QH_BITS;
static const unsigned VAL_RANGE = 1u << RANGE_BITS;

void arithmetic_encode_init_context(uint32_t* context)
{
  memset(context, 0, sizeof(uint32_t)*257);
}

void arithmetic_encode_chunk_callback(void* context, const uint8_t* src, int chunklen)
{
  uint32_t* histogram = static_cast<uint32_t*>(context);
  for (int i = 0; i < chunklen; ++i) {
    int c = src[i];
    if (c == 255) {
      // escape sequence for symbols 255 and 256
      c = int(src[i+1]) + 1;
      ++i;
    }
    ++histogram[c];
  }
}


static void quantize_histogram(uint8_t* __restrict qh, uint32_t *h, int len, uint32_t hTot)
{
  double invhLen = 2.0/hTot;
  double M_PI = 3.1415926535897932384626433832795;
  for (unsigned c = 0; c < len; ++c) {
    unsigned val = 0;
    unsigned cnt = h[c];
    if (cnt != 0) {
      val = round((asin(cnt*invhLen - 1.0)/M_PI + 0.5)*QH_SCALE);
      if (val < 1)
        val = 1;
      else if (val >= QH_SCALE)
        val = QH_SCALE-1;
    }
    qh[c] = val;
  }
}

// return value:
// maxC = the character with the biggest numeric value that appears in the source at least once
static int prepare1(
  uint32_t * __restrict histogram,
  uint8_t*   __restrict qh, // quantized histogram
  double* pInfo)
{
  // find highest-numbered character that occurred at least once
  unsigned maxC = 0;
  int32_t tot = 0;
  for (unsigned c = 0; c < 257; ++c) {
    tot += histogram[c];
    if (histogram[c] != 0)
      maxC = c;
  }

  if (pInfo) {
    // calculate source entropy of static part of histogram
    double entropy = log2(double(tot))*tot;
    for (unsigned c = 0; c <= maxC; ++c) {
      int32_t cnt = histogram[c];
      if (cnt)
        entropy -= log2(double(cnt))*cnt;
    }
    pInfo[0] = entropy;
    pInfo[3] = 0;
  }
  quantize_histogram(qh, histogram, maxC+1, tot);

  return maxC;
}

static void prepare2(
  uint16_t* __restrict c2low, unsigned maxC,
  const unsigned       h[257],
  const uint8_t        qh[257],
  double*              pQuantizedEntropy)
{
  quantized_histogram_to_range(c2low, maxC+1, qh, VAL_RANGE);

  // calculate entropy after quantization
  double entropy = 0;
  for (unsigned c = 0; c <= maxC; ++c) {
    unsigned cnt = h[c];
    if (cnt)
      entropy += log2(double(VAL_RANGE)/c2low[c])*cnt;
  }
  *pQuantizedEntropy = entropy;

  // c2low -> cumulative sums of ranges
  unsigned lo = 0;
  for (unsigned c = 0; c <= maxC; ++c) {
    unsigned range = c2low[c];
    // printf("%3u: %04x\n", c, range);
    c2low[c] = lo;
    lo += range;
  }
  c2low[maxC+1] = VAL_RANGE;
}

class CArithmeticEncoder {
public:
  uint8_t* put(uint64_t cScale, uint64_t cLo, uint64_t cRange, uint8_t* dst);
  void spillOverflow(uint8_t* dst);
  void init() {
    m_lo = 0;
    m_range = uint64_t(1) << 63;
  }
  uint64_t m_lo;    // scaled by 2**63
  uint64_t m_range; // scaled by 2**63
};

static uint8_t* encode_val(uint8_t* dst, unsigned val, unsigned offset)
{
  val += 2;
  do {
    *dst++ = (val % 2) + offset;
    val /= 2;
  } while (val > 1);
  return dst;
}

static unsigned encode_qh(uint8_t* dst0, const uint8_t* qh, unsigned len)
{
  uint8_t* dst = dst0;
  int zr = -1;
  for (unsigned c = 0; c < len; ++c) {
    unsigned val = qh[c];
    if (val > 0) {
      if (zr >= 0)
        dst = encode_val(dst, zr, 0);
      zr = 0;
      dst = encode_val(dst, val-1, 2);
    }
  }
  dst = encode_val(dst, 257-len, 0);
  return dst - dst0;
}

// return quantized code of h0/(h0+h1) in range [0..scale]
static unsigned quantize_histogram_pair(unsigned h0, unsigned h1, unsigned scale)
{
  if (h0 == 0)
    return 0;
  if (h1 == 0)
    return scale;
  unsigned hTot = h0 + h1;
  printf("%u+%u => %.2f\n", h0, h1, (log2(hTot)*hTot-log2(h0)*h0-log2(h1)*h1)/8);
  double M_PI = 3.1415926535897932384626433832795;
  unsigned val = round((asin(h0*2.0/hTot - 1.0)/M_PI + 0.5)*scale);
  if (val < 1)
    val = 1;
  else if (val >= scale)
    val = scale-1;
  return val;
}

// return the number of stored octets
static int store_model(uint8_t* dst, const uint8_t qh[257], unsigned maxC, double* pNbits, CArithmeticEncoder* pEnc)
{
  // encode qh
  uint8_t eqh[257*QH_BITS];
  unsigned qhlen = encode_qh(eqh, qh, maxC + 1);
  unsigned hist[4] = {0};
  for (unsigned c = 0; c < qhlen; ++c)
    ++hist[eqh[c]];

  unsigned hh02 = quantize_histogram_pair(hist[0]+hist[1], hist[2]+hist[3], 9); // range [1..8], because both sums > 0
  unsigned hh01 = quantize_histogram_pair(hist[0], hist[1], 8);                 // range [1..8], because hist[0] > 0
  unsigned hh23 = quantize_histogram_pair(hist[2], hist[3], 7);                 // range [0..7]

  unsigned range_tab[3][3];
  range_tab[0][1] = quantized_histogram_pair_to_range_qh_scale9(hh02, VAL_RANGE);
  range_tab[1][1] = quantized_histogram_pair_to_range_qh_scale8(hh01, VAL_RANGE);
  range_tab[2][1] = quantized_histogram_pair_to_range_qh_scale7(hh23, VAL_RANGE);
  for (int k = 0; k < 3; ++k) {
    range_tab[k][0] = 0;
    range_tab[k][2] = VAL_RANGE;
  }
  printf("%.5f %.5f\n", double(hist[0]+hist[1])/(hist[0]+hist[1]+hist[2]+hist[3]), double(range_tab[0][1])/VAL_RANGE);
  printf("%.5f %.5f\n", double(hist[0])/(hist[0]+hist[1]), double(range_tab[1][1])/VAL_RANGE);
  printf("%.5f %.5f\n", double(hist[2])/(hist[2]+hist[3]), double(range_tab[2][1])/VAL_RANGE);

  unsigned hhw = (hh02-1) + (hh01-1)*8 + hh23*64; // combine all hh in a single word

  uint8_t* p = dst;

  // store hhw
  p = pEnc->put(512, hhw, 1, p);
  // store eqh
  for (unsigned c = 0; c < qhlen; ++c)
  {
    unsigned val = eqh[c];
    unsigned val0 = val / 2;
    unsigned val1 = val % 2;
    unsigned lo0 = range_tab[0][val0];
    unsigned ra0 = range_tab[0][val0+1] - lo0;
    unsigned* range_tab1 = range_tab[val0+1];
    unsigned lo1 = range_tab1[val1];
    unsigned ra1 = range_tab1[val1+1] - lo1;
    p = pEnc->put(VAL_RANGE, lo0, ra0, p);
    if (ra1 != 0)
      p = pEnc->put(VAL_RANGE, lo1, ra1, p);
  }

  int len = p - dst;
  *pNbits = len*8.0 + 63 - log2(pEnc->m_range);
  return len;
}

// static uint8_t* inc_dst(uint8_t* dst0, uint8_t* dst) {
  // uint8_t val = 0;
  // uint8_t* inc;
  // for (inc = dst-1; val == 0 && inc != dst0; --inc)
    // *inc = (val = *inc + 1);
  // if (val == 0) {
    // memmove(dst0+1, dst0, dst-dst0);
    // dst0[0] = 1;
    // dst += 1;
  // }
  // return dst;
// }

static void inc_dst(uint8_t* dst) {
  uint8_t val;
  do {
    --dst;
    *dst = (val = *dst + 1);
  } while (val==0);
}

static int encode(uint8_t* dst, const uint8_t* src, unsigned srclen, const uint16_t c2low[257], CArithmeticEncoder* pEnc)
{
  const uint64_t MSB_MSK   = uint64_t(255) << 56;
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);
  uint64_t lo     = pEnc->m_lo << 1;              // scaled by 2**64
  uint64_t range  = pEnc->m_range >> (RANGE_BITS-1); // scaled by 2**50
  uint8_t* dst0   = dst;
  uint64_t prevLo = lo;
  // int dbg_i = 0;
  for (unsigned i = 0; i < srclen; ++i) {
    int c = src[i];
    if (c == 255) {
      // escape sequence for symbols 255 and 256
      c = int(src[i+1]) + 1;
      ++i;
    }
    // printf("[%3d]=%3d\n", dbg_i, c); ++dbg_i;
    uint64_t cLo = c2low[c+0];
    uint64_t cRa = c2low[c+1] - cLo;
    lo   += range * cLo;
    range = range * cRa;

    // at this point range is scaled by 2**64 - the same scale as lo
    uint64_t nxtRange = range >> RANGE_BITS;
    if (nxtRange <= MIN_RANGE) {
      // re-normalize
      if (lo < prevLo) // lo overflow
        inc_dst(dst); //dst = inc_dst(dst0, dst);
      dst[0] = lo >> (64-8*1);
      dst[1] = lo >> (64-8*2);
      dst[2] = lo >> (64-8*3);
      dst += 3;
      lo <<= 24;
      nxtRange = range << (24-RANGE_BITS);
      prevLo = lo;
    }
    range = nxtRange;
  }
  // output last bits
  if (lo < prevLo) // lo overflow
    inc_dst(dst); //dst = inc_dst(dst0, dst);
  uint64_t hi = lo + ((range<<RANGE_BITS)-1);
  if (hi > lo) {
    uint64_t dbits = lo ^ hi;
    while ((dbits & MSB_MSK)==0) {
      // lo and hi have the same upper octet
      *dst++ = uint8_t(hi>>56);
      hi    <<= 8;
      dbits <<= 8;
    }
    // put out last octet
    *dst++ = uint8_t(hi>>56);
  } else {
    inc_dst(dst); //dst = inc_dst(dst0, dst);
  }
  return dst - dst0;
}

uint8_t* CArithmeticEncoder::put(uint64_t cScale, uint64_t cLo, uint64_t cRange, uint8_t* dst)
{
  const uint64_t MIN_RANGE = uint64_t(1) << (33-1);
  const uint64_t BIT_63    = uint64_t(1) << 63;
  const uint64_t MSK_63    = BIT_63 - 1;
  uint64_t range = m_range / cScale;
  uint64_t lo = m_lo + range * cLo;
  range *= cRange;

  if (range <= MIN_RANGE) {
    // re-normalize
    if (lo >= BIT_63) // lo overflow
      inc_dst(dst);
    dst[0] = lo >> (63-8*1);
    dst[1] = lo >> (63-8*2);
    dst[2] = lo >> (63-8*3);
    dst += 3;
    lo = (lo << 24) & MSK_63;
    range <<= 24;
  }
  m_lo    = lo;
  m_range = range;
  return dst;
}

void CArithmeticEncoder::spillOverflow(uint8_t* dst)
{
  const uint64_t BIT_63 = uint64_t(1) << 63;
  const uint64_t MSK_63 = BIT_63 - 1;
  if (m_lo >= BIT_63) // lo overflow
    inc_dst(dst);
  m_lo &= MSK_63;
}

// return value:
//  0 - not compressible, because all input characters have approximately equal probability or because input is too short
// >0 - the length of compressed buffer
int arithmetic_encode(uint32_t* context, uint8_t* dst, const uint8_t* src, int srclen, int origlen, double* pInfo)
{
  uint8_t qh[257]; // quantized histogram
  int maxC = prepare1(context, qh, pInfo);

  CArithmeticEncoder enc;
  enc.init();
  double modelLenBits;
  unsigned modellen = store_model(dst, qh, maxC, &modelLenBits, &enc);
  if (pInfo)
    pInfo[1] = modelLenBits;

  uint16_t c2low[258];
  double quantizedEntropy;
  prepare2(c2low, maxC, context, qh, &quantizedEntropy);

  int lenEst = (modelLenBits+quantizedEntropy+7)/8;
  if (pInfo)
    pInfo[3] = quantizedEntropy;
  if (lenEst >= origlen)
    return 0; // not compressible

  int dstlen = encode(&dst[modellen], src, srclen, c2low, &enc);

  if (pInfo)
    pInfo[2] = dstlen*8.0;

  int reslen = modellen + dstlen;

  if (reslen >= origlen)
    return 0; // not compressible

  return reslen;
}
