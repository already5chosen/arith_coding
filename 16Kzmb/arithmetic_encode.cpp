// #include <cstdio>
#include <cmath>
#include <cstring>

#include "arithmetic_encode.h"


static const int RANGE_BITS = 14;
static const unsigned VAL_RANGE = 1u << RANGE_BITS;

void arithmetic_encode_calc_histograms(unsigned histograms[2][260], const uint8_t* src, int srclen)
{
  // calculated histograms of appearance of characters
  memset(histograms, 0, sizeof(histograms[0][0])*260*2);
  int c0 = -1;
  for (unsigned i = 0; i < srclen; ++i) {
    int c = src[i];
    ++histograms[1][c==c0 ? 0 : c+1]; // histogram that counts repetition as a special character
    ++histograms[0][c]; // simple histogram
    c0 = c;
  }
}

double arithmetic_encode_calc_entropy(const unsigned histogram[260], int srclen)
{
  double entropy = 0;
  int rem = srclen;
  double invSrclen = 1.0/srclen;
  for (int c = 0; rem > 0; ++c) {
    unsigned cnt = histogram[c];
    if (cnt) {
      entropy -= log2(cnt*invSrclen)*cnt;
    }
    rem -= cnt;
  }
  return entropy < 1.0 ? 0 : entropy;
}

// return value:
// -1  - source consists of repetition of the same character
// >=0 - maxC = the character with the biggest numeric value that appears in the source at least once
static int prepare1(const uint8_t* src, unsigned srclen, const unsigned histogram[260], uint16_t c2low[257], double* pQuantizedEntropy, double* pInfo)
{
  // find maximal frequency and highest-numbered character that occurred at least once
  unsigned maxC = 256;
  unsigned maxCnt = 0;
  for (unsigned c = 0; c < 257; ++c) {
    unsigned cnt = histogram[c];
    if (cnt != 0) {
      maxC = c;
      if (maxCnt < cnt)
        maxCnt = cnt;
    }
  }

  if (pInfo) {
    pInfo[3] = 0;
  }

  // it's responsibility of callee to assure that maxCnt < srclen

  // translate counts to ranges and store in c2low
  // 1st pass - translate characters with counts that are significantly lower than maxCnt
  unsigned thr = maxCnt - maxCnt/8;
  unsigned remCnt   = srclen;
  unsigned remRange = VAL_RANGE;
  for (unsigned c = 0; c <= maxC; ++c) {
    unsigned cnt = histogram[c];
    if (cnt < thr) {
      unsigned range = 0;
      if (cnt != 0) {
        // calculate range from full statistics
        range = (uint64_t(cnt)*(VAL_RANGE*2) + srclen)/(srclen*2);
        if (range == 0)
          range = 1;
        remCnt   -= cnt;
        remRange -= range;
      }
      c2low[c] = range;
    }
  }
  // 2nd pass - translate characters with higher counts
  for (unsigned c = 0; remCnt != 0; ++c) {
    unsigned cnt = histogram[c];
    if (cnt >= thr) {
      // Calculate range from the remaining range and count
      // It is non-ideal, but this way we distribute the worst rounding errors
      // relatively evenly among higher ranges, where it hase the smallest impact
      unsigned range = (uint64_t(cnt)*(remRange*2) + remCnt)/(remCnt*2);
      // (range < VAL_RANGE) is guaranteed , because we already handled the case of repetition of the same character
      remCnt   -= cnt;
      remRange -= range;
      c2low[c] = range;
    }
  }

  // calculate entropy after quantization
  double entropy = 0;
  for (unsigned c = 0; c <= maxC; ++c) {
    unsigned cnt = histogram[c];
    if (cnt)
      entropy += log2(double(VAL_RANGE)/c2low[c])*cnt;
  }
  *pQuantizedEntropy = entropy;
  if (pInfo)
    pInfo[3] = entropy;
  return maxC;
}

static void prepare2(uint16_t c2low[258], unsigned maxC)
{
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

static int inline floor_log2(unsigned x)
{
  return sizeof(unsigned)*8 - 1 - __builtin_clz(x); // floor(log2(x))
}

// return the number of stored octets
static int store_model(uint8_t* dst, const uint16_t c2range[257], unsigned maxC, double* pNbits, CArithmeticEncoder* pEnc, bool rep)
{
  int nRanges = 0;
  unsigned hist[RANGE_BITS+1] = {0};
  for (unsigned c = 0; c < maxC; ++c) {
    uint32_t range = c2range[c];
    if (range > 0) {
      ++nRanges;
      hist[1+floor_log2(range)] += 1;
    }
  }
  // No need to store a last range.
  // Decoder will calculate it as VAL_RANGE-sum of previous ranges
  hist[0] = maxC - nRanges; // # of zero ranges before maxC

  // for (int i = 0; i< RANGE_BITS+1; ++i)
    // printf("hist[%2d]=%d\n", i, hist[i]);
  // for (int i = 0; i< 256; ++i)
    // printf("range[%3d]=%5d\n", i, c2range[i]);

  uint8_t* p = dst;

  // store maxC or maxC-1
  p = pEnc->put(256, rep ? maxC-1 : maxC, 1, p);
  // store histogram of log2
  unsigned rem = maxC;
  for (int i = 0; i < RANGE_BITS && rem > 0; ++i) {
    p = pEnc->put(rem+1, hist[i], 1, p);
    rem -= hist[i];
  }

  // store c2range
  for (unsigned c = 0; c < maxC; ++c) {
    uint32_t range = c2range[c];
    int log2_i = range == 0 ? 0 : 1+floor_log2(range);
    int lo = 0;
    for (int i = 0; i < log2_i; ++i)
      lo += hist[i];
    p = pEnc->put(maxC, lo, hist[log2_i], p); // exp.range
    if (log2_i > 1) {
      uint32_t expRangeSz = uint32_t(1) << (log2_i-1);
      p = pEnc->put(expRangeSz, range-expRangeSz, 1, p); // offset within range
    }
  }
  pEnc->spillOverflow(p);

  int len = p - dst;
  *pNbits = len*8.0 + 63 - log2(pEnc->m_range);
  return len;
}

static void inc_dst(uint8_t* dst) {
  uint8_t val;
  do {
    --dst;
    *dst = (val = *dst + 1);
  } while (val==0);
}

static int encode(uint8_t* dst, const uint8_t* src, unsigned srclen, const uint16_t c2low[258], CArithmeticEncoder* pEnc)
{
  const uint64_t MSB_MSK   = uint64_t(255) << 56;
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);
  uint64_t lo     = pEnc->m_lo << 1;              // scaled by 2**64
  uint64_t range  = pEnc->m_range >> (RANGE_BITS-1); // scaled by 2**50
  uint8_t* dst0   = dst;
  uint64_t prevLo = lo;
  for (unsigned i = 0; i < srclen; ++i) {
    int c = src[i];
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


static int encode_rep(uint8_t* dst, const uint8_t* src, unsigned srclen, const uint16_t c2low[258], CArithmeticEncoder* pEnc)
{
  const uint64_t MSB_MSK   = uint64_t(255) << 56;
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);
  uint64_t lo     = pEnc->m_lo << 1;              // scaled by 2**64
  uint64_t range  = pEnc->m_range >> (RANGE_BITS-1); // scaled by 2**50
  uint8_t* dst0   = dst;
  uint64_t prevLo = lo;
  int c0 = -1;
  for (unsigned i = 0; i < srclen; ++i) {
    int c1 = src[i];
    int c = c0 == c1 ? 0 : c1 + 1;
    c0 = c1;
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
int arithmetic_encode(
  std::vector<uint8_t>* dst,
  const uint8_t* src, int srclen,
  const unsigned histogram[260], bool rep, double* pInfo)
{
  uint16_t c2low[258];
  double quantizedEntropy;
  int maxC = prepare1(src, srclen, histogram, c2low, &quantizedEntropy, pInfo);

  size_t sz0 = dst->size();
  dst->resize(sz0 + 640);
  CArithmeticEncoder enc;
  enc.init();
  double modelLenBits;
  unsigned modellen = store_model(&dst->at(sz0), c2low, maxC, &modelLenBits, &enc, rep);
  if (pInfo)
    pInfo[1] = modelLenBits;

  int lenEst = (modelLenBits+quantizedEntropy+7)/8;
  if (lenEst >= srclen)
    return 0; // not compressible

  // printf("ml=%u\n", modellen);
  dst->resize(sz0 + lenEst + 64);
  prepare2(c2low, maxC);
  int dstlen = rep ?
   encode_rep(&dst->at(sz0+modellen), src, srclen, c2low, &enc):
   encode    (&dst->at(sz0+modellen), src, srclen, c2low, &enc);

  if (pInfo)
    pInfo[2] = dstlen*8.0;

  int reslen = modellen + dstlen;

  if (reslen >= srclen)
    return 0; // not compressible

  dst->resize(sz0 + reslen);
  return reslen;
}
