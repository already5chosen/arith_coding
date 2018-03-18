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

static void resize_up(std::vector<uint32_t>* vec, size_t len)
{
  if (vec->size() < len)
    vec->resize(len);
}

struct context_chunk_t {
  uint32_t histogram[ARITH_CODER_N_DYNAMIC_SYMBOLS];
  bool     hasDictionary;
  uint8_t  qHistogram[ARITH_CODER_N_DYNAMIC_SYMBOLS];
  uint16_t ranges[ARITH_CODER_N_DYNAMIC_SYMBOLS];
};

static const int CONTEX_HDR_LEN = 1 + 258;
static const int CONTEX_CHUNK_LEN = (sizeof(context_chunk_t)-1)/sizeof(uint32_t) + 1;

void arithmetic_encode_init_context(std::vector<uint32_t>* context, int tilelen)
{
  int nChunks = tilelen/(ARITH_CODER_RUNS_PER_CHUNK*2) + 1; // estimate # of chunks
  resize_up(context, CONTEX_HDR_LEN+CONTEX_CHUNK_LEN*nChunks);
  memset(&context->at(0), 0, sizeof(uint32_t)*CONTEX_HDR_LEN);
}

int arithmetic_encode_chunk_callback(void* context, const uint8_t* chunk, int chunklen, int nRuns)
{
  if (chunk) {
    std::vector<uint32_t>* vec = static_cast<std::vector<uint32_t>*>(context);
    uint32_t chunk_i = (*vec)[0];
    resize_up(vec, CONTEX_HDR_LEN+CONTEX_CHUNK_LEN*(chunk_i+1));
    (*vec)[0] = chunk_i + 1;

    uint32_t  tmp[ARITH_CODER_N_DYNAMIC_SYMBOLS];
    uint32_t* dynSum = &vec->at(1);
    memcpy(tmp, dynSum, sizeof(tmp));
    memset(dynSum, 0, sizeof(tmp));

    uint32_t* histogram = &dynSum[1];
    for (int i = 0; i < chunklen; ++i) {
      int c = chunk[i];
      if (c == 255) {
        // escape sequence for symbols 255 and 256
        c = int(chunk[i+1]) + 1;
        ++i;
      }
      ++histogram[c];
    }

    context_chunk_t* chunk = reinterpret_cast<context_chunk_t*>(&vec->at(CONTEX_HDR_LEN+CONTEX_CHUNK_LEN*chunk_i));
    uint32_t dynLast = chunklen; // # of symbols with value >= ARITH_CODER_N_DYNAMIC_SYMBOLS-1
    for (int i = 0; i < ARITH_CODER_N_DYNAMIC_SYMBOLS-1; ++i) {
      uint32_t h = histogram[i];
      chunk->histogram[i] = h;
      dynSum[i] = tmp[i] + h;
      dynLast -= h;
    }
    dynSum[ARITH_CODER_N_DYNAMIC_SYMBOLS-1] = tmp[ARITH_CODER_N_DYNAMIC_SYMBOLS-1] + dynLast;
    chunk->histogram[ARITH_CODER_N_DYNAMIC_SYMBOLS-1] = dynLast;
  }
  return ARITH_CODER_RUNS_PER_CHUNK;
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
// maxC = the character with the bighest numeric value that appears in the source at least once
static int prepare1(
  uint32_t * __restrict context,
  uint8_t*   __restrict qh, // quantized histogram
  double* pInfo)
{
  memset(qh, 0, 258);

  // find highest-numbered character that occurred at least once
  unsigned maxC = 0;
  uint32_t *h = &context[1];
  for (unsigned c = 0; c < ARITH_CODER_N_DYNAMIC_SYMBOLS; ++c) {
    if (h[c] != 0)
      maxC = c;
  }

  double entropy = 0;
  int32_t tot_s = h[ARITH_CODER_N_DYNAMIC_SYMBOLS-1];
  if (tot_s != 0) {
    if (h[maxC] != tot_s) {
      if (pInfo) {
        // calculate source entropy of static part of histogram
        entropy = log2(double(tot_s))*tot_s;
        for (unsigned c = ARITH_CODER_N_DYNAMIC_SYMBOLS; c <= maxC; ++c) {
          int32_t cnt = h[c];
          if (cnt)
            entropy -= log2(double(cnt))*cnt;
        }
      }
      quantize_histogram(&qh[ARITH_CODER_N_DYNAMIC_SYMBOLS], &h[ARITH_CODER_N_DYNAMIC_SYMBOLS], maxC+1-ARITH_CODER_N_DYNAMIC_SYMBOLS, tot_s);
    } else {
      qh[maxC] = 255;  // static part of histogram consists of repetition of the same character
    }
  }

  int nChunks = context[0];
  for (int chunk_i = 0; chunk_i < nChunks; ++chunk_i) {
    context_chunk_t* chunk = reinterpret_cast<context_chunk_t*>(&context[CONTEX_HDR_LEN+CONTEX_CHUNK_LEN*chunk_i]);
    memset(chunk->qHistogram, 0, sizeof(chunk->qHistogram));

    int32_t tot_d = 0;
    unsigned maxC_d = 0;
    for (unsigned c = 0; c < ARITH_CODER_N_DYNAMIC_SYMBOLS; ++c) {
      int32_t cnt = chunk->histogram[c];
      tot_d += cnt;
      if (cnt)
        maxC_d = c;
    }
    if (tot_d) {
      if (chunk->histogram[maxC_d] != tot_d) {
        if (pInfo) {
          // calculate source entropy of dynamic chunk
          entropy += log2(double(tot_d))*tot_d;
          for (unsigned c = 0; c < ARITH_CODER_N_DYNAMIC_SYMBOLS; ++c) {
            int32_t cnt = chunk->histogram[c];
            if (cnt)
              entropy -= log2(double(cnt))*cnt;
          }
        }
        quantize_histogram(chunk->qHistogram, chunk->histogram, maxC_d+1, tot_d);
      } else {
        chunk->qHistogram[maxC_d] = 255; // chunk consists of repetition of the same character
      }
    }
  }

  if (pInfo) {
    pInfo[0] = entropy;
    pInfo[3] = 0;
  }

  return maxC;
}

static void prepare2(
  uint16_t* __restrict c2low, unsigned maxC,
  const unsigned       h[256],
  const uint8_t        qh[256],
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

static int inline floor_log2(unsigned x)
{
  return sizeof(unsigned)*8 - 1 - __builtin_clz(x); // floor(log2(x))
}

// return the number of stored octets
static int store_model(uint8_t* dst, const uint8_t qh[256], unsigned maxC, double* pNbits, CArithmeticEncoder* pEnc)
{
  int nRanges = 0;
  unsigned hist[QH_BITS+1] = {0};
  for (unsigned c = 0; c <= maxC; ++c) {
    uint32_t range = qh[c];
    if (range > 0) {
      ++nRanges;
      hist[1+floor_log2(range)] += 1;
    }
  }
  hist[0] = maxC + 1 - nRanges; // # of zero ranges before maxC

  // for (int i = 0; i< QH_BITS+1; ++i)
    // printf("hist[%2d]=%d\n", i, hist[i]);
  // for (int i = 0; i< 256; ++i)
    // printf("range[%3d]=%5d\n", i, c2range[i]);

  uint8_t* p = dst;

  // store maxC
  p = pEnc->put(256, maxC, 1, p);
  // store histogram of log2
  unsigned rem = maxC + 1;
  for (int i = 0; i < QH_BITS && rem > 0; ++i) {
    p = pEnc->put(rem+1, hist[i], 1, p);
    rem -= hist[i];
  }

  // store c2range
  for (unsigned c = 0; c <= maxC; ++c) {
    uint32_t range = qh[c];
    int log2_i = range == 0 ? 0 : 1+floor_log2(range);
    int lo = 0;
    for (int i = 0; i < log2_i; ++i)
      lo += hist[i];
    p = pEnc->put(maxC+1, lo, hist[log2_i], p); // exp.range
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
// -1 - source consists of repetition of the same character
//  0 - not compressible, because all input characters have approximately equal probability or because input is too short
// >0 - the length of compressed buffer
int arithmetic_encode(uint32_t* context, uint8_t* dst, const uint8_t* src, int nRuns, int origlen, double* pInfo)
{
  uint8_t qh[258]; // quantized histogram
  int maxC = prepare1(context, qh, pInfo);
#if 0
  unsigned h[256];  // histogram
  int maxC = prepare1(h, qh, src, srclen, pInfo);


  if (maxC == -1)
    return -1; // source consists of repetition of the same character

  size_t sz0 = dst->size();
  dst->resize(sz0 + 640);
  CArithmeticEncoder enc;
  enc.init();
  double modelLenBits;
  unsigned modellen = store_model(&dst->at(sz0), qh, maxC, &modelLenBits, &enc);
  if (pInfo)
    pInfo[1] = modelLenBits;

  uint16_t c2low[257];
  double quantizedEntropy;
  prepare2(c2low, maxC, h, qh, &quantizedEntropy);

  int lenEst = (modelLenBits+quantizedEntropy+7)/8;
  if (pInfo)
    pInfo[3] = quantizedEntropy;
  if (lenEst >= srclen)
    return 0; // not compressible

  // printf("ml=%u\n", modellen);
  dst->resize(sz0 + lenEst + 64);
  int dstlen = encode(&dst->at(sz0+modellen), src, srclen, c2low, &enc);

  if (pInfo)
    pInfo[2] = dstlen*8.0;

  int reslen = modellen + dstlen;

  if (reslen >= srclen)
    return 0; // not compressible

  dst->resize(sz0 + reslen);
  return reslen;
#endif
  return 0;
}
