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

static const int N_COMMON_SYMBOLS = 257 + 1 - ARITH_CODER_N_DYNAMIC_SYMBOLS;
struct context_common_c2low_t {
  uint16_t c2low[N_COMMON_SYMBOLS+1];
};

struct context_chunk_c2low_t {
  uint32_t srclen;
  uint16_t c2low[ARITH_CODER_N_DYNAMIC_SYMBOLS+1];
};

enum {
  // context header
  CONTEXT_HDR_NCHUNKS_I = 0,
  CONTEXT_HDR_LASTRUN_I,
  CONTEXT_HDR_COMMON_MAXC_I,
  CONTEXT_HDR_MAXMAXC_I,
  CONTEXT_HDR_COMMON_HISTOGRAM_I,
  CONTEXT_HDR_LEN = CONTEXT_HDR_COMMON_HISTOGRAM_I+N_COMMON_SYMBOLS,
  // context histogram chunk
  CONTEXT_CHK_MAXC_I = 0,
  CONTEXT_CHK_SRCLEN_I,
  CONTEXT_CHK_HISTOGRAM_I,
  CONTEXT_CHUNK_H_LEN = CONTEXT_CHK_HISTOGRAM_I + ARITH_CODER_N_DYNAMIC_SYMBOLS,
  // context quantized histogram
  CONTEXT_COMMON_QH_LEN = N_COMMON_SYMBOLS,
  CONTEXT_CHUNK_QH_LEN  = ARITH_CODER_N_DYNAMIC_SYMBOLS,
  // context quantized c2low
  CONTEXT_COMMON_C2LOW_SZ = (sizeof(context_common_c2low_t)-1)/sizeof(uint32_t) + 1,
  CONTEXT_CHUNK_C2LOW_SZ  = (sizeof(context_chunk_c2low_t)-1)/sizeof(uint32_t) + 1,
};

static inline size_t GetContextLength(int nChunks) {
  return
    CONTEXT_HDR_LEN+
    CONTEXT_CHUNK_H_LEN*nChunks+
    CONTEXT_COMMON_QH_LEN+
    (CONTEXT_CHUNK_QH_LEN*nChunks+sizeof(uint32_t)-1)/sizeof(uint32_t);
}


void arithmetic_encode_init_context(std::vector<uint32_t>* context, int tilelen)
{
  int nChunks = tilelen/(ARITH_CODER_RUNS_PER_CHUNK*2) + 1; // estimate # of chunks
  resize_up(context, GetContextLength(nChunks));
  (*context)[CONTEXT_HDR_NCHUNKS_I] = 0;
  (*context)[CONTEXT_HDR_LASTRUN_I] = 0;
}

int arithmetic_encode_chunk_callback(void* context, const uint8_t* src, int chunklen, int nRuns)
{
  if (src) {
    std::vector<uint32_t>* vec = static_cast<std::vector<uint32_t>*>(context);
    uint32_t chunk_i = (*vec)[CONTEXT_HDR_NCHUNKS_I];
    resize_up(vec, GetContextLength(chunk_i+1));
    (*vec)[CONTEXT_HDR_NCHUNKS_I] = chunk_i + 1;
    (*vec)[CONTEXT_HDR_LASTRUN_I] = nRuns;

    // calculate histogram
    uint32_t histogram[257]={0};
    for (int i = 0; i < chunklen; ++i) {
      int c = src[i];
      if (c == 255) {
        // escape sequence for symbols 255 and 256
        c = int(src[i+1]) + 1;
        ++i;
      }
      ++histogram[c];
    }

    // add common part of histogram to common histogram in the header
    uint32_t  commonTot=0;
    uint32_t* commonHistogram = &vec->at(CONTEXT_HDR_COMMON_HISTOGRAM_I);
    if (chunk_i == 0) {
      for (int i = 0; i < N_COMMON_SYMBOLS; ++i) {
        uint32_t val = histogram[ARITH_CODER_N_DYNAMIC_SYMBOLS-1+i];
        commonHistogram[i] = val;
        commonTot += val;
      }
    } else {
      for (int i = 0; i < N_COMMON_SYMBOLS; ++i) {
        uint32_t val = histogram[ARITH_CODER_N_DYNAMIC_SYMBOLS-1+i];
        commonHistogram[i] += val;
        commonTot += val;
      }
    }

    // store dynamic (per-chunk) part of histogram in its slot
    uint32_t* chunk = &vec->at(CONTEXT_HDR_LEN+CONTEXT_CHUNK_H_LEN*chunk_i);
    chunk[CONTEXT_CHK_SRCLEN_I] = chunklen;
    uint32_t* dynHistogram = &chunk[CONTEXT_CHK_HISTOGRAM_I];
    dynHistogram[ARITH_CODER_N_DYNAMIC_SYMBOLS-1] = commonTot;
    memcpy(dynHistogram, histogram, sizeof(*dynHistogram)*(ARITH_CODER_N_DYNAMIC_SYMBOLS-1));
  }
  return ARITH_CODER_RUNS_PER_CHUNK;
}

static void quantize_histogram(uint8_t* __restrict qh, const uint32_t *h, int len, uint32_t hTot)
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

static void prepare1(uint32_t * context, double* pInfo)
{
  int nChunks = context[CONTEXT_HDR_NCHUNKS_I];
  int runsInLastChunk = context[CONTEXT_HDR_LASTRUN_I];
  if (runsInLastChunk < ARITH_CODER_RUNS_PER_CHUNK/2 && nChunks > 1) {
    // add last chunk to before-last
    uint32_t* dstChunk = &context[CONTEXT_HDR_LEN+CONTEXT_CHUNK_H_LEN*(nChunks-2)];
    uint32_t* srcChunk = &context[CONTEXT_HDR_LEN+CONTEXT_CHUNK_H_LEN*(nChunks-1)];
    dstChunk[CONTEXT_CHK_SRCLEN_I] += srcChunk[CONTEXT_CHK_SRCLEN_I];
    for (int i = 0; i < ARITH_CODER_N_DYNAMIC_SYMBOLS; ++i) {
      dstChunk[CONTEXT_CHK_HISTOGRAM_I+i] += srcChunk[CONTEXT_CHK_HISTOGRAM_I+i];
    }
    nChunks -= 1;
    context[CONTEXT_HDR_NCHUNKS_I] = nChunks;
  }

  double entropy = 0;
  uint8_t* qHistogram = reinterpret_cast<uint8_t*>(&context[CONTEXT_HDR_LEN+CONTEXT_CHUNK_H_LEN*nChunks]);
  {
    uint32_t* histogram = &context[CONTEXT_HDR_COMMON_HISTOGRAM_I];
    // find highest-numbered character that occurred at least once
    int maxC = -1;
    int32_t tot = 0;
    for (int c = 0; c < CONTEXT_COMMON_QH_LEN; ++c) {
      tot += histogram[c];
      if (histogram[c] != 0)
        maxC = c;
    }
    context[CONTEXT_HDR_COMMON_MAXC_I] = maxC;
    if (maxC >= 0) {
      if (pInfo) {
        // calculate source entropy of static part of histogram
        entropy += log2(double(tot))*tot;
        for (int c = 0; c <= maxC; ++c) {
          int32_t cnt = histogram[c];
          if (cnt)
            entropy -= log2(double(cnt))*cnt;
        }
      }
      quantize_histogram(qHistogram, histogram, maxC+1, tot);
    }
  }
  qHistogram += CONTEXT_COMMON_QH_LEN;

  int maxMaxC = -1;
  for (int chunk_i = 0; chunk_i < nChunks; ++chunk_i) {
    uint32_t* chunk = &context[CONTEXT_HDR_LEN+CONTEXT_CHUNK_H_LEN*chunk_i];
    uint32_t* histogram = &chunk[CONTEXT_CHK_HISTOGRAM_I];
    // find highest-numbered character that occurred at least once
    int maxC = -1;
    int32_t tot = 0;
    for (int c = 0; c < CONTEXT_CHUNK_QH_LEN; ++c) {
      tot += histogram[c];
      if (histogram[c] != 0)
        maxC = c;
    }
    chunk[CONTEXT_CHK_MAXC_I] = maxC;
    if (maxC > maxMaxC)
      maxMaxC = maxC;

    if (maxC >= 0) {
      if (pInfo) {
        // calculate source entropy of static part of histogram
        entropy += log2(double(tot))*tot;
        for (int c = 0; c <= maxC; ++c) {
          int32_t cnt = histogram[c];
          if (cnt)
            entropy -= log2(double(cnt))*cnt;
        }
      }
      quantize_histogram(qHistogram, histogram, maxC+1, tot);
    }
    qHistogram += CONTEXT_CHUNK_QH_LEN;
  }
  context[CONTEXT_HDR_MAXMAXC_I] = maxMaxC;

  if (pInfo) {
    pInfo[0] = entropy;
    pInfo[3] = 0;
    pInfo[4] = nChunks;
  }
}


static void range2low(uint16_t* c2low, const uint16_t* c2range, unsigned len)
{
  unsigned lo = 0;
  for (unsigned c = 0; c < len; ++c) {
    c2low[c] = lo;
    lo += c2range[c];
  }
  c2low[len] = VAL_RANGE;
}

// return entropy estimate after quantization
static double prepare2(uint32_t * context)
{
  int nChunks = context[CONTEXT_HDR_NCHUNKS_I];
  double entropy = 0;
  uint8_t* qHistogram = reinterpret_cast<uint8_t*>(&context[CONTEXT_HDR_LEN+CONTEXT_CHUNK_H_LEN*nChunks]);
  uint32_t* dst = &context[CONTEXT_HDR_COMMON_HISTOGRAM_I];
  {
    int maxC = int32_t(context[CONTEXT_HDR_COMMON_MAXC_I]);
    if (maxC >= 0) {
      uint16_t ranges[CONTEXT_COMMON_QH_LEN];
      quantized_histogram_to_range(ranges, maxC+1, qHistogram, VAL_RANGE);
      // calculate entropy after quantization
      for (int c = 0; c <= maxC; ++c) {
        unsigned cnt = context[CONTEXT_HDR_COMMON_HISTOGRAM_I+c];
        if (cnt)
          entropy -= log2(ranges[c]/double(VAL_RANGE))*cnt;
      }
      context_common_c2low_t* dstCommon = reinterpret_cast<context_common_c2low_t*>(dst);
      range2low(dstCommon->c2low, ranges, maxC+1);
    }
  }
  qHistogram += CONTEXT_COMMON_QH_LEN;
  dst        += CONTEXT_COMMON_C2LOW_SZ;

  for (int chunk_i = 0; chunk_i < nChunks; ++chunk_i) {
    uint32_t* srcChunk = &context[CONTEXT_HDR_LEN+CONTEXT_CHUNK_H_LEN*chunk_i];
    int maxC = int32_t(srcChunk[CONTEXT_CHK_MAXC_I]);
    uint32_t srclen = srcChunk[CONTEXT_CHK_SRCLEN_I];
    context_chunk_c2low_t* dstChunk = reinterpret_cast<context_chunk_c2low_t*>(dst);
    if (maxC >= 0) {
      uint16_t ranges[CONTEXT_CHUNK_QH_LEN];
      quantized_histogram_to_range(ranges, maxC+1, qHistogram, VAL_RANGE);
      // calculate entropy after quantization
      for (int c = 0; c <= maxC; ++c) {
        unsigned cnt = srcChunk[CONTEXT_CHK_HISTOGRAM_I+c];
        if (cnt)
          entropy -= log2(ranges[c]/double(VAL_RANGE))*cnt;
      }
      // printf("%10d\n", int(entropy/8));
      range2low(dstChunk->c2low, ranges, maxC+1);
    }
    dstChunk->srclen = srclen;

    qHistogram += CONTEXT_CHUNK_QH_LEN;
    dst        += CONTEXT_CHUNK_C2LOW_SZ;
  }
  return entropy;
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
  prepare1(context, pInfo);

  CArithmeticEncoder enc;
  enc.init();
  double modelLenBits;
  unsigned modellen = store_model(dst, context, &modelLenBits, &enc);
  if (pInfo)
    pInfo[1] = modelLenBits;

  double quantizedEntropy = prepare2(context);

  int lenEst = (modelLenBits+quantizedEntropy+7)/8;
  if (pInfo)
    pInfo[3] = quantizedEntropy;
  if (lenEst >= origlen)
    return 0; // not compressible

  int dstlen = encode(&dst[modellen], src, context, &enc);

  if (pInfo)
    pInfo[2] = dstlen*8.0;

  int reslen = modellen + dstlen;
  // printf("%d+%d=%d(%d) <> %d\n", modellen, dstlen, reslen, lenEst, origlen);

  if (reslen >= origlen)
    return 0; // not compressible

  return reslen;
}
