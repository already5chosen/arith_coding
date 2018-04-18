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

static const int N_COMMON_SYMBOLS = 257 + 1 - ARITH_CODER_N_DYNAMIC_SYMBOLS;
struct context_common_c2low_t {
  uint16_t c2low[N_COMMON_SYMBOLS+1];
};

struct context_chunk_c2low_t {
  uint32_t nRanges;
  uint16_t c2low[ARITH_CODER_N_DYNAMIC_SYMBOLS+1];
};

enum {
  // context header
  CONTEXT_HDR_SRC_OFFSET_I=0,
  CONTEXT_HDR_HIST_OFFSET_I,
  CONTEXT_HDR_PLAIN1_LEN_I,
  CONTEXT_HDR_PLAIN2_LEN_I,
  CONTEXT_HDR_NCHUNKS_I,
  CONTEXT_HDR_COMMON_HLEN_I,
  CONTEXT_HDR_COMMON_NRANGES,
  CONTEXT_HDR_CHUNK_MAX_HLEN_I,
  CONTEXT_HDR_LEN,
  // context histogram chunk
  CONTEXT_CHK_HLEN_I = 0,
  CONTEXT_CHK_HISTOGRAM_I,
  CONTEXT_CHUNK_H_LEN = CONTEXT_CHK_HISTOGRAM_I + ARITH_CODER_N_DYNAMIC_SYMBOLS,
  // context quantized histogram
  CONTEXT_COMMON_QH_LEN = N_COMMON_SYMBOLS,
  CONTEXT_CHUNK_QH_LEN  = ARITH_CODER_N_DYNAMIC_SYMBOLS,
  // context quantized c2low
  CONTEXT_COMMON_C2LOW_SZ = (sizeof(context_common_c2low_t)-1)/sizeof(uint32_t) + 1,
  CONTEXT_CHUNK_C2LOW_SZ  = (sizeof(context_chunk_c2low_t)-1)/sizeof(uint32_t) + 1,
};

void arithmetic_encode_init_context(uint32_t* context, int tilelen)
{
  uint32_t srcOffset  = CONTEXT_HDR_LEN;
  uint32_t histOffset = srcOffset + (tilelen*2-1)/sizeof(uint32_t) + 1;
  context[CONTEXT_HDR_SRC_OFFSET_I]    = srcOffset;
  context[CONTEXT_HDR_HIST_OFFSET_I]   = histOffset;
  context[CONTEXT_HDR_PLAIN1_LEN_I]    = 0;
  context[CONTEXT_HDR_PLAIN2_LEN_I]    = 0;
  uint32_t* commonH  = &context[histOffset];
  memset(commonH, 0, N_COMMON_SYMBOLS*sizeof(uint32_t));
}

void arithmetic_encode_chunk_callback(void* context_ptr, const uint8_t* src, int srclen)
{
  uint32_t* context = static_cast<uint32_t*>(context_ptr);
  uint8_t*  dst = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_SRC_OFFSET_I]]);
  uint32_t* commonH  = &context[context[CONTEXT_HDR_HIST_OFFSET_I]];
  uint32_t* chunkH   = &commonH[N_COMMON_SYMBOLS];
  uint32_t plain1Len = context[CONTEXT_HDR_PLAIN1_LEN_I];
  uint32_t plain2Len = context[CONTEXT_HDR_PLAIN2_LEN_I];
  dst += plain1Len + plain2Len;

  uint32_t chunk_i = plain1Len / ARITH_CODER_SYMBOLS_PER_CHUNK;
  uint32_t sym_i   = plain1Len % ARITH_CODER_SYMBOLS_PER_CHUNK;
  chunkH += CONTEXT_CHUNK_H_LEN*chunk_i+CONTEXT_CHK_HISTOGRAM_I;

  // split source in two plains and update histograms
  for (int i = 0; i < srclen; ++i) {
    int c = src[i];
    if (c == 255) {
      // escape sequence for symbols 255 and 256
      c = int(src[i+1]) + 1;
      ++i;
    }

    dst[0] = c;
    if (c >= ARITH_CODER_N_DYNAMIC_SYMBOLS-1) {
      int c2 = c - (ARITH_CODER_N_DYNAMIC_SYMBOLS-1);
      ++commonH[c2];
      c = ARITH_CODER_N_DYNAMIC_SYMBOLS-1;
      dst[0] = c;
      dst[1] = c2;
      ++dst;
      ++plain2Len;
    }
    ++dst;

    if (sym_i == 0)
      memset(chunkH, 0, ARITH_CODER_N_DYNAMIC_SYMBOLS*sizeof(uint32_t)); // prepare new chunk of dynamic histogram

    ++plain1Len;
    ++chunkH[c];
    ++sym_i;
    if (sym_i == ARITH_CODER_SYMBOLS_PER_CHUNK) {
      sym_i = 0;
      chunkH += CONTEXT_CHUNK_H_LEN;
    }
  }

  context[CONTEXT_HDR_PLAIN1_LEN_I] = plain1Len;
  context[CONTEXT_HDR_PLAIN2_LEN_I] = plain2Len;
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

static void prepare1(uint32_t* context, double* pInfo)
{
  uint32_t* commonH  = &context[context[CONTEXT_HDR_HIST_OFFSET_I]];
  uint32_t* chunkH   = &commonH[N_COMMON_SYMBOLS];
  uint32_t plain1Len = context[CONTEXT_HDR_PLAIN1_LEN_I];
  uint32_t nChunks = plain1Len / ARITH_CODER_SYMBOLS_PER_CHUNK; // full chunks
  uint32_t symbolsInLastChunk = plain1Len % ARITH_CODER_SYMBOLS_PER_CHUNK;
  if (symbolsInLastChunk > 0) {
    // partial chunk exists
    if (symbolsInLastChunk < ARITH_CODER_SYMBOLS_PER_CHUNK/2 && nChunks > 0) {
      // add last chunk to before-last
      uint32_t* dstChunk = &chunkH[CONTEXT_CHUNK_H_LEN*(nChunks-1)];
      uint32_t* srcChunk = &chunkH[CONTEXT_CHUNK_H_LEN*(nChunks-0)];
      for (int i = 0; i < ARITH_CODER_N_DYNAMIC_SYMBOLS; ++i) {
        dstChunk[CONTEXT_CHK_HISTOGRAM_I+i] += srcChunk[CONTEXT_CHK_HISTOGRAM_I+i];
      }
    } else {
      nChunks += 1;
    }
  }
  context[CONTEXT_HDR_NCHUNKS_I] = nChunks;

  double entropy = 0;
  uint8_t* qHistogram = reinterpret_cast<uint8_t*>(&chunkH[CONTEXT_CHUNK_H_LEN*nChunks]);
  {
    uint32_t* histogram = commonH;
    // find highest-numbered character that occurred at least once
    int maxC = -1;
    int32_t tot = 0;
    for (int c = 0; c < CONTEXT_COMMON_QH_LEN; ++c) {
      tot += histogram[c];
      if (histogram[c] != 0)
        maxC = c;
    }
    context[CONTEXT_HDR_COMMON_HLEN_I] = maxC + 1;
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
    uint32_t* chunk = &chunkH[CONTEXT_CHUNK_H_LEN*chunk_i];
    uint32_t* histogram = &chunk[CONTEXT_CHK_HISTOGRAM_I];
    // find highest-numbered character that occurred at least once
    int maxC = -1;
    int32_t tot = 0;
    for (int c = 0; c < CONTEXT_CHUNK_QH_LEN; ++c) {
      tot += histogram[c];
      if (histogram[c] != 0)
        maxC = c;
    }
    chunk[CONTEXT_CHK_HLEN_I] = maxC + 1;
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
  context[CONTEXT_HDR_CHUNK_MAX_HLEN_I] = maxMaxC+1;

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
  uint32_t* commonH = &context[context[CONTEXT_HDR_HIST_OFFSET_I]];
  uint32_t* chunkH  = &commonH[N_COMMON_SYMBOLS];

  // calculate common c2low table
  double entropy = 0;
  uint8_t* qHistogram = reinterpret_cast<uint8_t*>(&chunkH[CONTEXT_CHUNK_H_LEN*nChunks]);
  uint32_t* dst = commonH;
  {
    int hlen = context[CONTEXT_HDR_COMMON_HLEN_I];
    int nRanges = 0;
    if (hlen > 0) {
      uint16_t ranges[CONTEXT_COMMON_QH_LEN];
      nRanges = quantized_histogram_to_range(ranges, hlen, qHistogram, VAL_RANGE);
      if (nRanges > 1) {
        // calculate entropy after quantization
        for (int c = 0; c < hlen; ++c) {
          unsigned cnt = commonH[c];
          if (cnt)
            entropy -= log2(ranges[c]/double(VAL_RANGE))*cnt;
        }
        context_common_c2low_t* dstCommon = reinterpret_cast<context_common_c2low_t*>(dst);
        range2low(dstCommon->c2low, ranges, hlen);
      }
    }
    context[CONTEXT_HDR_COMMON_NRANGES] = nRanges;
  }
  qHistogram += CONTEXT_COMMON_QH_LEN;
  dst        += CONTEXT_COMMON_C2LOW_SZ;

  // calculate c2low tables for dynamic chunks
  for (int chunk_i = 0; chunk_i < nChunks; ++chunk_i) {
    uint32_t* srcChunk = &chunkH[CONTEXT_CHUNK_H_LEN*chunk_i];
    int hlen = srcChunk[CONTEXT_CHK_HLEN_I];
    int nRanges = 0;
    context_chunk_c2low_t* dstChunk = reinterpret_cast<context_chunk_c2low_t*>(dst);
    if (hlen > 0) {
      uint16_t ranges[CONTEXT_CHUNK_QH_LEN];
      nRanges = quantized_histogram_to_range(ranges, hlen, qHistogram, VAL_RANGE);
      if (nRanges > 1) {
        // calculate entropy after quantization
        for (int c = 0; c < hlen; ++c) {
          unsigned cnt = srcChunk[CONTEXT_CHK_HISTOGRAM_I+c];
          if (cnt) {
            entropy -= log2(ranges[c]/double(VAL_RANGE))*cnt;
            ++nRanges;
          }
        }
        // printf("%10d\n", int(entropy/8));
        range2low(dstChunk->c2low, ranges, hlen);
      }
    }
    dstChunk->nRanges = nRanges;

    qHistogram += CONTEXT_CHUNK_QH_LEN;
    dst        += CONTEXT_CHUNK_C2LOW_SZ;
  }
  return entropy;
}

class CArithmeticEncoder {
public:
  uint8_t* put(uint64_t cScale, uint64_t cLo, uint64_t cRange, uint8_t* dst);
  void init() {
    m_lo = 0;
    m_range = uint64_t(1) << 63;
  }
  uint64_t m_lo;    // scaled by 2**64
  uint64_t m_range; // scaled by 2**63
};

// return quantized code of h0/(h0+h1) in range [0..scale]
static unsigned quantize_histogram_pair(unsigned h0, unsigned h1, unsigned scale)
{
  if (h0 == 0)
    return 0;
  if (h1 == 0)
    return scale;
  unsigned hTot = h0 + h1;
  // printf("%u+%u => %.2f\n", h0, h1, (log2(hTot)*hTot-log2(h0)*h0-log2(h1)*h1)/8);
  double M_PI = 3.1415926535897932384626433832795;
  unsigned val = round((asin(h0*2.0/hTot - 1.0)/M_PI + 0.5)*scale);
  if (val < 1)
    val = 1;
  else if (val >= scale)
    val = scale-1;
  return val;
}

static unsigned add_value_to_histogram(unsigned* h, unsigned val)
{
  val += 2;
  unsigned nbits = 0;
  do {
    ++nbits;
    ++h[val % 2];
    val /= 2;
  } while (val > 1);
  return nbits;
}

static uint8_t* encode_value(uint8_t* dst, unsigned val, const unsigned* range_tab0, const unsigned* range_tab1, CArithmeticEncoder* pEnc)
{
  unsigned lo0 = range_tab0[0];
  unsigned ra0 = range_tab0[1] - lo0;
  val += 2;
  do {
    dst = pEnc->put(VAL_RANGE, lo0, ra0, dst);
    unsigned bit = val % 2;
    unsigned lo1 = range_tab1[bit];
    unsigned ra1 = range_tab1[bit+1] - lo1;
    if (ra1 != 0 && ra1 != VAL_RANGE)
      dst = pEnc->put(VAL_RANGE, lo1, ra1, dst);
    val /= 2;
  } while (val > 1);
  return dst;
}

static uint8_t* store_model_store_nChunks(uint8_t* dst, int val, CArithmeticEncoder* pEnc)
{
  // use predefined value for sake of simplicity
  static const unsigned encTab[] = { 0, 4, 7, 8 };
  while (val > 1) {
    unsigned bit = val % 2;
    dst = pEnc->put(8, encTab[bit], encTab[bit+1]-encTab[bit], dst);
    val /= 2;
  }
  dst = pEnc->put(8, encTab[2], encTab[2+1]-encTab[2], dst);
  return dst;
}

static uint8_t* store_model_store_data(uint8_t* dst, const uint8_t* qh, int len, int nCol, int nChunks, CArithmeticEncoder* pEnc)
{
  if (len == 0)
    return pEnc->put(8*10*10*10+1, 0, 1, dst); // store special code instead of hhw

  // first pass - calculate histograms
  unsigned hist[8] = {0};
  int runlen = (nCol-len)*nChunks-2;
  unsigned prev_val = 0;
  for (int i = 0; i < len; ++i) {
    const uint8_t* col = &qh[len-1-i];
    for (int chunk_i = 0; chunk_i < nChunks; ++chunk_i, col += nCol) {
      unsigned val = *col;
      ++runlen;
      if (val != prev_val) {
        if (runlen >= 0)
          hist[0] += add_value_to_histogram(&hist[2], runlen);
        runlen = -1;
        unsigned diff;
        hist[1] += 1;
        if (val > prev_val) {
          hist[6] += (prev_val != 0); // '+'
          diff = val - prev_val - 1;
        } else {
          hist[7] += (prev_val != 255); // '-'
          diff = prev_val - val - 1;
        }
        if (diff != 0)
          hist[1] += add_value_to_histogram(&hist[4], diff-1);
        prev_val = val;
      }
    }
    prev_val = qh[len-1-i];
  }
  hist[0] += add_value_to_histogram(&hist[2], runlen+1);

#if 0
  double entr = 0;
  for (int i = 0; i < 4; ++i) {
    double h0 = hist[i*2+0], h1 = hist[i*2+1];
    if (h0 != 0 && h1 != 0) {
      entr += (h0+h1)*log2(h0+h1);
      entr -= h0*log2(h0);
      entr -= h1*log2(h1);
    }
  }
  printf("entr=%.2f\n", entr/8);
#endif

  unsigned hh01 = quantize_histogram_pair(hist[0], hist[1], 9); // range [1..8], because both sums > 0
  unsigned hh23 = quantize_histogram_pair(hist[2], hist[3], 9); // range [0..9]
  unsigned hh45 = quantize_histogram_pair(hist[4], hist[5], 9); // range [0..9]
  unsigned hh67 = quantize_histogram_pair(hist[6], hist[7], 9); // range [0..9]

  unsigned range_tab[4][3];
  range_tab[0][1] = quantized_histogram_pair_to_range_qh_scale9(hh01, VAL_RANGE);
  range_tab[1][1] = quantized_histogram_pair_to_range_qh_scale9(hh23, VAL_RANGE);
  range_tab[2][1] = quantized_histogram_pair_to_range_qh_scale9(hh45, VAL_RANGE);
  range_tab[3][1] = quantized_histogram_pair_to_range_qh_scale9(hh67, VAL_RANGE);
  for (int k = 0; k < 4; ++k) {
    range_tab[k][0] = 0;
    range_tab[k][2] = VAL_RANGE;
    // printf("range_tab[%d][1] = %5d\n", k, range_tab[k][1]);
  }
  // printf("%.5f %.5f\n", double(hist[0]+hist[1])/(hist[0]+hist[1]+hist[2]+hist[3]), double(range_tab[0][1])/VAL_RANGE);
  // printf("%.5f %.5f\n", double(hist[0])/(hist[0]+hist[1]), double(range_tab[1][1])/VAL_RANGE);
  // printf("%.5f %.5f\n", double(hist[2])/(hist[2]+hist[3]), double(range_tab[2][1])/VAL_RANGE);

  unsigned hhw = (hh01-1)+8*(hh23+(10*(hh45+10*hh67))); // combine all hh in a single word

  // store hhw
  dst = pEnc->put(8*10*10*10+1, hhw+1, 1, dst);

  // second pass - encode
  runlen = (nCol-len)*nChunks-2;
  prev_val = 0;
  for (int i = 0; i < len; ++i) {
    const uint8_t* col = &qh[len-1-i];
    for (int chunk_i = 0; chunk_i < nChunks; ++chunk_i, col += nCol) {
      unsigned val = *col;
      ++runlen;
      if (val != prev_val) {
        if (runlen >= 0)
          dst = encode_value(dst, runlen, &range_tab[0][0], range_tab[1], pEnc);
        runlen = -1;
        unsigned diff;
        dst = pEnc->put(VAL_RANGE, range_tab[0][1], VAL_RANGE-range_tab[0][1], dst);
        if (val > prev_val) {
          if (prev_val != 0 && range_tab[3][1] != 0)
            dst = pEnc->put(VAL_RANGE, 0, range_tab[3][1], dst); // '+'
          diff = val - prev_val - 1;
        } else {
          if (prev_val != 255 && range_tab[3][1] != VAL_RANGE)
            dst = pEnc->put(VAL_RANGE, range_tab[3][1], VAL_RANGE-range_tab[3][1], dst); // '-'
          diff = prev_val - val - 1;
        }
        if (diff != 0)
          dst = encode_value(dst, diff-1, &range_tab[0][1], range_tab[2], pEnc);
        prev_val = val;
      }
    }
    prev_val = qh[len-1-i];
  }
  dst = encode_value(dst, runlen+1, &range_tab[0][0], range_tab[1], pEnc);

  return dst;
}

// return the number of stored octets
static int store_model(uint8_t* dst, uint32_t * context, double* pNbits, CArithmeticEncoder* pEnc)
{
  uint8_t* dst0 = dst;
  int nChunks = context[CONTEXT_HDR_NCHUNKS_I];
  dst = store_model_store_nChunks(dst, nChunks, pEnc);

  uint32_t* commonH  = &context[context[CONTEXT_HDR_HIST_OFFSET_I]];
  uint32_t* chunkH   = &commonH[N_COMMON_SYMBOLS];
  uint8_t* qHistogram = reinterpret_cast<uint8_t*>(&chunkH[CONTEXT_CHUNK_H_LEN*nChunks]);
  int commonHlen = context[CONTEXT_HDR_COMMON_HLEN_I];
  dst = store_model_store_data(dst, qHistogram, commonHlen, CONTEXT_COMMON_QH_LEN, 1, pEnc);
  qHistogram += CONTEXT_COMMON_QH_LEN;

  int chunkMaxHlen = context[CONTEXT_HDR_CHUNK_MAX_HLEN_I];
  if (chunkMaxHlen > 0) {
    for (int chunk_i = 0; chunk_i < nChunks; ++chunk_i) {
      int hlen = chunkH[CONTEXT_CHUNK_H_LEN*chunk_i+CONTEXT_CHK_HLEN_I];
      if (hlen < chunkMaxHlen)
        memset(&qHistogram[CONTEXT_CHUNK_QH_LEN*chunk_i+hlen], 0, chunkMaxHlen-hlen);
    }
  }
  dst = store_model_store_data(dst, qHistogram, chunkMaxHlen, CONTEXT_CHUNK_QH_LEN, nChunks, pEnc);

  int len = dst - dst0;
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

static int encode(uint8_t* dst, const uint32_t* context, CArithmeticEncoder* pEnc)
{
  const uint8_t*  src = reinterpret_cast<const uint8_t*>(&context[context[CONTEXT_HDR_SRC_OFFSET_I]]);
  const uint32_t* commonH  = &context[context[CONTEXT_HDR_HIST_OFFSET_I]];
  uint32_t plain1Len = context[CONTEXT_HDR_PLAIN1_LEN_I];
  uint32_t nChunks = context[CONTEXT_HDR_NCHUNKS_I];

  const uint64_t MSB_MSK   = uint64_t(255) << 56;
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);
  uint64_t lo     = pEnc->m_lo;                      // scaled by 2**64
  uint64_t range  = pEnc->m_range >> (RANGE_BITS-1); // scaled by 2**(64-RANGE_BITS)
  uint8_t* dst0   = dst;
  uint64_t prevLo = lo;
  // int dbg_i = 0;
  const uint16_t* common_c2low = context[CONTEXT_HDR_COMMON_NRANGES] > 1 ?
    reinterpret_cast<const context_common_c2low_t*>(commonH)->c2low :
    0;
  for (uint32_t chunk_i = 0; chunk_i < nChunks; ++chunk_i) {
    const context_chunk_c2low_t* chunk = reinterpret_cast<const context_chunk_c2low_t*>(&commonH
      [CONTEXT_COMMON_C2LOW_SZ+CONTEXT_CHUNK_C2LOW_SZ*chunk_i]);
    const uint16_t* chunk_c2low = chunk->nRanges > 1 ? chunk->c2low : 0;
    unsigned srclen = (chunk_i == nChunks-1) ?
      plain1Len - chunk_i*ARITH_CODER_SYMBOLS_PER_CHUNK:
      ARITH_CODER_SYMBOLS_PER_CHUNK;
    // int dbg_i = -1;
    for (unsigned i = 0; i < srclen; ++i) {
      int common_c = -1;
      do {
        // ++dbg_i;
        const uint16_t* c2low = common_c2low;
        int c = common_c;
        common_c = -1;
        if (c < 0) {
          // process next input character
          c = *src++;
          c2low = chunk_c2low;
          if (c == ARITH_CODER_N_DYNAMIC_SYMBOLS-1) {
            common_c = *src++;
          }
        }
        // if (chunk_i==77 && dbg_i >= 0 && dbg_i < 500)
          // printf("%d:%3d %3d\n", dbg_i, common_c, c);
        // printf("[%3d]=%3d\n", dbg_i, c); ++dbg_i;
        if (c2low) {
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
      } while (common_c >= 0);
    }
    // printf(": %10d\n", dst-dst0);
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
  uint64_t range = m_range / cScale;
  uint64_t lo = m_lo + range * cLo * 2;
  range *= cRange;

  if (lo < m_lo) // lo overflow
    inc_dst(dst);

  if (range <= MIN_RANGE) {
    // re-normalize
    dst[0] = lo >> (64-8*1);
    dst[1] = lo >> (64-8*2);
    dst[2] = lo >> (64-8*3);
    dst += 3;
    lo    <<= 24;
    range <<= 24;
  }
  m_lo    = lo;
  m_range = range;
  return dst;
}

// return value:
//  0 - not compressible, because all input characters have approximately equal probability or because input is too short
// >0 - the length of compressed buffer
int arithmetic_encode(uint32_t* context, uint8_t* dst, int origlen, double* pInfo)
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

  int dstlen = encode(&dst[modellen], context, &enc);

  if (pInfo)
    pInfo[2] = dstlen*8.0;

  int reslen = modellen + dstlen;
  // printf("%d+%d=%d(%d) <> %d\n", modellen, dstlen, reslen, lenEst, origlen);

  if (reslen >= origlen)
    return 0; // not compressible

  return reslen;
}
