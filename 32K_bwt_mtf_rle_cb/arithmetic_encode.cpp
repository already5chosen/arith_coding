#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <algorithm>

#include "arithmetic_encode.h"
#include "arithmetic_coder_ut.h"
#include "arithmetic_coder_cfg.h"
#include "fast_log2.h"

static const int      QH_SCALE  = 1 << QH_BITS;
static const unsigned VAL_RANGE = 1u << RANGE_BITS;

struct context_plain_hdr_t {
  uint32_t headOffset;
  uint32_t tailOffset;
  uint32_t len;
  uint32_t nChunks;
  uint32_t maxHLen;
  uint32_t nSymbols;
  uint32_t pageSz;
  uint32_t contextSz;
};

struct context_plain_hdrs_t {
  context_plain_hdr_t a[9];
};

enum {
  // context header
  CONTEXT_HDR_SRC_OFFSET_I=0,
  CONTEXT_HDR_QH_OFFSET_I,
  CONTEXT_HDR_PLAIN_HDRS_I,
  CONTEXT_PLAINS_HDRS_SZ = (sizeof(context_plain_hdrs_t)-1)/sizeof(uint32_t) + 1,
  CONTEXT_HDR_LEN = CONTEXT_HDR_PLAIN_HDRS_I + CONTEXT_PLAINS_HDRS_SZ,
  // context plain histogram chunk
  CONTEXT_CHK_NEXT_I = 0,
  CONTEXT_CHK_PREV_I,
  CONTEXT_CHK_HLEN_I = CONTEXT_CHK_PREV_I, // PREV is used during adaptation, HLEN is used after adaptation
  CONTEXT_CHK_PG_PER_CHUNK_I,
  CONTEXT_CHK_HISTOGRAM_I,
  // context plain c2low chunk (overlaps histogram chunk)
  CONTEXT_CHK_N_RANGES_I = CONTEXT_CHK_HLEN_I,
  CONTEXT_CHK_SYMBOLS_PER_CHUNK_I,
  CONTEXT_CHK_C2LOW_I,
};

void arithmetic_encode_init_context(uint32_t* context, int tilelen)
{
  uint32_t srcOffset = CONTEXT_HDR_LEN;
  context[CONTEXT_HDR_SRC_OFFSET_I] = srcOffset;
  context[CONTEXT_HDR_QH_OFFSET_I]  = srcOffset + (tilelen*2-1)/sizeof(uint32_t) + 1;
  context_plain_hdrs_t* hdrs = reinterpret_cast<context_plain_hdrs_t*>(&context[CONTEXT_HDR_PLAIN_HDRS_I]);

  hdrs->a[0].nSymbols  = 8;
  hdrs->a[1].nSymbols  = 2;
  hdrs->a[2].nSymbols  = 2;
  hdrs->a[3].nSymbols  = 5;
  hdrs->a[4].nSymbols  = 8;
  hdrs->a[5].nSymbols  = 16;
  hdrs->a[6].nSymbols  = 32;
  hdrs->a[7].nSymbols  = 64;
  hdrs->a[8].nSymbols  = 128;

  hdrs->a[0].pageSz    = ARITH_CODER_L1_PAGE_SZ;
  hdrs->a[1].pageSz    = ARITH_CODER_L2_0_PAGE_SZ;
  hdrs->a[2].pageSz    = ARITH_CODER_L2_1_PAGE_SZ;
  hdrs->a[3].pageSz    = ARITH_CODER_L2_3_PAGE_SZ;
  hdrs->a[4].pageSz    = ARITH_CODER_L2_8_PAGE_SZ;
  hdrs->a[5].pageSz    = ARITH_CODER_L2_16_PAGE_SZ;
  hdrs->a[6].pageSz    = ARITH_CODER_L2_32_PAGE_SZ;
  hdrs->a[7].pageSz    = ARITH_CODER_L2_64_PAGE_SZ;
  hdrs->a[8].pageSz    = ARITH_CODER_L2_128_PAGE_SZ;

  for (int i = 0; i < 9; ++i) {
    hdrs->a[i].headOffset = 0;
    hdrs->a[i].tailOffset = 0;
    hdrs->a[i].len = 0;
    uint32_t nSymbols = hdrs->a[i].nSymbols;
    uint32_t contextSz = sizeof(uint32_t) > sizeof(uint16_t) ?
      CONTEXT_CHK_HISTOGRAM_I + nSymbols :
      CONTEXT_CHK_C2LOW_I + nSymbols + 1;
    hdrs->a[i].contextSz = contextSz;
  }
}

static uint32_t* allocate_chunk(uint32_t* context, int hdr_i)
{
  context_plain_hdrs_t* hdrs = reinterpret_cast<context_plain_hdrs_t*>(&context[CONTEXT_HDR_PLAIN_HDRS_I]);
  context_plain_hdr_t* hdr = &hdrs->a[hdr_i];
  uint32_t offs = context[CONTEXT_HDR_QH_OFFSET_I];
  context[CONTEXT_HDR_QH_OFFSET_I] = offs + hdr->contextSz;
  uint32_t* p = &context[offs];
  uint32_t tailOffset = hdr->tailOffset;
  hdr->tailOffset = offs;
  p[CONTEXT_CHK_PREV_I] = tailOffset;
  p[CONTEXT_CHK_NEXT_I] = 0;
  if (tailOffset != 0) {
    context[tailOffset+CONTEXT_CHK_NEXT_I] = offs;
  } else {
    hdr->headOffset = offs;
  }
  memset(&p[CONTEXT_CHK_HISTOGRAM_I], 0, hdr->nSymbols*sizeof(uint32_t)); // prepare new chunk of histogram
  return p;
}

void arithmetic_encode_chunk_callback(void* context_ptr, const uint8_t* src, int srclen)
{
  uint32_t* context = static_cast<uint32_t*>(context_ptr);
  uint8_t*  dst = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_SRC_OFFSET_I]]);
  context_plain_hdrs_t* hdrs = reinterpret_cast<context_plain_hdrs_t*>(&context[CONTEXT_HDR_PLAIN_HDRS_I]);

  unsigned lvl2_sym_i_arr[8];
  for (int i = 0; i < 8; ++i)
    lvl2_sym_i_arr[i] = hdrs->a[i+1].len % hdrs->a[i+1].pageSz;

  uint32_t lvl1Len = hdrs->a[0].len;
  dst += lvl1Len*2;

  uint32_t lvl1_sym_i = lvl1Len % ARITH_CODER_L1_PAGE_SZ;
  uint32_t* lvl1H = &context[hdrs->a[0].tailOffset+CONTEXT_CHK_HISTOGRAM_I];

  // split source in two two levels and update histograms
  for (int i = 0; i < srclen; ++i) {
    int c = src[i];
    if (c == 255) {
      // escape sequence for symbols 255 and 256
      c = int(src[i+1]) + 1;
      ++i;
    }

    static const uint8_t class_tab[32][2] = {
      {0, 0},  {0, 0},
      {1, 2},  {1, 2},
      {2, 4},  {2, 4}, {2, 4}, {2, 4}, {2, 4},
      {3, 9},  {3, 9}, {3, 9}, {3, 9},
      {3, 9},  {3, 9}, {3, 9}, {3, 9},
      {4, 17},
      {5, 33},  {5, 33},
      {6, 65},  {6, 65},  {6, 65},  {6, 65},
      {7, 129}, {7, 129}, {7, 129}, {7, 129},
      {7, 129}, {7, 129}, {7, 129}, {7, 129},
    };
    int tab_i = c < 17 ? c : (c + 17*15) >> 4;
    int c0 = class_tab[tab_i][0];
    int c1 = c - class_tab[tab_i][1];

    dst[0] = c0;
    dst[1] = c1;
    dst += 2;

    if (lvl1_sym_i == 0)
      lvl1H = allocate_chunk(context, 0) + CONTEXT_CHK_HISTOGRAM_I;

    ++lvl1Len;
    ++lvl1H[c0];
    ++lvl1_sym_i;
    if (lvl1_sym_i == ARITH_CODER_L1_PAGE_SZ)
      lvl1_sym_i = 0;

    context_plain_hdr_t* hdr = &hdrs->a[c0+1];
    unsigned lvl2_sym_i = lvl2_sym_i_arr[c0];
    uint32_t* lvl2H = &context[hdr->tailOffset];
    if (lvl2_sym_i == 0)
      lvl2H = allocate_chunk(context, c0+1);
    ++lvl2H[CONTEXT_CHK_HISTOGRAM_I+c1];
    ++hdr->len;
    ++lvl2_sym_i;
    if (lvl2_sym_i == hdr->pageSz)
      lvl2_sym_i = 0;
    lvl2_sym_i_arr[c0] = lvl2_sym_i;
  }
  hdrs->a[0].len = lvl1Len;
}

static void quantize_histogram(uint8_t* __restrict qh, const uint32_t *h, int len, uint32_t hTot)
{
  double invhLen = 2.0/hTot;
  double M_PI = 3.1415926535897932384626433832795;
  for (int c = 0; c < len; ++c) {
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

static void AdaptAdd(uint32_t* dst, const uint32_t* prev, unsigned len)
{
  for (unsigned c = 0; c < len; ++c)
    dst[c] += prev[c];
}

static double AdaptAddAndCalculateEntropy(uint32_t* hAcc, const uint32_t* h, unsigned  hlen)
{
  double entropy=0;
  uint32_t tot = 0;
  for (unsigned c = 0; c < hlen; ++c) {
    uint32_t val = h[c] + hAcc[c];
    hAcc[c] = val;
    if (val != 0) {
      tot += val;
      entropy -= fast_log2(val)*val;
    }
  }
  if (tot > 0) entropy += fast_log2(tot)*tot;
  return entropy;
}

static int EntrArrIdx(int iBeg, int iEnd) {
  int ri = iEnd - iBeg;
  int rowLen = 257-ri;
  return ((257+rowLen+1)*ri)/2 + (iBeg % rowLen);
}

static void Adapt(uint32_t* context, context_plain_hdr_t* hdr)
{
  uint32_t* chnk0 = &context[hdr->headOffset];
  const uint32_t plainLen = hdr->len;
  const uint32_t pageSz = hdr->pageSz;
  const int      nPages = (plainLen+pageSz-1)/pageSz;
  hdr->nChunks = nPages;
  if (nPages > 1) {
    chnk0[CONTEXT_CHK_PG_PER_CHUNK_I] = 0;
    unsigned nSymbols = hdr->nSymbols;
    int i0 = 0;
    const double DICT_LEN = (nSymbols+1)*5.5; // estimate for a size of the dictionary (model)

    float entrArr[(257*258)/2]; // entropy estimates for sections of pages
    // 2D matrix,
    // Each row contains entropies for given number of pages, starting with 1 page and up to 257
    // A column corresponds to particular starting index of the page, modulo # of columns in the row.
    // In order to minimize storage, only upper-left triangle+diagonal is stored

    uint32_t* chnk2 = chnk0;
    for (int i2 = 0; i2 < nPages; ++i2, chnk2 = &context[chnk2[CONTEXT_CHK_NEXT_I]]) {
      uint32_t accH[256]={0};
      // calculate and store entropy for all sections that start at [i0..i2] and end at page i2
      int entrArrRowBegI = 0;
      uint32_t* chnk1 = chnk2;
      for (int i1 = i2; i1 >= i0; --i1, chnk1 = &context[chnk1[CONTEXT_CHK_PREV_I]]) {
        double e = AdaptAddAndCalculateEntropy(accH, &chnk1[CONTEXT_CHK_HISTOGRAM_I], nSymbols);
        int rowLen = 257-(i2-i1);
        entrArr[entrArrRowBegI + (i1 % rowLen)] = static_cast<float>(e);
        entrArrRowBegI += rowLen;
      }

      chnk2[CONTEXT_CHK_PG_PER_CHUNK_I] = 0;
      while (i0 < i2) {
        double ee2_min  = 1e100;
        int i1_min = i0;
        // find the best split of pages [i0..i2]
        for (int i1 = i0; i1 < i2; ++i1) {
          double e1 = entrArr[EntrArrIdx(i0,i1)];   // entropy of segment[i0..i1]
          double e2 = entrArr[EntrArrIdx(i1+1,i2)]; // entropy of segment[i1+1..i2]
          double ee = e1 + e2;
          // printf("%d %d %d %d %f %f\n", i0, i1, i2, nSymbols, e1, e2);
          if (ee < ee2_min) {
            ee2_min = ee;
            i1_min  = i1;
          }
        }
        bool split = true;
        if (i2-i0 < 256) {
          double ee1 = entrArr[EntrArrIdx(i0,i2)];
          // printf("%d %d %d %d %f %f * %d\n", i0, i1_min, i2, nSymbols, ee1, ee2_min + DICT_LEN, nPages);
          split = ee2_min + DICT_LEN < ee1;
        }
        if (split) {
          // printf("%d %d %d %d * %d\n", i0, i1_min, i2, nSymbols, nPages);
          // segment ends at page i1_min
          while (i0 != i1_min) {
            chnk0 = &context[chnk0[CONTEXT_CHK_NEXT_I]];
            ++i0;
          }
          chnk0[CONTEXT_CHK_PG_PER_CHUNK_I] = 1; // mark
          chnk0 = &context[chnk0[CONTEXT_CHK_NEXT_I]];
          ++i0;
        } else {
          // the best split is not as good as no split
          break;
        }
      }
    }
    context[hdr->tailOffset+CONTEXT_CHK_PG_PER_CHUNK_I] = 1; // mark last page
    // printf("%d: %f\n", nSymbols, ee_sum/8);

    // database compaction
    int nSeg = 0;
    unsigned pgPerSeg = 0;
    uint32_t* pSrc = &context[hdr->headOffset];
    uint32_t* pDst = pSrc;
    for (int i = 0; i < nPages; ++i, pSrc = &context[pSrc[CONTEXT_CHK_NEXT_I]]) {
      if (pgPerSeg == 0) {
        if (pSrc != pDst)
          memcpy(&pDst[CONTEXT_CHK_HISTOGRAM_I], &pSrc[CONTEXT_CHK_HISTOGRAM_I], sizeof(pDst[0])*nSymbols);
      } else {
        AdaptAdd(&pDst[CONTEXT_CHK_HISTOGRAM_I], &pSrc[CONTEXT_CHK_HISTOGRAM_I], nSymbols);
      }
      ++pgPerSeg;
      if (pSrc[CONTEXT_CHK_PG_PER_CHUNK_I] != 0) {
        pDst[CONTEXT_CHK_PG_PER_CHUNK_I] = pgPerSeg;
        pDst = &context[pDst[CONTEXT_CHK_NEXT_I]];
        ++nSeg;
        pgPerSeg = 0;
      }
    }
    hdr->nChunks = nSeg;
  } else if (nPages > 0) {
    chnk0[CONTEXT_CHK_PG_PER_CHUNK_I] = 1;
  }
}

static double Prepare1Plain(uint32_t* context, double* pInfo, context_plain_hdr_t* hdr, uint8_t* qHistogram)
{
  double entropy = 0;
  int maxMaxC = -1;
  const uint32_t nChunks = hdr->nChunks;
  const int      nSymbols = hdr->nSymbols;
  uint32_t* chunk = &context[hdr->headOffset];
  for (unsigned chunk_i = 0; chunk_i < nChunks; ++chunk_i, chunk = &context[chunk[CONTEXT_CHK_NEXT_I]]) {
    qHistogram[0] = chunk[CONTEXT_CHK_PG_PER_CHUNK_I]-1;
    uint32_t* histogram = &chunk[CONTEXT_CHK_HISTOGRAM_I];
    // find highest-numbered character that occurred at least once
    int maxC = -1;
    int32_t tot = 0;
    for (int c = 0; c < nSymbols; ++c) {
      tot += histogram[c];
      if (histogram[c] != 0)
        maxC = c;
    }
    chunk[CONTEXT_CHK_HLEN_I] = maxC + 1;
    if (maxC > maxMaxC)
      maxMaxC = maxC;

    if (maxC >= 0) {
      if (pInfo) {
        // calculate source entropy of histogram
        entropy += log2(double(tot))*tot;
        for (int c = 0; c <= maxC; ++c) {
          int32_t cnt = histogram[c];
          if (cnt)
            entropy -= log2(double(cnt))*cnt;
        }
      }
      quantize_histogram(qHistogram+1, histogram, maxC+1, tot);
    }
    qHistogram += nSymbols+1;
  }
  hdr->maxHLen = maxMaxC+1;

  return entropy;
}

static void prepare1(uint32_t* context, double* pInfo)
{
  context_plain_hdrs_t* hdrs = reinterpret_cast<context_plain_hdrs_t*>(&context[CONTEXT_HDR_PLAIN_HDRS_I]);
  for (int i = 0; i < 9; ++i)
    Adapt(context, &hdrs->a[i]);
  uint32_t qhOffset = context[CONTEXT_HDR_QH_OFFSET_I];
  double entropy = 0;
  uint8_t* qHistogram = reinterpret_cast<uint8_t*>(&context[qhOffset]);
  for (int i = 0; i < 9; ++i) {
    entropy += Prepare1Plain(context, pInfo, &hdrs->a[i], qHistogram);
    qHistogram += hdrs->a[i].nChunks * (hdrs->a[i].nSymbols+1);
  }
  if (pInfo) {
    pInfo[0] = entropy;
    pInfo[3] = 0;
    pInfo[4] = hdrs->a[0].nChunks;
    pInfo[5] = hdrs->a[1].nChunks;
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
static double Prepare2Plain(uint32_t* context, context_plain_hdr_t* hdr, const uint8_t* qHistogram)
{
  double entropy = 0;
  // calculate c2low tables
  const uint32_t nChunks   = hdr->nChunks;
  const uint32_t nSymbols  = hdr->nSymbols;
  uint32_t* buf = &context[hdr->headOffset];
  for (uint32_t chunk_i = 0; chunk_i < nChunks; ++chunk_i) {
    int hlen       = buf[CONTEXT_CHK_HLEN_I];
    int pgPerChunk = buf[CONTEXT_CHK_PG_PER_CHUNK_I];
    int nRanges = 0;
    // context_plain1_c2low_t* dstChunk = reinterpret_cast<context_plain1_c2low_t*>(dst);
    if (hlen > 0) {
      uint16_t ranges[258];
      nRanges = quantized_histogram_to_range(ranges, hlen, qHistogram+1, VAL_RANGE);
      if (nRanges > 1) {
        // calculate entropy after quantization
        int32_t tot = 0;
        for (int c = 0; c < hlen; ++c) {
          unsigned cnt = buf[CONTEXT_CHK_HISTOGRAM_I+c];
          if (cnt)
            entropy -= fast_log2(ranges[c])*cnt;
          tot += cnt;
        }
        entropy += RANGE_BITS*tot;
        // printf("%10d\n", int(entropy/8));
        range2low(reinterpret_cast<uint16_t*>(&buf[CONTEXT_CHK_C2LOW_I]), ranges, hlen);
      }
    }
    buf[CONTEXT_CHK_N_RANGES_I] = nRanges;
    buf[CONTEXT_CHK_SYMBOLS_PER_CHUNK_I] = pgPerChunk*hdr->pageSz;

    buf = &context[buf[CONTEXT_CHK_NEXT_I]];
    qHistogram += nSymbols+1;
  }
  return entropy;
}

// return entropy estimate after quantization
static double prepare2(uint32_t * context)
{
  context_plain_hdrs_t* hdrs = reinterpret_cast<context_plain_hdrs_t*>(&context[CONTEXT_HDR_PLAIN_HDRS_I]);
  double entropy = 0;
  const uint8_t* qHistogram = reinterpret_cast<const uint8_t*>(&context[context[CONTEXT_HDR_QH_OFFSET_I]]);
  for (int i = 0; i < 9; ++i) {
    entropy += Prepare2Plain(context, &hdrs->a[i], qHistogram);
    qHistogram += hdrs->a[i].nChunks * (hdrs->a[i].nSymbols+1);
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
  context_plain_hdrs_t* hdrs = reinterpret_cast<context_plain_hdrs_t*>(&context[CONTEXT_HDR_PLAIN_HDRS_I]);

  dst = store_model_store_nChunks(dst, hdrs->a[0].nChunks, pEnc);
  for (int i = 1; i < 9; ++i)
    dst = store_model_store_nChunks(dst, hdrs->a[i].nChunks+1, pEnc);

  uint8_t* qHistogram = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_QH_OFFSET_I]]);
  for (int i = 0; i < 9; ++i) {
    context_plain_hdr_t* hdr = &hdrs->a[i];
    uint32_t nChunks  = hdr->nChunks;
    uint32_t nSymbols = hdr->nSymbols;
    uint32_t maxHlen  = hdr->maxHLen;
    if (nChunks > 0) {
      if (maxHlen > 0) {
        uint32_t* plainH = &context[hdr->headOffset];
        for (uint32_t chunk_i = 0; chunk_i < nChunks; ++chunk_i, plainH = &context[plainH[CONTEXT_CHK_NEXT_I]]) {
          uint32_t hlen = plainH[CONTEXT_CHK_HLEN_I];
          if (hlen < maxHlen)
            memset(&qHistogram[(nSymbols+1)*chunk_i+hlen+1], 0, maxHlen-hlen);
        }
      }
      dst = store_model_store_data(dst, qHistogram, maxHlen+1, nSymbols+1, nChunks, pEnc);
      qHistogram += (nSymbols+1)*nChunks;
    }
  }

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
  const context_plain_hdrs_t* hdrs = reinterpret_cast<const context_plain_hdrs_t*>(&context[CONTEXT_HDR_PLAIN_HDRS_I]);
  struct lvl2Rec_t {
    const uint16_t* c2low;
    uint32_t        lvl2_i;
    uint32_t        offset;
    uint32_t        srclen;
    uint32_t        nChunks;
  };
  lvl2Rec_t lvl2Records[8];
  for (int i = 0; i < 8; ++i) {
    lvl2Records[i].lvl2_i = 0;
    lvl2Records[i].offset = hdrs->a[i+1].headOffset;
    lvl2Records[i].nChunks = hdrs->a[i+1].nChunks;
  }

  uint32_t lvl1Len     = hdrs->a[0].len;
  uint32_t lvl1nChunks = hdrs->a[0].nChunks;

  const uint64_t MSB_MSK   = uint64_t(255) << 56;
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);
  uint64_t lo     = pEnc->m_lo;                      // scaled by 2**64
  uint64_t range  = pEnc->m_range >> (RANGE_BITS-1); // scaled by 2**(64-RANGE_BITS)
  uint8_t* dst0   = dst;
  uint64_t prevLo = lo;

  const uint32_t* p_p1_c2low = &context[hdrs->a[0].headOffset];
  for (uint32_t lvl1chunk_i = 0; lvl1chunk_i < lvl1nChunks; ++lvl1chunk_i, p_p1_c2low = &context[p_p1_c2low[CONTEXT_CHK_NEXT_I]]) {
    const uint16_t* lvl1_c2low = p_p1_c2low[CONTEXT_CHK_N_RANGES_I] > 1 ?
      reinterpret_cast<const uint16_t*>(&p_p1_c2low[CONTEXT_CHK_C2LOW_I]) : 0;
    unsigned srclen1 = (lvl1chunk_i == lvl1nChunks-1) ? lvl1Len : p_p1_c2low[CONTEXT_CHK_SYMBOLS_PER_CHUNK_I];
    lvl1Len -= srclen1;
    for (unsigned i = 0; i < srclen1; ++i) {
      const uint16_t* c2low = lvl1_c2low;
      lvl2Rec_t* lvl2Rec;
      uint32_t lvl2_i;
      for (int lvl = 0; lvl < 2; ++lvl) {
        int c = src[lvl];
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

        if (lvl==0) {
          lvl2Rec = &lvl2Records[c];
          c2low   = lvl2Rec->c2low;
          lvl2_i  = lvl2Rec->lvl2_i;
          if (lvl2_i == 0) {
            const uint32_t* c2lowRec = &context[lvl2Rec->offset];
            lvl2Rec->c2low =
            c2low = c2lowRec[CONTEXT_CHK_N_RANGES_I] > 1 ?
              reinterpret_cast<const uint16_t*>(&c2lowRec[CONTEXT_CHK_C2LOW_I]) : 0;
            lvl2Rec->srclen = (lvl2Rec->nChunks==1) ? UINT_MAX: c2lowRec[CONTEXT_CHK_SYMBOLS_PER_CHUNK_I];
          }
        }
      }
      src += 2;
      ++lvl2_i;
      if (lvl2_i == lvl2Rec->srclen) {
        --lvl2Rec->nChunks;
        lvl2Rec->offset = context[lvl2Rec->offset+CONTEXT_CHK_NEXT_I];
        lvl2_i = 0;
      }
      lvl2Rec->lvl2_i = lvl2_i;
    }
    // printf(": %10d\n", dst-dst0);
    // printf("%016I64x %d %d %d %d\n", range, p_p1_c2low->nRanges, plain2chunk_i, plain2_i, uu);
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
