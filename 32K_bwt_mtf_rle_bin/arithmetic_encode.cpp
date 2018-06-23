#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <algorithm>

#include "arithmetic_encode.h"
#include "arithmetic_coder_ut.h"
#include "arithmetic_coder_cfg.h"
#include "fast_log2.h"

// static const int      QH_SCALE  = 1 << QH_BITS;
static const unsigned VAL_RANGE = 1u << RANGE_BITS;

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

static uint8_t  h2qh_tab[ARITH_CODER_CNT_MAX+1];
static uint16_t qh2h_tab[ARITH_CODER_QH_SCALE+1];
static double   log2_tab[ARITH_CODER_CNT_MAX+1];
static double   entropy_tab[ARITH_CODER_CNT_MAX+1];
static double   quantized_entropy_tab[ARITH_CODER_CNT_MAX+1];

void arithmetic_encode_init_tables()
{
  for (int i = 1; i <= ARITH_CODER_CNT_MAX; ++i)
    log2_tab[i] = log2(i);

  for (int i = 1; i < ARITH_CODER_CNT_MAX; ++i)
    entropy_tab[i] =
      log2_tab[ARITH_CODER_CNT_MAX]*ARITH_CODER_CNT_MAX
      - log2_tab[i]*i
      - log2_tab[ARITH_CODER_CNT_MAX-i]*(ARITH_CODER_CNT_MAX-i);

  h2qh_tab[ARITH_CODER_CNT_MAX] = ARITH_CODER_QH_SCALE;
  for (int i = 1; i < ARITH_CODER_CNT_MAX; ++i)
    h2qh_tab[i] = quantize_histogram_pair(i, ARITH_CODER_CNT_MAX-i, ARITH_CODER_QH_SCALE);

  qh2h_tab[ARITH_CODER_QH_SCALE] = VAL_RANGE;
  double M_PI = 3.1415926535897932384626433832795;
  for (int i = 1; i < ARITH_CODER_QH_SCALE; ++i)
    qh2h_tab[i] = int(round((sin((i*2-ARITH_CODER_QH_SCALE)*(M_PI/2/ARITH_CODER_QH_SCALE))+1.0)*(VAL_RANGE/2)));

  for (int i = 1; i < ARITH_CODER_CNT_MAX; ++i) {
    int32_t ra = qh2h_tab[h2qh_tab[i]];
    quantized_entropy_tab[i] =
      RANGE_BITS*ARITH_CODER_CNT_MAX
      - log2_tab[ra]*i
      - log2_tab[VAL_RANGE-ra]*(ARITH_CODER_CNT_MAX-i);
  }
}

enum {
  // context header
  CONTEXT_HDR_SRC_OFFSET_I=0,
  CONTEXT_HDR_SRC_LEN_I =0,
  CONTEXT_HDR_H_OFFSET_I,
  CONTEXT_HDR_H_LEN_I,
  CONTEXT_HDR_PREV_C0_I,
  CONTEXT_HDR_CH0_CNT_I,
  CONTEXT_HDR_NC_I,
  CONTEXT_HDR_CNTRS_I = CONTEXT_HDR_NC_I + 256,
  CONTEXT_CNTRS_SZ = (256*2*sizeof(uint8_t))/sizeof(uint32_t),
  CONTEXT_HDR_LEN  = CONTEXT_HDR_CNTRS_I + CONTEXT_CNTRS_SZ,
};

void arithmetic_encode_init_context(uint32_t* context, int tilelen)
{
  uint32_t srcOffset = CONTEXT_HDR_LEN;
  context[CONTEXT_HDR_SRC_OFFSET_I] = srcOffset;
  context[CONTEXT_HDR_SRC_LEN_I]    = 0;
  context[CONTEXT_HDR_H_OFFSET_I]   = srcOffset + (tilelen*2-1)/sizeof(uint32_t) + 1;
  context[CONTEXT_HDR_H_LEN_I]      = 0;
  context[CONTEXT_HDR_PREV_C0_I]    = 1;
  context[CONTEXT_HDR_CH0_CNT_I]    = 0;
  memset(&context[CONTEXT_HDR_NC_I], 0, sizeof(uint32_t)*256);
  memset(&context[CONTEXT_HDR_CNTRS_I], 0, sizeof(uint32_t)*CONTEXT_CNTRS_SZ);
}

void arithmetic_encode_chunk_callback(void* context_ptr, const uint8_t* src, int srclen)
{
  uint32_t* context = static_cast<uint32_t*>(context_ptr);
  uint8_t*  cntrs = reinterpret_cast<uint8_t*>(&context[CONTEXT_HDR_CNTRS_I]);
  uint32_t  hlen = context[CONTEXT_HDR_H_LEN_I];
  uint8_t*  hbeg = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_H_OFFSET_I]]);
  uint8_t*  h = &hbeg[hlen*2];

  int  prevC = context[CONTEXT_HDR_PREV_C0_I];
  uint32_t ch0Cnt = 0;

  // split source in two two levels and update histograms
  for (int i = 0; i < srclen; ++i) {
    int cx = src[i];
    if (cx == 255) {
      // escape sequence for symbols 255 and 256
      cx = int(src[i+1]) + 1;
      ++i;
    }
    unsigned c = cx == 0 ? 0 : cx - 1;
    if (c == 0) {
      if (prevC != 0) {
        cntrs[0*2+0] += 1;
        cntrs[0*2+1] += cx;
        if (cntrs[0*2+0] == ARITH_CODER_CNT_MAX) {
          h[0] = 0;
          h[1] = cntrs[0*2+1];
          h += 2;
          cntrs[0*2+0] = cntrs[0*2+1] = 0;
          ++context[CONTEXT_HDR_NC_I+0];
        }
      } else {
        ++ch0Cnt;
      }
    }
    prevC = c;
    for (unsigned bit = 256; bit != 1; ) {
      unsigned pos = c & (0-bit);
      bit /= 2;
      pos += bit;
      unsigned b = (c & bit) != 0;
      cntrs[pos*2+0] += 1;
      cntrs[pos*2+1] += b;
      if (cntrs[pos*2+0] == ARITH_CODER_CNT_MAX) {
        h[0] = pos;
        h[1] = cntrs[pos*2+1];
        h += 2;
        cntrs[pos*2+0] = cntrs[pos*2+1] = 0;
        ++context[CONTEXT_HDR_NC_I+pos];
      }
    }
  }
  context[CONTEXT_HDR_PREV_C0_I] = prevC;
  context[CONTEXT_HDR_CH0_CNT_I] += ch0Cnt;

  hbeg = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_H_OFFSET_I]]);
  context[CONTEXT_HDR_H_LEN_I] = (h-hbeg)/2;

  uint8_t* dst = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_SRC_OFFSET_I]]);
  uint32_t dstlen = context[CONTEXT_HDR_SRC_LEN_I];
  context[CONTEXT_HDR_SRC_LEN_I] = dstlen + srclen;
  memcpy(&dst[dstlen], src, srclen);
}

// return estimate of result length
static int prepare(uint32_t* context, uint32_t qhOffsets[256], double* pInfo)
{
  uint8_t*  cntrs = reinterpret_cast<uint8_t*>(&context[CONTEXT_HDR_CNTRS_I]);
  uint32_t offset = 0;
  for (int i = 0; i < 256; ++i) {
    qhOffsets[i] = offset;
    offset += context[CONTEXT_HDR_NC_I+i] + (cntrs[2*i+0] != 0);
  }
  uint32_t qho[256];
  memcpy(qho, qhOffsets, sizeof(qho));

  const uint32_t hlen = context[CONTEXT_HDR_H_LEN_I];
  const uint8_t* h = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_H_OFFSET_I]]);
  uint8_t* qh = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_H_OFFSET_I]]) + hlen*2;
  double entropy = context[CONTEXT_HDR_CH0_CNT_I];
  double quantizedEntropy = entropy;

  // process full sections
  for (uint32_t src_i = 0; src_i < hlen; ++src_i) {
    int pos = h[src_i*2+0];
    int val = h[src_i*2+1];
    uint32_t offs = qho[pos];
    qh[offs] = h2qh_tab[val];
    qho[pos] = offs + 1;
    entropy += entropy_tab[val];
    quantizedEntropy += quantized_entropy_tab[val];
  }

  // process partial sections
  for (int i = 0; i < 256; ++i) {
    unsigned tot = cntrs[2*i+0];
    if (tot != 0) {
      unsigned h0 = cntrs[2*i+1];
      unsigned qhVal = quantize_histogram_pair(h0, tot-h0, ARITH_CODER_QH_SCALE);
      qh[qho[i]] = qhVal;
      entropy +=
          log2_tab[tot]*tot
        - log2_tab[h0]*h0
        - log2_tab[tot-h0]*(tot-h0);
      int32_t ra0 = qh2h_tab[qhVal];
      quantizedEntropy +=
          RANGE_BITS*tot
        - log2_tab[ra0]*h0
        - log2_tab[VAL_RANGE-ra0]*(tot-h0);
    }
  }

  double modelLen = 0;
  for (int i = 0; i < 256; ++i) {
    const uint32_t nc = context[CONTEXT_HDR_NC_I+i] + (cntrs[2*i+0] != 0);
    const uint8_t* src = &qh[qhOffsets[i]];
    int prev = 0;
    for (uint32_t c = 0; c < nc; ++c) {
      modelLen += 1;
      int val = src[c];
      if (val != prev) {
        if (val > prev) {
          modelLen += (prev != 0);
          if (prev < ARITH_CODER_QH_SCALE-1) {
            modelLen += 1;
            if (val > prev + 1)
              modelLen += log2_tab[ARITH_CODER_QH_SCALE-1-prev];
          }
        } else {
          modelLen += (prev != ARITH_CODER_QH_SCALE);
          if (prev > 1) {
            modelLen += 1;
            if (val < prev - 1)
              modelLen += log2_tab[prev-1];
          }
        }
        prev = val;
      }
    }
  }

  if (pInfo) {
    pInfo[0] = entropy;
    pInfo[1] = modelLen;
    // pInfo[2] = reslen*8.0 - modelLenBits;;
    pInfo[3] = quantizedEntropy;
  }

  return ceil(quantizedEntropy+modelLen);
}

#if 0
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

static void prepare1(uint32_t* context, double* pInfo)
{
  context_plain_hdrs_t* hdrs = reinterpret_cast<context_plain_hdrs_t*>(&context[CONTEXT_HDR_PLAIN_HDRS_I]);
  for (int i = 0; i < 9; ++i)
    Adapt(context, &hdrs->a[i]);
  uint32_t qhOffset = context[CONTEXT_HDR_QH_OFFSET_I];
  double entropy = 0;
  uint8_t* qHistogram = reinterpret_cast<uint8_t*>(&context[qhOffset]);
  double entrArr[9];
  for (int i = 0; i < 9; ++i) {
    entropy += entrArr[i] = Prepare1Plain(context, pInfo, &hdrs->a[i], qHistogram);
    qHistogram += hdrs->a[i].nChunks * (hdrs->a[i].nSymbols+1);
  }
  if (pInfo) {
    pInfo[0] = entropy+context[CONTEXT_HDR_CH0_CNT_I];
    pInfo[3] = 0;
    for (int i = 0; i < 9; ++i) {
      pInfo[4+9*0+i] = hdrs->a[i].nChunks;
      pInfo[4+9*1+i] = hdrs->a[i].len;
      pInfo[4+9*2+i] = entrArr[i];
    }
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
  const int      qHistogramNCol = nSymbols==2 ? 2 : nSymbols + 1;
  uint32_t* buf = &context[hdr->headOffset];
  for (uint32_t chunk_i = 0; chunk_i < nChunks; ++chunk_i) {
    int hlen       = buf[CONTEXT_CHK_HLEN_I];
    int pgPerChunk = buf[CONTEXT_CHK_PG_PER_CHUNK_I];
    int nRanges = 0;
    // context_plain1_c2low_t* dstChunk = reinterpret_cast<context_plain1_c2low_t*>(dst);
    if (hlen > 0) {
      uint16_t ranges[258];
      if (nSymbols > 2) {
        nRanges = quantized_histogram_to_range(ranges, hlen, qHistogram+1, VAL_RANGE);
      } else {
        nRanges = 1;
        unsigned qhVal = qHistogram[1];
        if (qhVal!=0 && qhVal!=255) {
          unsigned ra0 = quantized_histogram_pair_to_range_qh_scale255(qhVal, VAL_RANGE);
          ranges[0] = ra0;
          ranges[1] = VAL_RANGE-ra0;
          nRanges = 2;
        }
      }
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
    qHistogram += qHistogramNCol;
  }
  return entropy;
}

// return entropy estimate after quantization
static double prepare2(uint32_t * context)
{
  context_plain_hdrs_t* hdrs = reinterpret_cast<context_plain_hdrs_t*>(&context[CONTEXT_HDR_PLAIN_HDRS_I]);
  double entropy = context[CONTEXT_HDR_CH0_CNT_I];
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

template<class T>
void store_model_store_data_loop(T* obj,
  const uint8_t*             qh,
  const context_plain_hdr_t* hdr0,
  const context_plain_hdr_t* hdr1)
{
  for (const context_plain_hdr_t* hdr = hdr0; hdr != hdr1; ++hdr) {
    uint32_t nChunks  = hdr->nChunks;
    if (nChunks > 0) {
      const uint32_t nSymbols = hdr->nSymbols;
      const uint32_t maxHlen  = hdr->maxHLen;

      const int len  = nSymbols==2 ? 2 : maxHlen+1;
      const int nCol = nSymbols==2 ? 2 : nSymbols+1;

      int runlen = (nCol-len)*nChunks-2;
      unsigned prev_val = 0;
      for (int i = 0; i < len; ++i) {
        const uint8_t* col = &qh[len-1-i];
        for (uint32_t chunk_i = 0; chunk_i < nChunks; ++chunk_i, col += nCol) {
          unsigned val = *col;
          ++runlen;
          if (val != prev_val) {
            if (runlen >= 0)
              obj->add_zero_run(runlen);
            runlen = -1;
            unsigned diff;
            obj->add_nonzero();
            if (val > prev_val) {
              obj->add_plus(prev_val != 0); // '+'
              diff = val - prev_val - 1;
            } else {
              obj->add_minus(prev_val != 255); // '-'
              diff = prev_val - val - 1;
            }
            if (diff != 0)
              obj->add_difference(diff-1);
            prev_val = val;
          }
        }
        prev_val = qh[len-1-i];
      }
      obj->add_zero_run(runlen+1);
      qh += (nSymbols+1)*nChunks;
    }
  }
}

static uint8_t* store_model_store_data(
  uint8_t*                   dst,
  const uint8_t*             qh,
  const context_plain_hdr_t* hdr0,
  const context_plain_hdr_t* hdr1,
  CArithmeticEncoder*        pEnc)
{
  // first pass - calculate histograms
  struct histogram_pass_t {
    unsigned hist[8];
    void add_zero_run  (unsigned val) { hist[0] += add_value_to_histogram(&hist[2], val); }
    void add_difference(unsigned val) { hist[1] += add_value_to_histogram(&hist[4], val); }
    void add_nonzero()       { hist[1] += 1;   }
    void add_plus (bool ena) { hist[6] += ena; }
    void add_minus(bool ena) { hist[7] += ena; }
  };
  histogram_pass_t pass1 = {{0}};
  store_model_store_data_loop(&pass1, qh, hdr0, hdr1);

#if 0
  double entr = 0;
  for (int i = 0; i < 4; ++i) {
    double h0 = pass1.hist[i*2+0], h1 = pass1.hist[i*2+1];
    if (h0 != 0 && h1 != 0) {
      entr += (h0+h1)*log2(h0+h1);
      entr -= h0*log2(h0);
      entr -= h1*log2(h1);
    }
  }
  printf("entr=%.2f\n", entr/8);
#endif

  unsigned hh01 = quantize_histogram_pair(pass1.hist[0], pass1.hist[1], 9); // range [1..8], because both sums > 0
  unsigned hh23 = quantize_histogram_pair(pass1.hist[2], pass1.hist[3], 9); // range [0..9]
  unsigned hh45 = quantize_histogram_pair(pass1.hist[4], pass1.hist[5], 9); // range [0..9]
  unsigned hh67 = quantize_histogram_pair(pass1.hist[6], pass1.hist[7], 9); // range [0..9]

  // prepare second pass
  struct encode_pass_t {
    uint8_t*            dst;
    CArithmeticEncoder* pEnc;
    unsigned range_tab[4][3];
    void add_zero_run  (unsigned val) { dst = encode_value(dst, val, &range_tab[0][0], range_tab[1], pEnc); }
    void add_difference(unsigned val) { dst = encode_value(dst, val, &range_tab[0][1], range_tab[2], pEnc); }
    void add_nonzero()       { dst = pEnc->put(VAL_RANGE, range_tab[0][1], VAL_RANGE-range_tab[0][1], dst); }
    void add_plus (bool ena) {
      if (ena && range_tab[3][1] != 0)
        dst = pEnc->put(VAL_RANGE, 0, range_tab[3][1], dst); // '+'
    }
    void add_minus(bool ena) {
      if (ena && range_tab[3][1] != VAL_RANGE)
        dst = pEnc->put(VAL_RANGE, range_tab[3][1], VAL_RANGE-range_tab[3][1], dst); // '-'
    }
  };
  encode_pass_t pass2;
  pass2.range_tab[0][1] = quantized_histogram_pair_to_range_qh_scale9(hh01, VAL_RANGE);
  pass2.range_tab[1][1] = quantized_histogram_pair_to_range_qh_scale9(hh23, VAL_RANGE);
  pass2.range_tab[2][1] = quantized_histogram_pair_to_range_qh_scale9(hh45, VAL_RANGE);
  pass2.range_tab[3][1] = quantized_histogram_pair_to_range_qh_scale9(hh67, VAL_RANGE);
  for (int k = 0; k < 4; ++k) {
    pass2.range_tab[k][0] = 0;
    pass2.range_tab[k][2] = VAL_RANGE;
    // printf("pass2.range_tab[%d][1] = %5d\n", k, pass2.range_tab[k][1]);
  }
  // printf("%.5f %.5f\n", double(pass1.hist[0]+pass1.hist[1])/(pass1.hist[0]+pass1.hist[1]+pass1.hist[2]+pass1.hist[3]), double(pass2.range_tab[0][1])/VAL_RANGE);
  // printf("%.5f %.5f\n", double(pass1.hist[0])/(pass1.hist[0]+pass1.hist[1]), double(pass2.range_tab[1][1])/VAL_RANGE);
  // printf("%.5f %.5f\n", double(pass1.hist[2])/(pass1.hist[2]+pass1.hist[3]), double(pass2.range_tab[2][1])/VAL_RANGE);

  unsigned hhw = (hh01-1)+8*(hh23+(10*(hh45+10*hh67))); // combine all hh in a single word

  // store hhw
  dst = pEnc->put(8*10*10*10, hhw, 1, dst);

  // second pass - encode
  pass2.dst  = dst;
  pass2.pEnc = pEnc;
  store_model_store_data_loop(&pass2, qh, hdr0, hdr1);

  return pass2.dst;
}

// return the number of stored octets
static int store_model(uint8_t* dst, uint32_t * context, double* pNbits, CArithmeticEncoder* pEnc)
{
  uint8_t* dst0 = dst;
  context_plain_hdrs_t* hdrs = reinterpret_cast<context_plain_hdrs_t*>(&context[CONTEXT_HDR_PLAIN_HDRS_I]);

  dst = store_model_store_nChunks(dst, hdrs->a[0].nChunks, pEnc);
  for (int i = 1; i < 9; ++i)
    dst = store_model_store_nChunks(dst, hdrs->a[i].nChunks+1, pEnc);

  uint8_t* qHistogram0 = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_QH_OFFSET_I]]);
  uint8_t* qHistogram = qHistogram0;
  context_plain_hdr_t* hdr0 = hdrs->a;
  context_plain_hdr_t* hdr  = hdr0;
  for (int i = 0; i < 9; ++i) {
    uint32_t nChunks  = hdr->nChunks;
    uint32_t nSymbols = hdr->nSymbols;
    uint32_t maxHlen  = hdr->maxHLen;
    if (nChunks > 0 && nSymbols>2) {
      uint32_t* plainH = &context[hdr->headOffset];
      for (uint32_t chunk_i = 0; chunk_i < nChunks; ++chunk_i, plainH = &context[plainH[CONTEXT_CHK_NEXT_I]]) {
        uint32_t hlen = plainH[CONTEXT_CHK_HLEN_I];
        if (hlen < maxHlen)
          memset(&qHistogram[(nSymbols+1)*chunk_i+hlen+1], 0, maxHlen-hlen);
      }
    }
    qHistogram += (nSymbols+1)*nChunks;
    hdr        += 1;
    int nOctets = qHistogram - qHistogram0;
    if (nOctets > 0) {
      if (i == 8 || nOctets+(hdr->nSymbols+1)*hdr->nChunks > ARITH_CODER_MODEL_MAX_SEGMENT_SZ) {
        dst = store_model_store_data(dst, qHistogram0, hdr0, hdr, pEnc);
        qHistogram0 = qHistogram;
        hdr0        = hdr;
      }
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
  static const uint16_t c2low_half[3] = {0, VAL_RANGE/2, VAL_RANGE};
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
  int prevC0 = 1;
  for (uint32_t lvl1chunk_i = 0; lvl1chunk_i < lvl1nChunks; ++lvl1chunk_i, p_p1_c2low = &context[p_p1_c2low[CONTEXT_CHK_NEXT_I]]) {
    const uint16_t* lvl1_c2low = p_p1_c2low[CONTEXT_CHK_N_RANGES_I] > 1 ?
      reinterpret_cast<const uint16_t*>(&p_p1_c2low[CONTEXT_CHK_C2LOW_I]) : 0;
    unsigned srclen1 = (lvl1chunk_i == lvl1nChunks-1) ? lvl1Len : p_p1_c2low[CONTEXT_CHK_SYMBOLS_PER_CHUNK_I];
    lvl1Len -= srclen1;
    for (unsigned i = 0; i < srclen1; ++i) {
      const uint16_t* c2low = lvl1_c2low;
      lvl2Rec_t* lvl2Rec = 0;
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
          c2low   = c2low_half;
          if ((c | prevC0) != 0) {
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
          prevC0 = c;
        }
      }
      src += 2;
      if (lvl2Rec) {
        ++lvl2_i;
        if (lvl2_i == lvl2Rec->srclen) {
          --lvl2Rec->nChunks;
          lvl2Rec->offset = context[lvl2Rec->offset+CONTEXT_CHK_NEXT_I];
          lvl2_i = 0;
        }
        lvl2Rec->lvl2_i = lvl2_i;
      }
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
#endif

static int encode(uint8_t* dst, const uint32_t* context, uint32_t qhOffsets[256])
{
  return 0;
}


// return value:
//  0 - not compressible, because all input characters have approximately equal probability or because input is too short
// >0 - the length of compressed buffer
int arithmetic_encode(uint32_t* context, uint8_t* dst, int origlen, double* pInfo)
{
  uint32_t qhOffsets[256];
  int estlen = prepare(context, qhOffsets, pInfo);
  if (estlen >= origlen)
    return 0; // not compressible

  int reslen = encode(dst, context, qhOffsets);

  if (pInfo) {
    double modelLenBits = pInfo[1];
    pInfo[2] = reslen*8.0 - modelLenBits;
  }

  return reslen;
#if 0
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

  int reslen = modellen + dstlen;
  // printf("%d+%d=%d(%d) <> %d\n", modellen, dstlen, reslen, lenEst, origlen);

  if (pInfo)
    pInfo[2] = reslen*8.0 - modelLenBits;

  if (reslen >= origlen)
    return 0; // not compressible

  return reslen;
#endif
}

/*
    pInfo[0] = entropy;
    pInfo[1] = modelLenBits;
    pInfo[2] = reslen*8.0 - modelLenBits;;
    pInfo[3] = quantizedEntropy;
    pInfo[4] = hdrs->a[0].nChunks;
    pInfo[5] = hdrs->a[1].nChunks;
*/
