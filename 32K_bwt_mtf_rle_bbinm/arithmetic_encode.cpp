#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <algorithm>

#include "arithmetic_encode.h"
#include "arithmetic_coder_cfg.h"

static const unsigned VAL_RANGE = 1u << RANGE_BITS;
static const unsigned QQH_SCALE = 15;               // not necessarily the same as ARITH_CODER_QH_SCALE
static const unsigned MODEL_QH_LEN = 4;

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
static uint8_t  qhMtf0[ARITH_CODER_QH_SCALE+1];
static uint32_t model_qh_tail_tab[ARITH_CODER_QH_SCALE+1];
static uint32_t model_qh2h_tab[QQH_SCALE] = {
  0,
    46927670,   185659716,   410132882,   710536612,
  1073741824,  1483874706,  1923010482,  2371956814,
  2811092590,  3221225472,  3584430684,  3884834414,
  4109307580,  4248039626,
};

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
      - log2(ra)*i
      - log2(VAL_RANGE-ra)*(ARITH_CODER_CNT_MAX-i);
    // printf("%3d %2d %5d %.3f\n", i, h2qh_tab[i], ra, quantized_entropy_tab[i]);
  }

  uint64_t val = uint64_t(1) << 32;
  uint64_t sum = val;
  for (int i = MODEL_QH_LEN+1; i <= ARITH_CODER_QH_SCALE; ++i) {
    val -= val / ARITH_CODER_QH_DECEY;
    sum += val;
    model_qh_tail_tab[i] = (val << 32)/sum;
  }

  for (int i = 0; i <= ARITH_CODER_QH_SCALE/2; ++i)
    qhMtf0[i] = ARITH_CODER_QH_SCALE/2 - i;
  for (int i = ARITH_CODER_QH_SCALE/2+1; i <= ARITH_CODER_QH_SCALE; ++i)
    qhMtf0[i] = i;
}

struct encode_prm_t {
  uint8_t* qh;
  unsigned mid0;
  unsigned midFactors[2];
  uint8_t  bisectionTab[256][2];
  uint8_t  modelQh[MODEL_QH_LEN];
  uint16_t modelC2loTab[ARITH_CODER_QH_SCALE+2];
};

// Select MidFactor that generates the smallest number of bits
// The purpose is not the densest encoding, but the fastest decoding
static unsigned SelectBestMidFactor(const uint32_t* srcHistogram, unsigned lo, unsigned hi)
{
  uint32_t min_nBits = uint32_t(-1);
  uint32_t min_midFactor = 0;
  for (unsigned midFactor = 16; midFactor < 32; ++midFactor) {
    uint32_t nBits = 0;
    for (unsigned val = lo; val <= hi; ++val) {
      unsigned tLo = lo, tHi = hi;
      unsigned nb = 0;
      while (tLo != tHi) {
        unsigned tMid = (tLo*midFactor + tHi*(32-midFactor))/32;
        if (val > tMid)
          tLo = tMid + 1;
        else
          tHi = tMid;
        ++nb;
      }
      nBits += nb*srcHistogram[val];
    }
    if (nBits < min_nBits) {
      min_nBits = nBits;
      min_midFactor = midFactor;
    }
  }
  return min_midFactor;
}

static int BuildBisectionTab(encode_prm_t* prm, unsigned lo, unsigned hi, unsigned midFactor)
{
  if (lo != hi) {
    unsigned mid = (lo*midFactor + hi*(32-midFactor))/32;
    prm->bisectionTab[mid][0] = BuildBisectionTab(prm, lo, mid, midFactor);
    prm->bisectionTab[mid][1] = BuildBisectionTab(prm, mid+1, hi, midFactor);
    return mid;
  }
  return 1;
}

static int qh_mtf_encode(uint8_t t[], int c)
{ // move-to-front encoder
  int v1 = t[0];
  int k = 0;
  if (c != v1) {
    // c is not at front
    t[0] = c;
    int v0 = v1;
    for (k = 1; c != (v1=t[k]); ++k) {
      t[k] = v0;
      v0 = v1;
    }
    t[k] = v0;
  }
  return k;
}

static void build_qh_encode_table(uint16_t c2lo[ARITH_CODER_QH_SCALE+2], const uint8_t modelQh[MODEL_QH_LEN])
{
  unsigned ranges[ARITH_CODER_QH_SCALE+1]={0};
  unsigned remVal = VAL_RANGE;
  for (int i = 0; i < MODEL_QH_LEN; ++i) {
    uint32_t qhVal = modelQh[i];
    if (qhVal != QQH_SCALE) {
      unsigned ra = (uint64_t(model_qh2h_tab[qhVal]) * remVal + (uint32_t(1)<<31)) >> 32;
      if (remVal-ra < ARITH_CODER_QH_SCALE-i)
        ra = remVal - (ARITH_CODER_QH_SCALE-i);
      ranges[i] = ra;
      remVal   -= ra;
    } else {
      ranges[i] = remVal;
      remVal = 0;
      break;
    }
  }
  if (remVal > 0) {
    for (int i = ARITH_CODER_QH_SCALE; i > MODEL_QH_LEN; --i) {
      unsigned ra = (uint64_t(model_qh_tail_tab[i]) * remVal + (uint32_t(1)<<31)) >> 32;
      if (ra == 0) ra = 1;
      ranges[i] = ra;
      remVal -= ra;
    }
    ranges[MODEL_QH_LEN] = remVal;
  }

  unsigned acc = 0;
  for (int i = 0; i < ARITH_CODER_QH_SCALE+1; ++i) {
    c2lo[i] = acc;
    acc += ranges[i];
  }
  c2lo[ARITH_CODER_QH_SCALE+1] = VAL_RANGE;
}

// return estimate of model length
static double prepare_qh_encode(encode_prm_t* prm, const uint32_t histogram[ARITH_CODER_QH_SCALE+1])
{
  uint32_t tot = 0;
  for (int i = 0; i < ARITH_CODER_QH_SCALE+1; ++i)
    tot += histogram[i];

  memset(prm->modelQh, 0, sizeof(prm->modelQh));
  for (int i = 0; i < MODEL_QH_LEN && tot != 0; ++i) {
    uint32_t hVal = histogram[i];
    prm->modelQh[i] = quantize_histogram_pair(hVal, tot-hVal, QQH_SCALE);
    tot -= hVal;
  }

  build_qh_encode_table(prm->modelC2loTab, prm->modelQh);

  double entr = 0;
  for (int i = 0; i < ARITH_CODER_QH_SCALE+1; ++i) {
    uint32_t h = histogram[i];
    if (h > 0) {
      unsigned ra = prm->modelC2loTab[i+1] - prm->modelC2loTab[i+0];
      entr += (RANGE_BITS-log2(ra))*h;
    }
  }

  return entr;
}

// return estimate of result length
static int prepare(encode_prm_t* prm, const uint8_t* src, int srclen, uint32_t srcHistogram[258], double* pInfo)
{
  // calculate total # of non-zero mtf characters
  uint32_t srcNzSymbols = 0;
  for (unsigned i = 1; i < 256; ++i)
    srcNzSymbols += srcHistogram[i];
  // Find a midpoint of the non-zero histogram
  uint32_t midNzSymbols = 0;
  unsigned mid0 = 0;
  for (unsigned i = 1; ; ++i) {
    midNzSymbols += srcHistogram[i];
    if (midNzSymbols*2 >= srcNzSymbols) {
      mid0 = i;
      break;
    }
  }
  if (mid0 == 255) mid0 = 254;

  prm->mid0 = mid0;
  prm->midFactors[0] = SelectBestMidFactor(&srcHistogram[0], 1, mid0);
  prm->midFactors[1] = SelectBestMidFactor(&srcHistogram[0], mid0+1, 255);

  prm->bisectionTab[0][0] = 1;
  prm->bisectionTab[0][1] = 1;
  prm->bisectionTab[1][0] = 0;
  prm->bisectionTab[1][1] = mid0+1;
  prm->bisectionTab[mid0+1][0] = BuildBisectionTab(prm, 2,   mid0+1, prm->midFactors[0]);
  prm->bisectionTab[mid0+1][1] = BuildBisectionTab(prm, mid0+2, 256, prm->midFactors[1]);

  uint8_t cntrs[258][2] = {{0}};
  int nextTabOffset = 256;
  double entropy = 0;
  double quantizedEntropy = entropy;
  uint32_t qhOffsets[258], qhOffset = 0;
  uint8_t qhMtfTab[258][ARITH_CODER_QH_SCALE+1];
  for (int i = 0; i < 258; ++i)
    memcpy(qhMtfTab[i], qhMtf0, sizeof(qhMtf0));
  uint32_t qhHistogram[ARITH_CODER_QH_SCALE+1]={0};
  for (int src_i = 0; src_i < srclen; ++src_i) {
    unsigned c = src[src_i];
    if (c == 255) {
      // escape sequence for symbols 255 and 256
      c = unsigned(src[src_i+1]) + 1;
      ++src_i;
    }
    int tMid = 1;
    int tabOffset = nextTabOffset;
    nextTabOffset = c < 2 ? 0 : 256;
    do {
      unsigned b = c > tMid;
      int idx = tabOffset + tMid;
      tabOffset = b ? 0 : tabOffset;
      int cntr0 = cntrs[idx][0];
      if (cntr0 == 0) {
        qhOffsets[idx] = qhOffset;
        ++qhOffset;
      }
      ++cntr0;
      cntrs[idx][0] = cntr0;
      cntrs[idx][1] += b;
      if (cntr0 == ARITH_CODER_CNT_MAX) {
        unsigned hVal = cntrs[idx][1];
        entropy += entropy_tab[hVal];
        quantizedEntropy += quantized_entropy_tab[hVal];
        unsigned qhVal = h2qh_tab[hVal]; // quantized '1'-to-total ratio
        prm->qh[qhOffsets[idx]] = qhVal;
        cntrs[idx][0] = cntrs[idx][1] = 0;

        int mtf_k = qh_mtf_encode(qhMtfTab[idx], qhVal);  // mtf encoder
        ++qhHistogram[mtf_k];
      }
      tMid = prm->bisectionTab[tMid][b];
    } while (tMid != 1);
  }

  // process partial sections
  for (int idx = 0; idx < 258; ++idx) {
    unsigned tot = cntrs[idx][0];
    if (tot != 0) {
      unsigned h0 = cntrs[idx][1];
      unsigned qhVal = quantize_histogram_pair(h0, tot-h0, ARITH_CODER_QH_SCALE);
      prm->qh[qhOffsets[idx]] = qhVal;
      entropy +=
          log2_tab[tot]*tot
        - log2_tab[h0]*h0
        - log2_tab[tot-h0]*(tot-h0);
      if (h0 != 0 && h0 != tot) {
        int32_t ra0 = qh2h_tab[qhVal];
        quantizedEntropy +=
            RANGE_BITS*tot
          - log2(ra0)*h0
          - log2(VAL_RANGE-ra0)*(tot-h0);
      }

      int mtf_k = qh_mtf_encode(qhMtfTab[idx], qhVal);  // mtf encoder
      ++qhHistogram[mtf_k];
    }
  }

  double modelLen = 4*8 + prepare_qh_encode(prm, qhHistogram);
  if (pInfo) {
    pInfo[0] = entropy;
    pInfo[1] = modelLen;
    pInfo[3] = quantizedEntropy;
  }

  return ceil((quantizedEntropy+modelLen)/8);
}

static void inc_dst(uint8_t* dst) {
  uint8_t val;
  do {
    --dst;
    *dst = (val = *dst + 1);
  } while (val==0);
}

static int encode(uint8_t* dst, encode_prm_t* prm, const uint8_t* src, int srclen)
{
  uint8_t cntrs[258] = {0};
  uint16_t currH[258];
  uint8_t qhMtfTab[258][ARITH_CODER_QH_SCALE+1];
  for (int i = 0; i < 258; ++i)
    memcpy(qhMtfTab[i], qhMtf0, sizeof(qhMtf0));

  const uint64_t MSB_MSK   = uint64_t(255) << 56;
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);
  uint64_t lo     = 0;                              // scaled by 2**64
  uint64_t range  = uint64_t(1) << (64-RANGE_BITS); // scaled by 2**(64-RANGE_BITS)
  uint8_t* dst0   = dst;
  uint64_t prevLo = lo;

  int nextTabOffset = 256;
  const uint8_t* qh = prm->qh;
  for (int src_i = 0; src_i < srclen; ++src_i) {
    uint16_t bitsBuf[256*2];
    uint16_t* wrBits = bitsBuf;

    unsigned c = src[src_i];
    if (c == 255) {
      // escape sequence for symbols 255 and 256
      c = unsigned(src[src_i+1]) + 1;
      ++src_i;
    }

    int tabOffset = nextTabOffset;
    nextTabOffset = c < 2 ? 0 : 256;
    int tMid = 1;
    do {
      unsigned idx = tabOffset + tMid;
      int cntr = cntrs[idx];
      if (cntr == 0) {
        unsigned qhVal = *qh++;
        int mtf_k = qh_mtf_encode(qhMtfTab[idx], qhVal);  // mtf encoder
        unsigned cLo = prm->modelC2loTab[mtf_k+0];
        unsigned cHi = prm->modelC2loTab[mtf_k+1];
        unsigned cRa = cHi - cLo;
        wrBits[0] = cLo;
        wrBits[1] = cRa;
        if (cRa != VAL_RANGE)
          wrBits += 2;
        currH[idx] = qhVal != ARITH_CODER_QH_SCALE ? qh2h_tab[qhVal] : 0;
      }
      ++cntr;
      if (cntr == ARITH_CODER_CNT_MAX)
        cntr = 0;
      cntrs[idx] = cntr;

      unsigned cRa1 = currH[idx];
      unsigned cRa0 = VAL_RANGE - cRa1;
      unsigned b = c > tMid;
      tabOffset = b == 0 ? tabOffset : 0;
      wrBits[0] = b == 0 ? 0         : cRa0;
      wrBits[1] = b == 0 ? cRa0      : cRa1;
      if (cRa1 != 0)
        wrBits += 2;
      tMid = prm->bisectionTab[tMid][b];
    } while (tMid != 1);

    for (uint16_t* rdBits = bitsBuf; rdBits != wrBits; rdBits += 2) {
      uint64_t cLo = rdBits[0];
      uint64_t cRa = rdBits[1];
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

// return value:
//  0 - not compressible, because all input characters have approximately equal probability or because input is too short
// >0 - the length of compressed buffer
int arithmetic_encode(uint8_t* dst, const uint8_t* src, int srclen, uint32_t srcHistogram[258], int origlen, uint8_t* tmp, double* pInfo)
{
  encode_prm_t prm;
  prm.qh = tmp;
  int estlen = prepare(&prm, src, srclen, srcHistogram, pInfo);
  if (estlen >= origlen)
    return 0; // not compressible

  dst[0] = prm.mid0; // [1..254]
  dst[1] = (prm.midFactors[0]-16)*16 + (prm.midFactors[1]-16);
  dst[2] = prm.modelQh[0]*16 + prm.modelQh[1];
  dst[3] = prm.modelQh[2]*16 + prm.modelQh[3];

  int reslen = encode(&dst[4], &prm, src, srclen) + 4;

  if (pInfo) {
    double modelLenBits = pInfo[1];
    pInfo[2] = reslen*8.0 - modelLenBits;
    pInfo[4] = prm.mid0;
    pInfo[5] = prm.midFactors[0];
    pInfo[6] = prm.midFactors[1];
  }

  if (reslen >= origlen)
    return 0; // not compressible

  return reslen;
}


/*
    pInfo[0] = entropy;
    pInfo[1] = modelLenBits;
    pInfo[2] = reslen*8.0 - modelLenBits;;
    pInfo[3] = quantizedEntropy;
    pInfo[4] = hdrs->a[0].nChunks;
    pInfo[5] = hdrs->a[1].nChunks;
*/
