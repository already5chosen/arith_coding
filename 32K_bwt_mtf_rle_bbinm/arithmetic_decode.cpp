#include <cstdint>
#include <cstring>
#include <cstdio>
#include <cmath>

#include "arithmetic_decode.h"
#include "arithmetic_coder_cfg.h"

static const unsigned VAL_RANGE = 1u << RANGE_BITS;
static const unsigned QQH_SCALE = 15;               // not necessarily the same as ARITH_CODER_QH_SCALE
static const unsigned MODEL_QH_LEN = 4;

static uint16_t qh2h_tab[ARITH_CODER_QH_SCALE+1];
static uint8_t  qhMtf0[ARITH_CODER_QH_SCALE+1];
static uint32_t model_qh_tail_tab[ARITH_CODER_QH_SCALE+1];
static uint32_t model_qh2h_tab[QQH_SCALE] = {
  0,
    46927670,   185659716,   410132882,   710536612,
  1073741824,  1483874706,  1923010482,  2371956814,
  2811092590,  3221225472,  3584430684,  3884834414,
  4109307580,  4248039626,
};

void arithmetic_decode_init_tables()
{
  qh2h_tab[ARITH_CODER_QH_SCALE] = VAL_RANGE;
  double M_PI = 3.1415926535897932384626433832795;
  for (int i = 1; i < ARITH_CODER_QH_SCALE; ++i)
    qh2h_tab[i] = int(round((sin((i*2-ARITH_CODER_QH_SCALE)*(M_PI/2/ARITH_CODER_QH_SCALE))+1.0)*(VAL_RANGE/2)));

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

static int qh_mtf_decode(uint8_t mtf_t[], int mtfC)
{ // move-to-front decoder
  int dstC = mtf_t[mtfC];
  if (mtfC != 0) {
    // update move-to-front decoder table
    memmove(&mtf_t[1], &mtf_t[0], mtfC);
    mtf_t[0] = dstC;
  }
  return dstC;
}

static int decode(
  uint8_t*       dst0,
  int            dstlen,
  int32_t        histogram[256],
  const uint8_t* src,
  int            srclen,
  const uint8_t  bisectionTab[256][2],
  const uint16_t modelC2loTab[ARITH_CODER_QH_SCALE+2])
{
  const uint16_t* modelC2loBeg = modelC2loTab;
  while (modelC2loBeg[1] == 0)
    ++modelC2loBeg;
  int modelOffset = modelC2loBeg - modelC2loTab;

  uint8_t  cntrs[258] = {0};
  uint16_t currH[258];
  uint8_t qhMtfTab[258][ARITH_CODER_QH_SCALE+1];
  for (int i = 0; i < 258; ++i)
    memcpy(qhMtfTab[i], qhMtf0, sizeof(qhMtf0));

  // initialize move-to-front decoder table
  uint8_t mtf_t[256];
  for (int i = 0; i < 256; ++i)
    mtf_t[i] = i;

  // initialize RLE decoder
  uint32_t rlInc = 1;

  // initialize arithmetic decoder
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);

  uint8_t tmpbuf[16] = {0};
  bool useTmpbuf = false;
  uint64_t value = 0; // scaled by 2**64
  for (int k = 0; k < 8; ++k) {
    value <<= 8;
    if (srclen > 0)
      value |= *src++;
    --srclen;
  }
  uint64_t range = uint64_t(1) << (64-RANGE_BITS);   // scaled by 2**(64-rb).  Maintained in [2**(33-rb)+1..2*(64-rb)]

  int nextTabOffset = 256;
  uint8_t* dst = dst0;
  while (dstlen > 0) {
    int tabOffset = nextTabOffset;
    int tMid = 1;
    unsigned c;
    do {
      unsigned idx = tabOffset + tMid;
      int hVal = VAL_RANGE-modelC2loBeg[1];
      int cntr = cntrs[idx];
      unsigned loadQh = 1;
      if (cntr != 0) {
        ++cntr;
        if (cntr == ARITH_CODER_CNT_MAX)
          cntr = 0;
        cntrs[idx] = cntr;
        hVal = currH[idx];
        loadQh = 0;
      }

      unsigned b;
      for (;;) {
        b = 0;
        if (hVal != 0) {
          b = 1;
          if (hVal != VAL_RANGE) {
            b = (value >= range*(VAL_RANGE-hVal));
            // keep decoder in sync with encoder
            uint64_t cLo = b == 0 ? 0                : VAL_RANGE - hVal;
            uint64_t cHi = b == 0 ? VAL_RANGE - hVal : VAL_RANGE;
            if ((loadQh & b) != 0) {
              uint64_t pr;
              do {
                ++b;
              } while ((pr=modelC2loBeg[b]*range-1) < value-1);
              if (pr > value-1)
                b -= 1;
              if (b > ARITH_CODER_QH_SCALE) {
                return -104; // should not happen
              }
              cLo = modelC2loBeg[b];
              cHi = modelC2loBeg[b+1];
            }
            value -= range * cLo;
            range *= (cHi-cLo);  // at this point range is scaled by 2**64 - the same scale as value
            uint64_t nxtRange = range >> RANGE_BITS;

            if (nxtRange <= MIN_RANGE) {
              #ifdef ENABLE_PERF_COUNT
              ++pModel->m_renormalizationCount;
              #endif
              if (srclen < 8) {
                if (!useTmpbuf && srclen > 0) {
                  memcpy(tmpbuf, src, srclen);
                  src = tmpbuf;
                  useTmpbuf = true;
                }
              }

              uint32_t threeOctets =
                (uint32_t(src[2])       ) |
                (uint32_t(src[1]) << 1*8) |
                (uint32_t(src[0]) << 2*8);
              value = (value << 24) + threeOctets;
              src    += 3;
              srclen -= 3;
              nxtRange = range << (24 - RANGE_BITS);
              if (srclen < -7)
                return dst-dst0;
              if (value > ((nxtRange<<RANGE_BITS)-1)) {
                return -103; // should not happen
              }
            }
            range = nxtRange;
          }
        }
        if (loadQh == 0)
          break;

        unsigned qhVal = qh_mtf_decode(qhMtfTab[idx], b+modelOffset); // QH mtf decode
        currH[idx] = hVal = qh2h_tab[qhVal];
        cntrs[idx] = 1;
        loadQh = 0;
      }
      tabOffset = b ? 0 : tabOffset;
      c = tMid + b;
      tMid = bisectionTab[tMid][b];
    } while (tMid != 1);
    // arithmetic decode of one character done
    nextTabOffset = c < 2 ? 0 : 256;

    // RLE decode
    if (c < 2) {
      // zero run
      uint32_t delta_rl = c+rlInc;
      rlInc += delta_rl;

      if (delta_rl > uint32_t(dstlen))
        return -26; // zero run too long (A)

      // insert zero run
      int cm0 = mtf_t[0];
      histogram[cm0] += delta_rl;
      memset(dst, cm0, delta_rl);
      dst    += delta_rl;
      dstlen -= delta_rl;
    } else {
      rlInc = 1;
      // MTF decode
      int mtfC = c - 1;
      int dstC = mtf_t[mtfC];
      histogram[dstC] += 1;
      *dst++ = dstC;
      dstlen -= 1;
      // update move-to-front decoder table
      memmove(&mtf_t[1], &mtf_t[0], mtfC);
      mtf_t[0] = dstC;
    }
  }

  return dst - dst0;
}

static int BuildBisectionTab(uint8_t tab[][2], unsigned lo, unsigned hi, unsigned midFactor)
{
  if (lo != hi) {
    unsigned mid = (lo*midFactor + hi*(32-midFactor))/32;
    tab[mid][0] = BuildBisectionTab(tab, lo, mid, midFactor);
    tab[mid][1] = BuildBisectionTab(tab, mid+1, hi, midFactor);
    return mid;
  }
  return 1;
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

// Arithmetic decode followed by RLE and MTF decode
int arithmetic_decode(
  uint8_t*       dst,
  int            dstlen,
  int32_t        histogram[256],
  const uint8_t* src,
  int            srclen,
  int*           pInfo)
{
  if (pInfo) {
    pInfo[0] = 0;
    pInfo[1] = 0;
    pInfo[2] = 0;
    pInfo[3] = 0;
    pInfo[4] = 0;
  }
  if (srclen < 5)
    return -1;

  unsigned mid0 = src[0];
  if (mid0 < 1 || mid0 > 254)
    return -2; // mid0 out of range
  unsigned midFactors = src[1];
  unsigned midFactor0 = midFactors/16 + 16;
  unsigned midFactor1 = midFactors%16 + 16;

  uint8_t bisectionTab[256][2];
  bisectionTab[0][0] = 1;
  bisectionTab[0][1] = 1;
  bisectionTab[1][0] = 0;
  bisectionTab[1][1] = mid0+1;
  bisectionTab[mid0+1][0] = BuildBisectionTab(bisectionTab, 2,   mid0+1, midFactor0);
  bisectionTab[mid0+1][1] = BuildBisectionTab(bisectionTab, mid0+2, 256, midFactor1);

  uint8_t modelQh[MODEL_QH_LEN];
  modelQh[0] = src[2] / 16;
  modelQh[1] = src[2] % 16;
  modelQh[2] = src[3] / 16;
  modelQh[3] = src[3] % 16;
  uint16_t modelC2loTab[ARITH_CODER_QH_SCALE+2];
  build_qh_encode_table(modelC2loTab, modelQh);

  memset(histogram, 0, sizeof(histogram[0])*256);
  int textlen = decode(dst, dstlen, histogram, &src[4], srclen-4, bisectionTab, modelC2loTab);
  return textlen;
}
