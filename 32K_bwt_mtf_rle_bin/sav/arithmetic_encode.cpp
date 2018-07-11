#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <algorithm>

#include "arithmetic_encode.h"
#include "arithmetic_coder_cfg.h"

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
      - log2(ra)*i
      - log2(VAL_RANGE-ra)*(ARITH_CODER_CNT_MAX-i);
    // printf("%3d %2d %5d %.3f\n", i, h2qh_tab[i], ra, quantized_entropy_tab[i]);
  }
}

enum {
  // context header
  CONTEXT_HDR_SRC_OFFSET_I=0,
  CONTEXT_HDR_SRC_LEN_I,
  CONTEXT_HDR_QH_OFFSET_I,
  CONTEXT_HDR_QH_LEN_I,
  CONTEXT_HDR_PREV_C0_I,
  CONTEXT_HDR_MODEL_LEN_I,
  CONTEXT_HDR_S_ENTROPY_I,
  CONTEXT_HDR_Q_ENTROPY_I  = CONTEXT_HDR_S_ENTROPY_I + sizeof(double)/sizeof(uint32_t),
  CONTEXT_HDR_Q1_ENTROPY_I = CONTEXT_HDR_Q_ENTROPY_I + sizeof(double)/sizeof(uint32_t),
  CONTEXT_HDR_QH_OFFSETS_I = CONTEXT_HDR_Q1_ENTROPY_I + sizeof(double)/sizeof(uint32_t),
  CONTEXT_HDR_CNTRS_I   = CONTEXT_HDR_QH_OFFSETS_I + 258,
  CONTEXT_CNTRS_SZ      = (258*2*sizeof(uint8_t)-1)/sizeof(uint32_t) + 1,
  CONTEXT_HDR_PREV_QH_I = CONTEXT_HDR_CNTRS_I + CONTEXT_CNTRS_SZ,
  CONTEXT_PREV_QH_SZ    = (258*sizeof(uint8_t)-1)/sizeof(uint32_t) + 1,
  CONTEXT_HDR_PREV_H_I  = CONTEXT_HDR_PREV_QH_I + CONTEXT_PREV_QH_SZ,
  CONTEXT_PREV_H_SZ     = (258*sizeof(uint8_t)-1)/sizeof(uint32_t) + 1,
  CONTEXT_HDR_LEN  = CONTEXT_HDR_PREV_H_I + CONTEXT_PREV_H_SZ,
};

void arithmetic_encode_init_context(uint32_t* context, int tilelen)
{
  uint32_t srcOffset = CONTEXT_HDR_LEN;
  context[CONTEXT_HDR_SRC_OFFSET_I] = srcOffset;
  context[CONTEXT_HDR_SRC_LEN_I]    = 0;
  context[CONTEXT_HDR_QH_OFFSET_I]  = srcOffset + (tilelen*2-1)/sizeof(uint32_t) + 1;
  context[CONTEXT_HDR_QH_LEN_I]     = 0;
  context[CONTEXT_HDR_PREV_C0_I]    = 0;
  context[CONTEXT_HDR_MODEL_LEN_I]  = 0;
  memset(&context[CONTEXT_HDR_S_ENTROPY_I], 0, sizeof(double));
  memset(&context[CONTEXT_HDR_Q_ENTROPY_I], 0, sizeof(double));
  memset(&context[CONTEXT_HDR_Q1_ENTROPY_I], 0, sizeof(double));
  memset(&context[CONTEXT_HDR_QH_OFFSETS_I], 0, sizeof(uint32_t)*258);
  memset(&context[CONTEXT_HDR_CNTRS_I], 0, sizeof(uint32_t)*CONTEXT_CNTRS_SZ);
  memset(&context[CONTEXT_HDR_PREV_QH_I], ARITH_CODER_QH_SCALE/2, sizeof(uint32_t)*CONTEXT_PREV_QH_SZ);
  memset(&context[CONTEXT_HDR_PREV_H_I], ARITH_CODER_CNT_MAX/2, sizeof(uint32_t)*CONTEXT_PREV_H_SZ);
}

static int calculate_model_length(unsigned prev, unsigned val)
{
  int modelLen = 1; // isEqual ?
  if (prev != val) {
    unsigned ra = prev, diff = val - prev;
    if (val > prev) {
      modelLen += (prev != 0); // sign
      ra   = ARITH_CODER_QH_SCALE - prev;
    } else {
      modelLen += (prev != ARITH_CODER_QH_SCALE); // sign
      diff = prev - val;
    }
    prev = val;
    if (ra > 1) {
      modelLen += 1; // is +/-1
      if (diff > 1) {
        ra   -= 1;
        diff -= 2;
        while (ra > 1) {
          modelLen += 1;
          unsigned mid = ra / 2;
          if (diff < mid) {
            ra = mid;
          } else {
            ra   -= mid;
            diff -= mid;
          }
        }
      }
    }
  }
  return modelLen;
}

static void add_double(uint32_t* dst, double val) {
  double x;
  memcpy(&x, dst, sizeof(double));
  x += val;
  memcpy(dst, &x, sizeof(double));
}

void arithmetic_encode_chunk_callback(void* context_ptr, const uint8_t* src, int srclen)
{
  uint32_t* context = static_cast<uint32_t*>(context_ptr);
  uint8_t*  cntrs   = reinterpret_cast<uint8_t*>(&context[CONTEXT_HDR_CNTRS_I]);
  uint8_t*  prevQh  = reinterpret_cast<uint8_t*>(&context[CONTEXT_HDR_PREV_QH_I]);
  uint8_t*  prevH   = reinterpret_cast<uint8_t*>(&context[CONTEXT_HDR_PREV_H_I]);
  uint8_t*  qh      = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_QH_OFFSET_I]]);
  uint32_t  qhlen   = context[CONTEXT_HDR_QH_LEN_I];

  int  prevC0 = context[CONTEXT_HDR_PREV_C0_I];

  // update histograms
  double entropy = 0;
  double quantizedEntropy = 0;
  double quantizedEntropy1 = 0;
  for (int i = 0; i < srclen; ++i) {
    unsigned c = src[i];
    if (c == 255) {
      // escape sequence for symbols 255 and 256
      c = unsigned(src[i+1]) + 1;
      ++i;
    }
    unsigned tLo = 0, tHi = 256;
    do {
      unsigned tMid = (tLo*3 + tHi)/4;
      unsigned b = c > tMid;
      tLo = (b == 0) ? tLo  : tMid + 1;
      tHi = (b == 0) ? tMid : tHi;
      unsigned idx = (tMid < 2)  & prevC0 ? tMid : tMid + 2; // separate statistics for non-MS characters of zero run
      unsigned cntr0 = cntrs[idx*2+0];
      if (cntr0 == 0) {
        context[CONTEXT_HDR_QH_OFFSETS_I+idx] = qhlen;
        ++qhlen;
      }
      cntr0 += 1;
      cntrs[idx*2+0] = cntr0;
      cntrs[idx*2+1] += b;
      if (cntr0 == ARITH_CODER_CNT_MAX) {
        int hVal  = cntrs[idx*2+1];
        cntrs[idx*2+0] = cntrs[idx*2+1] = 0;

        entropy          += entropy_tab[hVal];
        {
        quantizedEntropy += quantized_entropy_tab[hVal];
        unsigned qhVal = h2qh_tab[hVal]; // quantized '1'-to-total ratio
        qh[context[CONTEXT_HDR_QH_OFFSETS_I+idx]] = qhVal;

        // estimate model length
        context[CONTEXT_HDR_MODEL_LEN_I] += calculate_model_length(prevQh[idx], qhVal);
        prevQh[idx] = qhVal;
        }

        {
        int prevHVal = prevH[idx];
        prevH[idx]   = hVal;
        int newQhVal = h2qh_tab[hVal];
        int oldQhVal = h2qh_tab[prevHVal];
        double oldQhEntropy = 0;
        if (hVal != 0 && hVal != ARITH_CODER_CNT_MAX) {
          unsigned ra = (uint32_t(prevHVal)*VAL_RANGE*2 + ARITH_CODER_CNT_MAX)/(ARITH_CODER_CNT_MAX*2);
          oldQhEntropy =
            RANGE_BITS*ARITH_CODER_CNT_MAX
              - log2(ra)*hVal
              - log2(VAL_RANGE-ra)*(ARITH_CODER_CNT_MAX-hVal);
        }
        int nxtQhVal = oldQhVal;
        if (newQhVal != oldQhVal) {
          double newQhEntropy = quantized_entropy_tab[hVal]; // case of change of qh
          if (newQhEntropy+1 < oldQhEntropy) {
            oldQhEntropy = newQhEntropy;
            nxtQhVal     = newQhVal;
          }
        }
        quantizedEntropy1 += oldQhEntropy;
        quantizedEntropy1 += calculate_model_length(oldQhVal, nxtQhVal);
        }
      }
    } while (tLo != tHi);
    prevC0 = (c < 2);
  }
  context[CONTEXT_HDR_PREV_C0_I] = prevC0;
  context[CONTEXT_HDR_QH_LEN_I]  = qhlen;
  add_double(&context[CONTEXT_HDR_S_ENTROPY_I], entropy);
  add_double(&context[CONTEXT_HDR_Q_ENTROPY_I], quantizedEntropy);
  add_double(&context[CONTEXT_HDR_Q1_ENTROPY_I], quantizedEntropy1);

  uint8_t* dst = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_SRC_OFFSET_I]]);
  uint32_t dstlen = context[CONTEXT_HDR_SRC_LEN_I];
  context[CONTEXT_HDR_SRC_LEN_I] = dstlen + srclen;
  memcpy(&dst[dstlen], src, srclen);
}

// return estimate of result length
static int prepare(uint32_t* context, double* pInfo)
{
  uint8_t*  cntrs  = reinterpret_cast<uint8_t*>(&context[CONTEXT_HDR_CNTRS_I]);
  uint8_t*  prevQh = reinterpret_cast<uint8_t*>(&context[CONTEXT_HDR_PREV_QH_I]);
  uint8_t*  prevH  = reinterpret_cast<uint8_t*>(&context[CONTEXT_HDR_PREV_H_I]);
  uint8_t*  qh     = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_QH_OFFSET_I]]);

  double entropy, quantizedEntropy, quantizedEntropy1;
  memcpy(&entropy, &context[CONTEXT_HDR_S_ENTROPY_I], sizeof(double));
  memcpy(&quantizedEntropy, &context[CONTEXT_HDR_Q_ENTROPY_I], sizeof(double));
  memcpy(&quantizedEntropy1, &context[CONTEXT_HDR_Q1_ENTROPY_I], sizeof(double));
  int32_t modelLen = context[CONTEXT_HDR_MODEL_LEN_I];

  // process partial sections
  for (int idx = 0; idx < 258; ++idx) {
    unsigned tot = cntrs[idx*2+0];
    if (tot != 0) {
      unsigned h0 = cntrs[idx*2+1];
      unsigned qhVal = quantize_histogram_pair(h0, tot-h0, ARITH_CODER_QH_SCALE);
      qh[context[CONTEXT_HDR_QH_OFFSETS_I+idx]] = qhVal;
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

        {
        int prevHVal = prevH[idx];
        int newQhVal = qhVal;
        int oldQhVal = h2qh_tab[prevHVal];
        double oldQhEntropy = 0;
        int ra = (uint32_t(prevHVal)*VAL_RANGE*2 + ARITH_CODER_CNT_MAX)/(ARITH_CODER_CNT_MAX*2);
        oldQhEntropy =
            RANGE_BITS*tot
          - log2(ra)*h0
          - log2(VAL_RANGE-ra)*(tot-h0);
        int nxtQhVal = oldQhVal;
        if (newQhVal != oldQhVal) {
          int32_t ra0 = qh2h_tab[qhVal];
          double newQhEntropy =
              RANGE_BITS*tot
            - log2(ra0)*h0
            - log2(VAL_RANGE-ra0)*(tot-h0);
          if (newQhEntropy < oldQhEntropy) {
            oldQhEntropy = newQhEntropy;
            nxtQhVal     = newQhVal;
          }
        }
        quantizedEntropy1 += oldQhEntropy;
        quantizedEntropy1 += calculate_model_length(oldQhVal, nxtQhVal);
        }
      }

      // estimate model length
      modelLen += calculate_model_length(prevQh[idx], qhVal);
    }
  }

  if (pInfo) {
    pInfo[0] = entropy;
    pInfo[1] = modelLen;
    pInfo[3] = quantizedEntropy;
    pInfo[4] = quantizedEntropy1/8;
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

static uint16_t* encodeQh(uint16_t* wrBits, unsigned val, unsigned prev)
{
  int isNotEqual = (val != prev);
  wrBits[0] = isNotEqual;
  wrBits[1] = VAL_RANGE/2;
  wrBits += 2;
  if (isNotEqual) {
    int up = val > prev;
    wrBits[0] = up;
    wrBits[1] = VAL_RANGE/2;
    unsigned ra = prev, diff = val - prev;
    if (up) {
      if (prev != 0)
        wrBits += 2;
      ra = ARITH_CODER_QH_SCALE - prev;
    } else {
      if (prev != ARITH_CODER_QH_SCALE)
        wrBits += 2;
      diff = prev - val;
    }
    prev = val;
    if (ra > 1) {
      int gtThanOne = diff > 1;
      wrBits[0] = gtThanOne;
      wrBits[1] = VAL_RANGE/2;
      wrBits += 2;
      if (gtThanOne) {
        ra   -= 1;
        diff -= 2;
        while (ra > 1) {
          unsigned mid = ra / 2;
          int ge = diff >= mid;
          wrBits[0] = ge;
          wrBits[1] = VAL_RANGE/2;
          wrBits += 2;
          if (ge) {
            ra   -= mid;
            diff -= mid;
          } else {
            ra = mid;
          }
        }
      }
    }
  }
  return wrBits;
}

static int encode(uint8_t* dst, const uint32_t* context)
{
  uint8_t cntrs[258] = {0};
  uint8_t prevQh[258];
  uint16_t currH[258];
  memset(prevQh, ARITH_CODER_QH_SCALE/2, sizeof(prevQh));

  const uint64_t MSB_MSK   = uint64_t(255) << 56;
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);
  uint64_t lo     = 0;                              // scaled by 2**64
  uint64_t range  = uint64_t(1) << (64-RANGE_BITS); // scaled by 2**(64-RANGE_BITS)
  uint8_t* dst0   = dst;
  uint64_t prevLo = lo;

  const uint8_t* qh  = reinterpret_cast<const uint8_t*>(&context[context[CONTEXT_HDR_QH_OFFSET_I]]);
  const uint8_t* src = reinterpret_cast<const uint8_t*>(&context[context[CONTEXT_HDR_SRC_OFFSET_I]]);
  const int srclen = context[CONTEXT_HDR_SRC_LEN_I];
  unsigned prevC0 = 0;
  for (int src_i = 0; src_i < srclen; ++src_i) {
    uint16_t bitsBuf[9*9*2];
    uint16_t* wrBits = bitsBuf;

    unsigned c = src[src_i];
    if (c == 255) {
      // escape sequence for symbols 255 and 256
      c = unsigned(src[src_i+1]) + 1;
      ++src_i;
    }
    unsigned tLo = 0, tHi = 256;
    do {
      unsigned tMid = (tLo*3 + tHi)/4;
      unsigned idx = (tMid < 2)  & prevC0 ? tMid : tMid + 2; // separate statistics for non-MS characters of zero run
      int hVal = VAL_RANGE/2;
      int cntr = cntrs[idx];
      if (cntr == 0) {
        unsigned qhVal = *qh++;
        unsigned prevQhVal = prevQh[idx];
        prevQh[idx] = qhVal;
        wrBits = encodeQh(wrBits, qhVal, prevQhVal);
        currH[idx] = qhVal != ARITH_CODER_QH_SCALE ? qh2h_tab[qhVal] : 0;
      }
      ++cntr;
      if (cntr == ARITH_CODER_CNT_MAX)
        cntr = 0;
      cntrs[idx] = cntr;
      hVal = currH[idx];

      unsigned b = c > tMid;
      wrBits[0] = b;
      wrBits[1] = hVal;
      if (hVal != 0)
        wrBits += 2;
      tLo = (b == 0) ? tLo  : tMid + 1;
      tHi = (b == 0) ? tMid : tHi;
    } while (tLo != tHi);
    prevC0 = (c < 2);

    for (uint16_t* rdBits = bitsBuf; rdBits != wrBits; rdBits += 2) {
      unsigned b    = rdBits[0];
      unsigned hVal = rdBits[1];
      uint64_t cLo = b == 0 ? 0                : VAL_RANGE - hVal;
      uint64_t cRa = b == 0 ? VAL_RANGE - hVal : hVal;
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
int arithmetic_encode(uint32_t* context, uint8_t* dst, int origlen, double* pInfo)
{
  int estlen = prepare(context, pInfo);
  if (estlen >= origlen)
    return 0; // not compressible

  int reslen = encode(dst, context);

  if (pInfo) {
    double modelLenBits = pInfo[1];
    pInfo[2] = reslen*8.0 - modelLenBits;
    // pInfo[4] = 0;
    pInfo[5] = 0;
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
