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
  CONTEXT_HDR_H_OFFSET_I,
  CONTEXT_HDR_H_LEN_I,
  CONTEXT_HDR_PREV_C0_I,
  CONTEXT_HDR_NC_I,
  CONTEXT_HDR_CNTRS_I = CONTEXT_HDR_NC_I + 258,
  CONTEXT_CNTRS_SZ = (258*2*sizeof(uint8_t))/sizeof(uint32_t),
  CONTEXT_HDR_LEN  = CONTEXT_HDR_CNTRS_I + CONTEXT_CNTRS_SZ,
};

void arithmetic_encode_init_context(uint32_t* context, int tilelen)
{
  uint32_t srcOffset = CONTEXT_HDR_LEN;
  context[CONTEXT_HDR_SRC_OFFSET_I] = srcOffset;
  context[CONTEXT_HDR_SRC_LEN_I]    = 0;
  context[CONTEXT_HDR_H_OFFSET_I]   = srcOffset + (tilelen*2-1)/sizeof(uint32_t) + 1;
  context[CONTEXT_HDR_H_LEN_I]      = 0;
  context[CONTEXT_HDR_PREV_C0_I]    = 0;
  memset(&context[CONTEXT_HDR_NC_I], 0, sizeof(uint32_t)*258);
  memset(&context[CONTEXT_HDR_CNTRS_I], 0, sizeof(uint32_t)*CONTEXT_CNTRS_SZ);
}

void arithmetic_encode_chunk_callback(void* context_ptr, const uint8_t* src, int srclen)
{
  uint32_t* context = static_cast<uint32_t*>(context_ptr);
  uint8_t*  cntrs = reinterpret_cast<uint8_t*>(&context[CONTEXT_HDR_CNTRS_I]);
  uint32_t  hlen = context[CONTEXT_HDR_H_LEN_I];
  uint8_t*  hbeg = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_H_OFFSET_I]]);
  uint8_t*  h = &hbeg[hlen];

  int  prevC0 = context[CONTEXT_HDR_PREV_C0_I];

  // update histograms
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
      cntrs[idx*2+0] += 1;
      cntrs[idx*2+1] += b;
      if (cntrs[idx*2+0] == ARITH_CODER_CNT_MAX) {
        unsigned bSum = cntrs[idx*2+1];
        cntrs[idx*2+0] = cntrs[idx*2+1] = 0;
        ++context[CONTEXT_HDR_NC_I+idx];
        if (idx >= 255) {
          *h++ = 255;
          idx -= 255;
        }
        h[0] = idx;
        h[1] = bSum;
        h += 2;
      }
    } while (tLo != tHi);
    prevC0 = (c < 2);
  }
  context[CONTEXT_HDR_PREV_C0_I] = prevC0;

  hbeg = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_H_OFFSET_I]]);
  context[CONTEXT_HDR_H_LEN_I] = (h-hbeg);

  uint8_t* dst = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_SRC_OFFSET_I]]);
  uint32_t dstlen = context[CONTEXT_HDR_SRC_LEN_I];
  context[CONTEXT_HDR_SRC_LEN_I] = dstlen + srclen;
  memcpy(&dst[dstlen], src, srclen);
}

// return estimate of result length
static int prepare(uint32_t* context, uint32_t qhOffsets[258], double* pInfo)
{
  uint8_t*  cntrs = reinterpret_cast<uint8_t*>(&context[CONTEXT_HDR_CNTRS_I]);
  uint32_t offset = 0;
  for (int i = 0; i < 258; ++i) {
    qhOffsets[i] = offset;
    offset += context[CONTEXT_HDR_NC_I+i] + (cntrs[2*i+0] != 0);
  }
  uint32_t qho[258];
  memcpy(qho, qhOffsets, sizeof(qho));

  const uint32_t hlen = context[CONTEXT_HDR_H_LEN_I];
  const uint8_t* h = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_H_OFFSET_I]]);
  uint8_t* qh = reinterpret_cast<uint8_t*>(&context[context[CONTEXT_HDR_H_OFFSET_I]]) + hlen;
  double entropy = 0;
  double quantizedEntropy = entropy;

  // process full sections
  for (uint32_t src_i = 0; src_i < hlen; src_i += 2) {
    int pos = h[src_i+0];
    int val = h[src_i+1];
    if (pos == 255) {
      pos += val;
      val = h[src_i+2];
      ++src_i;
    }
    uint32_t offs = qho[pos];
    qh[offs] = h2qh_tab[val];
    qho[pos] = offs + 1;
    entropy += entropy_tab[val];
    quantizedEntropy += quantized_entropy_tab[val];
  }

  // process partial sections
  for (int i = 0; i < 258; ++i) {
    unsigned tot = cntrs[2*i+0];
    if (tot != 0) {
      unsigned h0 = cntrs[2*i+1];
      unsigned qhVal = quantize_histogram_pair(h0, tot-h0, ARITH_CODER_QH_SCALE);
      qh[qho[i]] = qhVal;
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
    }
  }

  int32_t modelLen = 0;
  for (int i = 0; i < 258; ++i) {
    const uint32_t nc = context[CONTEXT_HDR_NC_I+i] + (cntrs[2*i+0] != 0);
    const uint8_t* src = &qh[qhOffsets[i]];
    unsigned prev = 0;
    for (uint32_t c = 0; c < nc; ++c) {
      unsigned val = src[c];
      modelLen += 1; // isEqual ?
      if (val != prev) {
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
    }
  }

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

static int encode(uint8_t* dst, const uint32_t* context, uint32_t qhOffsets[258])
{
  uint8_t cntrs[258] = {0};
  uint8_t prevQh[258] = {0};
  uint16_t currH[258];

  const uint64_t MSB_MSK   = uint64_t(255) << 56;
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);
  uint64_t lo     = 0;                              // scaled by 2**64
  uint64_t range  = uint64_t(1) << (64-RANGE_BITS); // scaled by 2**(64-RANGE_BITS)
  uint8_t* dst0   = dst;
  uint64_t prevLo = lo;


  const uint8_t* qh  = reinterpret_cast<const uint8_t*>(&context[context[CONTEXT_HDR_H_OFFSET_I]]) + context[CONTEXT_HDR_H_LEN_I];
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
        uint32_t qhOffset = qhOffsets[idx];
        unsigned qhVal = qh[qhOffset];
        unsigned prevQhVal = prevQh[idx];
        qhOffsets[idx] = qhOffset + 1;
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
  uint32_t qhOffsets[258];
  int estlen = prepare(context, qhOffsets, pInfo);
  if (estlen >= origlen)
    return 0; // not compressible

  int reslen = encode(dst, context, qhOffsets);

  if (pInfo) {
    double modelLenBits = pInfo[1];
    pInfo[2] = reslen*8.0 - modelLenBits;
    pInfo[4] = 0;
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
