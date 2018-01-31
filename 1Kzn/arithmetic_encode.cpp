#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <algorithm>

#include "arithmetic_encode.h"


static const int RANGE_BITS = 10;
static const unsigned VAL_RANGE = 1u << RANGE_BITS;

typedef struct {
  unsigned firstSmallI; // index of the first small range
  unsigned smallSubrangeVal;
  struct {
    uint16_t val;
    uint8_t  alias; // 0..254 = index small range group, 255 - normal range
  } ar[257];
} c2range_t;

// return # of small ranges
static unsigned histogram_to_range(c2range_t* c2range, unsigned maxC, const unsigned* h, unsigned srclen)
{
  // translate counts to ranges and store in c2range
  // 1st pass - find small ranges, those that should be rounded up to reach 1
  unsigned nSmallRanges   = 0;
  unsigned smallRangesTot = 0;
  for (unsigned c = 0; c <= maxC; ++c) {
    unsigned cnt = h[c];
    unsigned alias = 255;
    if (cnt != 0) {
      if (uint64_t(cnt)*VAL_RANGE < srclen) {
        smallRangesTot += cnt;
        alias = nSmallRanges;
        if (nSmallRanges == 0)
          c2range->firstSmallI = c; // save index of the first small range
        ++nSmallRanges;
      }
    }
    c2range->ar[c].alias = alias;
  }


  // 2nd pass - translate remaining characters while rounding toward zero
  int remRange = VAL_RANGE;
  int n1cnt = 0; // count # of character with exact range >= 1
  int n2cnt = 0; // count # of character with exact range >= 2
  for (unsigned c = 0; c <= maxC; ++c) {
    unsigned alias = c2range->ar[c].alias;
    unsigned cnt = alias==255 ? h[c] : (alias==0 ? smallRangesTot : 0);
    unsigned range = 0;
    if (cnt != 0) {
      n1cnt += 1;
      range = 1;
      if (uint64_t(cnt)*VAL_RANGE >= srclen) {
        range = (uint64_t(cnt)*VAL_RANGE)/srclen;
        n2cnt += (range > 1);
      }
      remRange -= range;
    }
    c2range->ar[c].val = range;
  }

  if (remRange < 0) {
    while (-remRange >= n2cnt) {
      // very unlikely case when all ranges >1 have to be decremented
      n2cnt = 0;
      remRange = VAL_RANGE;
      for (unsigned c = 0; c <= maxC; ++c) {
        unsigned range = c2range->ar[c].val;
        if (range > 1) {
          c2range->ar[c].val = (range -= 1);
          n2cnt += (range > 1);
        }
        remRange -= range;
      }
    }
  } else if (remRange > 0) {
    while (remRange >= n1cnt) {
      // very unlikely case when all ranges > 0 have to be incremented
      remRange = VAL_RANGE;
      for (unsigned c = 0; c <= maxC; ++c) {
        unsigned range = c2range->ar[c].val;
        if (range > 0)
          c2range->ar[c].val = (range += 1);
        remRange -= range;
      }
    }
  }

  if (remRange != 0) {
    // calculate effect of increment/decrement of range on entropy
    double de[256], denz[256];
    int incrDecr = remRange > 0 ? 1 : -1;
    int denzLen = 0;
    for (unsigned c = 0; c <= maxC; ++c) {
      double deltaE = 0;
      unsigned alias = c2range->ar[c].alias;
      unsigned cnt = alias==255 ? h[c] : (alias==0 ? smallRangesTot : 0);
      if (cnt != 0) {
        int range = c2range->ar[c].val;
        int nxtRange = range+incrDecr;
        if (nxtRange != 0) {
          denz[denzLen] = deltaE = log2(double(range)/nxtRange)*cnt;
          ++denzLen;
        }
      }
      de[c] = deltaE;
    }

    int indx = remRange > 0 ? remRange - 1 : -remRange - 1;
    std::nth_element(&denz[0], &denz[indx], &denz[denzLen]);
    double thr = denz[indx];

    // when remRange > 0 increment ranges that will have maximal effect on entropy
    // when remRange > 0 decrement ranges that will have minimal effect on entropy
    for (unsigned c = 0; c <= maxC; ++c) {
      double deltaE = de[c];
      if (deltaE != 0 && deltaE < thr) {
        c2range->ar[c].val += incrDecr;
        remRange   -= incrDecr;
      }
    }
    for (unsigned c = 0; remRange != 0; ++c) {
      double deltaE = de[c];
      if (deltaE == thr) {
        c2range->ar[c].val += incrDecr;
        remRange   -= incrDecr;
      }
    }
  }
  return nSmallRanges;
}

// return value:
// -1  - source consists of repetition of the same character
// >=0 - maxC = the character with the biggest numeric value that appears in the source at least once
static int prepare1(const uint8_t* src, unsigned srclen, c2range_t* c2low, double* pQuantizedEntropy, double* pInfo)
{
  // calculated statistics of appearance
  unsigned stat[256]={0};
  for (unsigned i = 0; i < srclen; ++i)
    ++stat[src[i]];

  // find maximal frequency and highest-numbered character that occurred at least once
  unsigned maxC = 255;
  unsigned maxCnt = 0;
  for (unsigned c = 0; c < 256; ++c) {
    unsigned cnt = stat[c];
    if (cnt != 0) {
      maxC = c;
      if (maxCnt < cnt)
        maxCnt = cnt;
    }
  }

  if (pInfo) {
    // calculate source entropy
    double entropy = 0;
    for (unsigned c = 0; c <= maxC; ++c) {
      unsigned cnt = stat[c];
      if (cnt)
        entropy += log2(double(srclen)/cnt)*cnt;
    }
    pInfo[0] = entropy;
    pInfo[3] = 0;
  }

  if (maxCnt==srclen)
    return -1; // source consists of repetition of the same character

  unsigned nSmallRanges = histogram_to_range(c2low, maxC, stat, srclen);
  unsigned smallSubrangeVal = 0;
  double smallRangeNbits = 0;
  if (nSmallRanges > 0) {
    smallSubrangeVal = VAL_RANGE/nSmallRanges;
    smallRangeNbits  = log2(double(VAL_RANGE)/c2low->ar[c2low->firstSmallI].val*double(VAL_RANGE)/smallSubrangeVal);
  }
  c2low->smallSubrangeVal = smallSubrangeVal;

  // calculate entropy after quantization
  double entropy = 0;
  for (unsigned c = 0; c <= maxC; ++c) {
    unsigned cnt = stat[c];
    if (cnt) {
      double nbits = c2low->ar[c].alias == 255 ? log2(double(VAL_RANGE)/c2low->ar[c].val) : smallRangeNbits;
      entropy += nbits*cnt;
    }
  }
  *pQuantizedEntropy = entropy;
  if (pInfo)
    pInfo[3] = entropy;
  return maxC;
}

static void prepare2(c2range_t* c2low, unsigned maxC)
{
  // c2low -> cumulative sums of ranges
  unsigned lo = 0;
  for (unsigned c = 0; c <= maxC; ++c) {
    unsigned range = c2low->ar[c].val;
    // printf("%3u: %04x\n", c, range);
    c2low->ar[c].val = lo;
    lo += range;
  }
  c2low->ar[maxC+1].val = VAL_RANGE;
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
static int store_model(uint8_t* dst, const c2range_t* c2range, unsigned maxC, double* pNbits, CArithmeticEncoder* pEnc)
{
  int nRanges = 0;
  unsigned hist[RANGE_BITS+1] = {0};
  for (unsigned c = 0; c < maxC; ++c) {
    uint32_t range = c2range->ar[c].val;
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

  // store maxC
  p = pEnc->put(256, maxC, 1, p);
  // store histogram of log2
  unsigned rem = maxC;
  for (int i = 0; i < RANGE_BITS && rem > 0; ++i) {
    p = pEnc->put(rem+1, hist[i], 1, p);
    rem -= hist[i];
  }

  // store c2range
  for (unsigned c = 0; c < maxC; ++c) {
    uint32_t range = c2range->ar[c].val;
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

static int encode(uint8_t* dst, const uint8_t* src, unsigned srclen, const c2range_t* c2low, CArithmeticEncoder* pEnc)
{
  const uint64_t MSB_MSK   = uint64_t(255) << 56;
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);
  uint64_t lo     = pEnc->m_lo << 1;              // scaled by 2**64
  uint64_t range  = pEnc->m_range >> (RANGE_BITS-1); // scaled by 2**50
  uint8_t* dst0   = dst;
  uint64_t prevLo = lo;
  int alias = 255;
  for (unsigned i = 0; i < srclen; ) {
    int c;
    uint64_t cLo, cRa;
    if (alias == 255) {
      // process next input character
      int c0 = c = src[i];
      alias = c2low->ar[c0].alias;
      if (alias != 255)
        c = c2low->firstSmallI; // remap all 'small' characters to the first small character
      ++i;
      cLo = c2low->ar[c+0].val;
      cRa = c2low->ar[c+1].val - cLo;
    } else {
      // complete processing of 'small' character
      cLo = c2low->smallSubrangeVal*alias;
      cRa = c2low->smallSubrangeVal;
      alias = 255;
    }
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
int arithmetic_encode(std::vector<uint8_t>* dst, const uint8_t* src, int srclen, double* pInfo)
{
  c2range_t c2low;
  double quantizedEntropy;
  int maxC = prepare1(src, srclen, &c2low, &quantizedEntropy, pInfo);

  if (maxC == -1)
    return -1; // source consists of repetition of the same character

  size_t sz0 = dst->size();
  dst->resize(sz0 + 640);
  CArithmeticEncoder enc;
  enc.init();
  double modelLenBits;
  unsigned modellen = store_model(&dst->at(sz0), &c2low, maxC, &modelLenBits, &enc);
  if (pInfo)
    pInfo[1] = modelLenBits;

  int lenEst = (modelLenBits+quantizedEntropy+7)/8;
  if (lenEst >= srclen)
    return 0; // not compressible

  // printf("ml=%u\n", modellen);
  dst->resize(sz0 + lenEst + 64);
  prepare2(&c2low, maxC);
  int dstlen = encode(&dst->at(sz0+modellen), src, srclen, &c2low, &enc);

  if (pInfo)
    pInfo[2] = dstlen*8.0;

  int reslen = modellen + dstlen;

  if (reslen >= srclen)
    return 0; // not compressible

  dst->resize(sz0 + reslen);
  return reslen;
}
