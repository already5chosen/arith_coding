// #include <immintrin.h>
#include <cstdio>
#include <cfenv>
#include <cstring>
#include <cmath>

#include "arithmetic_encode.h"

#define _LIKELY(cond)   __builtin_expect((cond), 1)
#define _UNLIKELY(cond) __builtin_expect((cond), 0)

static const int      RANGE_BITS = 14;
static const unsigned VAL_RANGE = 1u << RANGE_BITS;
static const float    INV_VAL_RANGE = 1.0f / VAL_RANGE;
static const double   RANGE_LEAKAGE_FACTOR = 1.0 - 1.0/(int64_t(1)<<29);

typedef union {
  uint32_t u;
  float    f;
} c2low_t;

// return value:
// -1  - source consists of repetition of the same character
// >=0 - maxC = the character with the biggest numeric value that appears in the source at least once
static int prepare1(const uint8_t* src, unsigned srclen, c2low_t c2low[257], double* pQuantizedEntropy, double* pInfo)
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

  // translate counts to ranges and store in c2low
  // 1st pass - translate characters with counts that are significantly lower than maxCnt
  unsigned thr = maxCnt - maxCnt/8;
  unsigned remCnt   = srclen;
  unsigned remRange = VAL_RANGE;
  for (unsigned c = 0; c <= maxC; ++c) {
    unsigned cnt = stat[c];
    if (cnt < thr) {
      unsigned range = 0;
      if (cnt != 0) {
        // calculate range from full statistics
        range = (uint64_t(cnt)*(VAL_RANGE*2) + srclen)/(srclen*2);
        if (range == 0)
          range = 1;
        remCnt   -= cnt;
        remRange -= range;
      }
      c2low[c].u = range;
    }
  }
  // 2nd pass - translate characters with higher counts
  for (unsigned c = 0; remCnt != 0; ++c) {
    unsigned cnt = stat[c];
    if (cnt >= thr) {
      // Calculate range from the remaining range and count
      // It is non-ideal, but this way we distribute the worst rounding errors
      // relatively evenly among higher ranges, where it hase the smallest impact
      unsigned range = (uint64_t(cnt)*(remRange*2) + remCnt)/(remCnt*2);
      // (range < VAL_RANGE) is guaranteed , because we already handled the case of repetition of the same character
      remCnt   -= cnt;
      remRange -= range;
      c2low[c].u = range;
    }
  }

  // calculate entropy after quantization
  double entropy = 0;
  for (unsigned c = 0; c <= maxC; ++c) {
    unsigned cnt = stat[c];
    if (cnt)
      entropy += log2(double(VAL_RANGE)/c2low[c].u)*cnt;
  }
  *pQuantizedEntropy = entropy;
  if (pInfo)
    pInfo[3] = entropy;

  return maxC;
}

static void prepare2(c2low_t c2low[257], unsigned maxC)
{
  // c2low -> cumulative sums of ranges
  unsigned lo = 0;
  for (unsigned c = 0; c <= maxC; ++c) {
    unsigned range = c2low[c].u;
    // printf("%3u: %04x\n", c, range);
    c2low[c].f = lo*INV_VAL_RANGE;
    lo += range;
  }
  c2low[maxC+1].f = 1.0f;
}

static int inline floor_log4(unsigned x) {
  return unsigned(sizeof(unsigned)*8 - 1 - __builtin_clz(x))/2; // floor(log4(x))
}

// return the number of stored octets
static int store_model(uint8_t* dst, const c2low_t c2low[256], unsigned maxC, double* pInfo)
{
  uint16_t c2range[256];
  // copy and count non-zero ranges
  int nRanges = 0;
  unsigned hist[8] = {0};
  for (unsigned c = 0; c < maxC; ++c) {
    uint32_t range = c2low[c].u;
    c2range[c] = range;
    if (range > 0) {
      ++nRanges;
      hist[1+floor_log4(range)] += 1;
    }
  }
  // No need to store a last range.
  // Decoder will calculated is as VAL_RANGE-sum of previous ranges
  hist[0] = maxC - nRanges; // # of zero ranges before maxC

  uint8_t* p = dst;
  *p++ = maxC;
  unsigned rem = maxC;
  for (int i = 0; i < 7 && rem > 0; ++i) {
    rem -= hist[i];
    *p++ = hist[i];
  }

  // arithmetic encoder
  const uint64_t MSK31_0  = (uint64_t(1) << 32)-(uint64_t(1) << 0);
  const uint64_t MSK39_0  = (uint64_t(1) << 40)-(uint64_t(1) << 0);
  const uint64_t MSK47_40 = (uint64_t(1) << 48)-(uint64_t(1) << 32);
  uint64_t lo    = 0;                  // 48 bits
  uint64_t range = uint64_t(-1) >> 16; // in fact, range-1
  int pending_bytes = 0;
  int ix = 0;
  int phase0 = 1;
  for (unsigned c = 0; c < maxC; c += phase0) {
    // prevent range from becoming too small
    while (range < (uint64_t(1) << 28)) {
      // squeeze out all ones in bits[39..32]
      lo = (lo & MSK47_40) | ((lo & MSK31_0) << 8);
      range = (range << 8) | 255;
      ++pending_bytes;
    }

    uint32_t val = c2range[c];
    if (phase0) {
      uint32_t cnt  = 0;
      ix = 0;
      // printf("%3u lo=%010llx ra=%010llx <= %04x\n", c, lo, range, val);
      if (val != 0) {
        ix = 1+floor_log4(val);
        for (int i = 0; i < ix; ++i)
          cnt += hist[i];
      }

      // insert code of ix
      uint32_t cnt1 = hist[ix];
      hist[ix] -= 1;

      uint32_t l_num = cnt;    // up to 255
      uint32_t r_num = cnt1;   // up to 255
      uint32_t den = (maxC-c); // up to 255

      lo   += ((range+1) * l_num + den - 1)/den; // ceil
      range = ((range+1) * r_num)/den - 1;       // floor

      if (val != 0)
        phase0 = 0; // stage insertion of code of specific value within hist range
    } else {
      // insert code of specific value within hist range
      int lshift = (ix-1)*2;
      uint32_t valNum = val - (1u << lshift); // up to 2^14-2^12-1
      uint32_t valDen = 3u << lshift;         // up to 2^14-2^12

      lo   += ((range+1) * valNum + valDen - 1)/valDen; // ceil
      range = (range+1)/valDen - 1;                     // floor
      phase0 = 1;
    }

    uint64_t hi = lo + range;
    uint64_t dbits = lo ^ hi;
    while ((dbits >> 40)==0) {
      // lo and hi have 8 common MS bits
      *p++ = uint8_t(lo >> 40);
      if (pending_bytes) {
        uint8_t pending_byte = 0-((lo>>39) & 1);
        do {
          *p++ = pending_byte;
          --pending_bytes;
        } while (pending_bytes);
      }
      lo    = (lo & MSK39_0) << 8;
      range = (range << 8) | 255;
      dbits <<= 8;
    }
  }

  // last octet(s)
  // We want a shift register at decoder to be in range [lo..hi]
  // regardless of alien characters that a possibly resides after
  // our code stream
  uint64_t lastOctet = ((lo + range) >> 40)-1;
  *p++ = uint8_t(lastOctet);
  while (lastOctet == (lo>>40)) {
    lastOctet = 255;
    *p++ = uint8_t(lastOctet);
    lo = (lo & MSK39_0) << 8;
  }

  int len = p - dst;
  if (pInfo)
    pInfo[1] = len*8;
  return len;
}

// static uint64_t to_u(double x){
  // uint64_t u;
  // memcpy(&u, &x, sizeof(u));
  // return u;
// }

static void inc_dst(uint8_t* dst) {
  uint8_t val;
  do {
    --dst;
    *dst = (val = *dst + 1);
  } while (val==0);
}

int storeLastOctets(uint8_t* dst, double lo, double hi) {
  const uint64_t MSB_MSK = uint64_t(-1) << 40; // 16 overhead bits + 8 most significant bits
  uint64_t uLo; memcpy(&uLo, &lo, sizeof(uint64_t));
  uint64_t uHi; memcpy(&uHi, &hi, sizeof(uint64_t));
  dst[5] = uHi >> (0*8);
  dst[4] = uHi >> (1*8);
  dst[3] = uHi >> (2*8);
  dst[2] = uHi >> (3*8);
  dst[1] = uHi >> (4*8);
  dst[0] = uHi >> (5*8);
  uint64_t dbits = uLo ^ uHi;
  if (dbits==0)
    return 6;
  int ret = 0;
  while ((dbits & MSB_MSK)==0) {
    // lo and hi have the same upper octet
    ++ret;
    dbits <<= 8;
  }
  return ret + 1;
}

static int encode(uint8_t* dst, const uint8_t* src, unsigned srclen, const c2low_t c2low[257], int maxC)
{
  const double MIN_RANGE      = 1.0 /(uint64_t(1) << 49); // 2**(-49)
  const double LO_INCR_OFFSET = MIN_RANGE;                // 2**(-49)
  const double MAX_LO_L       = MIN_RANGE*0.5;            // 2**(-50)
  double lo_h  = 1.0;   // upper part of accumulator + offset
  double lo_l  = 0.0;   // lower part of accumulator
  double range = 1.0;
  uint8_t* dst0 = dst;
  int64_t prevDstW = -1;
  for (unsigned i = 0; i < srclen; ++i) {
    if (lo_l >= MAX_LO_L) {
      double nxtLo = lo_h + lo_l; lo_l -= nxtLo - lo_h; lo_h = nxtLo;
    }

    range *= RANGE_LEAKAGE_FACTOR;
    int c = src[i];
    float cLo = c2low[c+0].f;
    float cRa = c2low[c+1].f - cLo;
    double loIncr = range * cLo;
    // Reduce precision of loIncr in controlled manner, in order to prevent uncontrolled leakage of LS bits
    // After reduction loIncr still contains no less than 39 significant bits, so efficiency of compression does not suffer.
    loIncr += LO_INCR_OFFSET; loIncr -= LO_INCR_OFFSET;

    // if (i < 8420) {
      // fesetround(FE_TONEAREST);
      // printf("[%d]=%3d: lo_h %.20e lo_l %.20e range %.20e. Incr %.20e\n", i, c, lo_h, lo_l, range, loIncr);
      // fesetround(FE_TOWARDZERO);
    // }

    // dual-double style addition
    double nxtLo = lo_h + loIncr;
    loIncr -= (nxtLo-lo_h);
    lo_l += loIncr;
    lo_h = nxtLo;
    range *= cRa;

    if (range <= MIN_RANGE) {
      // re-normalize
      {double nxtLo = lo_h + lo_l; lo_l -= nxtLo - lo_h; lo_h = nxtLo;}
      lo_h  *= double(uint64_t(1) << 48);
      int64_t dstW = lo_h;
      lo_h -= dstW-1;
      lo_l  *= double(uint64_t(1) << 48);
      range *= double(uint64_t(1) << 48);
      {double nxtLo = lo_h + lo_l; lo_l -= nxtLo - lo_h; lo_h = nxtLo;}
      dstW -= (int64_t(1)<<48);
      if (_LIKELY(prevDstW >= 0)) {
        prevDstW += (dstW >> 48);
        if (_UNLIKELY(prevDstW >= (int64_t(1)<<48)))
          inc_dst(dst); // very unlikely
        dstW &= uint64_t(-1) >> (64-48);
        // store 6 octets
        dst[5] = prevDstW >> (0*8);
        dst[4] = prevDstW >> (1*8);
        dst[3] = prevDstW >> (2*8);
        dst[2] = prevDstW >> (3*8);
        dst[1] = prevDstW >> (4*8);
        dst[0] = prevDstW >> (5*8);
        dst += 6;
      }
      prevDstW = dstW;
    }
  }

  if (prevDstW >= 0) {
    // store 6 octets
    dst[5] = prevDstW >> (0*8);
    dst[4] = prevDstW >> (1*8);
    dst[3] = prevDstW >> (2*8);
    dst[2] = prevDstW >> (3*8);
    dst[1] = prevDstW >> (4*8);
    dst[0] = prevDstW >> (5*8);
    dst += 6;
  }

  // output last bits
  lo_h -= 1.0;
  {double nxtLo = lo_h + lo_l; lo_l -= nxtLo - lo_h; lo_h = nxtLo;}
  if (lo_h >= 1.0) {
    inc_dst(dst);
    lo_h -= 1.0;
    {double nxtLo = lo_h + lo_l; lo_l -= nxtLo - lo_h; lo_h = nxtLo;}
  }

  while (range <= (1.0/256)) {
    // octet is common in lo and lo+range-ULP
    lo_h  *= 256.0;
    lo_l  *= 256.0;
    int oct = lo_h;
    *dst++ = oct;
    lo_h -= oct;
    {double nxtLo = lo_h + lo_l; lo_l -= nxtLo - lo_h; lo_h = nxtLo;}
    range *= 256.0;
  }
  // last octet
  int oct = lo_h*256.0 + 1.0;
  if (oct < 256)
    *dst++ = oct;
  else
    inc_dst(dst);

  return dst - dst0;
}

// return value:
// -1 - source consists of repetition of the same character
//  0 - not compressible, because all input characters have approximately equal probability or because input is too short
// >0 - the length of compressed buffer
int arithmetic_encode(std::vector<uint8_t>* dst, const uint8_t* src, int srclen, double* pInfo)
{
  c2low_t c2low[257];
  double quantizedEntropy;
  int maxC = prepare1(src, srclen, c2low, &quantizedEntropy, pInfo);

  if (maxC == -1)
    return -1; // source consists of repetition of the same character

  size_t sz0 = dst->size();
  dst->resize(sz0 + 640);
  unsigned modellen = store_model(&dst->at(sz0), c2low, maxC, pInfo);

  if ((quantizedEntropy+7)/8 + modellen >= srclen)
    return 0; // not compressible

  int rdir = fegetround();
  fesetround(FE_TOWARDZERO);
  prepare2(c2low, maxC);
  // printf("ml=%u\n", modellen);
  dst->resize(sz0 + modellen + int(quantizedEntropy/8)+64);
  fesetround(FE_TOWARDZERO);
  int dstlen = encode(&dst->at(sz0+modellen), src, srclen, c2low, maxC);
  fesetround(rdir);

  if (pInfo)
    pInfo[2] = dstlen*8.0;

  int reslen = modellen + dstlen;

  if (reslen >= srclen)
    return 0; // not compressible

  dst->resize(sz0 + reslen);
  return reslen;
}
