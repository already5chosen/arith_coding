// #include <cstdio>
#include <cstdint>
#include <cstring>
#include <cfenv>
#include <cmath>

#include "arithmetic_decode.h"

#define _LIKELY(cond)   __builtin_expect((cond), 1)
#define _UNLIKELY(cond) __builtin_expect((cond), 0)

static const int RANGE_BITS = 14;
static const unsigned VAL_RANGE = 1u << RANGE_BITS;

// load_ranges
// return  value:
//   on success the # of processed source octets,
//   on parsing error negative error code
static
int load_ranges(uint16_t* ranges, const uint8_t* src, int srclen, int* pInfo)
{
  if (srclen == 0) return -1;
  int maxC = src[0];
  int srcI = 1;
  const int HIST_LEN = (RANGE_BITS+1)/2 + 1;
  int hist[HIST_LEN] = {0};
  int rem = maxC;
  for (int i = 0; i < HIST_LEN-1 && rem > 0; ++i) {
    if (srcI == srclen) return -1;
    hist[i] = src[srcI++];
    rem -= hist[i];
  }
  if (rem < 0) return -2; // inconsistent histogram table
  if (rem > 0)
    hist[HIST_LEN-1] = rem;

  // arithmetic decoder
  const uint64_t MSK31_0  = (uint64_t(1) << 32)-(uint64_t(1) << 0);
  const uint64_t MSK39_0  = (uint64_t(1) << 40)-(uint64_t(1) << 0);
  const uint64_t MSK47_40 = (uint64_t(1) << 48)-(uint64_t(1) << 32);
  uint64_t lo    = 0;                  // 48 bits
  uint64_t range = uint64_t(-1) >> 16; // in fact, range-1

  uint64_t value = 0;  // 40 bits
  for (int i = 0; i < 6; ++i, ++srcI) {
    uint8_t srcO = srcI < srclen ? src[srcI] : 0;
    value = (value << 8) | srcO;
  }

  uint32_t sum    = 0;
  int ix = 0;
  int phase0 = 1;
  for (int c = 0; c < maxC; c += phase0) {
    // prevent range from becoming too small
    while (range < (uint64_t(1) << 28)) {
      // squeeze out all ones in bits[39..32]
      lo = (lo & MSK47_40) | ((lo & MSK31_0) << 8);
      range = (range << 8) | 255;
      value = (value & MSK47_40) | ((value & MSK31_0) << 8);
      value |= (srcI < srclen) ? src[srcI] : 0;
      ++srcI;
      if (value < lo || value > lo + range) return -201; // should not happen
    }

    uint64_t deltaV = value - lo;
    if (phase0) {
      // printf("%3u lo=%010llx ra=%010llx va=%010llx\n", c, lo, range, value); fflush(stdout);
      // course search
      uint32_t cnt = 0;
      for (ix = 0; deltaV*(maxC-c) >= (cnt+hist[ix])*(range+1) ; ++ix)
        cnt += hist[ix];

      // insert code of ix
      uint32_t cnt1 = hist[ix];
      hist[ix] -= 1;

      uint32_t l_num = cnt;  // up to 255
      uint32_t r_num = cnt1; // up to 255
      uint32_t den = (maxC-c); // up to 255

      lo   += ((range+1) * l_num + den - 1)/den; // ceil
      range = ((range+1) * r_num)/den - 1;       // floor

      ranges[c] = 0;
      if (ix != 0)
        phase0 = 0; // stage look up for specific value within hist[ix] range
    } else {
      // fine search
      int lshift = (ix-1)*2;
      uint32_t val0   = 1u << lshift; // up to 2^12
      uint32_t valDen = 3u << lshift; // up to 2^14-2^12

      uint32_t valNum = (deltaV*valDen)/(range+1); // floor
      // insert code of specific value within hist range
      lo   += ((range+1) * valNum + valDen - 1)/valDen; // ceil
      range = (range+1)/valDen - 1;                     // floor
      phase0 = 1;

      uint32_t val = val0 + valNum;
      ranges[c] = val;
      sum += val;
    }

    uint64_t hi = lo + range;
    uint64_t dbits = lo ^ hi;
    while ((dbits >> 40)==0) {
      // lo and hi have 8 common MS bits
      lo    = (lo & MSK39_0) << 8;
      range = (range << 8) | 255;
      value = (value & MSK39_0) << 8;
      value |= (srcI < srclen) ? src[srcI] : 0;
      ++srcI;
      if (srcI > srclen+5) return -3;
      if (value < lo || value > lo + range) return -202; // should not happen
      dbits <<= 8;
    }
  }

  srcI -= (6-1);
  // Imitate encoder's logic for last octets in order to find out
  // an exact length of encoded stream.
  uint64_t lastOctet = ((lo + range) >> 40)-1;
  while (lastOctet == (lo>>40)) {
    lastOctet = 255;
    ++srcI;
    lo = (lo & MSK39_0) << 8;
  }
  if (srcI > srclen) return -4;

  if (sum >= VAL_RANGE)
    return -5; // parsing error

  ranges[maxC] = VAL_RANGE - sum;
  for (unsigned c = maxC+1; c < 256; ++c)
    ranges[c] = 0;

  // for (int c = 0; c <= maxC; ++c)
    // printf("%3u: %04x\n", c, ranges[c]);
  // printf("len=%d, sum=%04x\n", srcI, sum);

  if (pInfo)
    pInfo[0] = srcI*8;

  return srcI; // success
}

static inline uint64_t umulh(uint64_t a, uint64_t b) {
  return ((unsigned __int128)a * b) >> 64;
}

namespace {

struct arithmetic_decode_model_t {
  #ifdef ENABLE_PERF_COUNT
  int m_extLookupCount;
  int m_longLookupCount;
  int m_renormalizationCount;
  #endif
  uint16_t m_ci2lo[257];
  double   m_ci2invRange[256]; // 1/range
  uint8_t  m_ci2c[256];

  int  load_and_prepare(const uint8_t* src, int srclen, int* pInfo);
  int  decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen);
private:
  static const int LO2CI_NBITS = 9;
  static const int LO2CI_SZ = 1 << LO2CI_NBITS;
  uint8_t  m_lo2ci[LO2CI_SZ+1];
  unsigned m_maxCi;

  void prepare(const uint16_t c2range[256]);
  int val2ci(double valRatio)
  #ifndef ENABLE_PERF_COUNT
  const
  #endif
  {
    int lo = valRatio;
    unsigned ri = lo >> (RANGE_BITS-LO2CI_NBITS);
    unsigned ci = m_lo2ci[ri]; // ci is the biggest character index for which m_ci2lo[ci] <= (lo/32)*32
    if (m_ci2lo[ci+1] <= lo) {
      do {
        #ifdef ENABLE_PERF_COUNT
        ++m_extLookupCount;
        #endif
        ++ci;
      } while (m_ci2lo[ci+1] <= lo);
    }
    return ci;
  }
};


int arithmetic_decode_model_t::load_and_prepare(const uint8_t* src, int srclen, int* pInfo)
{
  #ifdef ENABLE_PERF_COUNT
  m_extLookupCount = 0;
  m_longLookupCount = 0;
  m_renormalizationCount = 0;
  #endif
  uint16_t c2range[256];
  int ret = load_ranges(c2range, src, srclen, pInfo);
  if (ret >= 0)
    prepare(c2range);
  return ret;
}

void arithmetic_decode_model_t::prepare(const uint16_t c2range[256])
{
  // m_c2low -> cumulative sums of ranges
  unsigned lo   = 0;
  unsigned invI = 0;
  unsigned ci   = 0;
  for (int c = 0; c < 256; ++c) {
    unsigned range = c2range[c];
    // printf("%3u: %04x %04x\n", c, range, lo);
    if (range != 0) {
      m_ci2lo[ci] = lo;
      m_ci2c[ci]  = c;
      m_ci2invRange[ci] = double(VAL_RANGE)/range;
      // build inverse index m_lo2ci
      lo += range;
      for (; invI <= ((lo-1) >> (RANGE_BITS-LO2CI_NBITS)); ++invI) {
        m_lo2ci[invI] = ci;
      }
      ++ci;
    }
  }
  for (unsigned i = ci; i < 257; ++i)
    m_ci2lo[i] = VAL_RANGE;
  for (; invI < LO2CI_SZ+1; ++invI)
    m_lo2ci[invI] = ci-1;
  m_maxCi = ci-1;
}

// static uint64_t to_u(double x){
  // uint64_t u;
  // memcpy(&u, &x, sizeof(u));
  // return u;
// }

int arithmetic_decode_model_t::decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen)
{
// const uint8_t* src0 = src;
  const double MIN_RANGE  = 1.0 /(uint64_t(1) << 17); // 2**(-17)
  const double TWO_POW16  = int64_t(1) << 16;
  const double TWO_POWn16 = 1.0/TWO_POW16;
  const double TWO_POWn52 = 1.0/(int64_t(1) << 52);

  if (srclen < 2)
    return -11;

  uint8_t tmpbuf[16] = {0};
  bool useTmpbuf = false;
  if (srclen < 7) {
    memcpy(tmpbuf, src, srclen);
    src = tmpbuf;
    useTmpbuf = true;
  }

  double val = 0.0;
  uint64_t val56 = 0;
  for (int k = 0; k < 7; ++k)
    val56 = val56*256+src[k];
  val = (val56 >> 4)*TWO_POWn52;
  src    += 6;
  srclen -= 6;
  double range    = 1.0;
  double invRange = (VAL_RANGE);
  for (int i = 0; ; ) {

    if ((i & 31)==0) {
      // invRange = 1.0/range;
      invRange *= 2.0-range*invRange*(1.0/VAL_RANGE);
    }

    if (val >= range) return -12;
    // That is an input error, rather than internal error of decoder.
    // Due to way that encoder works, not any bit stream is possible as its output
    // That's the case of illegal code stream.
    // The case is extremely unlikely, but not impossible.

    unsigned ci = val2ci(val*invRange); // can be off by -1

    // fesetround(FE_TONEAREST);
    // printf("[%2d]=%3d: val %.20e range %.20e. [%5d..%5d) %23.17f  %23.17f %23.17f %.20e\n"
    // , i, c, val, range, m_c2low[c+0], m_c2low[c+1]
    // , val*invRange
    // , val/range*VAL_RANGE
    // , val/range*VAL_RANGE
    // , range * invRange
    // );
    // fesetround(FE_TOWARDZERO);

    //
    // Keep decoder in sync with encoder
    //
    // Reduce precision of range in controlled manner, in order to prevent uncontrolled leakage of LS bits
    // After reduction range still contains no less than 49 significant bits, so efficiency of compression does not suffer.
    // range = range*(1.0+1.0/VAL_RANGE) - range;
    range *= (1.0/VAL_RANGE);
    range += 1.0;
    range -= 1.0;
    int64_t cLo = m_ci2lo[ci+0];
    int64_t cHi = m_ci2lo[ci+1];
    double valDecr = range * cLo;
    val   -= valDecr;
    double nxtRange = range * (cHi-cLo);
    if (_UNLIKELY(val >= nxtRange)) {
      #ifdef ENABLE_PERF_COUNT
      ++m_longLookupCount;
      #endif
      // fesetround(FE_TONEAREST);
      // printf("[%2d]=%3d: val %.20e nxtRange %.20e [%5d..%5d..%5d)\n", i, c, val, nxtRange, m_c2low[c+1], m_c2low[c+2], m_c2low[c+3]);
      // fesetround(FE_TOWARDZERO);
      val -= nxtRange;
      cLo = cHi;
      cHi = m_ci2lo[ci+2];
      ++ci;
      nxtRange = range * (cHi-cLo);
      if (_UNLIKELY(val >= nxtRange)) {
        return -102; // should never happen
      }
    }
    range = nxtRange;
    // update invRange
    invRange *= m_ci2invRange[ci];

    dst[i] = m_ci2c[ci];
    ++i;
    if (i == dstlen)
      break; // success

    if (range <= MIN_RANGE) {
      // printf("%d: renorm\n", i);
      #ifdef ENABLE_PERF_COUNT
      ++m_renormalizationCount;
      #endif
      if (srclen < 3) {
        if (!useTmpbuf) {
          if (srclen > 0)
            memcpy(tmpbuf, src, srclen);
          src = tmpbuf;
          useTmpbuf = true;
        } else {
          if (src > &tmpbuf[16-3])
            return i;
        }
      }
      int32_t threeOctets =
        (int32_t(src[2])       ) |
        (int32_t(src[1]) << 1*8) |
        (int32_t(src[0]) << 2*8);
      // if (i < 20) {
        // fesetround(FE_TONEAREST);
        // printf("%5d: %06x\n", src-src0, threeOctets);
        // fesetround(FE_TOWARDZERO);
      // }
      src    += 2;
      srclen -= 2;
      const int32_t MASK16 = (int32_t(1)<<16) - 1;
      double valIncr = int64_t((threeOctets>>4) & MASK16) * TWO_POWn52;
      val      *= TWO_POW16;
      val      += valIncr;
      invRange *= TWO_POWn16;
      range    *= TWO_POW16;
    }

  }
  return dstlen;
}

}

int arithmetic_decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen, int* pInfo)
{
  arithmetic_decode_model_t model;
  int rdir = fegetround();
  fesetround(FE_TOWARDZERO);
  int modellen = model.load_and_prepare(&src[0], srclen, pInfo);
  if (modellen >= 0) {
    if (pInfo)
      pInfo[1] = srclen-modellen;
  fesetround(FE_TOWARDZERO);
    int textlen = model.decode(dst, dstlen, &src[modellen], srclen-modellen);
    fesetround(rdir);
    #ifdef ENABLE_PERF_COUNT
    if (pInfo) {
      pInfo[2] = model.m_extLookupCount;
      pInfo[3] = model.m_longLookupCount;
      pInfo[4] = model.m_renormalizationCount;
    }
    #endif
    return textlen;
    // return textlen <= 0 ? textlen : dstlen;
  }
  fesetround(rdir);
  return modellen;
}

