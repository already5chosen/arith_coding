// #include <immintrin.h>
// #include <cstdio>
#include <cstdint>
#include <cstring>
#include <cfenv>
#include <cmath>
// #include <cctype>

#include "arithmetic_decode.h"

#define _LIKELY(cond)   __builtin_expect((cond), 1)
#define _UNLIKELY(cond) __builtin_expect((cond), 0)

static const int RANGE_BITS = 14;
static const unsigned VAL_RANGE = 1u << RANGE_BITS;
static const double VAL_RATIO_ADJ = 1.0+1.0/(uint64_t(1) << 32);

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
  uint16_t m_c2low[258];
  double m_c2invRange[256]; // 1/range

  int  load_and_prepare(const uint8_t* src, int srclen, int* pInfo);
  int  decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen);
private:
  static const int RANGE2C_NBITS = 9;
  static const int RANGE2C_SZ = 1 << RANGE2C_NBITS;
  uint8_t  m_range2c[RANGE2C_SZ+1];
  unsigned m_maxC;

  void prepare(uint16_t c2range[256]);
  int val2c(double valRatio) {
    int lo = valRatio;
    unsigned ri = lo >> (RANGE_BITS-RANGE2C_NBITS);
    unsigned c = m_range2c[ri]; // c is the biggest character for which m_c2low[c] <= (val/32)*32
    if (m_c2low[c+1] <= lo) {
      do {
        #ifdef ENABLE_PERF_COUNT
        ++m_extLookupCount;
        #endif
        ++c;
      } while (m_c2low[c+1] <= lo);
    }
    return c;
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

void arithmetic_decode_model_t::prepare(uint16_t c2range[256])
{
  // m_c2low -> cumulative sums of ranges
  unsigned maxC = 255;
  unsigned lo = 0;
  unsigned invI = 0;
  for (int c = 0; c < 256; ++c) {
    unsigned range = c2range[c];
    // printf("%3u: %04x %04x\n", c, range, lo);
    m_c2low[c] = lo;
    if (range != 0) {
      m_c2invRange[c] = double(VAL_RANGE)/range;
      // build inverse index m_range2c
      maxC = c;
      lo += range;
      for (; invI <= ((lo-1) >> (RANGE_BITS-RANGE2C_NBITS)); ++invI) {
        m_range2c[invI] = maxC;
      }
    }
  }
  m_c2low[maxC + 1] = VAL_RANGE;
  for (unsigned c = maxC + 2; c < 258; ++c)
    m_c2low[c] = VAL_RANGE+c;
  for (; invI < RANGE2C_SZ+1; ++invI) {
    m_range2c[invI] = maxC;
  }
  m_maxC = maxC;
}

// static uint64_t to_u(double x){
  // uint64_t u;
  // memcpy(&u, &x, sizeof(u));
  // return u;
// }

// static int64_t load48bits(const uint8_t* src) {
  // int64_t value = 0;
  // for (int k = 0; k < 6; ++k)
    // value = (value << 8) | src[k];
  // return value;
// }

static inline double and_sd(double x, uint64_t mask) {
  uint64_t u;
  memcpy(&u, &x, sizeof(u));
  u &= mask;
  memcpy(&x, &u, sizeof(u));
  return x;
}

int arithmetic_decode_model_t::decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen)
{
  const double MIN_RANGE      = 1.0 /(uint64_t(1) << 49); // 2**(-49)
  const double TWO_POW48 = int64_t(1) << 48;
  const double TWO_POWn48 = 1.0/TWO_POW48;
  const double TWO_POWn56 = 1.0/(int64_t(1) << 56);
  const double TWO_POWn104 = TWO_POWn48*TWO_POWn56;

  if (srclen < 2)
    return -11;

  uint8_t tmpbuf[32] = {0};
  bool useTmpbuf = false;
  if (srclen < 13) {
    memcpy(tmpbuf, src, srclen);
    src = tmpbuf;
    useTmpbuf = true;
  }

  double val_h = 0.0, val_l = 0;
  for (int k = 0; k < 13; ++k) {
    double x = src[12-k];
    {double s = val_h + x; double d1 = s - x; x = val_h - d1; val_h = s; }
    val_l += x;
    { double nxtVal = val_h + val_l; val_l -= nxtVal - val_h; val_h = nxtVal; }
    val_h *= 1.0/256;
    val_l *= 1.0/256;
  }
  src    += 13;
  srclen -= 13;
  double range    = 1.0;
  double invRange = (VAL_RANGE*VAL_RATIO_ADJ);
  for (int i = 0; ; ) {

    if ((i & 63)==0) {
      // invRange = 1.0/range;
      invRange *= 2.0-range*invRange*(1.0/VAL_RATIO_ADJ/VAL_RANGE);
    }

    if (val_h >= range) return -12;
    // That is an input error, rather than internal error of decoder.
    // Due to way that encoder works, not any bit stream is possible as its output
    // That's the case of illegal code stream.
    // The case is extremely unlikely, but not impossible.

    unsigned c = val2c(val_h*invRange); // can be off by +1

    // if (i >= 0*4600 && i < 4650) {
      // double ra = range * (1.0/VAL_RANGE);
      // int64_t cLo = m_c2low[c+0];
      // double valDecr = ra * cLo;
      // double valDecrI = (long double)ra * cLo - valDecr;
      // fesetround(FE_TONEAREST);
      // printf("[%2d]=%3d: val_h %.20e val_l %.20e range %.20e. [%5d..%5d) %23.17f  %23.17f %23.17f.  %.20e %.20e\n"
      // , i, c, val_h, val_l, range
      // , m_c2low[c+0], m_c2low[c+1]
      // , val_h*invRange
      // , val_h/range*VAL_RANGE*VAL_RATIO_ADJ
      // , val_h/range*VAL_RANGE
      // , valDecr
      // , valDecrI
      // );
      // fesetround(FE_TOWARDZERO);
    // }

    // keep decoder in sync with encoder
    range *= (1.0/VAL_RANGE);
    int64_t cLo = m_c2low[c+0];
    double valDecr = range * cLo;
    if (_UNLIKELY(valDecr > val_h)) {
      cLo = m_c2low[c-1];
      --c;
      valDecr = range * cLo;
      if (valDecr > val_h) {
        return -102; // should never happen
      }
    }
    int64_t cRa = m_c2low[c+1] - cLo;
    range *= cRa;
    // Reduce precision of valDecr in controlled manner, in order to prevent uncontrolled leakage of LS bits
    range = and_sd(range, -int64_t(VAL_RANGE));
    dst[i] = c;

    // dual-double style subtraction
    double nxtVal = val_h - valDecr; valDecr += (nxtVal - val_h); val_l -= valDecr; val_h = nxtVal;
    nxtVal = val_h + val_l; val_l -= (nxtVal - val_h); val_h = nxtVal;

    ++i;
    if (i == dstlen)
      break; // success

    // update invRange
    invRange *= m_c2invRange[c];

    // printf("was here A2. %016llx %016llx %016llx\n",  lo, value, hi); fflush(stdout);
    if (range <= MIN_RANGE) {
      // printf("%d: renorm\n", i);
      #ifdef ENABLE_PERF_COUNT
      ++m_renormalizationCount;
      #endif
      if (srclen < 6) {
        if (!useTmpbuf) {
          if (srclen > 0)
            memcpy(tmpbuf, src, srclen);
          src = tmpbuf;
          useTmpbuf = true;
        } else {
          if (src > &tmpbuf[32-6])
            return i;
        }
      }
      int64_t sixOctets =
        (int64_t(src[5])       ) |
        (int64_t(src[4]) << 1*8) |
        (int64_t(src[3]) << 2*8) |
        (int64_t(src[2]) << 3*8) |
        (int64_t(src[1]) << 4*8) |
        (int64_t(src[0]) << 5*8);
      src    += 6;
      srclen -= 6;
      // if (i > 4600 && i < 4650) {
        // fesetround(FE_TONEAREST);
        // printf("%012llx\n", sixOctets);
        // fesetround(FE_TOWARDZERO);
      // }
      double valIncr = sixOctets*TWO_POWn104;
      { double nxtVal = val_h + val_l; val_l -= nxtVal - val_h; val_h = nxtVal; }
      val_h *= TWO_POW48;
      val_l *= TWO_POW48;
      // dual-double style addition
      {double nxtVal = val_h + valIncr; valIncr -= nxtVal - val_h; val_h = nxtVal;}
      val_l += valIncr;
      invRange *= TWO_POWn48;
      range    *= TWO_POW48;
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
  }
  fesetround(rdir);
  return modellen;
}

