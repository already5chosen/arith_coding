#include <cstdint>
#include <cstring>
#include <cstdio>
// #include <cmath>
#include <cctype>

#include "arithmetic_decode.h"


static const int RANGE_BITS = 14;
static const unsigned VAL_RANGE = 1u << RANGE_BITS;


// load_ranges
// return  value:
//   on success the # of processed source octets,
//   on parsing error negative error code
//static
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
  if (hist[0]==maxC) return -2; // inconsistent histogram table
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

namespace {

struct arithmetic_decode_model_t {
  #ifdef ENABLE_PERF_COUNT
  int m_longLookupCount;
  int m_renormalizationCount;
  #endif
  uint16_t m_c2low[257];
  float    m_c2invRange[256]; // (VAL_RANGE/range) rounded down

  int  load_and_prepare(const uint8_t* src, int srclen, int* pInfo);
  int  decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen);
private:
  static const int RANGE2C_NBITS = 9;
  static const int RANGE2C_SZ = 1 << RANGE2C_NBITS;
  uint8_t  m_range2c[RANGE2C_SZ+1];
  unsigned m_maxC;

  void prepare();
  int val2c(uint64_t value, uint64_t range, double invRange) {
    #if 1
    int ri = int(int64_t(value<<RANGE2C_NBITS)*invRange);
    #else
    // this trick relies on knowledge of IEEE-754 binary64 format
    // It works, but the speedup is not significant
    double dRi = int64_t(value)*invRange + 1.0;
    uint64_t u64ri;
    memcpy(&u64ri, &dRi, sizeof(u64ri));
    int ri = (u64ri >> (52 - RANGE2C_NBITS)) & (RANGE2C_SZ-1);
    #endif
    // unsigned ri = (value*RANGE2C_SZ)/range;
    unsigned c = m_range2c[ri]; // c is the biggest character for which m_c2low[c] <= (val/32)*32
    if (c != m_range2c[ri+1]) {
      #ifdef ENABLE_PERF_COUNT
      ++m_longLookupCount;
      #endif
      while (((m_c2low[c+1]*range+VAL_RANGE-1) >> RANGE_BITS) <= value)
        ++c;
    }
    return c;
  }
};


int arithmetic_decode_model_t::load_and_prepare(const uint8_t* src, int srclen, int* pInfo)
{
  #ifdef ENABLE_PERF_COUNT
  m_longLookupCount = 0;
  m_renormalizationCount = 0;
  #endif
  int ret = load_ranges(m_c2low, src, srclen, pInfo);
  if (ret >= 0)
    prepare();
  return ret;
}

void arithmetic_decode_model_t::prepare()
{
  // m_c2low -> cumulative sums of ranges
  unsigned maxC = 255;
  unsigned lo = 0;
  unsigned invI = 0;
  for (int c = 0; c < 256; ++c) {
    unsigned range = m_c2low[c];
    // printf("%3u: %04x %04x\n", c, range, lo);
    m_c2low[c] = lo;
    if (range != 0) {
      const float VAL_RANGE_MINUS_EPS = float(0.999999940395355*VAL_RANGE);
      m_c2invRange[c] = VAL_RANGE_MINUS_EPS/range; // (VAL_RANGE/range) rounded down
      // build inverse index m_range2c
      maxC = c;
      lo += range;
      for (; invI <= ((lo-1) >> (RANGE_BITS-RANGE2C_NBITS)); ++invI) {
        m_range2c[invI] = maxC;
      }
    }
  }
  m_c2low[256] = VAL_RANGE;
  for (; invI <= RANGE2C_SZ; ++invI) {
    m_range2c[invI] = maxC;
  }
  m_maxC = maxC;
}

int arithmetic_decode_model_t::decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen)
{
  uint64_t VAL_MSK  = uint64_t(-1) >> 16;
  uint64_t MSK31_0  = (uint64_t(1) << 32)-(uint64_t(1) << 0);
  uint64_t MSK47_40 = (uint64_t(1) << 48)-(uint64_t(1) << 40);

  if (srclen < 2)
    return -11;

  uint8_t tmpbuf[16] = {0};
  bool useTmpbuf = false;
  if (srclen < 6) {
    memcpy(tmpbuf, src, srclen);
    src = tmpbuf;
    useTmpbuf = true;
  }

  uint64_t value = 0;
  for (int k = 0; k < 6; ++k)
    value = (value << 8) | src[k];
  src    += 6;
  srclen -= 6;

  uint64_t lo = 0;
  int64_t range  = int64_t(1) << 48; // scaled to 2**48. Maintained in (2**31..2**48]
  double invRange = 1.0/range;       // approximation of (1/range) rounded down. Maintained in [2**48..2**-31)
  for (int i = 0; ; ) {
    unsigned c = val2c(value - lo, range, invRange);
    // if (i % 1000 < 10 || i > dstlen - 20000) {
      // printf(": %6d [%012llx %012llx] %012llx %012llx. %08x c=%3d '%c' [%04x..%04x) => [%012llx %012llx]\n"
       // , i, lo, hi, value, range
       // , unsigned(round(double(value - lo)*0x100000000/range))
       // , c, isprint(c) ? c : '.', m_c2low[c+0], m_c2low[c+1]
       // , lo + ((range * m_c2low[c+0])>>15)
       // , lo + ((range * m_c2low[c+1])>>15) - 1);
      // fflush(stdout);
    // }
    dst[i] = c;
    ++i;
    if (i == dstlen)
      break; // success

    // keep decoder in sync with encoder
    uint32_t cLo = m_c2low[c+0];
    uint32_t cHi = m_c2low[c+1];
    lo   += (range * cLo + VAL_RANGE-1) >> RANGE_BITS;
    range = (range * (cHi-cLo)) >> RANGE_BITS;

    if (value - lo >= uint64_t(range)) return -102; // should not happen

    // update invRange
    invRange = invRange*m_c2invRange[c];

    if (range <= (1u << 31)) {
      #ifdef ENABLE_PERF_COUNT
      ++m_renormalizationCount;
      #endif
      if (srclen < 6) {
        if (!useTmpbuf && srclen > 0) {
          memcpy(tmpbuf, src, srclen);
          src = tmpbuf;
          useTmpbuf = true;
        }
      }
      uint64_t hi = lo + range -1;
      uint64_t dbits = lo ^ hi;
      int octetsShifted = 0;
      for (; (dbits >> 40)==0; ++octetsShifted) {
        // lo and hi have the same upper octet
        value = (value << 8) | src[octetsShifted];
        dbits <<= 8;
      }
      value &= VAL_MSK;
      range <<= octetsShifted*8;
      lo = (lo << octetsShifted*8) & VAL_MSK;
      for (; range <= (1u << 31); ++octetsShifted) {
        // squeeze out bits[39..32]
        lo    = (lo    & MSK47_40) | ((lo    & MSK31_0) << 8);
        value = (value & MSK47_40) | ((value & MSK31_0) << 8);
        value |= src[octetsShifted];
        range <<= 8;
      }
      static const double invRangeDivTab[7] = {
        1./(int64_t(1)<<(8*0)), 1./(int64_t(1)<<(8*1)),
        1./(int64_t(1)<<(8*2)), 1./(int64_t(1)<<(8*3)),
        1./(int64_t(1)<<(8*4)), 1./(int64_t(1)<<(8*5)),
        1./(int64_t(1)<<(8*6)),
      };
      invRange *= invRangeDivTab[octetsShifted];
      // do one NR iteration to increase precision of invRange and assure that invRange*range <= 2**64
      invRange *= 2.0 - invRange * range;

      src    += octetsShifted;
      srclen -= octetsShifted;
      if (srclen < -5)
        return i;
      if (value - lo >= uint64_t(range)) return -103; // should not happen
    }
  }
  return dstlen;
}

}

int arithmetic_decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen, int* pInfo)
{
  arithmetic_decode_model_t model;
  int modellen = model.load_and_prepare(&src[0], srclen, pInfo);
  if (modellen >= 0) {
    if (pInfo)
      pInfo[1] = srclen-modellen;
    int textlen = model.decode(dst, dstlen, &src[modellen], srclen-modellen);
    #ifdef ENABLE_PERF_COUNT
    if (pInfo) {
      pInfo[2] = model.m_longLookupCount;
      pInfo[3] = model.m_renormalizationCount;
    }
    #endif
    return textlen;
  }
  return modellen;
}

