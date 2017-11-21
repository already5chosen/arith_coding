#include <cstdint>
#include <cstring>
#include <cfenv>
// #include <cstdio>
// #include <cmath>
// #include <cctype>

#include "arithmetic_decode.h"


static const int RANGE_BITS = 14;
static const unsigned VAL_RANGE = 1u << RANGE_BITS;
static const float INV_VAL_RANGE = 1.0f / VAL_RANGE;

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
  float m_c2low[257];
  float m_c2invRange[256]; // 1/range

  int  load_and_prepare(const uint8_t* src, int srclen, int* pInfo);
  int  decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen);
private:
  static const int RANGE2C_NBITS = 9;
  static const int RANGE2C_SZ = 1 << RANGE2C_NBITS;
  uint8_t  m_range2c[RANGE2C_SZ+2];
  unsigned m_maxC;

  void prepare(uint16_t c2range[256]);
  int val2c(double value, double lo, double range, float invRange) {
    // uint32_t loEst31b = umulh(value,invRange);
    // unsigned ri = loEst31b>>(31-RANGE2C_NBITS);
    unsigned ri = int(float(value-lo)*invRange*512.0f);
    unsigned c = m_range2c[ri]; // c is the biggest character for which m_c2low[c] <= (val/32)*32
    if (c != m_range2c[ri+2]) {
      #ifdef ENABLE_PERF_COUNT
      ++m_extLookupCount;
      #endif
      // unsigned loEst = loEst31b>>(31-RANGE_BITS);
      // while (m_c2low[c+1] <= loEst) {
        // ++c;
      // }
      while (m_c2low[c+1]*range+lo <= value) {
        #ifdef ENABLE_PERF_COUNT
        ++m_longLookupCount;
        #endif
        ++c;
      }
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
    m_c2low[c] = lo*INV_VAL_RANGE;
    if (range != 0) {
      m_c2invRange[c] = float(VAL_RANGE)/range;
      // build inverse index m_range2c
      maxC = c;
      lo += range;
      for (; invI <= ((lo-1) >> (RANGE_BITS-RANGE2C_NBITS)); ++invI) {
        m_range2c[invI] = maxC;
      }
    }
  }
  m_c2low[256] = 1.0f;
  for (; invI < RANGE2C_SZ+2; ++invI) {
    m_range2c[invI] = maxC;
  }
  m_maxC = maxC;
}

// static uint64_t to_u(double x){
  // uint64_t u;
  // memcpy(&u, &x, sizeof(u));
  // return u;
// }

int arithmetic_decode_model_t::decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen)
{
  const uint64_t MSB_MSK = uint64_t(-1) << 40; // 16 overhead bits + 8 most significant bits
  const uint64_t VAL_MSK = uint64_t(-1) >> 16;
  const double MIN_RANGE = 1.0 /(int32_t(1) << 26); // 2**(-26)

  if (srclen < 2)
    return -11;

  uint8_t tmpbuf[16] = {0};
  bool useTmpbuf = false;
  if (srclen < 8) {
    memcpy(tmpbuf, src, srclen);
    src = tmpbuf;
    useTmpbuf = true;
  }

  double value = 0;
  for (int k = 0; k < 6; ++k)
    value = (value+src[5-k]) * (1.0/256);
  src    += 6;
  srclen -= 6;
  value  += 16.0;

  double lo    = 16.0;    // offset by 16 to push significant bits into 48 LS bits
  double range = 1.0;
  float  invRange = 1.0f; // approximation of 1/range from below

  for (int i = 0; ; ) {
    if ((i & 63)==0) {
      invRange = float(1.0/range);
    }

    if (value - lo >= range) return -12;
    // That is an input error, rather than internal error of decoder.
    // Due to way that encoder works, not any bit stream is possible as its output
    // That's the case of illegal code stream.
    // The case is extremely unlikely, but not impossible.

    unsigned c = val2c(value, lo, range, invRange);
    dst[i] = c;
    ++i;
    if (i == dstlen)
      break; // success
    // if (i > 5880 && i < 5890) {
      // unsigned ri = float(dValue)*invRange*512.0f;
      // printf("%d:%3d %8.1f %8.1f %8.1f    %18.11f %18.11f   %20.18f  %20.18f %d %d\n"//"=%016llx %016llx\n"
     // , i
     // , c, m_c2low[c+0]*VAL_RANGE, m_c2low[c+1]*VAL_RANGE, c < 254 ? m_c2low[c+2]*VAL_RANGE : m_c2low[c+1]*VAL_RANGE
     // , dValue/range*VAL_RANGE, float(dValue)*invRange*VAL_RANGE
     // ,m_c2low[c+1]*range+lo, value
     // ,m_range2c[ri+0]
     // ,m_range2c[ri+1]
     // //, to_u(value), to_u(range+16)
     // );
    // }

    // keep decoder in sync with encoder
    float cLo = m_c2low[c+0];
    float cRa = m_c2low[c+1] - cLo;
    lo   += range * cLo;
    range = (16.0 + range * cRa) - 16.0;
    // update invRange
    invRange *= m_c2invRange[c];

    if (value-lo >= range) return -102; // should never happen

    // printf("was here A2. %016llx %016llx %016llx\n",  lo, value, hi); fflush(stdout);

    if (range <= MIN_RANGE) {
      #ifdef ENABLE_PERF_COUNT
      ++m_renormalizationCount;
      #endif
      if (srclen < 8) {
        if (!useTmpbuf && srclen > 0) {
          memcpy(tmpbuf, src, srclen);
          src = tmpbuf;
          useTmpbuf = true;
        }
      }
      double hi = lo + range;
      uint64_t uLo; memcpy(&uLo, &lo, sizeof(uint64_t));
      uint64_t uHi; memcpy(&uHi, &hi, sizeof(uint64_t));
      uHi -= 1;
      uint64_t dbits = uLo ^ uHi;
      dbits |= (1<<15); // we don't want to renormalize by more than 32 bits
      unsigned bitsShifted = (__builtin_clzll(dbits)-16) & (-8); // # of common MS bits in lo and hi
      uint64_t uVa; memcpy(&uVa, &value, sizeof(uint64_t));
      if (bitsShifted != 0) {
        uLo = (uLo & ~VAL_MSK) | ((uLo << bitsShifted) & VAL_MSK); // shift significant bits, keep upper bits intact
        uVa = (uVa & ~VAL_MSK) | ((uVa << bitsShifted) & VAL_MSK); // shift significant bits, keep upper bits intact
        uint32_t fourOctets =
          (uint32_t(src[3])       ) |
          (uint32_t(src[2]) << 1*8) |
          (uint32_t(src[1]) << 2*8) |
          (uint32_t(src[0]) << 3*8);
        uVa    |= fourOctets >> (32-bitsShifted);
        unsigned octetsShifted = bitsShifted/8;                    // # of common MS octets in lo and hi
        src    += octetsShifted;
        srclen -= octetsShifted;
      }
      if (__builtin_expect(bitsShifted < 16, 0)) {
        do {
          // squeeze out bits[56..49]
          uLo = (uLo & MSB_MSK) | ((uLo << 8) & (~MSB_MSK));
          uVa = (uVa & MSB_MSK) | ((uVa << 8) & (~MSB_MSK));
          uVa |= *src++;
          srclen--;
          bitsShifted += 8;
        } while (bitsShifted < 16);
      }
      memcpy(&lo,    &uLo, sizeof(uint64_t));
      memcpy(&value, &uVa, sizeof(uint64_t));
      unsigned octetsShifted = bitsShifted/8;                    // # of common MS octets in lo and hi
      const double pow256_tab[5]    = { 1., 256., 256.*256., 256.*256*256., 256.*256.*256*256};
      const float  invPow256_tab[5] = { 1., 1.0f/256., 1.0f/(256.*256.), 1.0f/(256.*256*256.), 1.0f/(256.*256.*256*256)};
      range    *= pow256_tab[octetsShifted];
      invRange *= invPow256_tab[octetsShifted];
      if (srclen < -7)
        return i;
      if (value < lo || value - lo >= range) {
        // printf(
         // "%016llx < %016llx [<=  %016llx] || %016llx > %016llx\n"
         // "%016llx %016llx %016llx %d\n"
         // , lo, value, lo+((nxtRange<<RANGE_BITS)-1)
         // , value - lo, (nxtRange<<RANGE_BITS)-1
         // , hi ^ dbits, hi, dbits, octetsShifted
         // );
        return -103; // should not happen
      }
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

