#include <cstdint>
#include <cstdio>
// #include <cmath>
// #include <cctype>

#include "arithmetic_decode.h"


static const unsigned VAL_RANGE = 1u << 16;


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
  int hist[9] = {0};
  int rem = maxC;
  for (int i = 0; i < 8 && rem > 0; ++i) {
    if (srcI == srclen) return -1;
    hist[i] = src[srcI++];
    rem -= hist[i];
  }
  if (rem < 0) return -2; // inconsistent histogram table
  if (rem > 0)
    hist[8] = rem;

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
      uint32_t val0   = 1u << lshift; // up to 2^14
      uint32_t valDen = 3u << lshift; // up to 2^16-2^14

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

  // for (unsigned c = 0; c <= maxC; ++c)
    // printf("%3u: %04x\n", c, ranges[c]);
  // printf("len=%d, sum=%04x\n", int(p - src), sum);
  if (sum >= VAL_RANGE)
    return -5; // parsing error

  ranges[maxC] = VAL_RANGE - sum;
  for (unsigned c = maxC+1; c < 256; ++c)
    ranges[c] = 0;

  if (pInfo)
    pInfo[0] = srcI*8;

  return srcI; // success
}

namespace {

struct arithmetic_decode_model_t {
  uint16_t m_c2low[257];

  int  load_and_prepare(const uint8_t* src, int srclen, int* pInfo);
  int  decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen);
private:
  uint8_t  m_range2c[513];
  unsigned m_maxC;

  void prepare();
  int val2c(uint64_t value, uint64_t range) {
    unsigned ri = (value*512)/range;
    unsigned c = m_range2c[ri]; // c is the biggest character for which m_c2low[c] <= (val/128)*128
    if (c != m_range2c[ri+1]) {
      while (((m_c2low[c+1]*range) >> 16) <= value && c < m_maxC)
        ++c;
    }
    return c;
  }

  // const uint16_t* c2val(int c) const
  // {
    // return &m_i2low[m_c2i[c]];
  // }
  // uint16_t m_val2c[512];
  // const uint16_t* findC(int val)
  // {
    // int c = m_val2c[val]; // m_c2low[c] <= val
    // while (c < 255) {
      // if (m_c2low[c+1] > val)
        // break;
      // ++c;
    // }
    // return &m_c2low[c];
  // }
};


int arithmetic_decode_model_t::load_and_prepare(const uint8_t* src, int srclen, int* pInfo)
{
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
    m_c2low[c] = lo;
    // printf("%3d : %04x\n", c, lo);
    if (range != 0) {
      // build inverse index m_range2c
      maxC = c;
      lo += range;
      for (; invI <= ((lo-1) >> 7); ++invI) {
        m_range2c[invI] = maxC;
        // printf("%04x => %3d\n", invI << 7, maxC);
      }
    }
  }
  for (; invI < 513; ++invI) {
    m_range2c[invI] = maxC;
    // printf("%04x => %3d\n", invI << 7, maxC);
  }
  m_maxC = maxC;
}

int arithmetic_decode_model_t::decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen)
{
  uint64_t VAL_MSK  = uint64_t(-1) >> 16;
  uint64_t MSK31_0  = (uint64_t(1) << 32)-(uint64_t(1) << 0);
  uint64_t MSK47_40 = (uint64_t(1) << 48)-(uint64_t(1) << 40);
  uint64_t lo = 0;
  uint64_t hi = VAL_MSK;
  if (srclen < 2)
    return -11;
  uint64_t value = 0;
  for (int k = 0; k < 6; ++k) {
    value += uint64_t(*src) << ((5-k)*8);
    src++;
    srclen--;
    if (srclen == 0)
      break;
  }

  for (int i = 0; ; ) {
    uint64_t range;
    // keep decoder in sync with encoder (a)
    while ((range = hi - lo) < (1u << 28)) {
      // squeeze out bits[39..32]
      lo = (lo & MSK47_40) | ((lo & MSK31_0) << 8);
      hi = (hi & MSK47_40) | ((hi & MSK31_0) << 8);
      value = (value & MSK47_40) | ((value & MSK31_0) << 8);
      if (srclen > 0)
        value += *src++;
      --srclen;
      if (srclen < -5)
        return i;
      if (value < lo || value > hi) return -101; // should not happen
    }
    range += 1;

    unsigned c = val2c(value - lo, range);
    // if (i % 1000 < 10 || i > dstlen - 20000) {
      // printf(": %6d [%012llx %012llx] %012llx %012llx. %08x c=%3d '%c' [%04x..%04x) => [%012llx %012llx]\n"
       // , i, lo, hi, value, range
       // , unsigned(round(double(value - lo)*0x100000000/range))
       // , c, isprint(c) ? c : '.', m_c2low[c+0], m_c2low[c+1]
       // , lo + ((range * m_c2low[c+0])>>16)
       // , (c < m_maxC) ? lo + ((range * m_c2low[c+1])>>16) - 1 : hi);
      // fflush(stdout);
    // }
    dst[i] = c;
    ++i;
    if (i == dstlen)
      break; // success

    // keep decoder in sync with encoder (b)
    if (c < m_maxC)
      hi = lo + ((range * m_c2low[c+1])>>16) - 1;
    lo   = lo + ((range * m_c2low[c+0])>>16);

    if (value < lo || value > hi) return -102; // should not happen

    while (((lo ^ hi) >> 40)==0) {
      // lo and hi have the same upper octet
      if (((lo ^ value) >> 40)!=0)
        return -12;
      lo = (lo << 8) & VAL_MSK;
      hi = (hi << 8) & VAL_MSK;
      value = (value << 8) & VAL_MSK;
      if (srclen > 0)
        value += *src++;
      --srclen;
      if (srclen < -5)
        return i;
      if (value < lo || value > hi) return -103; // should not happen
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
    return textlen;
  }
  return modellen;
}

