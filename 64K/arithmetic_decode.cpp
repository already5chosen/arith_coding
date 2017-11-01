#include <cstdint>
#include <cstdio>

#include "arithmetic_decode.h"


static const unsigned VAL_RANGE = 1u << 16;


// load_ranges
// return  value:
//   on succes the # of processed source octets,
//   on parsing error negative error code
//static
int load_ranges(uint16_t* ranges, const uint8_t* src, unsigned srclen, int* pInfo)
{
  unsigned maxC = src[0];
  unsigned c = 0;
  unsigned i = 8;
  uint32_t sum = 0;
  while (c <= maxC) {
    unsigned di0 = (i+0) / 8;
    unsigned di1 = (i+7) / 8;
    unsigned ri = i % 8;
    if (di1 >= srclen)
      return -1; // parsing error

    unsigned b0 = src[di0];
    unsigned bx = src[di1];
    uint32_t val = (bx * 256u + b0) >> ri;
    uint32_t range  = 0;
    unsigned runlen = 1;
    // printf(">> %02x\n", val);
    if ((val & 1)==0) {
      // '0'   - 7-bit range or single zero, 8 bits toatal
      range = (val >> 1) & 127;
      i += 8;
    } else if ((val & 3)==1) {
      // '01'  - 9-bit range, 11 bits toatal
      di1 = (i+10) / 8;
      if (di1 >= srclen)
        return -2; // parsing error

      unsigned b1 = src[di0+1];
      bx = src[di1];
      val = ((bx*256u+b1)*256u+b0) >> ri;
      range = ((val >> 2) & 511) + 128;
      i += 11;
    } else if ((val & 7)==3) {
      // '011' - 16-bit range, 19 bits toatal
      di1 = (i+18) / 8;
      if (di1 >= srclen)
        return -3; // parsing error

      unsigned b1 = src[di0+1];
      unsigned b2 = src[di0+2];
      bx = src[di1];
      val = (((bx*256u+b2)*256u+b1)*256u+b0) >> ri;
      range = ((val >> 3) & 0xFFFF) + 128;
      i += 19;
    } else {
      // '111' - 5-bit zero run (length 2 to 33), 8 bits toatal
      runlen = ((val >> 3) & 31) + 2;
      i += 8;
    }
    sum += range;
    if (runlen + c > maxC+1) {
      // printf("%u: %u + %u > %u\n", i, runlen, c, maxC);
      return -4; // parsing error
    }
    do {
      // printf("%3u %04x\n", c, range);
      ranges[c] = range;
      ++c;
      --runlen;
    } while (runlen > 0);
  }

  // printf("sum = %u (%05x) %.3f %u\n", sum, sum, i / 8.0, maxC);
  if (sum < VAL_RANGE-1 || sum > VAL_RANGE)
    return -5; // parsing error

  for (; c < 256; ++c)
    ranges[c] = 0;

  if (pInfo)
    pInfo[0] = i;

  return (i + 7) / 8; // success
}

namespace {

struct arithmetic_decode_model_t {
  uint16_t m_c2low[257];

  int  load_and_prepare(const uint8_t* src, unsigned srclen, int* pInfo);
  int  decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen);
private:
  uint8_t  m_range2c[512];
  void prepare();
  int val2c(unsigned val) {
    unsigned c = m_range2c[val>>7]; // c is the biggest character for which m_c2low[c] <= (val/128)*128
    for (;m_c2low[c+1] < val && m_c2low[c+1] != 0;++c) ;
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


int arithmetic_decode_model_t::load_and_prepare(const uint8_t* src, unsigned srclen, int* pInfo)
{
  int ret = load_ranges(m_c2low, src, srclen, pInfo);
  if (ret >= 0)
    prepare();
  return ret;
}

void arithmetic_decode_model_t::prepare()
{
  // cummulative sum
  unsigned lo = 0;
  for (int c = 0; c < 256; ++c) {
    unsigned range = m_c2low[c];
    m_c2low[c] = lo;
    unsigned hi = lo + range;
    // build inverse index m_range2c
    for (unsigned i = (lo >> 7); i < (hi >> 7); ++i)
      m_range2c[i] = c;
    lo = hi;
  }
  m_c2low[256] = lo;
  for (int i = lo >> 7; i < 512; ++i)
    m_range2c[i] = 255;
}

int arithmetic_decode_model_t::decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen)
{
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

