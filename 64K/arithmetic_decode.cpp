#include <cstdint>
//#include <cstring>
// #include <algorithm>
// #include <stdio.h>

//#include "arithmetic_encode.h"

static const unsigned VAL_RANGE = 1u << 16;


// load_ranges
// return  value:
//   on succes the # of processed source octets,
//   on parsing error -1
//static
int load_ranges(uint16_t* ranges, const uint8_t* src, unsigned srclen)
{
  unsigned c = 0;
  unsigned i = 0;
  uint32_t sum = 0;
  while (c < 256) {
    unsigned di0 = (i+0) / 8;
    unsigned di1 = (i+7) / 8;
    unsigned ri = i / 8;
    if (di1 >= srclen)
      return -1; // parsing error

    unsigned b0 = src[di0];
    unsigned bx = src[di1];
    uint32_t val = (bx * 256u + b0) >> ri;
    uint32_t range  = 0;
    unsigned runlen = 1;
    if ((val & 1)==0) {
      // '0'   - 7-bit range or single zero, 8 bits toatal
      range = (val >> 1) & 127;
      i += 8;
    } else if ((val & 3)==1) {
      // '01'  - 9-bit range, 11 bits toatal
      di1 = (i+10) / 8;
      if (di1 >= srclen)
        return -1; // parsing error

      unsigned b1 = src[di0+1];
      bx = src[di1];
      val = ((bx*256u+b1)*256u+b0) >> ri;
      range = ((val >> 2) & 511) + 128;
      i += 11;
    } else if ((val & 7)==3) {
      // '011' - 16-bit range, 19 bits toatal
      di1 = (i+18) / 8;
      if (di1 >= srclen)
        return -1; // parsing error

      unsigned b1 = src[di0+1];
      unsigned b2 = src[di0+2];
      bx = src[di1];
      val = (((bx*256u+b2)*256u+b1)*256u+b0) >> ri;
      range = ((val >> 3) & 0xFFFF) + 128;
      ++c;
      i += 19;
    } else {
      // '111' - 5-bit zero run (length 2 to 33), 8 bits toatal
      runlen = ((val >> 3) & 31) + 2;
      i += 8;
    }
    sum += range;
    if (runlen + c > 256)
      return -1; // parsing error
    do {
      ranges[c] = range;
      ++c;
      --runlen;
    } while (runlen > 0);
  }

  if (sum < VAL_RANGE-1 || sum > VAL_RANGE)
    return -1; // parsing error

  return (i + 7) / 8; // success
}

namespace {
struct arithmetic_decode_model_t {
  uint16_t m_c2low[257];

  int  load_and_prepare(const uint8_t* src, unsigned srclen);
  int  decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen);
private:
  uint8_t  m_range2c[256];
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


int arithmetic_decode_model_t::load_and_prepare(const uint8_t* src, unsigned srclen)
{
  int ret = load_ranges(m_c2low, src, srclen);
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
  return srclen;
}

#if 0

// int arithmetic_decode_model_t::load(uint16_t* ranges, const uint8_t* src, unsigned srclen)
// {
  // unsigned sum = 0;
  // unsigned c = 0;
  // unsigned i = 0;
  // while (c < 256) {
    // unsigned remlen = srclen - i;
    // if (remlen == 0)
      // return -1; // parsing error
    // unsigned b0 = src[i+0];
    // if (b0 >= 128) {
      // if (b0 == 255) {
        // // big range
        // if (remlen < 3)
          // return -1; // parsing error
        // ranges[c] = src[i+1]*256+src[i+2];
        // c += 1;
        // i += 3;
      // } else {
        // // zero run
        // unsigned runlen = b - 127;
        // if (c + runlen > 256)
          // return -1; // parsing error
        // while (runlen > 0) {
          // ranges[c] = 0;
          // ++c;
          // --runlen;
        // }
        // i += 1;
      // }
    // } else {
      // // normal range
      // if (remlen < 2)
        // return -1; // parsing error
      // ranges[c] = b0*256+src[i+1];
      // c += 1;
      // i += 2;
    // }
  // }
  // return i;
// }

void arithmetic_decode_model_t::prepare(uint16_t* ranges)
{
  // calculated statistics of appearance
  unsigned stat[256]={0};
  for (unsigned i = 0; i < srclen; ++i)
    ++stat[src[i]];

  struct stat_and_c_t {
    unsigned cnt, c;
    bool operator< (const stat_and_c_t& b) const {
      return cnt < b.cnt;
    }
  };

  // sort statistics in ascending order
  stat_and_c_t statAndC[256];
  for (unsigned i = 0; i < 256; ++i) {
    statAndC[i].cnt = stat[i];
    statAndC[i].c   = i;
  }
  std::sort(&statAndC[0], &statAndC[256]);

  // generate m_i2low and m_c2i
  unsigned i0 = 0;
  while (statAndC[i0].cnt == 0) {
    m_c2i[statAndC[i0].c] = 255;
    ++i0;
  }
  m_i2low[255] = m_i2low[256] = 0;

  unsigned val = 0;
  for (unsigned i = i0; i < 256; ++i) {
    m_c2i[statAndC[i].c] = i-i0;
    m_i2low[i-i0] = val;
    unsigned cnt = statAndC[i].cnt;
    unsigned range = (uint64_t(cnt)*((VAL_RANGE-val)*2) + srclen)/(srclen*2);
    if (range == 0)
      range = 1;
    else if (range == VAL_RANGE)
      range = VAL_RANGE-1;
    val    += range;
    srclen -= cnt;
  }
  m_i2low[256-i0] = val;

  // generate inverse index m_val2c
}

void arithmetic_encode_model_t::store_model(std::vector<uint8_t>* dst) const
{
  int zero_run = 0;
  for (unsigned c = 0; c < 256; ++c) {
    const uint16_t* pVal = c2val(c);
    unsigned range = (uint32_t(pVal[1]) - uint32_t(pVal[0])) & (VAL_RANGE-1);
    // printf("%3u %5u %5u %u\n", c, pVal[0], pVal[1], range);
    if (range == 0) {
      ++zero_run;
      if (zero_run == 127) {
        dst->push_back(127+127);
        zero_run = 0;
      }
    } else {
      if (zero_run) {
        dst->push_back(127+zero_run);
        zero_run = 0;
      }
      unsigned b0 = range >> 8;
      if (b0 >= 128) {
        // it happens at most twice
        dst->push_back(255);
      }
      dst->push_back(uint8_t(b0));
      dst->push_back(uint8_t(range));
    }
  }
  if (zero_run)
    dst->push_back(127+zero_run);
}

void arithmetic_encode_model_t::encode(std::vector<uint8_t>* dst, const uint8_t* src, unsigned srclen) const
{
  uint64_t VAL_MSK  = uint64_t(-1) >> 16;
  uint64_t MSK31_0  = (uint64_t(1) << 32)-(uint64_t(1) << 0);
  uint64_t MSK47_40 = (uint64_t(1) << 48)-(uint64_t(1) << 40);
  uint64_t lo = 0;
  uint64_t hi = VAL_MSK;
  int pending_bytes = 0;
  for (unsigned i = 0; i < srclen; ++i) {
    uint64_t range;
    while ((range = hi - lo) < (1u << 28)) {
      // sqweeze out bits[39..32]
      lo = (lo & MSK47_40) | ((lo & MSK31_0) << 8);
      hi = (hi & MSK47_40) | ((hi & MSK31_0) << 8);
    }
    range += 1;

    const uint16_t* pVal = c2val(src[i]);
    uint32_t cLo = pVal[0];
    uint32_t cHi = pVal[1];
    if (cHi==0)
      cHi = (1u << 16);
    hi = lo + ((range * cHi)>>16);
    lo = lo + ((range * cLo)>>16);

    while (((lo ^ hi) >> 40)==0) {
      // lo and hi have the same upper octet
      dst->push_back(uint8_t(lo>>40));
      lo = (lo << 8) & VAL_MSK;
      hi = (hi << 8) & VAL_MSK;
      while (pending_bytes) {
        uint8_t pending_byte = 0-((lo>>47) & 1);
        dst->push_back(pending_byte);
        --pending_bytes;
      }
    }
  }
  // put out last octets, at least 2
  lo = (lo+hi)>>32;
  dst->push_back(uint8_t(lo>>8));
  while (pending_bytes) {
    uint8_t pending_byte = 0-((lo>>7) & 1);
    dst->push_back(pending_byte);
    --pending_bytes;
  }
  dst->push_back(uint8_t(lo));
}

}

#endif

}

bool arithmetic_decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen)
{
  arithmetic_decode_model_t model;
  int modellen = model.load_and_prepare(&src[0], srclen);
  if (modellen < 0)
    return false;
  int codelen = model.decode(dst, dstlen, &src[modellen], srclen-modellen);
  return codelen==srclen-modellen;
}

