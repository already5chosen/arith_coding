#include <algorithm>
#include <cstdio>
#include <cmath>

#include "arithmetic_encode.h"


static const unsigned VAL_RANGE = 1u << 16;

static void prepare(const uint8_t* src, unsigned srclen, uint16_t c2low[257], double* pInfo)
{
  // calculated statistics of appearance
  unsigned stat[256]={0};
  for (unsigned i = 0; i < srclen; ++i)
    ++stat[src[i]];

  if (pInfo) {
    // calculate source entropy
    double entropy = 0;
    for (unsigned c = 0; c < 256; ++c) {
      unsigned cnt = stat[c];
      if (cnt)
        entropy += log2(double(srclen)/cnt)*cnt;
    }
    pInfo[0] = entropy;
  }

  // sort statistics in ascending order
  struct stat_and_c_t {
    unsigned cnt, c;
    bool operator< (const stat_and_c_t& b) const {
      return cnt < b.cnt;
    }
  };
  stat_and_c_t statAndC[256];
  for (unsigned c = 0; c < 256; ++c) {
    statAndC[c].cnt = stat[c];
    statAndC[c].c   = c;
  }
  std::sort(&statAndC[0], &statAndC[256]);

  // translate counts to ranges and store in c2low
  unsigned i = 0;
  while (statAndC[i].cnt == 0) {
    c2low[statAndC[i].c] = 0;
    ++i;
  }
  unsigned remRange = VAL_RANGE;
  for ( ; srclen > 0; ++i) {
    unsigned cnt = statAndC[i].cnt;
    unsigned range = (uint64_t(cnt)*(remRange*2) + srclen)/(srclen*2);
    if (range == 0)
      range = 1;
    else if (range == VAL_RANGE)
      range = VAL_RANGE-1;
    c2low[statAndC[i].c] = range;
    remRange -= range;
    srclen   -= cnt;
  }

  // c2low -> cumulative sums of ranges
  unsigned lo = 0;
  for (unsigned c = 0; c < 256; ++c) {
    unsigned range = c2low[c];
    c2low[c] = lo;
    lo += range;
  }
  c2low[256] = 0;
}

static unsigned insert_bits(uint8_t* dst, unsigned i, uint32_t val, int nbits)
{
  // printf("<< %2d %05x\n", nbits, val);
  int di = i / 8;
  int ri = i % 8;
  if (ri != 0) {
    val = (val << ri) | dst[di]; // it's o.k, because nbits <= 25
    nbits += ri;
  }
  while (nbits > 0) {
    dst[di] = uint8_t(val);
    val  >>= 8;
    di    += 1;
    nbits -= 8;
  }
  return di*8 + nbits;
}

// model tags:
// '0'   - 7-bit range or single zero
// '01'  - 9-bit range
// '011' - 16-bit range
// '111' - 5-bit zero run (length 2 to 33)
static int insert_number(uint8_t* dst, unsigned i, unsigned val)
{
  if (val < 128)
    return insert_bits(dst, i, val*2 + 0, 8);
  val -= 128;
  if (val < 512)
    return insert_bits(dst, i, val*4 + 1, 11);
  else
    return insert_bits(dst, i, val*8 + 3, 19);
}
static int insert_zeros(uint8_t* dst, unsigned i, unsigned runlen)
{
  if (runlen == 1)
    return insert_bits(dst, i, 0*2 + 0, 8);
  else
    return insert_bits(dst, i, (runlen-2)*8 + 7, 8);
}

// return the number of stored octets
static int store_model(uint8_t* dst, const uint16_t c2low[257], double* pInfo)
{
  unsigned zero_run = 0;
  unsigned i = 0;
  for (unsigned c = 0; c < 256; ++c) {
    unsigned range = (uint32_t(c2low[c+1]) - uint32_t(c2low[c])) & (VAL_RANGE-1);
    // printf("%3u %04x\n", c, range);
    if (range == 0) {
      ++zero_run;
      if (zero_run == 33) {
        i = insert_zeros(dst, i, zero_run);
        zero_run = 0;
      }
    } else {
      if (zero_run) {
        i = insert_zeros(dst, i, zero_run);
        zero_run = 0;
      }
      i = insert_number(dst, i, range);
    }
  }
  if (zero_run)
    i = insert_zeros(dst, i, zero_run);
  
  if (pInfo)
    pInfo[1] = i;
  return (i + 7) / 8;
}

static void encode(std::vector<uint8_t>* dst, const uint8_t* src, unsigned srclen, const uint16_t c2low[257])
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

    int c = src[i];
    uint32_t cLo = c2low[c+0];
    uint32_t cHi = c2low[c+1];
    if (cHi==0)
      cHi = VAL_RANGE;
    hi = lo + ((range * cHi)>>16) - 1;
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


void arithmetic_encode(std::vector<uint8_t>* dst, const uint8_t* src, int srclen, double* pInfo)
{
  uint16_t c2low[257];
  prepare(src, srclen, c2low, pInfo);
  size_t sz0 = dst->size();
  dst->resize(sz0 + 640);
  unsigned modellen = store_model(&dst->at(sz0), c2low, pInfo);
  // printf("ml=%u\n", modellen);
  dst->resize(sz0 + modellen);
  encode(dst, src, srclen, c2low);
  if (pInfo)
    pInfo[2] = (dst->size()-sz0-modellen)*8.0;
}
