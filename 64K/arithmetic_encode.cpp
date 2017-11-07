#include <algorithm>
//#include <cstdio>
#include <cmath>

#include "arithmetic_encode.h"


static const unsigned VAL_RANGE = 1u << 16;

// return value:
// -1  - source consists of repetition of the same character
// >=0 - maxC = the character with the biggest numeric value that appears in the source at least once
static int prepare(const uint8_t* src, unsigned srclen, uint16_t c2low[257], double* pQuantizedEntropy, double* pInfo)
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
    pInfo[3] = 0;
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

  if (statAndC[254].cnt==0)
    return -1; // source consists of repetition of the same character

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
    // here range < VAL_RANGE, because we already handled the case of repetition of the same character
    c2low[statAndC[i].c] = range;
    remRange -= range;
    srclen   -= cnt;
  }

  // calculate entropy after quanization
  double entropy = 0;
  for (unsigned c = 0; c < 256; ++c) {
    unsigned cnt = stat[c];
    if (cnt)
      entropy += log2(double(VAL_RANGE)/c2low[c])*cnt;
  }
  *pQuantizedEntropy = entropy;
  if (pInfo)
    pInfo[3] = entropy;

  // c2low -> cumulative sums of ranges
  unsigned lo = 0;
  unsigned maxC = 255;
  for (unsigned c = 0; c < 256; ++c) {
    unsigned range = c2low[c];
    c2low[c] = lo;
    lo += range;
    if (range > 0)
      maxC = c;
  }
  return maxC;
}

// return the number of stored octets
static int store_model(uint8_t* dst, const uint16_t c2low[256], unsigned maxC, double* pInfo)
{
  uint16_t c2range[256];
  uint32_t ranges[256];
  // translate cummulative sums o ranges to individual ranges
  int nRanges = 0;
  for (unsigned c = 0; c < maxC; ++c) {
    uint32_t range = c2low[c+1] - c2low[c];
    c2range[c] = range;
    if (range > 0)
      ranges[nRanges++] = (range << 8) | c;
  }
  uint32_t lastRange = VAL_RANGE - c2low[maxC];
  c2range[maxC] = lastRange;
  ranges[nRanges++] = (lastRange << 8) | maxC;

  uint8_t* p = dst;
  *p++ = 256 - nRanges; // # of zero ranges
  if (nRanges < 8) {
    // first variant of encoder - for small number of ranges
    // Arithmetic encode could have been used here too, but the savings are too small
    for (int i = 0; i < nRanges; ++i) {
      uint32_t rangeAndC = ranges[i];
      *p++ = uint8_t(rangeAndC);
      *p++ = uint8_t(rangeAndC >> 8);
      *p++ = uint8_t(rangeAndC >> 16);
    }
  } else {
    // main variant of encoder
    std::sort(ranges, &ranges[nRanges]);
    struct enc_rec_t {
      uint32_t cnt;
      uint32_t valx;
    } enc_tab[5];

    enc_tab[0].cnt = maxC + 1 - nRanges;
    enc_tab[0].valx = 256;
    uint32_t prevCnt = 0;
    *p++ = enc_tab[0].cnt;
    for (int i = 1; i < 5; ++i) {
      uint32_t cnt = (i*nRanges)/4;
      uint32_t valx = ranges[cnt-1];
      enc_tab[i].cnt = cnt-prevCnt;
      enc_tab[i].valx = valx;
      valx >>= 8;
      *p++ = uint8_t(valx>>0);
      *p++ = uint8_t(valx>>8);
      prevCnt = cnt;
    }
    // for (int i = 0; i < 5; ++i)
      // printf(" %u:%u.%u", enc_tab[i].cnt, enc_tab[i].valx>>8, enc_tab[i].valx & 255);
    // printf("\n");
    // for (unsigned c = 0; c <= maxC; ++c)
      // printf("%3u: %04x\n", c, c2range[c]);

    const uint64_t MSK23_0  = (uint64_t(1) << 24)-(uint64_t(1) << 0);
    const uint64_t MSK31_0  = (uint64_t(1) << 32)-(uint64_t(1) << 0);
    const uint64_t MSK39_32 = (uint64_t(1) << 40)-(uint64_t(1) << 32);
    uint64_t lo    = 0;                  // 40 bits
    uint64_t range = uint64_t(-1) >> 24; // in fact, range-1
    int pending_bytes = 0;
    unsigned nc = maxC + 1;
    for (unsigned c = 0; c < nc; ++c) {
      // prevent range from becoming too small
      while (range < (uint64_t(1) << 20)) {
        // squeeze out all ones in bits[39..32]
        lo = (lo & MSK39_32) | ((lo & MSK23_0) << 8);
        range = (range << 8) | 255;
        ++pending_bytes;
      }

      uint32_t val = c2range[c];
      // printf("%3u lo=%010llx ra=%010llx <= %04x\n", c, lo, range, val);
      int ix = 0;
      uint32_t val0 = 0;
      uint32_t val1 = 1;
      uint32_t cnt  = 0;
      if (val != 0) {
        uint32_t valx = val * 256 + c;
        cnt = enc_tab[0].cnt;
        for (ix = 1; enc_tab[ix].valx < valx; ++ix) {
          cnt += enc_tab[ix].cnt;
        }
        val0 = enc_tab[ix-1].valx >> 8;
        val1 = (enc_tab[ix].valx >> 8) + 1;
      }
      uint32_t cnt1 = enc_tab[ix].cnt;
      enc_tab[ix].cnt -= 1;
      uint32_t tc = val  - val0;
      uint32_t td = val1 - val0;
      uint32_t l_num = cnt*td + cnt1*tc; // up to 2^24-1
      uint32_t r_num = cnt1;             // up to 255
      uint32_t den = (nc-c)*td;

      lo   += ((range+1) * l_num + den - 1)/den; // ceil
      range = ((range+1) * r_num)/den - 1;       // floor
      // printf("%3u*lo=%010llx ra=%010llx %u/%u %u:%u:%u %u+%u<=%u\n"
       // , c, lo, range, l_num, den
       // , val0, val, val1
       // , cnt, cnt1, nc-c);

      uint64_t hi = lo + range;
      uint64_t dbits = lo ^ hi;
      while ((dbits >> 32)==0) {
        // lo and hi have 8 common MS bits
        *p++ = uint8_t(lo >> 32);
        if (pending_bytes) {
          uint8_t pending_byte = 0-((lo>>31) & 1);
          do {
            *p++ = pending_byte;
            --pending_bytes;
          } while (pending_bytes);
        }
        lo    = (lo & MSK31_0) << 8;
        range = (range << 8) | 255;
        dbits <<= 8;
      }
    }
    // last octet(s)
    // We want a shift register at decoder to be in range [lo..hi]
    // regardless of alien characters that a possibly resides after
    // our code stream
    uint64_t lastOctet = ((lo + range) >> 32)-1;
    *p++ = uint8_t(lastOctet);
    while (lastOctet == (lo>>32)) {
      lastOctet = 255;
      *p++ = uint8_t(lastOctet);
      lo = (lo & MSK31_0) << 8;
    }
  }

  int len = p - dst;
  if (pInfo)
    pInfo[1] = len*8;
  return len;
}

static int encode(uint8_t* dst, const uint8_t* src, unsigned srclen, const uint16_t c2low[257], int maxC)
{
  uint64_t VAL_MSK  = uint64_t(-1) >> 16;
  uint64_t MSK31_0  = (uint64_t(1) << 32)-(uint64_t(1) << 0);
  uint64_t MSK47_40 = (uint64_t(1) << 48)-(uint64_t(1) << 40);
  uint64_t lo = 0;
  uint64_t hi = VAL_MSK;
  int pending_bytes = 0;
  uint8_t* dst0 = dst;
  for (unsigned i = 0; i < srclen; ++i) {
    uint64_t range;
    while ((range = hi - lo) < (1u << 28)) {
      // squeeze out bits[39..32]
      lo = (lo & MSK47_40) | ((lo & MSK31_0) << 8);
      hi = (hi & MSK47_40) | ((hi & MSK31_0) << 8);
      ++pending_bytes;
    }
    range += 1;

    int c = src[i];
    if (c < maxC)
      hi = lo + ((range * c2low[c+1])>>16) - 1;
    lo   = lo + ((range * c2low[c+0])>>16);

    while (((lo ^ hi) >> 40)==0) {
      // lo and hi have the same upper octet
      *dst++ = uint8_t(lo>>40);
      lo = (lo << 8) & VAL_MSK;
      hi = (hi << 8) & VAL_MSK;
      while (pending_bytes) {
        uint8_t pending_byte = 0-((lo>>47) & 1);
        *dst++ = pending_byte;
        --pending_bytes;
      }
    }
  }
  // put out last octet
  *dst++ = uint8_t(hi>>40);
  return dst - dst0;
}

// return value:
// -1 - source consists of repetition of the same character
//  0 - not compressible, because all input characters have approximately equal probability or because input is too short
// >0 - the length of compressed buffer
int arithmetic_encode(std::vector<uint8_t>* dst, const uint8_t* src, int srclen, double* pInfo)
{
  uint16_t c2low[257];
  double quantizedEntropy;
  int maxC = prepare(src, srclen, c2low, &quantizedEntropy, pInfo);

  if (maxC == -1)
    return -1; // source consists of repetition of the same character

  size_t sz0 = dst->size();
  dst->resize(sz0 + 640);
  unsigned modellen = store_model(&dst->at(sz0), c2low, maxC, pInfo);

  if ((quantizedEntropy+7)/8 + modellen >= srclen)
    return 0; // not compressible

  // printf("ml=%u\n", modellen);
  dst->resize(sz0 + modellen + int(quantizedEntropy/8)+64);
  int dstlen = encode(&dst->at(sz0+modellen), src, srclen, c2low, maxC);

  if (pInfo)
    pInfo[2] = dstlen*8.0;

  int reslen = modellen + dstlen;

  if (reslen >= srclen)
    return 0; // not compressible

  dst->resize(sz0 + reslen);
  return reslen;
}
