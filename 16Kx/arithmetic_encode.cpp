//#include <cstdio>
#include <cmath>

#include "arithmetic_encode.h"


static const int RANGE_BITS = 14;
static const unsigned VAL_RANGE = 1u << RANGE_BITS;

// return value:
// -1  - source consists of repetition of the same character
// >=0 - maxC = the character with the biggest numeric value that appears in the source at least once
static int prepare(const uint8_t* src, unsigned srclen, uint16_t c2low[257], double* pQuantizedEntropy, double* pInfo)
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
      c2low[c] = range;
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
      c2low[c] = range;
    }
  }

  // calculate entropy after quantization
  double entropy = 0;
  for (unsigned c = 0; c <= maxC; ++c) {
    unsigned cnt = stat[c];
    if (cnt)
      entropy += log2(double(VAL_RANGE)/c2low[c])*cnt;
  }
  *pQuantizedEntropy = entropy;
  if (pInfo)
    pInfo[3] = entropy;

  // c2low -> cumulative sums of ranges
  unsigned lo = 0;
  for (unsigned c = 0; c <= maxC; ++c) {
    unsigned range = c2low[c];
    // printf("%3u: %04x\n", c, range);
    c2low[c] = lo;
    lo += range;
  }
  c2low[maxC+1] = VAL_RANGE;
  return maxC;
}

static int inline floor_log4(unsigned x) {
  return unsigned(sizeof(unsigned)*8 - 1 - __builtin_clz(x))/2; // floor(log4(x))
}

// return the number of stored octets
static int store_model(uint8_t* dst, const uint16_t c2low[256], unsigned maxC, double* pInfo)
{
  uint16_t c2range[256];
  // translate cumulative sums of ranges to individual ranges
  int nRanges = 0;
  unsigned hist[8] = {0};
  for (unsigned c = 0; c < maxC; ++c) {
    uint32_t range = c2low[c+1] - c2low[c];
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

static inline uint64_t umulh(uint64_t a, uint64_t b) {
  return ((unsigned __int128)a * b) >> 64;
}

static int encode(uint8_t* dst, const uint8_t* src, unsigned srclen, const uint16_t c2low[257], int maxC)
{
  uint64_t MSB_MSK  = uint64_t(255) << 55;
  uint64_t ACC_MSK  = uint64_t(-1) >> 1;
  uint64_t lo    = 0;                  // scaled by 2**63
  uint64_t range = uint64_t(-1) << 63; // scaled by 2**63
  int pending_bytes = 0;
  uint8_t* dst0 = dst;
  for (unsigned i = 0; i < srclen; ++i) {

    int c = src[i];
    uint64_t cLo = c2low[c+0];
    uint64_t cHi = c2low[c+1];
    uint64_t incLo = umulh(range, cLo << (63-RANGE_BITS));
    uint64_t incHi = umulh(range, cHi << (63-RANGE_BITS));
    lo   += incLo*2;
    range = (incHi-incLo)*2;

    if (range <= (1u << 30)) {
      uint64_t hi = lo + range -1;
      uint64_t dbits = lo ^ hi;
      while ((dbits & MSB_MSK)==0) {
        // lo and hi have the same upper octet
        *dst++ = uint8_t(lo>>55);
        while (pending_bytes) {
          uint8_t pending_byte = 0-((lo>>54) & 1);
          *dst++ = pending_byte;
          --pending_bytes;
        }
        lo    <<= 8;
        range <<= 8;
        dbits <<= 8;
      }
      while (range <= (1u << 30)) {
        // squeeze out bits[55..48]
        lo = (lo & MSB_MSK) | (lo << 8);
        range <<= 8;
        ++pending_bytes;
      }
      lo &= ACC_MSK;
    }
  }
  uint64_t hi = lo + range -1;
  uint64_t dbits = lo ^ hi;
  while ((dbits & MSB_MSK)==0) {
    // lo and hi have the same upper octet
    *dst++ = uint8_t(hi>>55);
    hi    <<= 8;
    dbits <<= 8;
  }
  // put out last octet
  *dst++ = uint8_t(hi>>55);
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
