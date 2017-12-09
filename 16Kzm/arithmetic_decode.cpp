#include <cstdint>
#include <cstring>
// #include <cstdio>
// #include <cmath>
// #include <cctype>

#include "arithmetic_decode.h"


static const int RANGE_BITS = 14;
static const unsigned VAL_RANGE = 1u << RANGE_BITS;

class CArithmeticDecoder {
public:
  uint64_t get(uint64_t cScale)  {
    m_range /= cScale;
    return (m_val >> 1)/m_range;
  }
  void init(const uint8_t* src, int srclen);
  int put(uint64_t cLo, uint64_t cRange); // return 0 on success, negative number on error
  uint64_t       m_val;   // scaled by 2**64
  uint64_t       m_range; // scaled by 2**63
  const uint8_t* m_src;
  int            m_srclen;
};

void CArithmeticDecoder::init(const uint8_t* src, int srclen)
{
  uint64_t value = 0; // scaled by 2**64
  for (int k = 0; k < 8; ++k) {
    value <<= 8;
    if (srclen > 0)
      value |= *src++;
    --srclen;
  }
  m_val    = value;
  m_range  = uint64_t(1) << 63;
  m_src    = src;
  m_srclen = srclen;
}

// return 0 on success,
// return negative number on error
int CArithmeticDecoder::put(uint64_t cLo, uint64_t cRange)
{
  const uint64_t MIN_RANGE = uint64_t(1) << (33-1);
  uint64_t value = m_val;      // scaled by 2**64
  uint64_t range = m_range;    // scaled by 2**63

  // keep decoder in sync with encoder
  value -= range * cLo * 2;
  range *= cRange;

  if (range <= MIN_RANGE) {
    uint32_t threeOctets = 0;
    for (int k = 0; k < 3; ++k) {
      threeOctets <<= 8;
      if (m_srclen > 0)
        threeOctets |= *m_src++;
      --m_srclen;
    }
    value = (value << 24) + threeOctets;
    range <<= 24;
    if (m_srclen < -7)
      return -14;
  }

  if ((value>>1) >= range)
    return -13;
  // That is an input error, rather than internal error of decoder.
  // Due to way that encoder works, not any bit stream is possible as its output
  // That's the case of illegal code stream.
  // The case is extremely unlikely, but not impossible.

  m_val   = value;
  m_range = range;
  return 0;
}

// load_ranges
// return  value:
//   on success 0
//   on parsing error negative error code
static
int load_ranges(uint16_t* ranges, CArithmeticDecoder* pDec)
{
  // load histogram of log2
  unsigned hist[RANGE_BITS+1] = {0};
  hist[0] = pDec->get(255); // hist[0] in range [0..254]
  int err = pDec->put(hist[0], 1);
  if (err)
    return err;
  unsigned rem = 256 - hist[0];
  for (int i = 1; i < RANGE_BITS && rem > 0; ++i) {
    hist[i] = pDec->get(rem+1); // hist[i] in range [0..rem]
    err = pDec->put(hist[i], 1);
    if (err)
      return err;
    rem -= hist[i];
  }
  if (rem > 0)
    hist[RANGE_BITS] = rem;

  // load c2range
  unsigned sum = 0;
  int nRanges = 0; // count non-zero ranges
  unsigned maxC = 255;
  for (unsigned c = 0; sum < VAL_RANGE && c < 255; ++c) {
    // extract exp.range
    unsigned ix = pDec->get(256-c);
    int log2_i;
    unsigned lo = 0;
    for (log2_i = 0; lo+hist[log2_i] <= ix; ++log2_i)
      lo += hist[log2_i];
    err = pDec->put(lo, hist[log2_i]);
    if (err)
      return err;
    hist[log2_i] -= 1; // update histogram

    unsigned range = 0;
    if (log2_i != 0) {
      ++nRanges;
      range = uint32_t(1) << (log2_i-1);
      if (log2_i > 1) {
        // extract offset within range
        unsigned offset = pDec->get(range);
        err = pDec->put(offset, 1);
        if (err)
          return err;
        range += offset;
      }
      sum += range;
    }
    ranges[c] = range;
    maxC = c;
  }

  if (sum < VAL_RANGE) {
    ranges[255] = VAL_RANGE - sum; // ranges[255] implied
    ++nRanges;
  } else if (sum == VAL_RANGE) {
    for (unsigned c = maxC+1; c < 256; ++c)
      ranges[c] = 0;
  } else { // sum > VAL_RANGE
    return -5; // parsing error
  }

  return 0; // success
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
  uint16_t m_c2low[257];
  uint32_t m_c2invRange[256]; // floor(2^31/range)

  int load_and_prepare(CArithmeticDecoder* pDec);
  int val2c_estimate(uint64_t value, uint64_t invRange) {
    uint64_t loEst33b = umulh(value,invRange);
    unsigned ri = loEst33b>>(33-RANGE2C_NBITS);
    unsigned lo = loEst33b>>(33-RANGE_BITS);
    unsigned c = m_range2c[ri]; // c is the biggest character for which m_c2low[c] <= (lo/32)*32
    if (__builtin_expect(m_c2low[c+1] <= lo, 0)) {
      do {
        #ifdef ENABLE_PERF_COUNT
        ++m_extLookupCount;
        #endif
        ++c;
      } while (m_c2low[c+1] <= lo);
    }
    return c;
    #if 0
    int ret;
    do {
      #ifdef ENABLE_PERF_COUNT
      ++m_extLookupCount;
      #endif
      ret = c;
      ++c;
    } while (m_c2low[c] <= lo);
    return ret;
    #endif
  }
private:
  static const int RANGE2C_NBITS = 9;
  static const int RANGE2C_SZ = 1 << RANGE2C_NBITS;
  uint8_t  m_range2c[RANGE2C_SZ+1];
  unsigned m_maxC;

  void prepare();
};


int arithmetic_decode_model_t::load_and_prepare(CArithmeticDecoder* pDec)
{
  #ifdef ENABLE_PERF_COUNT
  m_extLookupCount = 0;
  m_longLookupCount = 0;
  m_renormalizationCount = 0;
  #endif
  int ret = load_ranges(m_c2low, pDec);
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
      m_c2invRange[c] = (1u << 31)/range; // floor(2^31/range)
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

int decode(arithmetic_decode_model_t* pModel, uint8_t* dst, int dstlen, CArithmeticDecoder* pDec)
{
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);
  const double TWO_POW83  = double(int64_t(1) << 40) * (int64_t(1) << 43);

  const uint8_t* src = pDec->m_src;
  int srclen = pDec->m_srclen;

  uint8_t tmpbuf[16] = {0};
  bool useTmpbuf = false;
  if (srclen <= 0) {
    src = tmpbuf;
    useTmpbuf = true;
  }

  uint64_t value = pDec->m_val;                     // scaled by 2**64
  uint64_t range = pDec->m_range >> (RANGE_BITS-1); // scaled by 2**50.  Maintained in [2**19+1..2*50]
  uint64_t invRange = int64_t(TWO_POW83/range);     // approximation of floor(2**83/range). Maintained in [2**33..2*64)

  // uint64_t mxProd = 0, mnProd = uint64_t(-1);
  for (int i = 0; ; ) {
    #if 0
    uint64_t prod = umulh(invRange, range);
    if (prod > mxProd || prod < mnProd) {
      const int64_t PROD_ONE = int64_t(1) << 19;
      if (prod > mxProd) mxProd=prod;
      if (prod < mnProd) mnProd=prod;
      printf("%016llx*%016llx=%3lld [%3lld..%3lld] %d (%d)\n"
        ,(unsigned long long)range, (unsigned long long)invRange
        ,(long long)(prod-PROD_ONE), (long long)(mnProd-PROD_ONE), (long long)(mxProd-PROD_ONE)
        ,i, i & 15
        );
    }
    #endif

    if (value > (range << RANGE_BITS)-1) return -12;
    // That is an input error, rather than internal error of decoder.
    // Due to way that encoder works, not any bit stream is possible as its output
    // That's the case of illegal code stream.
    // The case is extremely unlikely, but not impossible.

    unsigned c = pModel->val2c_estimate(value, invRange); // can be off by -1, much less likely by -2
    // keep decoder in sync with encoder
    uint64_t cLo = pModel->m_c2low[c+0];
    uint64_t cHi = pModel->m_c2low[c+1];
    value -= range * cLo;
    uint64_t nxtRange = range * (cHi-cLo);
    // at this point range is scaled by 2**64 - the same scale as value
    while (__builtin_expect(value >= nxtRange, 0)) {
      #ifdef ENABLE_PERF_COUNT
      ++m_longLookupCount;
      #endif
      value -= nxtRange;
      cLo = cHi;
      cHi = pModel->m_c2low[c+2];
      nxtRange = range * (cHi-cLo);
      ++c;
    }
    range = nxtRange;
    dst[i] = c;
    ++i;
    if (i == dstlen)
      break; // success

    nxtRange = range >> RANGE_BITS;
    if (nxtRange <= MIN_RANGE) {
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

      uint32_t threeOctets =
        (uint32_t(src[2])       ) |
        (uint32_t(src[1]) << 1*8) |
        (uint32_t(src[0]) << 2*8);
      value = (value << 24) + threeOctets;
      src    += 3;
      srclen -= 3;
      nxtRange = range << (24 - RANGE_BITS);
      invRange >>= 24;
      if (srclen < -7)
        return i;
      if (value > ((nxtRange<<RANGE_BITS)-1)) {
        return -103; // should not happen
      }
    }
    range = nxtRange;

    // update invRange
    invRange = umulh(invRange, uint64_t(pModel->m_c2invRange[c]) << (63-31)) << (RANGE_BITS+1);

    if ((i & 15)==0) {
      // do one NR iteration to increase precision of invRange and assure that invRange*range <= 2**83
      const uint64_t PROD_ONE = uint64_t(1) << 19;
      uint64_t prod = umulh(invRange, range); // scaled to 2^20
      invRange = umulh(invRange, (PROD_ONE*2-1-prod)<<44)<<1;
    }
  }
  return dstlen;
}

}

int arithmetic_decode(uint8_t* dst, int dstlen, const uint8_t* src, int srclen, int* pInfo)
{
  if (pInfo) {
    pInfo[0] = 0;
    pInfo[1] = 0;
    pInfo[2] = 0;
    pInfo[3] = 0;
    pInfo[4] = 0;
  }
  if (srclen < 1)
    return -1;

  CArithmeticDecoder dec;
  dec.init(src, srclen);
  arithmetic_decode_model_t model;
  int err = model.load_and_prepare(&dec);
  if (err >= 0) {
    int modellen = srclen - dec.m_srclen;
    if (pInfo) {
      pInfo[0] = modellen*8;
      pInfo[1] = dec.m_srclen;
    }
    int textlen = decode(&model, dst, dstlen, &dec);
    #ifdef ENABLE_PERF_COUNT
    if (pInfo) {
      pInfo[2] = model.m_extLookupCount;
      pInfo[3] = model.m_longLookupCount;
      pInfo[4] = model.m_renormalizationCount;
    }
    #endif
    return textlen;
  }
  return err;
}
