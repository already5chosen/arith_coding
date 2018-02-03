#include <cstdint>
#include <cstring>
#include <climits>
// #include <cstdio>
// #include <cmath>
// #include <cctype>

#include "arithmetic_decode.h"


static const int RANGE_BITS = 10;
static const unsigned VAL_RANGE = 1u << RANGE_BITS;
static const uint32_t SUBRANGE_SCALE = 1u << (31-RANGE_BITS);

typedef struct {
  unsigned nSmallRanges;
  uint16_t val[257];     // ranges/lo
  uint8_t  small2c[256]; // map index of small range to character

  unsigned GetFirstSmallI() const {  // index of the first small range
    return nSmallRanges != 0 ? small2c[0] : UINT_MAX;
  }
} c2range_t;

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
int load_ranges(c2range_t* ranges, CArithmeticDecoder* pDec)
{
  // load maxValC
  unsigned maxValC = pDec->get(256);
  int err = pDec->put(maxValC, 1);
  if (err)
    return err;

  // load histogram of log2
  unsigned hist[RANGE_BITS+1] = {0};
  unsigned rem = maxValC;
  for (int i = 0; i < RANGE_BITS && rem > 0; ++i) {
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
  for (unsigned c = 0; c < maxValC; ++c) {
    // extract exp.range
    unsigned ix = pDec->get(maxValC);
    int log2_i;
    unsigned lo = 0;
    for (log2_i = 0; lo+hist[log2_i] <= ix; ++log2_i)
      lo += hist[log2_i];
    err = pDec->put(lo, hist[log2_i]);
    if (err)
      return err;

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
    ranges->val[c] = range;
  }

  if (sum >= VAL_RANGE) {
    return -5; // parsing error
  }

  ranges->val[maxValC] = VAL_RANGE - sum; // ranges[maxValC] implied
  for (unsigned c = maxValC+1; c < 256; ++c)
    ranges->val[c] = 0;

  // load information about small ranges
  ranges->nSmallRanges = 0;
  unsigned nZeros = hist[0] + 255 - maxValC; // total number of zero ranges
  if (nZeros > 0) {
    // load nSmallZeroRanges
    unsigned nSmallZeroRanges = pDec->get(nZeros+1);
    int err = pDec->put(nSmallZeroRanges, 1);
    if (err)
      return err;
    if (nSmallZeroRanges > 0) {
      // small ranges present
      ranges->nSmallRanges = nSmallZeroRanges + 1;

      // load position of 'master' small range
      unsigned masterSmallRangeI = pDec->get(nRanges+1);
      int err = pDec->put(masterSmallRangeI, 1);
      if (err)
        return err;

      // find firstSmallI
      unsigned firstSmallI = 0;
      for (unsigned nnz = 0; ; ++firstSmallI) {
        nnz += (ranges->val[firstSmallI] != 0);
        if (nnz > masterSmallRangeI)
          break;
      }
      ranges->small2c[0] = firstSmallI;

      // calculate # of zeros between firstSmallI and maxValC
      unsigned nZerosA = 0;
      for (unsigned c = firstSmallI+1; c < maxValC; ++c)
        nZerosA += (ranges->val[c] == 0);

      unsigned aliasI = 1;
      if (nZerosA > 0) {
        // load # of small ranges between firstSmallI and maxValC
        unsigned nSmallZeroRangesA = pDec->get(nZerosA+1);
        int err = pDec->put(nSmallZeroRangesA, 1);
        if (err)
          return err;
        if (nSmallZeroRangesA > 0) {
          if (nSmallZeroRangesA < nZerosA) {
            // load alias flags
            unsigned rem = nSmallZeroRangesA;
            for (unsigned c = firstSmallI+1; rem > 0; ++c) {
              if (ranges->val[c] == 0) {
                unsigned cVal = pDec->get(nZerosA);
                unsigned cLo    = 0; // regular zero range
                unsigned cRange = nZerosA - nSmallZeroRangesA;
                if (cVal >= cRange) {
                  ranges->small2c[aliasI] = c;
                  ++aliasI;
                  cLo    = cRange; // alias range
                  cRange = nSmallZeroRangesA;
                  --rem;
                }
                int err = pDec->put(cLo, cRange);
                if (err)
                  return err;
              }
            }
          } else {
            // mark all zeros as aliases
            for (unsigned c = firstSmallI+1; c < maxValC; ++c) {
              if (ranges->val[c] == 0) {
                ranges->small2c[aliasI] = c;
                ++aliasI;
              }
            }
          }
          nSmallZeroRanges -= nSmallZeroRangesA;
        }
      }

      if (nSmallZeroRanges > 0) {
        // map small ranges that follow last non-zero range
        unsigned nZerosB = 255 - maxValC;
        if (nSmallZeroRanges < nZerosB) {
          // load alias flags
          unsigned rem = nSmallZeroRanges;
          for (unsigned c = maxValC+1; rem > 0; ++c) {
            unsigned cVal = pDec->get(nZerosB);
            unsigned cLo    = 0; // regular zero range
            unsigned cRange = nZerosB - nSmallZeroRanges;
            if (cVal >= cRange) {
              ranges->small2c[aliasI] = c;
              ++aliasI;
              cLo    = cRange; // alias range
              cRange = nSmallZeroRanges;
              --rem;
            }
            int err = pDec->put(cLo, cRange);
            if (err)
              return err;
          }
        } else {
          // mark all zeros as aliases
          for (unsigned c = maxValC+1; c < maxValC+nSmallZeroRanges; ++c) {
            ranges->small2c[aliasI] = c;
            ++aliasI;
          }
        }
      }
    }
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
  uint32_t  m_smallSubrangeScale;
  uint32_t  m_invSmallSubrange[2]; // indexed my (cRa % 2)
  c2range_t m_c2low;
  uint32_t  m_c2invRange[256]; // floor(2^31/range)

  int load_and_prepare(CArithmeticDecoder* pDec);
  int val2c(uint64_t value, uint64_t invRange, uint64_t range) {
    uint64_t loEst33b = umulh(value,invRange);
    unsigned lo = loEst33b>>(33-RANGE_BITS);
    lo += (value-lo*range >= range);
    return m_lo2c[lo]; // c is the biggest character for which m_c2low[c] <= lo
  }
  int val2smallI(uint64_t value, uint64_t invRange, uint64_t range) {
    uint64_t loEst33b = umulh(value,invRange);
    unsigned idx = (loEst33b*m_c2low.nSmallRanges) >> 33;
    unsigned lo  = (idx*m_smallSubrangeScale)/SUBRANGE_SCALE;
    idx -= (lo*range > value);
    return idx;
  }
private:
  uint8_t  m_lo2c[VAL_RANGE+1];
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
  int ret = load_ranges(&m_c2low, pDec);
  if (ret >= 0)
    prepare();
  return ret;
}

void arithmetic_decode_model_t::prepare()
{
  // m_c2low.val[] -> cumulative sums of ranges
  unsigned maxC = 255;
  unsigned lo = 0;
  unsigned invI = 0;
  for (int c = 0; c < 256; ++c) {
    unsigned range = m_c2low.val[c];
    m_c2low.val[c] = lo;
    if (range != 0) {
      m_c2invRange[c] = (1u << 31)/range; // floor(2^31/range)
      // build inverse index m_lo2c
      maxC = c;
      lo += range;
      for (; invI < lo; ++invI) {
        m_lo2c[invI] = maxC;
      }
    }
  }
  m_c2low.val[256] = VAL_RANGE;
  for (; invI <= VAL_RANGE; ++invI) {
    m_lo2c[invI] = maxC;
  }
  m_maxC = maxC;

  if (m_c2low.nSmallRanges) {
    m_smallSubrangeScale = (VAL_RANGE*SUBRANGE_SCALE)/m_c2low.nSmallRanges;
    unsigned cRa = VAL_RANGE/m_c2low.nSmallRanges;
    for (unsigned i = 0; i < 2; ++i)
      m_invSmallSubrange[(cRa+i) % 2] = (1u << 31)/(cRa+i); // floor(2^31/range)
  }
}

int decode(arithmetic_decode_model_t* pModel, uint8_t* dst, int dstlen, CArithmeticDecoder* pDec)
{
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);
  const double TWO_POW87   = double(int64_t(1) << 40) * (int64_t(1) << 47);

  const uint8_t* src = pDec->m_src;
  int srclen = pDec->m_srclen;

  uint8_t tmpbuf[16] = {0};
  bool useTmpbuf = false;
  if (srclen <= 0) {
    src = tmpbuf;
    useTmpbuf = true;
  }

  uint64_t value = pDec->m_val;                     // scaled by 2**64
  uint64_t range = pDec->m_range >> (RANGE_BITS-1); // scaled by 2**54.  Maintained in [2**23+1..2*54]
  uint64_t invRange = uint64_t(TWO_POW87/range);    // approximation of floor(2**87/range). Maintained in [2**33..2*64)
  unsigned firstSmallI = pModel->m_c2low.GetFirstSmallI();

  // uint64_t mxProd = 0, mnProd = uint64_t(-1);
  bool alias = false;
  for (int i = 0; ; ) {
    #if 0
    uint64_t prod = umulh(invRange, range);
    if (prod > mxProd || prod < mnProd) {
      const int64_t PROD_ONE = int64_t(1) << 23;
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

    unsigned c;
    uint64_t nxtRange;
    uint32_t cInvRange;

    if (__builtin_expect(!alias, 1)) {
      c = pModel->val2c(value, invRange, range); // can be off by -1
      // keep decoder in sync with encoder
      uint64_t cLo = pModel->m_c2low.val[c+0];
      uint64_t cHi = pModel->m_c2low.val[c+1];
      value -= range * cLo;
      nxtRange = range * (cHi-cLo);
      // at this point range is scaled by 2**64 - the same scale as value
      if (__builtin_expect(value >= nxtRange, 0)) {
        #ifdef ENABLE_PERF_COUNT
        ++pModel->m_longLookupCount;
        #endif
        value -= nxtRange;
        cLo = cHi;
        cHi = pModel->m_c2low.val[c+2];
        nxtRange = range * (cHi-cLo);
        ++c;
        if (__builtin_expect(value >= nxtRange, 0)) return -104;
      }
      // printf("%7d: %3u%c\n", i, c, alias ? '*' : '.');
      cInvRange = pModel->m_c2invRange[c];
      alias = (c == firstSmallI);
    } else {
      alias = false;
      unsigned idx = pModel->val2smallI(value, invRange, range); // can be off by -1
      // keep decoder in sync with encoder
      uint64_t cLo_x = pModel->m_smallSubrangeScale * idx;
      uint64_t cHi_x = cLo_x + pModel->m_smallSubrangeScale;
      uint64_t cLo = cLo_x / SUBRANGE_SCALE;
      uint64_t cHi = cHi_x / SUBRANGE_SCALE;
      if (idx+1==pModel->m_c2low.nSmallRanges)
        cHi = VAL_RANGE;
      value -= range * cLo;
      nxtRange = range * (cHi-cLo);
      // at this point range is scaled by 2**64 - the same scale as value
      if (__builtin_expect(value >= nxtRange, 0)) {
        #ifdef ENABLE_PERF_COUNT
        ++pModel->m_longLookupCount;
        #endif
        ++idx;
        value -= nxtRange;
        cLo = cHi;
        cHi = (cHi_x + pModel->m_smallSubrangeScale) / SUBRANGE_SCALE;
        if (idx+1==pModel->m_c2low.nSmallRanges)
          cHi = VAL_RANGE;
        nxtRange = range * (cHi-cLo);
        if (__builtin_expect(value >= nxtRange, 0)) return -105;
      }
      c = pModel->m_c2low.small2c[idx];
      // printf("%7d: %3u(%3d) %lld %lld !\n", i, c, idx, cLo, cHi-cLo);
      cInvRange = pModel->m_invSmallSubrange[(cHi-cLo)%2];
    }

    range = nxtRange;
    if (!alias) {
      dst[i] = c;
      ++i;
    }

    if (i == dstlen)
      break; // success

    nxtRange = range >> RANGE_BITS;
    if (nxtRange <= MIN_RANGE) {
      #ifdef ENABLE_PERF_COUNT
      ++pModel->m_renormalizationCount;
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
    invRange = umulh(invRange, uint64_t(cInvRange) << (63-31)) << (RANGE_BITS+1);

    if ((i & 15)==0) {
      // do one NR iteration to increase precision of invRange and assure that invRange*range <= 2**87
      const uint64_t PROD_ONE = uint64_t(1) << 23;
      uint64_t prod = umulh(invRange, range); // scaled to 2^23
      invRange = umulh(invRange, (PROD_ONE*2-1-prod)<<40)<<1;
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
