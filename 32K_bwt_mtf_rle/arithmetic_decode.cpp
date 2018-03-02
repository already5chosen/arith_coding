#include <cstdint>
#include <cstring>
#include <cstdio>
// #include <cmath>
// #include <cctype>

#include "arithmetic_decode.h"
#include "arithmetic_coder_ut.h"
#include "arithmetic_coder_cfg.h"

static const unsigned VAL_RANGE = 1u << RANGE_BITS;

class CArithmeticDecoder {
public:
  uint64_t get(uint64_t cScale)  {
    m_range /= cScale;
    return (m_val >> 1)/m_range;
  }
  void init(const uint8_t* src, int srclen);
  int put(uint64_t cLo, uint64_t cRange); // return 0 on success, negative number on error
  int extract_1bit(const unsigned range_tab[3], int32_t* pRes); // return 0 on success, negative number on error
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

// return 0 on success,
// return negative number on error
int CArithmeticDecoder::extract_1bit(const unsigned range_tab[3], int32_t* pRes)
{
  unsigned mid = range_tab[1];
  if (0==mid)
    return 1;
  if (VAL_RANGE==mid)
    return 0;
  unsigned val = get(VAL_RANGE);
  int res = (val >= mid);
  *pRes = res;
  unsigned lo = range_tab[res];
  unsigned ra = range_tab[res+1]-lo;
  return put(lo, ra);
}

// load_quantized_histogram
// return  value:
//   on success maxC >= 0
//   on parsing error negative error code
static
int load_quantized_histogram(uint8_t* qh, CArithmeticDecoder* pDec)
{
  // load hhw - combined hh word
  unsigned hhw = pDec->get(8*9*10*9);
  int err = pDec->put(hhw, 1);
  if (err)
    return err;

  // split hhw into individual hh words
  unsigned hh02 = (hhw % 8) + 1; hhw /= 8;
  unsigned hh01 = (hhw % 9) + 1; hhw /= 9;
  unsigned hh23 = (hhw %10) + 0; hhw /= 10;
  unsigned hh45 = hhw + 1;

  unsigned range_tab[4][3];
  range_tab[0][1] = quantized_histogram_pair_to_range_qh_scale9(hh02, VAL_RANGE);
  range_tab[1][1] = quantized_histogram_pair_to_range_qh_scale9(hh01, VAL_RANGE);
  range_tab[2][1] = quantized_histogram_pair_to_range_qh_scale9(hh23, VAL_RANGE);
  range_tab[3][1] = quantized_histogram_pair_to_range_qh_scale9(hh45, VAL_RANGE);
  for (int k = 0; k < 4; ++k) {
    range_tab[k][0] = 0;
    range_tab[k][2] = VAL_RANGE;
  }

  unsigned val = 0;
  int32_t prev_msb = 0;
  uint32_t rlAcc = 0;
  uint32_t rlMsb = 1;
  int down = 0;
  bool first = true;
  int maxC = -1;
  for (unsigned i = 0; ; )
  {
    int32_t msb = 0;
    if (!first) {
      err = pDec->extract_1bit(range_tab[0], &msb);
      if (err)
        return err;
    }
    first = false;

    if (msb != prev_msb) {
      // end of run
      int32_t rl = rlAcc + rlMsb - 2;
      rlMsb = 1;
      rlAcc = 0;
      if (prev_msb == 0) {
        // end of zero run
        i += rl;
        if (rl > 0)
          memset(&qh[257-i], val, rl);
      } else {
        rl += 1;
        // end of diff run
        if (down) {
          if (rl + 1 > val)
            return -15;
          val -= rl + 1;
        } else {
          val += rl + 1;
          if (val > 255)
            return -16;
        }
        qh[256-i] = val;
        if (maxC < 0)
          maxC = 256-i;
        i += 1;
      }
    }

    unsigned tabId = (msb==0) ? 0 : 2 - prev_msb;
    int32_t lsb;
    err = pDec->extract_1bit(range_tab[tabId+1], &lsb);
    if (err)
      return err;

    if (tabId != 2) {
      // zero run or difference
      rlAcc |= (-lsb) & rlMsb;
      rlMsb += rlMsb;
      uint32_t rl = rlAcc + rlMsb - 2;
      if (msb == 0) {
        if (rl + i >= 257) {
          if (rl + i == 257) {
            if (rl > 0)
              memset(&qh[0], val, rl);
            break; // done
          } else {
            return -17;
          }
        }
      } else {
        if (rl >= 255)
          return -18;
      }
    } else {
      down = lsb; // sign of difference
    }
    prev_msb = msb;
  }

  // for (int i = 0; i <= maxC; ++i)
    // printf("range[%3d]=%5d\n", i, qh[i]);

  return maxC; // success
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
  uint16_t m_c2low[258];
  uint32_t m_c2invRange[257]; // floor(2^31/range)

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
  uint8_t m_range2c[RANGE2C_SZ+1];
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
  uint8_t qh[257];
  int maxC = load_quantized_histogram(qh, pDec);
  if (maxC >= 0) {
    quantized_histogram_to_range(m_c2low, maxC+1, qh, VAL_RANGE);
    for (unsigned c = maxC+1; c < 257; ++c)
      m_c2low[c] = 0;
    prepare();
  }
  return maxC;
}

void arithmetic_decode_model_t::prepare()
{
  // m_c2low -> cumulative sums of ranges
  unsigned maxC = 0;
  unsigned lo = 0;
  unsigned invI = 0;
  for (int c = 0; c < 257; ++c) {
    unsigned range = m_c2low[c];
    // printf("%3u: %04x %04x\n", c, range, lo);
    m_c2low[c] = lo;
    if (range != 0) {
      m_c2invRange[c] = (1u << 31)/range; // floor(2^31/range)
      // build inverse index m_range2c
      maxC = c;
      lo += range;
      unsigned r2c = c <= 255 ? c : 255;
      for (; invI <= ((lo-1) >> (RANGE_BITS-RANGE2C_NBITS)); ++invI) {
        m_range2c[invI] = r2c;
      }
    }
  }
  m_c2low[257] = VAL_RANGE;
  unsigned r2c = maxC <= 255 ? maxC : 255;
  for (; invI <= RANGE2C_SZ; ++invI) {
    m_range2c[invI] = r2c;
  }
  m_maxC = maxC;
}

int decode(
  arithmetic_decode_model_t* pModel,
  uint8_t*                   dst0,
  int                        dstlen,
  int32_t                    histogram[256],
  CArithmeticDecoder*        pDec)
{
  // initialize move-to-front decoder table
  uint8_t mtf_t[256];
  for (int i = 0; i < 256; ++i)
    mtf_t[i] = i;

  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);
  const double INV_RANGE_SCALE = double(int64_t(1) << (54-RANGE_BITS)) * (int64_t(1) << 43); // 2**(97-rb)

  const uint8_t* src = pDec->m_src;
  int srclen = pDec->m_srclen;

  uint8_t tmpbuf[16] = {0};
  bool useTmpbuf = false;
  if (srclen <= 0) {
    src = tmpbuf;
    useTmpbuf = true;
  }

  uint64_t value = pDec->m_val;                       // scaled by 2**64
  uint64_t range = pDec->m_range >> (RANGE_BITS-1);   // scaled by 2**(64-rb).  Maintained in [2**(33-rb)+1..2*(64-rb)]
  uint64_t invRange = int64_t(INV_RANGE_SCALE/range); // approximation of floor(2**(97-rb)/range). Maintained in [2**33..2*64)

  // uint64_t mxProd = 0, mnProd = uint64_t(-1);
  uint32_t rlAcc = 0;
  uint32_t rlMsb = 1;
  uint8_t* dst = dst0;
  // int dbg_i = 0;
  for (int i = 0; ; ) {
    #if 0
    uint64_t prod = umulh(invRange, range);
    if (prod > mxProd || prod < mnProd) {
      const int64_t PROD_ONE = int64_t(1) << (33-RANGE_BITS);
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

    // {
    // uint64_t loEst33b = umulh(value,invRange);
    // unsigned ri = loEst33b>>(33-9);
    // unsigned lo = loEst33b>>(33-RANGE_BITS);
    // }
    int c = pModel->val2c_estimate(value, invRange); // can be off by -1, much less likely by -2
    // keep decoder in sync with encoder
    uint64_t cLo = pModel->m_c2low[c+0];
    uint64_t cHi = pModel->m_c2low[c+1];
    value -= range * cLo;
    uint64_t nxtRange = range * (cHi-cLo);
    // at this point range is scaled by 2**64 - the same scale as value
    while (__builtin_expect(value >= nxtRange, 0)) {
      #ifdef ENABLE_PERF_COUNT
      ++pModel->m_longLookupCount;
      #endif
      value -= nxtRange;
      cLo = cHi;
      cHi = pModel->m_c2low[c+2];
      nxtRange = range * (cHi-cLo);
      ++c;
    }
    range = nxtRange;

    // RLE and MTF decode
    // printf("[%3d]=%3d\n", dbg_i, c); ++dbg_i;
    if (c < 2) {
      // zero run
      rlAcc |= (-c) & rlMsb;
      rlMsb += rlMsb;
      uint32_t rl = rlAcc + rlMsb - 1;
      if (rl >= uint32_t(dstlen)) {
        if (rl == uint32_t(dstlen)) {
          // last run
          int c0 = mtf_t[0];
          // for (int ii=0; ii < rl; ++ii)
            // {printf("[%3d]=%3d\n", dbg_i, c0); ++dbg_i;}
          histogram[c0] += rl;
          memset(dst, c0, rl);
          dst += rl;
          break;
        }
        return -21; // zero run too long (A)
      }
    } else {
      if (rlMsb > 1) {
        // insert zero run
        uint32_t rl = rlAcc + rlMsb - 1;
        // for (int ii=0; ii < rl; ++ii)
          // {printf("[%3d]=%3d\n", dbg_i, 0); ++dbg_i;}
        rlMsb = 1;
        rlAcc = 0;
        if (rl >= uint32_t(dstlen))
          return -22; // zero run too long (B)
        int c0 = mtf_t[0];
        // for (int ii=0; ii < rl; ++ii)
          // {printf("[%3d]=%3d\n", dbg_i, c0); ++dbg_i;}
        histogram[c0] += rl;
        memset(dst, c0, rl);
        dst += rl;
        dstlen -= rl;
      }

      if (dstlen == 0)
        return -23; // decoded section is longer than expected
      dstlen -= 1;

      int mtfC = c - 1;
      // printf("[%3d]=%3d\n", dbg_i, mtfC); ++dbg_i;
      int dstC = mtf_t[mtfC];
      // printf("[%3d]=%3d\n", dbg_i, dstC); ++dbg_i;
      histogram[dstC] += 1;
      *dst++ = dstC;

      if (dstlen == 0)
        break; // success

      // update move-to-front decoder table
      memmove(&mtf_t[1], &mtf_t[0], mtfC);
      mtf_t[0] = dstC;
    }

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
    invRange = umulh(invRange, uint64_t(pModel->m_c2invRange[c]) << (63-31)) << (RANGE_BITS+1);

    if ((i & 15)==0) {
      // do one NR iteration to increase precision of invRange and assure that invRange*range <= 2**(97-rb)
      const uint64_t PROD_ONE = uint64_t(1) << (33-RANGE_BITS);
      uint64_t prod = umulh(invRange, range); // scaled to 2^(33-RANGE_BITS)
      invRange = umulh(invRange, (PROD_ONE*2-1-prod)<<(RANGE_BITS+30))<<1;
    }
  }
  return dst - dst0;
}

}

// Arithmetic decode followed by RLE and MTF decode
int arithmetic_decode(
  uint8_t*       dst,
  int            dstlen,
  int32_t        histogram[256],
  const uint8_t* src,
  int            srclen,
  int*           pInfo)
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
    memset(histogram, 0, sizeof(histogram[0])*256);
    int textlen = decode(&model, dst, dstlen, histogram, &dec);
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
