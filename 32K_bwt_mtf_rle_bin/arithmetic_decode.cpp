#include <cstdint>
#include <cstring>
#include <cstdio>
#include <climits>
#include <vector>
#include <cmath>
// #include <cctype>

#include "arithmetic_decode.h"
#include "arithmetic_coder_ut.h"
#include "arithmetic_coder_cfg.h"

static const unsigned VAL_RANGE = 1u << RANGE_BITS;

static uint16_t qh2h_tab[ARITH_CODER_QH_SCALE+1];

void arithmetic_decode_init_tables()
{
  qh2h_tab[ARITH_CODER_QH_SCALE] = VAL_RANGE;
  double M_PI = 3.1415926535897932384626433832795;
  for (int i = 1; i < ARITH_CODER_QH_SCALE; ++i)
    qh2h_tab[i] = int(round((sin((i*2-ARITH_CODER_QH_SCALE)*(M_PI/2/ARITH_CODER_QH_SCALE))+1.0)*(VAL_RANGE/2)));
}

static const uint8_t nSymbolsTab[9] = { 8, 2, 2, 5, 8, 16, 32, 64, 128 };
static const uint8_t lvl2BaseTab[8] = { 0, 1, 3, 8, 16, 32, 64, 128 };
static const int pageSzTab[9] = {
  ARITH_CODER_L1_PAGE_SZ    ,
  ARITH_CODER_L2_0_PAGE_SZ  ,
  ARITH_CODER_L2_1_PAGE_SZ  ,
  ARITH_CODER_L2_3_PAGE_SZ  ,
  ARITH_CODER_L2_8_PAGE_SZ  ,
  ARITH_CODER_L2_16_PAGE_SZ ,
  ARITH_CODER_L2_32_PAGE_SZ ,
  ARITH_CODER_L2_64_PAGE_SZ ,
  ARITH_CODER_L2_128_PAGE_SZ,
};

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
  if (0==mid || VAL_RANGE==mid) {
    *pRes = (0==mid);
    return 0;
  }
  unsigned val = get(VAL_RANGE);
  int res = (val >= mid);
  *pRes = res;
  unsigned lo = range_tab[res];
  unsigned ra = range_tab[res+1]-lo;
  return put(lo, ra);
}

// load_nChunks
// return  value:
//   on success nChunks > 0
//   on parsing error negative error code
static
int load_nChunks(CArithmeticDecoder* pDec)
{
  // use predefined value for sake of simplicity
  static const unsigned encTab[] = { 0, 4, 7, 8 };
  static const unsigned decTab[] = { 0, 0, 0, 0, 1, 1, 1, 2 };
  uint32_t val = 0;
  uint32_t msb = 1;
  for (int i = 0; i < 30; ++i) {
    // int symb = decTab[pDec->get(8)];
    int dVal = pDec->get(8);
    int symb = decTab[dVal];
    // printf("%d: [%d]=%d\n", i, dVal, symb);
    int ret = pDec->put(encTab[symb], encTab[symb+1]-encTab[symb]);
    if (ret < 0)
      return ret;
    if (symb == 2)
      return val+msb;
    val += (-symb) & msb;
    msb += msb;
  }
  return -15;
}

// load_quantized_histogram_tables
// return  value:
//   on success  0
//   on parsing error negative error code
static
int load_quantized_histogram_range_tab(unsigned range_tab[4][3], CArithmeticDecoder* pDec)
{
  // load hhw - combined hh word
  unsigned hhw = pDec->get(8*10*10*10);
  int err = pDec->put(hhw, 1);
  if (err)
    return err;

  // split hhw into individual hh words
  unsigned hh02 = (hhw %  8) + 1; hhw /= 8;
  unsigned hh01 = (hhw % 10) + 0; hhw /= 10;
  unsigned hh23 = (hhw % 10) + 0; hhw /= 10;
  unsigned hh45 = hhw + 0;

  ;
  range_tab[0][1] = quantized_histogram_pair_to_range_qh_scale9(hh02, VAL_RANGE);
  range_tab[1][1] = quantized_histogram_pair_to_range_qh_scale9(hh01, VAL_RANGE);
  range_tab[2][1] = quantized_histogram_pair_to_range_qh_scale9(hh23, VAL_RANGE);
  range_tab[3][1] = quantized_histogram_pair_to_range_qh_scale9(hh45, VAL_RANGE);
  for (int k = 0; k < 4; ++k) {
    range_tab[k][0] = 0;
    range_tab[k][2] = VAL_RANGE;
  }
  return 0;
}

// load_quantized_histogram
// return  value:
//   on success hlen >= 0
//   on parsing error negative error code
static
int load_quantized_histogram(uint8_t* qh, unsigned nCol, unsigned nChunks, CArithmeticDecoder* pDec, const unsigned range_tab[4][3])
{
  int val = 0;
  int32_t  prev_msb = 0;
  uint32_t rlAcc = 0;
  uint32_t rlMsb = 1;
  int  down  = 0;
  bool first = 1;
  int  maxMaxC  = -1;
  for (unsigned xi = 0, yi = 0; ; )
  {
    int32_t msb;
    int err = pDec->extract_1bit(range_tab[0], &msb);
    if (err)
      return err;

    if (msb != prev_msb) {
      // end of run
      int32_t rl = rlAcc + rlMsb - 2 + first;
      rlMsb = 1;
      rlAcc = 0;
      if (prev_msb == 0) {
        // end of zero run
        // printf("zero run %d at %d.%d\n", rl, yi, xi); fflush(stdout);
        while (rl > 0) {
          qh[(yi+1)*nCol-1-xi] = val;
          ++yi;
          if (yi == nChunks) {
            val = qh[nCol-1-xi];
            yi = 0;
            ++xi;
            if (xi == nCol)
              return -21;
          }
          --rl;
        }
      } else {
        rl += 1;
        // end of diff run

        if (down) {
          if (rl + 1 > val)
            return -17;
          val -= rl + 1;
        } else {
          val += rl + 1;
          if (val > 255)
            return -18;
        }
        // printf("val %d at %d.%d\n", val, yi, xi); fflush(stdout);
        qh[(yi+1)*nCol-1-xi] = val;
        if (maxMaxC < 0)
          maxMaxC = nCol-1-xi;
        ++yi;
        if (yi == nChunks) {
          val = qh[nCol-1-xi];
          yi = 0;
          ++xi;
        }
      }
      first = 0;
    }

    unsigned tabId = (msb==0) ? 0 : 2 - prev_msb;
    if (tabId != 2 || (val != 0 && val != 255)) {
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
          unsigned nxt_i = rl + xi*nChunks + yi;
          if (nxt_i >= nCol*nChunks) {
            if (nxt_i == nCol*nChunks) {
              while (rl > 0) {
                qh[(yi+1)*nCol-1-xi] = val;
                ++yi;
                if (yi == nChunks) {
                  val = qh[nCol-1-xi];
                  yi = 0;
                  ++xi;
                }
                --rl;
              }
              break; // done
            } else {
              return -19;
            }
          }
        } else {
          if (rl >= 255)
            return -20;
        }
      } else {
        down = lsb; // sign of difference
      }
    } else {
      down = (val != 0);
    }
    prev_msb = msb;
  }

  // for (int yi = 0; yi < nChunks; ++yi)
    // for (int xi = 0; xi <= maxMaxC; ++xi)
      // printf("qh[%3d][%3d]=%5d\n", yi, xi, qh[yi*257+xi]);

  return maxMaxC + 1; // success
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
  unsigned  m_nSymbols;
  uint16_t* m_c2low;      // [m_nSymbols+1]
  uint32_t* m_c2invRange; // [m_nSymbols];  floor(2^31/range)

  void init(unsigned nSymbols) {
    #ifdef ENABLE_PERF_COUNT
    m_extLookupCount       = 0;
    m_longLookupCount      = 0;
    m_renormalizationCount = 0;
    #endif
    m_c2low      = new uint16_t[nSymbols+1];
    m_c2invRange = new uint32_t[nSymbols];
    m_nSymbols   = nSymbols;
  }

  ~arithmetic_decode_model_t() {
    delete m_c2low;
    delete m_c2invRange;
  }

  int dequantize_and_prepare(const uint8_t* qh, unsigned hlen);
  int dequantize_and_prepare_pair(unsigned qh0);
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

int arithmetic_decode_model_t::dequantize_and_prepare(const uint8_t* qh, unsigned hlen)
{
  if (hlen > 0) {
    int nRanges = quantized_histogram_to_range(m_c2low, hlen, qh, VAL_RANGE);
    if (nRanges > 1) {
      for (unsigned c = hlen; c < m_nSymbols; ++c)
        m_c2low[c] = 0;
      prepare();
    }
    return nRanges;
  }
  return 0;
}

int arithmetic_decode_model_t::dequantize_and_prepare_pair(unsigned qh0)
{
  if (qh0 != 0) {
    if (qh0 != 255) {
      unsigned ra0 = quantized_histogram_pair_to_range_qh_scale255(qh0, VAL_RANGE);
      m_c2low[0] = ra0;
      m_c2low[1] = VAL_RANGE-ra0;
      prepare();
      return 2;
    } else {
      m_c2low[0] = 0;
    }
  } else {
    m_c2low[0] = 1;
  }
  return 1;
}

void arithmetic_decode_model_t::prepare()
{
  // m_c2low -> cumulative sums of ranges
  unsigned maxC = 0;
  unsigned lo = 0;
  unsigned invI = 0;
  for (unsigned c = 0; c < m_nSymbols; ++c) {
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
  m_c2low[m_nSymbols] = VAL_RANGE;
  unsigned r2c = maxC <= 255 ? maxC : 255;
  for (; invI <= RANGE2C_SZ; ++invI) {
    m_range2c[invI] = r2c;
  }
  m_maxC = maxC;
}

struct qh_descr_t {
  int32_t                    nChunks;
  int                        maxHLen;
  const uint8_t*             qh;
  arithmetic_decode_model_t* pModel;
  int                        theOnlyC;
  int                        sym_i;
};
int decode(
  uint8_t*                   dst0,
  int                        dstlen,
  int32_t                    histogram[256],
  CArithmeticDecoder*        pDec,
  qh_descr_t                 qhda[9])
{
  arithmetic_decode_model_t models[9];
  for (int i = 0; i < 9; ++i) {
    const int nSymbols = nSymbolsTab[i];
    qhda[i].sym_i = 0;
    models[i].init(nSymbols);
  }
  arithmetic_decode_model_t modelHalf;
  modelHalf.init(2);
  static const uint8_t qh_half[2] = {128,128};
  modelHalf.dequantize_and_prepare(qh_half, 2);

  // initialize move-to-front decoder table
  uint8_t mtf_t[256];
  for (int i = 0; i < 256; ++i)
    mtf_t[i] = i;
  // initialize RLE decoder
  uint32_t rlInc = 1;

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

  uint8_t* dst = dst0;
  int dst_i = 0;
  int prevC0 = 1;
  for (int lvl1chunk_i = 0; lvl1chunk_i < qhda[0].nChunks; ++lvl1chunk_i) {
    const uint8_t* lvl1Qh = qhda[0].qh;
    qhda[0].qh = lvl1Qh + nSymbolsTab[0] + 1;
    int nLvl1Ranges = models[0].dequantize_and_prepare(&lvl1Qh[1], qhda[0].maxHLen);
    const int32_t symbolsPerLvl1Chunk = ARITH_CODER_L1_PAGE_SZ*(int(lvl1Qh[0]) + 1);
    arithmetic_decode_model_t* pLvl1Model = &models[0];
    int theOnlyLvl1C = 0;
    if (nLvl1Ranges <= 1) {
      if (__builtin_expect(nLvl1Ranges == 0, 0))
        return -22; // chunk has to contain something
      theOnlyLvl1C = models[0].m_c2low[0];
      pLvl1Model = 0;
    }

    for (int lvl1sym_i = 0; lvl1sym_i < symbolsPerLvl1Chunk && dstlen != 0; ++lvl1sym_i) {
      int c0, c1;
      int c = theOnlyLvl1C;
      arithmetic_decode_model_t* pModel = pLvl1Model;
      for (int lvl = 0; lvl < 2; ++lvl) {
        if (pModel) {
          // begin arithmetic decode of one character
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
          c = pModel->val2c_estimate(value, invRange); // can be off by -1, much less likely by -2
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
              return dst-dst0;
            if (value > ((nxtRange<<RANGE_BITS)-1)) {
              return -103; // should not happen
            }
          }
          range = nxtRange;

          // update invRange
          invRange = umulh(invRange, uint64_t(pModel->m_c2invRange[c]) << (63-31)) << (RANGE_BITS+1);

          if ((dst_i & 15)==0) {
            // do one NR iteration to increase precision of invRange and assure that invRange*range <= 2**(97-rb)
            const uint64_t PROD_ONE = uint64_t(1) << (33-RANGE_BITS);
            uint64_t prod = umulh(invRange, range); // scaled to 2^(33-RANGE_BITS)
            invRange = umulh(invRange, (PROD_ONE*2-1-prod)<<(RANGE_BITS+30))<<1;
          }
          ++dst_i;
        }
        // arithmetic decode of one character done

        c1 = c;
        if (lvl == 0) {
          c0 = c;
          pModel = &modelHalf;
          if ((c0 | prevC0) != 0) {
            qh_descr_t* qhd = &qhda[c+1];
            pModel = qhd->pModel;
            c      = qhd->theOnlyC;
            int lvl2sym_i = qhd->sym_i;
            if (lvl2sym_i == 0) {
              int32_t lvl2nChunks = qhd->nChunks;
              if (__builtin_expect(lvl2nChunks == 0, 0))
                return -24; // chunk has to exist
              qhd->nChunks = lvl2nChunks - 1;
              const int nSymbols = nSymbolsTab[c0+1];
              const uint8_t* lvl2qh = qhd->qh;
              unsigned nCol = nSymbols==2 ? 2 : nSymbols+1;
              qhd->qh = lvl2qh + nCol;
              lvl2sym_i = pageSzTab[c0+1]*int32_t(lvl2qh[0]+1);
              pModel = &models[c0+1];
              int nLvl2Ranges = (nSymbols != 2) ?
                pModel->dequantize_and_prepare(&lvl2qh[1], qhd->maxHLen):
                pModel->dequantize_and_prepare_pair(lvl2qh[1]);
              qhd->theOnlyC = c = pModel->m_c2low[0];
              if (nLvl2Ranges <= 1) {
                if (__builtin_expect(nLvl2Ranges == 0, 0)) {
                  return -25; // chunk has to contain something
                }
                pModel = 0;
              }
              qhd->pModel = pModel;
            }
            qhd->sym_i = lvl2sym_i - 1;
          }
          prevC0 = c0;
        }
      }

      // RLE decode
      if (c0 == 0) {
        // zero run
        uint32_t delta_rl = c+rlInc;
        rlInc += delta_rl;

        if (delta_rl > uint32_t(dstlen))
          return -26; // zero run too long (A)

        // insert zero run
        int cm0 = mtf_t[0];
        histogram[cm0] += delta_rl;
        memset(dst, cm0, delta_rl);
        dst    += delta_rl;
        dstlen -= delta_rl;
      } else {
        rlInc = 1;
        // MTF decode
        int mtfC = lvl2BaseTab[c0] + c1;
        int dstC = mtf_t[mtfC];
        histogram[dstC] += 1;
        *dst++ = dstC;
        dstlen -= 1;
        // update move-to-front decoder table
        memmove(&mtf_t[1], &mtf_t[0], mtfC);
        mtf_t[0] = dstC;
      }
    }
    // printf("%016I64x %d %d %d %d\n", range, nLvl1Ranges, plain2chunk_i, plain2sym_i, theOnlyPlain2C);
  }
  return dst - dst0;
}

int decode(
  uint8_t*       dst0,
  int            dstlen,
  int32_t        histogram[256],
  const uint8_t* src,
  int            srclen)
{
  uint8_t  cntrs[256] = {0};
  uint8_t  prevQh[256] = {0};
  uint16_t currH[256];

  // initialize move-to-front decoder table
  uint8_t mtf_t[256];
  for (int i = 0; i < 256; ++i)
    mtf_t[i] = i;

  // initialize RLE decoder
  uint32_t rlInc = 1;

  // initialize arithmetic decoder
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);

  uint8_t tmpbuf[16] = {0};
  bool useTmpbuf = false;
  uint64_t value = 0; // scaled by 2**64
  for (int k = 0; k < 8; ++k) {
    value <<= 8;
    if (srclen > 0)
      value |= *src++;
    --srclen;
  }
  uint64_t range = uint64_t(1) << (64-RANGE_BITS);   // scaled by 2**(64-rb).  Maintained in [2**(33-rb)+1..2*(64-rb)]

  unsigned prevC0 = 0;
  uint8_t* dst = dst0;
  while (dstlen > 0) {
    unsigned tLo = 0, tHi = 256;
    do {
      unsigned tMid = (tLo + tHi)/2;
      int rleNonMsb = (tMid == 0) & prevC0; // Not a MS bit of RUNA/RUNB, assume equal probability of RUNA/RUNB
      int hVal = VAL_RANGE/2;
      if (!rleNonMsb) {
        int cntr = cntrs[tMid];
        if (cntr == 0) {
          unsigned prevQhVal = prevQh[tMid];
          unsigned qhVal = 42;//TODO;
          prevQh[tMid] = qhVal;
          currH[tMid] = qh2h_tab[qhVal];
        }
        ++cntr;
        if (cntr == ARITH_CODER_CNT_MAX)
          cntr = 0;
        cntrs[tMid] = cntr;
        hVal = currH[tMid];
      }
      unsigned b = 0;
      if (hVal != 0) {
        b = 1;
        if (hVal != VAL_RANGE) {
          b = (value >= range*(VAL_RANGE-hVal));
          // keep decoder in sync with encoder
          uint64_t cLo = b == 0 ? 0                : VAL_RANGE - hVal;
          uint64_t cHi = b == 0 ? VAL_RANGE - hVal : VAL_RANGE;
          value -= range * cLo;
          range *= (cHi-cLo);  // at this point range is scaled by 2**64 - the same scale as value
          uint64_t nxtRange = range >> RANGE_BITS;

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
            if (srclen < -7)
              return dst-dst0;
            if (value > ((nxtRange<<RANGE_BITS)-1)) {
              return -103; // should not happen
            }
          }
          range = nxtRange;
        }
      }
      tLo = (b == 0) ? tLo  : tMid + 1;
      tHi = (b == 0) ? tMid : tHi;
    } while (tLo != tHi);
    unsigned c = tLo;
    // arithmetic decode of one character done
    prevC0 = (c < 2);

    // RLE decode
    if (c < 2) {
      // zero run
      uint32_t delta_rl = c+rlInc;
      rlInc += delta_rl;

      if (delta_rl > uint32_t(dstlen))
        return -26; // zero run too long (A)

      // insert zero run
      int cm0 = mtf_t[0];
      histogram[cm0] += delta_rl;
      memset(dst, cm0, delta_rl);
      dst    += delta_rl;
      dstlen -= delta_rl;
    } else {
      rlInc = 1;
      // MTF decode
      int mtfC = c - 1;
      int dstC = mtf_t[mtfC];
      histogram[dstC] += 1;
      *dst++ = dstC;
      dstlen -= 1;
      // update move-to-front decoder table
      memmove(&mtf_t[1], &mtf_t[0], mtfC);
      mtf_t[0] = dstC;
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

  qh_descr_t qhda[9];
  int qhLen = 0;
  for (int i = 0; i < 9; ++i) {
    int nChunks = load_nChunks(&dec);
    if (nChunks < 0)
      return nChunks;
    nChunks -= (i > 0);
    if (nChunks > 1) {
      if (nChunks > INT_MAX/ARITH_CODER_L2_128_PAGE_SZ/2)
        return -31;
      if ((nChunks-1) * pageSzTab[i] > dstlen)
        return -32;
    }
    qhda[i].nChunks = nChunks;
    qhLen += nChunks* (nSymbolsTab[i]+1);
  }

  if (qhda[0].nChunks == 0)
    return -33;

  std::vector<uint8_t> qhVec(qhLen);
  uint8_t* qh = &qhVec.at(0);
  uint32_t qhSegSz = ARITH_CODER_MODEL_MAX_SEGMENT_SZ;
  unsigned qh_range_tab[4][3];
  for (int i = 0; i < 9; ++i) {
    int maxHLen = 0;
    int nChunks = qhda[i].nChunks;
    qhda[i].qh = qh;
    if (nChunks > 0) {
      uint32_t qhSz = nChunks* (nSymbolsTab[i]+1);
      qhSegSz += qhSz;
      if (qhSegSz > ARITH_CODER_MODEL_MAX_SEGMENT_SZ) {
        qhSegSz = qhSz;
        int ret = load_quantized_histogram_range_tab(qh_range_tab, &dec);
        if (ret != 0)
          return ret;
      }
      unsigned nCol = nSymbolsTab[i]==2 ? 2 : nSymbolsTab[i]+1;
      maxHLen = load_quantized_histogram(qh, nCol, nChunks, &dec, qh_range_tab);
      if (maxHLen < 0)
        return maxHLen;
      qh += qhSz;
    }
    qhda[i].maxHLen = maxHLen-1;
  }

  int modellen = srclen - dec.m_srclen;
  if (pInfo) {
    pInfo[0] = modellen*8;
    pInfo[1] = dec.m_srclen;
  }
  memset(histogram, 0, sizeof(histogram[0])*256);
  int textlen = decode(dst, dstlen, histogram, &dec, qhda);
  #ifdef ENABLE_PERF_COUNT
  if (pInfo) {
    pInfo[2] = model.m_extLookupCount;
    pInfo[3] = model.m_longLookupCount;
    pInfo[4] = model.m_renormalizationCount;
  }
  #endif
  return textlen;
}
