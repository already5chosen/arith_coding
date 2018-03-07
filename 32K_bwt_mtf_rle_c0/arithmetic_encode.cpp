#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <algorithm>

#include "arithmetic_encode.h"
#include "arithmetic_coder_ut.h"
#include "arithmetic_coder_cfg.h"

static const int      QH_SCALE  = 1 << QH_BITS;
static const unsigned VAL_RANGE = 1u << RANGE_BITS;

static void resize_up(std::vector<uint32_t>* vec, size_t len)
{
  if (vec->size() < len)
    vec->resize(len);
}

struct context_chunk_t {
  int      maxC;
  int      srclen;
  uint32_t histogram[257];
  uint8_t  qHistogram[257];
  uint16_t ranges[257];
};

static const int CONTEX_HDR_LEN = 3;
static const int CONTEX_CHUNK_LEN = (sizeof(context_chunk_t)-1)/sizeof(uint32_t) + 1;

void arithmetic_encode_init_context(std::vector<uint32_t>* context, int tilelen)
{
  int nChunks = tilelen/(ARITH_CODER_RUNS_PER_CHUNK*2) + 1; // estimate # of chunks
  resize_up(context, CONTEX_HDR_LEN+CONTEX_CHUNK_LEN*nChunks);
  (*context)[0] = 0;
  (*context)[1] = 0;
}

int arithmetic_encode_chunk_callback(void* context, const uint8_t* src, int chunklen, int nRuns)
{
  if (src) {
    std::vector<uint32_t>* vec = static_cast<std::vector<uint32_t>*>(context);
    uint32_t chunk_i = (*vec)[0];
    resize_up(vec, CONTEX_HDR_LEN+CONTEX_CHUNK_LEN*(chunk_i+1)+((chunk_i+1)*258)/sizeof(uint32_t)+1);
    (*vec)[0] = chunk_i + 1;
    (*vec)[1] = nRuns;

    context_chunk_t* chunk = reinterpret_cast<context_chunk_t*>(&vec->at(CONTEX_HDR_LEN+CONTEX_CHUNK_LEN*chunk_i));
    chunk->srclen = chunklen;
    uint32_t* histogram = chunk->histogram;
    memset(histogram, 0, sizeof(chunk->histogram));
    for (int i = 0; i < chunklen; ++i) {
      int c = src[i];
      if (c == 255) {
        // escape sequence for symbols 255 and 256
        c = int(src[i+1]) + 1;
        ++i;
      }
      ++histogram[c];
    }
  }
  return ARITH_CODER_RUNS_PER_CHUNK;
}


static void quantize_histogram(uint8_t* __restrict qh, uint32_t *h, int len, uint32_t hTot)
{
  double invhLen = 2.0/hTot;
  double M_PI = 3.1415926535897932384626433832795;
  for (unsigned c = 0; c < len; ++c) {
    unsigned val = 0;
    unsigned cnt = h[c];
    if (cnt != 0) {
      val = round((asin(cnt*invhLen - 1.0)/M_PI + 0.5)*QH_SCALE);
      if (val < 1)
        val = 1;
      else if (val >= QH_SCALE)
        val = QH_SCALE-1;
    }
    qh[c] = val;
  }
}

static void prepare1(uint32_t * context, double* pInfo)
{
  int nChunks = context[0];
  int runsInLastChunk = context[1];
  if (runsInLastChunk < ARITH_CODER_RUNS_PER_CHUNK/2 && nChunks > 1) {
    // add last chunk to before-last
    context_chunk_t* dstChunk = reinterpret_cast<context_chunk_t*>(&context[CONTEX_HDR_LEN+CONTEX_CHUNK_LEN*(nChunks-2)]);
    context_chunk_t* srcChunk = reinterpret_cast<context_chunk_t*>(&context[CONTEX_HDR_LEN+CONTEX_CHUNK_LEN*(nChunks-1)]);
    dstChunk->srclen += srcChunk->srclen;
    for (int i = 0; i < 257; ++i) {
      dstChunk->histogram[i] += srcChunk->histogram[i];
    }
    nChunks -= 1;
    context[0] = nChunks;
  }

  double entropy = 0;
  unsigned maxMaxC = 0;
  for (int chunk_i = 0; chunk_i < nChunks; ++chunk_i) {
    context_chunk_t* chunk = reinterpret_cast<context_chunk_t*>(&context[CONTEX_HDR_LEN+CONTEX_CHUNK_LEN*chunk_i]);
    uint32_t * histogram = chunk->histogram;
    // find highest-numbered character that occurred at least once
    unsigned maxC = 0;
    int32_t tot = 0;
    for (unsigned c = 0; c < 257; ++c) {
      tot += histogram[c];
      if (histogram[c] != 0)
        maxC = c;
    }
    chunk->maxC = maxC;
    if (maxC > maxMaxC)
      maxMaxC = maxC;

    if (pInfo) {
      // calculate source entropy of static part of histogram
      entropy += log2(double(tot))*tot;
      for (unsigned c = 0; c <= maxC; ++c) {
        int32_t cnt = histogram[c];
        if (cnt)
          entropy -= log2(double(cnt))*cnt;
      }
    }
    quantize_histogram(chunk->qHistogram, histogram, maxC+1, tot);
  }
  context[2] = maxMaxC;

  if (pInfo) {
    pInfo[0] = entropy;
    pInfo[3] = 0;
  }
}

// return entropy estimate after quantization
static double prepare2(uint32_t * context)
{
  int nChunks = context[0];
  double entropy = 0;
  for (int chunk_i = 0; chunk_i < nChunks; ++chunk_i) {
    context_chunk_t* chunk = reinterpret_cast<context_chunk_t*>(&context[CONTEX_HDR_LEN+CONTEX_CHUNK_LEN*chunk_i]);
    unsigned maxC = chunk->maxC;
    quantized_histogram_to_range(chunk->ranges, maxC+1, chunk->qHistogram, VAL_RANGE);
    // calculate entropy after quantization
    for (unsigned c = 0; c <= maxC; ++c) {
      unsigned cnt = chunk->histogram[c];
      if (cnt)
        entropy -= log2(chunk->ranges[c]/double(VAL_RANGE))*cnt;
    }
    // printf("%10d\n", int(entropy/8));
  }
  return entropy;
}

class CArithmeticEncoder {
public:
  uint8_t* put(uint64_t cScale, uint64_t cLo, uint64_t cRange, uint8_t* dst);
  void spillOverflow(uint8_t* dst);
  void init() {
    m_lo = 0;
    m_range = uint64_t(1) << 63;
  }
  uint64_t m_lo;    // scaled by 2**63
  uint64_t m_range; // scaled by 2**63
};

// return quantized code of h0/(h0+h1) in range [0..scale]
static unsigned quantize_histogram_pair(unsigned h0, unsigned h1, unsigned scale)
{
  if (h0 == 0)
    return 0;
  if (h1 == 0)
    return scale;
  unsigned hTot = h0 + h1;
  // printf("%u+%u => %.2f\n", h0, h1, (log2(hTot)*hTot-log2(h0)*h0-log2(h1)*h1)/8);
  double M_PI = 3.1415926535897932384626433832795;
  unsigned val = round((asin(h0*2.0/hTot - 1.0)/M_PI + 0.5)*scale);
  if (val < 1)
    val = 1;
  else if (val >= scale)
    val = scale-1;
  return val;
}

static unsigned add_value_to_histogram(unsigned* h, unsigned val)
{
  val += 2;
  unsigned nbits = 0;
  do {
    ++nbits;
    ++h[val % 2];
    val /= 2;
  } while (val > 1);
  return nbits;
}

static uint8_t* encode_value(uint8_t* dst, unsigned val, const unsigned* range_tab0, const unsigned* range_tab1, CArithmeticEncoder* pEnc)
{
  unsigned lo0 = range_tab0[0];
  unsigned ra0 = range_tab0[1] - lo0;
  val += 2;
  do {
    dst = pEnc->put(VAL_RANGE, lo0, ra0, dst);
    unsigned bit = val % 2;
    unsigned lo1 = range_tab1[bit];
    unsigned ra1 = range_tab1[bit+1] - lo1;
    if (ra1 != 0)
      dst = pEnc->put(VAL_RANGE, lo1, ra1, dst);
    val /= 2;
  } while (val > 1);
  return dst;
}

static uint8_t* store_model_store_base_data(uint8_t* dst, const uint8_t* qh, int len, int nChunks, CArithmeticEncoder* pEnc)
{
  // encode backward
  // first pass - calculate histograms
  unsigned hist[8] = {0};
  int runlen = 257-len-2;
  unsigned prev_val = 0;
  for (unsigned c = 0; c < len; ++c) {
    ++runlen;
    unsigned val = qh[len-1-c];
    if (val != prev_val) {
      if (runlen >= 0)
        hist[0] += add_value_to_histogram(&hist[2], runlen);
      runlen = -1;
      unsigned diff;
      hist[1] += 1;
      if (val > prev_val) {
        hist[6] += (prev_val != 0); // '+'
        diff = val - prev_val - 1;
      } else {
        hist[7] += (prev_val != 255); // '-'
        diff = prev_val - val - 1;
      }
      prev_val = val;
      if (diff != 0)
        hist[1] += add_value_to_histogram(&hist[4], diff-1);
    }
  }
  if (nChunks == 1) {
    hist[0] += add_value_to_histogram(&hist[2], runlen+2);
  } else {
    hist[0] += add_value_to_histogram(&hist[2], runlen+1);
    hist[1] += add_value_to_histogram(&hist[4], nChunks-2);
    hist[0] += add_value_to_histogram(&hist[2], 0);
  }

  unsigned hh01 = quantize_histogram_pair(hist[0], hist[1], 9); // range [1..8], because both sums > 0
  unsigned hh23 = quantize_histogram_pair(hist[2], hist[3], 9); // range [0..9]
  unsigned hh45 = quantize_histogram_pair(hist[4], hist[5], 9); // range [0..9]
  unsigned hh67 = quantize_histogram_pair(hist[6], hist[7], 9); // range [0..9]

  unsigned range_tab[4][3];
  range_tab[0][1] = quantized_histogram_pair_to_range_qh_scale9(hh01, VAL_RANGE);
  range_tab[1][1] = quantized_histogram_pair_to_range_qh_scale9(hh23, VAL_RANGE);
  range_tab[2][1] = quantized_histogram_pair_to_range_qh_scale9(hh45, VAL_RANGE);
  range_tab[3][1] = quantized_histogram_pair_to_range_qh_scale9(hh67, VAL_RANGE);
  for (int k = 0; k < 4; ++k) {
    range_tab[k][0] = 0;
    range_tab[k][2] = VAL_RANGE;
    // printf("range_tab[%d][1] = %5d\n", k, range_tab[k][1]);
  }
  // printf("%.5f %.5f\n", double(hist[0]+hist[1])/(hist[0]+hist[1]+hist[2]+hist[3]), double(range_tab[0][1])/VAL_RANGE);
  // printf("%.5f %.5f\n", double(hist[0])/(hist[0]+hist[1]), double(range_tab[1][1])/VAL_RANGE);
  // printf("%.5f %.5f\n", double(hist[2])/(hist[2]+hist[3]), double(range_tab[2][1])/VAL_RANGE);

  unsigned hhw = (hh01-1)+8*(hh23+(10*(hh45+10*hh67))); // combine all hh in a single word

  // store hhw
  dst = pEnc->put(8*10*10*10, hhw, 1, dst);

  // second pass - encode qh and nChunks
  runlen = 257-len-2;
  prev_val = 0;
  for (unsigned c = 0; c < len; ++c) {
    ++runlen;
    unsigned val = qh[len-1-c];
    if (val != prev_val) {
      if (runlen >= 0)
        dst = encode_value(dst, runlen, &range_tab[0][0], range_tab[1], pEnc);
      runlen = -1;
      unsigned diff;
      dst = pEnc->put(VAL_RANGE, range_tab[0][1], VAL_RANGE-range_tab[0][1], dst);
      if (val > prev_val) {
        if (prev_val != 0 && range_tab[3][1] != 0)
          dst = pEnc->put(VAL_RANGE, 0, range_tab[3][1], dst); // '+'
        diff = val - prev_val - 1;
      } else {
        if (prev_val != 255 && range_tab[3][1] != VAL_RANGE)
          dst = pEnc->put(VAL_RANGE, range_tab[3][1], VAL_RANGE-range_tab[3][1], dst); // '-'
        diff = prev_val - val - 1;
      }
      prev_val = val;
      if (diff != 0)
        dst = encode_value(dst, diff-1, &range_tab[0][1], range_tab[2], pEnc);
    }
  }
  if (nChunks == 1) {
    dst = encode_value(dst, runlen+2, &range_tab[0][0], range_tab[1], pEnc);
  } else {
    dst = encode_value(dst, runlen+1,  &range_tab[0][0], range_tab[1], pEnc);
    dst = encode_value(dst, nChunks-2, &range_tab[0][1], range_tab[2], pEnc);
    dst = encode_value(dst, 0,         &range_tab[0][0], range_tab[1], pEnc);
  }

  return dst;
}

static uint8_t* store_model_store_diff_data(uint8_t* dst, const uint8_t* baseQh, const uint8_t* qhT, int len, int nChunks, CArithmeticEncoder* pEnc)
{
  // first pass - calculate histograms
  unsigned hist[8] = {0};
  int runlen = -2;
  for (int i = 0; i < len; ++i) {
    const uint8_t* row = &qhT[(len-1-i)*nChunks];
    unsigned baseVal = baseQh[len-1-i];
    for (int chunk_i = 0; chunk_i < nChunks; ++chunk_i) {
      unsigned val = row[chunk_i];
      ++runlen;
      if (val != baseVal) {
        if (runlen >= 0)
          hist[0] += add_value_to_histogram(&hist[2], runlen);
        runlen = -1;
        unsigned diff;
        hist[1] += 1;
        if (val > baseVal) {
          hist[6] += (baseVal != 0); // '+'
          diff = val - baseVal - 1;
        } else {
          hist[7] += (baseVal != 255); // '-'
          diff = baseVal - val - 1;
        }
        if (diff != 0)
          hist[1] += add_value_to_histogram(&hist[4], diff-1);
      }
    }
  }
  hist[0] += add_value_to_histogram(&hist[2], runlen+1);


  unsigned hh01 = quantize_histogram_pair(hist[0], hist[1], 9); // range [1..8], because both sums > 0
  unsigned hh23 = quantize_histogram_pair(hist[2], hist[3], 9); // range [0..9]
  unsigned hh45 = quantize_histogram_pair(hist[4], hist[5], 9); // range [0..9]
  unsigned hh67 = quantize_histogram_pair(hist[6], hist[7], 9); // range [0..9]

  unsigned range_tab[4][3];
  range_tab[0][1] = quantized_histogram_pair_to_range_qh_scale9(hh01, VAL_RANGE);
  range_tab[1][1] = quantized_histogram_pair_to_range_qh_scale9(hh23, VAL_RANGE);
  range_tab[2][1] = quantized_histogram_pair_to_range_qh_scale9(hh45, VAL_RANGE);
  range_tab[3][1] = quantized_histogram_pair_to_range_qh_scale9(hh67, VAL_RANGE);
  for (int k = 0; k < 4; ++k) {
    range_tab[k][0] = 0;
    range_tab[k][2] = VAL_RANGE;
    // printf("range_tab[%d][1] = %5d\n", k, range_tab[k][1]);
  }
  // printf("%.5f %.5f\n", double(hist[0]+hist[1])/(hist[0]+hist[1]+hist[2]+hist[3]), double(range_tab[0][1])/VAL_RANGE);
  // printf("%.5f %.5f\n", double(hist[0])/(hist[0]+hist[1]), double(range_tab[1][1])/VAL_RANGE);
  // printf("%.5f %.5f\n", double(hist[2])/(hist[2]+hist[3]), double(range_tab[2][1])/VAL_RANGE);

  unsigned hhw = (hh01-1)+8*(hh23+(10*(hh45+10*hh67))); // combine all hh in a single word

  // store hhw
  dst = pEnc->put(8*10*10*10, hhw, 1, dst);

  // second pass - encode
  runlen = -2;
  for (int i = 0; i < len; ++i) {
    const uint8_t* row = &qhT[(len-1-i)*nChunks];
    unsigned baseVal = baseQh[len-1-i];
    for (int chunk_i = 0; chunk_i < nChunks; ++chunk_i) {
      unsigned val = row[chunk_i];
      ++runlen;
      if (val != baseVal) {
        if (runlen >= 0)
          dst = encode_value(dst, runlen, &range_tab[0][0], range_tab[1], pEnc);
        runlen = -1;
        unsigned diff;
        dst = pEnc->put(VAL_RANGE, range_tab[0][1], VAL_RANGE-range_tab[0][1], dst);
        if (val > baseVal) {
          if (baseVal != 0 && range_tab[3][1] != 0)
            dst = pEnc->put(VAL_RANGE, 0, range_tab[3][1], dst); // '+'
          diff = val - baseVal - 1;
        } else {
          if (baseVal != 255 && range_tab[3][1] != VAL_RANGE)
            dst = pEnc->put(VAL_RANGE, range_tab[3][1], VAL_RANGE-range_tab[3][1], dst); // '-'
          diff = baseVal - val - 1;
        }
        if (diff != 0)
          dst = encode_value(dst, diff-1, &range_tab[0][1], range_tab[2], pEnc);
      }
    }
  }
  dst = encode_value(dst, runlen+1, &range_tab[0][0], range_tab[1], pEnc);

  return dst;
}

// return the number of stored octets
static int store_model(uint8_t* dst, uint32_t * context, double* pNbits, CArithmeticEncoder* pEnc)
{
  int nChunks = context[0];
  unsigned maxMaxC = context[2];
  uint8_t* dst0 = dst;

  if (nChunks == 1) {
    context_chunk_t* chunk = reinterpret_cast<context_chunk_t*>(&context[CONTEX_HDR_LEN]);
    dst = store_model_store_base_data(dst, chunk->qHistogram, maxMaxC+1, nChunks, pEnc);
  } else {
    // transpose qHistogram[][]
    uint8_t* qhT = reinterpret_cast<uint8_t*>(&context[CONTEX_HDR_LEN+CONTEX_CHUNK_LEN*nChunks]);
    for (int chunk_i = 0; chunk_i < nChunks; ++chunk_i) {
      context_chunk_t* chunk = reinterpret_cast<context_chunk_t*>(&context[CONTEX_HDR_LEN+CONTEX_CHUNK_LEN*chunk_i]);
      int maxC = chunk->maxC;
      uint8_t* dst = &qhT[chunk_i];
      for (int i = 0; i <= maxC; ++i, dst += nChunks)
        *dst = chunk->qHistogram[i];
      for (int i = maxC+1; i <= maxMaxC; ++i, dst += nChunks)
        *dst = 0;
    }

    // fill qh with median values of rows of qhT[][]
    uint8_t qh[257];
    uint8_t* row = &qhT[(maxMaxC+1)*nChunks];
    for (int i = 0; i <= maxMaxC; ++i) {
      memcpy(row, &qhT[i*nChunks], nChunks);
      std::nth_element(&row[0], &row[nChunks/2], &row[nChunks]);
      qh[i] = row[nChunks/2];
    }

    // store qh and nChunks
    dst = store_model_store_base_data(dst, qh, maxMaxC+1, nChunks, pEnc);
    // store qhT-qh
    dst = store_model_store_diff_data(dst, qh, qhT, maxMaxC+1, nChunks, pEnc);
  }

  int len = dst - dst0;
  *pNbits = len*8.0 + 63 - log2(pEnc->m_range);
  return len;
}

// static uint8_t* inc_dst(uint8_t* dst0, uint8_t* dst) {
  // uint8_t val = 0;
  // uint8_t* inc;
  // for (inc = dst-1; val == 0 && inc != dst0; --inc)
    // *inc = (val = *inc + 1);
  // if (val == 0) {
    // memmove(dst0+1, dst0, dst-dst0);
    // dst0[0] = 1;
    // dst += 1;
  // }
  // return dst;
// }

static void inc_dst(uint8_t* dst) {
  uint8_t val;
  do {
    --dst;
    *dst = (val = *dst + 1);
  } while (val==0);
}

static void range2low(uint16_t* c2low, const uint16_t* c2range, unsigned len)
{
  unsigned lo = 0;
  for (unsigned c = 0; c < len; ++c) {
    c2low[c] = lo;
    lo += c2range[c];
  }
  c2low[len] = VAL_RANGE;
}

static int encode(uint8_t* dst, const uint8_t* src, uint32_t* context, CArithmeticEncoder* pEnc)
{
  const uint64_t MSB_MSK   = uint64_t(255) << 56;
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);
  uint64_t lo     = pEnc->m_lo << 1;              // scaled by 2**64
  uint64_t range  = pEnc->m_range >> (RANGE_BITS-1); // scaled by 2**50
  uint8_t* dst0   = dst;
  uint64_t prevLo = lo;
  // int dbg_i = 0;
  int nChunks = context[0];
  for (int chunk_i = 0; chunk_i < nChunks; ++chunk_i) {
    context_chunk_t* chunk = reinterpret_cast<context_chunk_t*>(&context[CONTEX_HDR_LEN+CONTEX_CHUNK_LEN*chunk_i]);
    uint16_t c2low[258]; // c2low -> cumulative sums of ranges
    range2low(c2low, chunk->ranges, chunk->maxC+1);

    unsigned srclen = chunk->srclen;
    for (unsigned i = 0; i < srclen; ++i) {
      int c = src[i];
      if (c == 255) {
        // escape sequence for symbols 255 and 256
        c = int(src[i+1]) + 1;
        ++i;
      }
      // printf("[%3d]=%3d\n", dbg_i, c); ++dbg_i;
      uint64_t cLo = c2low[c+0];
      uint64_t cRa = c2low[c+1] - cLo;
      lo   += range * cLo;
      range = range * cRa;

      // at this point range is scaled by 2**64 - the same scale as lo
      uint64_t nxtRange = range >> RANGE_BITS;
      if (nxtRange <= MIN_RANGE) {
        // re-normalize
        if (lo < prevLo) // lo overflow
          inc_dst(dst); //dst = inc_dst(dst0, dst);
        dst[0] = lo >> (64-8*1);
        dst[1] = lo >> (64-8*2);
        dst[2] = lo >> (64-8*3);
        dst += 3;
        lo <<= 24;
        nxtRange = range << (24-RANGE_BITS);
        prevLo = lo;
      }
      range = nxtRange;
    }
    src += srclen;
    // printf(": %10d\n", dst-dst0);
  }

  // output last bits
  if (lo < prevLo) // lo overflow
    inc_dst(dst); //dst = inc_dst(dst0, dst);
  uint64_t hi = lo + ((range<<RANGE_BITS)-1);
  if (hi > lo) {
    uint64_t dbits = lo ^ hi;
    while ((dbits & MSB_MSK)==0) {
      // lo and hi have the same upper octet
      *dst++ = uint8_t(hi>>56);
      hi    <<= 8;
      dbits <<= 8;
    }
    // put out last octet
    *dst++ = uint8_t(hi>>56);
  } else {
    inc_dst(dst); //dst = inc_dst(dst0, dst);
  }
  return dst - dst0;
}

uint8_t* CArithmeticEncoder::put(uint64_t cScale, uint64_t cLo, uint64_t cRange, uint8_t* dst)
{
  const uint64_t MIN_RANGE = uint64_t(1) << (33-1);
  const uint64_t BIT_63    = uint64_t(1) << 63;
  const uint64_t MSK_63    = BIT_63 - 1;
  uint64_t range = m_range / cScale;
  uint64_t lo = m_lo + range * cLo;
  range *= cRange;

  if (range <= MIN_RANGE) {
    // re-normalize
    if (lo >= BIT_63) // lo overflow
      inc_dst(dst);
    dst[0] = lo >> (63-8*1);
    dst[1] = lo >> (63-8*2);
    dst[2] = lo >> (63-8*3);
    dst += 3;
    lo = (lo << 24) & MSK_63;
    range <<= 24;
  }
  m_lo    = lo;
  m_range = range;
  return dst;
}

void CArithmeticEncoder::spillOverflow(uint8_t* dst)
{
  const uint64_t BIT_63 = uint64_t(1) << 63;
  const uint64_t MSK_63 = BIT_63 - 1;
  if (m_lo >= BIT_63) // lo overflow
    inc_dst(dst);
  m_lo &= MSK_63;
}

// return value:
//  0 - not compressible, because all input characters have approximately equal probability or because input is too short
// >0 - the length of compressed buffer
int arithmetic_encode(uint32_t* context, uint8_t* dst, const uint8_t* src, int srclen, int origlen, double* pInfo)
{
  prepare1(context, pInfo);

  CArithmeticEncoder enc;
  enc.init();
  double modelLenBits;
  unsigned modellen = store_model(dst, context, &modelLenBits, &enc);
  if (pInfo)
    pInfo[1] = modelLenBits;

  double quantizedEntropy = prepare2(context);

  int lenEst = (modelLenBits+quantizedEntropy+7)/8;
  if (pInfo)
    pInfo[3] = quantizedEntropy;
  if (lenEst >= origlen)
    return 0; // not compressible

  int dstlen = encode(&dst[modellen], src, context, &enc);

  if (pInfo)
    pInfo[2] = dstlen*8.0;

  int reslen = modellen + dstlen;
  // printf("%d+%d=%d(%d) <> %d\n", modellen, dstlen, reslen, lenEst, origlen);

  if (reslen >= origlen)
    return 0; // not compressible

  return reslen;
}
