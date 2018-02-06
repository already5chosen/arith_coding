// #include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>

#include "arithmetic_coder_ut.h"

static const uint32_t dequantization_tab[256] =
{ // round((sin((double(qhVal)/256 - 0.5)*M_PI) + 1.0)*2**24);
        0,      632,     2526,     5684,
    10104,    15786,    22729,    30932,
    40393,    51112,    63086,    76314,
    90794,   106523,   123500,   141721,
   161185,   181887,   203825,   226996,
   251396,   277021,   303868,   331933,
   361211,   391697,   423388,   456279,
   490364,   525638,   562097,   599734,
   638545,   678522,   719661,   761954,
   805396,   849980,   895699,   942547,
   990516,  1039599,  1089789,  1141078,
  1193459,  1246923,  1301463,  1357070,
  1413735,  1471452,  1530209,  1590000,
  1650815,  1712644,  1775479,  1839309,
  1904126,  1969920,  2036680,  2104396,
  2173059,  2242659,  2313183,  2384623,
  2456966,  2530203,  2604323,  2679313,
  2755163,  2831862,  2909397,  2987758,
  3066931,  3146907,  3227671,  3309213,
  3391520,  3474579,  3558379,  3642906,
  3728147,  3814090,  3900723,  3988031,
  4076002,  4164622,  4253878,  4343757,
  4434246,  4525329,  4616995,  4709228,
  4802016,  4895344,  4989197,  5083563,
  5178427,  5273774,  5369590,  5465860,
  5562571,  5659707,  5757255,  5855198,
  5953524,  6052216,  6151259,  6250640,
  6350343,  6450352,  6550654,  6651232,
  6752072,  6853158,  6954476,  7056009,
  7157744,  7259663,  7361753,  7463997,
  7566381,  7668888,  7771504,  7874212,
  7976999,  8079847,  8182741,  8285667,
  8388608,  8491549,  8594475,  8697369,
  8800217,  8903004,  9005712,  9108328,
  9210835,  9313219,  9415463,  9517553,
  9619472,  9721207,  9822740,  9924058,
 10025144, 10125984, 10226562, 10326864,
 10426873, 10526576, 10625957, 10725000,
 10823692, 10922018, 11019961, 11117509,
 11214645, 11311356, 11407626, 11503442,
 11598789, 11693653, 11788019, 11881872,
 11975200, 12067988, 12160221, 12251887,
 12342970, 12433459, 12523338, 12612594,
 12701214, 12789185, 12876493, 12963126,
 13049069, 13134310, 13218837, 13302637,
 13385696, 13468003, 13549545, 13630309,
 13710285, 13789458, 13867819, 13945354,
 14022053, 14097903, 14172893, 14247013,
 14320250, 14392593, 14464033, 14534557,
 14604157, 14672820, 14740536, 14807296,
 14873090, 14937907, 15001737, 15064572,
 15126401, 15187216, 15247007, 15305764,
 15363481, 15420146, 15475753, 15530293,
 15583757, 15636138, 15687427, 15737617,
 15786700, 15834669, 15881517, 15927236,
 15971820, 16015262, 16057555, 16098694,
 16138671, 16177482, 16215119, 16251578,
 16286852, 16320937, 16353828, 16385519,
 16416005, 16445283, 16473348, 16500195,
 16525820, 16550220, 16573391, 16595329,
 16616031, 16635495, 16653716, 16670693,
 16686422, 16700902, 16714130, 16726104,
 16736823, 16746284, 16754487, 16761430,
 16767112, 16771532, 16774690, 16776584,
};

static void histogram_to_range(uint16_t* ranges, unsigned maxC, const uint32_t* h, uint32_t hTot, unsigned range_scale)
{
  // translate counts to ranges and store in ranges

  memset(ranges, 0, sizeof(*ranges)*(maxC+1));

  // 1st pass - translate characters that are rounded range==1 from below
  uint32_t remCnt   = hTot;
  unsigned remRange = range_scale;
  unsigned nRanges = 0;
  do {
    uint32_t cntSum  = 0;
    nRanges = 0;
    for (unsigned c = 0; c <= maxC; ++c) {
      uint32_t cnt = h[c];
      if (cnt != 0 && ranges[c]==0) {
        if (uint64_t(cnt)*remRange < remCnt) {
          cntSum    += cnt;
          nRanges   += 1;
          ranges[c] = 1;
        }
      }
    }
    remCnt   -= cntSum;
    remRange -= nRanges;
  } while (nRanges != 0);

  // 2nd pass - translate remaining characters while rounding toward zero
  unsigned rangeSum = 0;
  for (unsigned c = 0; c <= maxC; ++c) {
    uint32_t cnt = h[c];
    if (cnt != 0) {
      unsigned range = ranges[c];
      if (range == 0)
        ranges[c] = range = (uint64_t(cnt)*remRange)/remCnt;
      rangeSum += range;
    }
  }

  if (rangeSum < range_scale) {
    remRange = range_scale - rangeSum;
    // calculate effect of increment of range on entropy
    double de[256], denz[256];
    int denzLen = 0;
    for (unsigned c = 0; c <= maxC; ++c) {
      double deltaE = 0;
      unsigned cnt = h[c];
      if (cnt != 0) {
        int range = ranges[c];
        int nxtRange = range+1;
        if (nxtRange != 0) {
          denz[denzLen] = deltaE = log2(double(range)/nxtRange)*cnt;
          ++denzLen;
        }
      }
      de[c] = deltaE;
    }

    int indx = remRange - 1;
    std::nth_element(&denz[0], &denz[indx], &denz[denzLen]);
    double thr = denz[indx];

    // increment ranges that will have maximal effect on entropy
    for (unsigned c = 0; c <= maxC; ++c) {
      double deltaE = de[c];
      if (deltaE != 0 && deltaE < thr) {
        ranges[c] += 1;
        remRange   -= 1;
      }
    }
    for (unsigned c = 0; remRange != 0; ++c) {
      double deltaE = de[c];
      if (deltaE == thr) {
        ranges[c] += 1;
        remRange   -= 1;
      }
    }
  }
}

void quantized_histogram_to_range(uint16_t* ranges, unsigned maxC, const uint8_t* qh, unsigned range_scale)
{
  // translate quantized histogram from asin domain to linear domain
  uint32_t qhr[256];
  uint32_t qhrTot = 0;
  for (unsigned c = 0; c <= maxC; ++c)
    qhrTot += (qhr[c] = dequantization_tab[qh[c]]);
 
  histogram_to_range(ranges, maxC, qhr, qhrTot, range_scale);
}

