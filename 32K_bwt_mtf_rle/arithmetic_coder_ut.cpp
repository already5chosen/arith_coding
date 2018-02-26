// #include <cstdio>
// #include <cmath>
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

static const uint32_t logEst_tab[256] =
{ // round(log((i+2)/(i+1))*(2**31-1))
 2977044471, 1741459379, 1235585093,  958394255,
  783065124,  662072128,  573512964,  505874286,
  452519969,  409354105,  373711018,  343780812,
  318291317,  296322127,  277190838,  260380768,
  245493518,  232216947,  220303022,  209552159,
  199801946,  190918866,  182792152,  175329131,
  168451680,  162093474,  156197842,  150716071,
  145606056,  140831217,  136359621,  132163268,
  128217500,  124500523,  120992995,  117677698,
  114539249,  111563865,  108739157,  106053964,
  103498196,  101062717,   98739229,   96520181,
   94398686,   92368448,   90423704,   88559164,
   86769967,   85051637,   83400044,   81811374,
   80282100,   78808950,   77388892,   76019105,
   74696966,   73420032,   72186023,   70992811,
   69838405,   68720943,   67638678,   66589974,
   65573293,   64587191,   63630309,   62701366,
   61799157,   60922543,   60070452,   59241867,
   58435831,   57651433,   56887816,   56144163,
   55419702,   54713699,   54025458,   53354317,
   52699646,   52060847,   51437349,   50828609,
   50234108,   49653354,   49085875,   48531220,
   47988961,   47458685,   46940001,   46432531,
   45935917,   45449813,   44973890,   44507831,
   44051333,   43604103,   43165864,   42736346,
   42315291,   41902452,   41497591,   41100479,
   40710895,   40328628,   39953472,   39585232,
   39223718,   38868748,   38520144,   38177738,
   37841366,   37510870,   37186096,   36866898,
   36553134,   36244665,   35941359,   35643087,
   35349725,   35061152,   34777253,   34497915,
   34223028,   33952487,   33686191,   33424039,
   33165935,   32911788,   32661506,   32415001,
   32172190,   31932989,   31697319,   31465103,
   31236263,   31010729,   30788428,   30569291,
   30353252,   30140245,   29930207,   29723076,
   29518792,   29317297,   29118534,   28922448,
   28728985,   28538094,   28349722,   28163821,
   27980342,   27799238,   27620464,   27443974,
   27269725,   27097675,   26927783,   26760007,
   26594310,   26430651,   26268995,   26109304,
   25951543,   25795677,   25641672,   25489495,
   25339114,   25190496,   25043612,   24898431,
   24754923,   24613060,   24472814,   24334157,
   24197063,   24061504,   23927456,   23794894,
   23663792,   23534126,   23405874,   23279013,
   23153519,   23029370,   22906547,   22785026,
   22664788,   22545812,   22428079,   22311569,
   22196263,   22082143,   21969190,   21857387,
   21746716,   21637161,   21528703,   21421328,
   21315018,   21209758,   21105533,   21002327,
   20900125,   20798914,   20698678,   20599403,
   20501076,   20403683,   20307212,   20211648,
   20116980,   20023194,   19930278,   19838221,
   19747011,   19656635,   19567083,   19478343,
   19390404,   19303256,   19216888,   19131289,
   19046449,   18962359,   18879008,   18796386,
   18714484,   18633293,   18552803,   18473006,
   18393892,   18315453,   18237680,   18160565,
   18084100,   18008275,   17933084,   17858518,
   17784569,   17711230,   17638494,   17566353,
   17494799,   17423826,   17353427,   17283594,
   17214321,   17145601,   17077427,   17009794,
   16942694,   16876121,   16810070,   16744533,
};

// log_estimate - approximate round(log((x+1)/x)*(2**31-1))
// where x in range [1..2^16]
static uint32_t log_estimate(uint32_t x)
{
  if (x <= 256)
    return logEst_tab[x-1];
  else
    return uint32_t(-1)/x;
}

static void histogram_to_range(uint16_t* ranges, unsigned len, const uint32_t* h, uint32_t hTot, unsigned range_scale)
{
  // translate counts to ranges and store in ranges

  memset(ranges, 0, sizeof(*ranges)*len);

  // 1st pass - translate characters that are rounded range==1 from below
  uint32_t remCnt   = hTot;
  unsigned remRange = range_scale;
  unsigned nRanges = 0;
  do {
    uint32_t cntSum  = 0;
    nRanges = 0;
    for (unsigned c = 0; c < len; ++c) {
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
  for (unsigned c = 0; c < len; ++c) {
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
    int64_t de[260], denz[260];
    int denzLen = 0;
    for (unsigned c = 0; c < len; ++c) {
      int64_t deltaE = 0;
      unsigned cnt = h[c];
      if (cnt != 0) {
        int range = ranges[c];
        int nxtRange = range+1;
        if (nxtRange != 0) {
          // use integer approximation of log(), bit-precise because repeatability is more important than absolute precision
          denz[denzLen] = deltaE = -int64_t(log_estimate(range))*cnt;
          ++denzLen;
        }
      }
      de[c] = deltaE;
    }

    int indx = remRange - 1;
    std::nth_element(&denz[0], &denz[indx], &denz[denzLen]);
    int64_t thr = denz[indx];

    // increment ranges that will have maximal effect on entropy
    for (unsigned c = 0; c < len; ++c) {
      int64_t deltaE = de[c];
      if (deltaE != 0 && deltaE < thr) {
        ranges[c] += 1;
        remRange   -= 1;
      }
    }
    for (unsigned c = 0; remRange != 0; ++c) {
      int64_t deltaE = de[c];
      if (deltaE == thr) {
        ranges[c] += 1;
        remRange   -= 1;
      }
    }
  }
}

void quantized_histogram_to_range(uint16_t* ranges, unsigned len, const uint8_t* qh, unsigned range_scale)
{
  // translate quantized histogram from asin domain to linear domain
  uint32_t qhr[260];
  uint32_t qhrTot = 0;
  for (unsigned c = 0; c < len; ++c)
    qhrTot += (qhr[c] = dequantization_tab[qh[c]]);

  histogram_to_range(ranges, len, qhr, qhrTot, range_scale);
}

static const uint32_t dequantization_tab7[6] =
{ // round((sin((double(i+1)/7 - 0.5)*M_PI) + 1.0)*2**24);
 1661467,  6316793,  13043934,  20510498,  27237639, 31892965,
};
static const uint32_t dequantization_tab8[7] =
{ // round((sin((double(i+1)/8 - 0.5)*M_PI) + 1.0)*2**24);
 1277090, 4913933, 10356853, 16777216, 23197579, 28640499, 32277342,
};
static const uint32_t dequantization_tab9[8] =
{ // round((sin((double(i+1)/9 - 0.5)*M_PI) + 1.0)*2**24);
  1011790, 3925123, 8388608, 13863883, 19690549, 25165824, 29629309, 32542642,
};


// qh_scale in range [7..9]
// return range[0]
static unsigned quantized_histogram_pair_to_range(unsigned qh, unsigned range_scale, unsigned qh_scale, const uint32_t* tab)
{
  if (qh <= 0)
    return 0;
  if (qh >= qh_scale)
    return range_scale;
  return (uint64_t(tab[qh-1])*range_scale + (uint64_t(1) << 24)) >> 25;
}

unsigned quantized_histogram_pair_to_range_qh_scale7(unsigned qh, unsigned range_scale) {
  return quantized_histogram_pair_to_range(qh, range_scale, 7, dequantization_tab7);
}
unsigned quantized_histogram_pair_to_range_qh_scale8(unsigned qh, unsigned range_scale) {
  return quantized_histogram_pair_to_range(qh, range_scale, 8, dequantization_tab8);
}
unsigned quantized_histogram_pair_to_range_qh_scale9(unsigned qh, unsigned range_scale) {
  return quantized_histogram_pair_to_range(qh, range_scale, 9, dequantization_tab9);
}

