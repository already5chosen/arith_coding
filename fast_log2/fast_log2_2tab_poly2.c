#include <stdint.h>
#include <stdio.h>
#include <math.h>

static const uint8_t mx1_tab[64] = {
 255, 248, 240, 233,
 225, 218, 212, 205,
 199, 192, 186, 180,
 175, 169, 164, 158,
 153, 148, 143, 138,
 134, 129, 125, 120,
 116, 112, 108, 104,
 100, 96 , 92 , 88 ,
 85 , 81 , 78 , 74 ,
 71 , 68 , 65 , 62 ,
 59 , 56 , 53 , 50 ,
 47 , 44 , 41 , 39 ,
 36 , 33 , 31 , 28 ,
 26 , 24 , 21 , 19 ,
 17 , 14 , 12 , 10 ,
 8  , 6  , 4  , 2  ,
};

static const uint32_t lg1_tab[64] = {
 536870912  , 682398916  , 937091716  , 1103362014 ,
 1382486880 , 1570799997 , 1664072623 , 1872638678 ,
 1983759008 , 2213741774 , 2343750303 , 2482852922 ,
 2516421264 , 2672857097 , 2721242383 , 2895923748 ,
 2959916243 , 3031281479 , 3110202637 , 3196869810 ,
 3164538868 , 3265651333 , 3245145759 , 3361388697 ,
 3353276123 , 3350864557 , 3354278275 , 3363645664 ,
 3379099398 , 3400776639 , 3428819239 , 3463373953 ,
 3359437008 , 3405756636 , 3310825865 , 3369405046 ,
 3283863116 , 3202468219 , 3125297270 , 3052429348 ,
 2983945771 , 2919930186 , 2860468654 , 2805649749 ,
 2755564656 , 2710307269 , 2669974306 , 2466913573 ,
 2435010941 , 2408305329 , 2214485088 , 2196683046 ,
 2008902480 , 1823597751 , 1818795892 , 1639837179 ,
 1463481486 , 1472347728 , 1302668303 , 1135729139 ,
 971571443  , 810237358  , 651769995  , 496213457  ,
};

static const uint8_t mx2_tab[71] = {
 137 ,
 136 ,
 134 ,
 132 ,
 130 ,
 128 ,
 126 ,
 124 ,
 122 ,
 120 ,
 118 ,
 116 ,
 114 ,
 112 ,
 110 ,
 108 ,
 106 ,
 104 ,
 102 ,
 100 ,
  98 ,
  96 ,
  94 ,
  92 ,
  90 ,
  88 ,
  86 ,
  84 ,
  82 ,
  80 ,
  78 ,
  76 ,
  74 ,
  72 ,
  70 ,
  68 ,
  66 ,
  64 ,
  62 ,
  60 ,
  58 ,
  56 ,
  54 ,
  52 ,
  50 ,
  48 ,
  47 ,
  45 ,
  43 ,
  41 ,
  39 ,
  37 ,
  35 ,
  33 ,
  31 ,
  29 ,
  27 ,
  25 ,
  23 ,
  21 ,
  19 ,
  17 ,
  15 ,
  13 ,
  11 ,
  10 ,
   8 ,
   6 ,
   4 ,
   2 ,
   0 ,
};

static const uint32_t lg2_tab[71] = {
          0,
   24205884,
   72626518,
  121058980,
  169503274,
  217959408,
  266427385,
  314907214,
  363398898,
  411902444,
  460417859,
  508945146,
  557484313,
  606035366,
  654598309,
  703173149,
  751759892,
  800358543,
  848969109,
  897591595,
  946226007,
  994872351,
 1043530633,
 1092200858,
 1140883032,
 1189577163,
 1238283254,
 1287001312,
 1335731344,
 1384473354,
 1433227349,
 1481993334,
 1530771317,
 1579561301,
 1628363294,
 1677177302,
 1726003329,
 1774841383,
 1823691468,
 1872553592,
 1921427759,
 1970313977,
 2019212250,
 2068122585,
 2117044987,
 2165979464,
 2190451231,
 2239403829,
 2288368515,
 2337345296,
 2386334178,
 2435335166,
 2484348266,
 2533373485,
 2582410828,
 2631460302,
 2680521912,
 2729595665,
 2778681566,
 2827779622,
 2876889838,
 2926012220,
 2975146776,
 3024293509,
 3073452428,
 3098036458,
 3147213665,
 3196403071,
 3245604684,
 3294818508,
 3344044550,
};

static const double poly_tab[] = { -0.72110585339188815190, -2.0006702299034562209 };
double fast_log2(uint32_t x)
{
  static const double INV_POW2_37 = 1.0/ ((int64_t)1 << 37);
  static const double INV_POW2_53 = 1.0/ ((int64_t)1 << 53);
  int lz = __builtin_clz(x);
  int64_t iRes = (uint64_t)(31 - lz) << 37; // integer part of result, scaled by 2^37
  x <<= lz;
  uint32_t idx1 = (x >> 25) & 63;
  uint64_t x2   = (uint64_t)x * (mx1_tab[idx1] + 257);
  uint32_t idx2 = (x2 >> 28) & 127;
  int64_t  x3   = x2 * (mx2_tab[idx2] + 8055);

  iRes += (uint64_t)idx1 << 31;
  iRes -= (uint64_t)1    << 31;
  iRes += (uint64_t)lg1_tab[idx1] << 2;
  iRes += lg2_tab[idx2];
  double res = iRes*INV_POW2_37;
  double dx  = (x3 - ((int64_t)1 << 53))*INV_POW2_53;
  double polyval0 = dx * poly_tab[0];
  double polyval1 = (dx + poly_tab[1]);
  double polyval = polyval0*polyval1;
  res += polyval;

  return res;
}