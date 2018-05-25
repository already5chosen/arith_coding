#include <stdint.h>
#include <stdio.h>
#include <math.h>

static const uint8_t tab1i[64] = {
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

static const uint32_t tab1[64] = {
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

// Polynomial [4.6920775424952743116e-1 -7.2123175429236729972e-1 1.4426947140791750446 0]
static const double poly0_tab[] = { 0.46920775424952743116 };
static const double poly2_tab[] = { -1.5371266731214592305, 3.0747461034327696616 };
double fast_log2(uint32_t x)
{
  static const double INV_POW2_35 = 1.0/ ((int64_t)1 << 35);
  static const double INV_POW2_40 = 1.0/ ((int64_t)1 << 40);
  int lz = __builtin_clz(x);
  int64_t iRes = (uint64_t)(31 - lz) << 35; // integer part of result, scaled by 2^35
  x <<= lz;
  uint32_t idx1 = (x >> 25) & 63;
  int64_t x2 = (uint64_t)x * (tab1i[idx1] + 257);

  iRes += (uint64_t)idx1 << 29;
  iRes -= (uint64_t)1    << 29;
  iRes += tab1[idx1];
  double res = iRes*INV_POW2_35;
  double dx  = (x2 - ((int64_t)1 << 40))*INV_POW2_40;
  double polyval0 = dx * poly0_tab[0];
  double polyval2 = (dx + poly2_tab[0]) * dx + poly2_tab[1];
  double polyval = polyval0*polyval2;
  res += polyval;

  return res;
}