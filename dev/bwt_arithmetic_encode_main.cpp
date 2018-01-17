#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <x86intrin.h>

#include "bwt_e.h"
#include "arithmetic_encode.h"

static void tst1(const uint8_t* src, int srclen);
static void tst2(const uint8_t* src, int srclen);
static void tst3(const uint8_t* src, int srclen);
static void tst4(const uint8_t* src, int srclen);
static void tst5(const uint8_t* src, int srclen);
static void tst6(const uint8_t* src, int srclen);
static void tst7(const uint8_t* src, int srclen);
static void tst8(const uint8_t* src, int srclen);
static void tst9(const uint8_t* src, int srclen);
static void tst10(const uint8_t* src, int srclen);
int main(int argz, char** argv)
{
  if (argz < 3) {
    printf(
      "Usage:\n"
      "bwt_ari_encode input-file output-file [-v]\n"
      );
    return 1;
  }
  char* inpfilename = argv[1];
  char* outfilename = argv[2];
  bool vFlag = (argz > 3) && (strcmp("-v", argv[3])==0);

  int ret = 1;
  FILE* fpinp = fopen(inpfilename, "rb");
  if (fpinp) {
    FILE* fpout = fopen(outfilename, "wb");
    if (fpout) {
      const size_t TILE_SIZE = 1024*1024;
      uint8_t* inptile = new uint8_t[TILE_SIZE];
      std::vector<uint8_t> bwtOut;
      std::vector<uint8_t> dst;
      ret = 0;
      bool done = false;
      while (!done) {
        size_t tilelen = fread(inptile, 1, TILE_SIZE, fpinp);
        if (tilelen < TILE_SIZE) {
          if (!feof(fpinp)) {
            perror(inpfilename);
            ret = 1;
            break;
          }
          done = true;
        }
        if (tilelen > 0) {
          uint8_t hdr[9];
          hdr[0] = uint8_t(tilelen >> 0);
          hdr[1] = uint8_t(tilelen >> 8);
          hdr[2] = uint8_t(tilelen >> 16);
          uint8_t* pRes = 0;
          double info[8]={0};
          int ressz = -1;

          uint64_t t0 = __rdtsc();
          uint64_t t1 = t0;
          unsigned srcHistograms[2][260];
          arithmetic_encode_calc_histograms(srcHistograms, inptile, tilelen);
          double entropy[3];
          entropy[0] = arithmetic_encode_calc_entropy(srcHistograms[0], tilelen);
          int srcI = 0;
          if (entropy[0] != 0) {
            // input does not consists of repetition of the same character
            entropy[1] = arithmetic_encode_calc_entropy(srcHistograms[1], tilelen);
            bwtOut.resize(tilelen+1);
            int primary_i = bwt(&bwtOut.at(0), inptile, tilelen);
            hdr[6+0] = (primary_i >> 0) % 256;
            hdr[6+1] = (primary_i >> 8) % 256;
            hdr[6+2] = (primary_i >>16) % 256;
            t1 = __rdtsc();
            bwtOut[tilelen] = bwtOut[tilelen-1] + 1;
            tst1(&bwtOut.at(0), tilelen);
            tst2(&bwtOut.at(0), tilelen);
            tst3(&bwtOut.at(0), tilelen);
            tst4(&bwtOut.at(0), tilelen);
            tst5(&bwtOut.at(0), tilelen);
            tst6(&bwtOut.at(0), tilelen);
            tst7(&bwtOut.at(0), tilelen);
            tst8(&bwtOut.at(0), tilelen);
            tst9(&bwtOut.at(0), tilelen);
            tst10(&bwtOut.at(0), tilelen);
            unsigned bwtHistograms[2][260];
            arithmetic_encode_calc_histograms(bwtHistograms, &bwtOut.at(0), tilelen);
            entropy[2] = arithmetic_encode_calc_entropy(bwtHistograms[1], tilelen);

            double minSrcEntropy = entropy[0];
            if (entropy[1] < minSrcEntropy) {
              srcI = 1;
              minSrcEntropy = entropy[1];
            }
            double adjBwtEntropy = entropy[2]*1.0025+24;
            if (adjBwtEntropy < minSrcEntropy) {
              minSrcEntropy = adjBwtEntropy;
              srcI = 2;
            }

            dst.clear();
            info[0] = entropy[srcI];
            if (srcI == 0) {
              ressz = arithmetic_encode(&dst, inptile, tilelen, srcHistograms[0], false, vFlag ? info : 0);
            } else if (srcI == 1) {
              ressz = arithmetic_encode(&dst, inptile, tilelen, srcHistograms[1], true,  vFlag ? info : 0);
            } else {
              ressz = arithmetic_encode(&dst, &bwtOut.at(0), tilelen, bwtHistograms[1], true,  vFlag ? info : 0);
              if (ressz+3 >= tilelen)
                ressz = 0;
            }
          }
          uint64_t t2 = __rdtsc();

          if (vFlag) {
            printf("%7u -> %7u. %c Model %7.3f. Coded %10.0f. Entropy %11.3f (%11.3f). %10.0f clocks. %6.1f + %4.1f = %6.1f clocks/char\n"
              ,unsigned(tilelen)
              ,ressz < 0 ? 0 : (ressz == 0 ? unsigned(tilelen) : unsigned(ressz))
              , ressz > 0 ? 'A' + srcI : (ressz == 0 ? 'N' : 'R')
              ,info[1]/8
              ,info[2]/8
              ,info[0]/8
              ,info[3]/8
              ,double(t1-t0)
              ,double(t1-t0)/tilelen
              ,double(t2-t1)/tilelen
              ,double(t2-t0)/tilelen
              );
          }

          int hdrlen = 6;
          if (ressz > 0) {
            // normal compression
            hdr[3] = uint8_t(ressz >> 0);
            hdr[4] = uint8_t(ressz >> 8);
            hdr[5] = uint8_t(ressz >> 16) | (srcI << 6);
            pRes = &dst.at(0);
            if (srcI == 2)
              hdrlen += 3;
          } else {
            // special cases
            if (ressz == 0) {
              // not compressible
              hdr[3] = 0;
              hdr[4] = 0;
              pRes   = inptile;
              ressz  = tilelen;
            } else {
              // input consists of repetition of the same character
              hdr[3] = 1;
              hdr[4] = inptile[0];
            }
            hdr[5] = 255;
          }

          size_t wrlen = fwrite(hdr, 1, hdrlen, fpout);
          if (wrlen != hdrlen) {
            perror(outfilename);
            ret = 1;
            break;
          }
          if (pRes) {
            wrlen = fwrite(pRes, 1, ressz, fpout);
            if (wrlen != size_t(ressz)) {
              perror(outfilename);
              ret = 1;
              break;
            }
          }
        }
      }
      delete [] inptile;
      fclose(fpout);
      if (ret != 0)
        remove(outfilename);
    } else {
      perror(outfilename);
    }
    fclose(fpinp);
  } else {
    perror(inpfilename);
  }
  return ret;
}

static double calc_entr(int* h, int hlen)
{
  int tot = 0;
  double res = 0;
  for (int i = 0; i < hlen; ++i) {
    int v = h[i];
    if (v > 0) {
      res -= log2(v)*v;
      tot += v;
    }
  }
  res += log2(tot)*tot;
  return res;
}

static void tst1(const uint8_t* src, int srclen)
{
  int h[256] = {0};
  for (int i = 0; i < srclen; ++i)
    ++h[src[i]];
  double e0 = calc_entr(h, 256)/8;
  double e1 = 0/8;
  printf("%11.3f + %11.3f = %11.3f - plain\n", e0, e1, e0+e1);
}

static void tst2(const uint8_t* src, int srclen)
{
  int h[2][256] = {{0}};
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++h[0][c];
    int i0 = i;
    do ++i; while (src[i] == c);
    int rl = i - i0 - 1;
    ++h[1][rl<256 ? rl : 255];
  }
  double e0 = calc_entr(h[0], 256)/8;
  double e1 = calc_entr(h[1], 256)/8;
  printf("%11.3f + %11.3f = %11.3f - simple rle\n", e0, e1, e0+e1);
}

static void tst3(const uint8_t* src, int srclen)
{
  int h[512] = {0};
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++h[c];
    int i0 = i;
    do ++i; while (src[i] == c);
    int rl = i - i0;
    if (rl > 1) {
      int rl_code = 254 + rl;
      ++h[rl_code < 512 ? rl_code : 511];
    }
  }
  double e0 = calc_entr(h, 512)/8;
  double e1 = 0;
  printf("%11.3f + %11.3f = %11.3f - intermixed rle\n", e0, e1, e0+e1);
}


static void tst4(const uint8_t* src, int srclen)
{
  int h[2][257] = {{0}};
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++h[0][c];
    int i0 = i;
    do ++i; while (src[i] == c);
    int rl = i - i0;
    if (rl > 1) {
      ++h[0][256];
      int rl_code = rl-2;
      ++h[1][rl_code<257 ? rl_code : 256];
    }
  }
  double e0 = calc_entr(h[0], 257)/8;
  double e1 = calc_entr(h[1], 257)/8;
  printf("%11.3f + %11.3f = %11.3f - char&rle_marker+rle\n", e0, e1, e0+e1);
}
static void tst5(const uint8_t* src, int srclen)
{
  int h[2][256] = {{0}};
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++h[0][c];
    int i0 = i;
    do ++i; while (src[i] == c);
    int rl = i - i0;
    if (rl > 1) {
      ++h[0][c];
      int rl_code = rl-2;
      ++h[1][rl_code<256 ? rl_code : 255];
    }
  }
  double e0 = calc_entr(h[0], 256)/8;
  double e1 = calc_entr(h[1], 256)/8;
  printf("%11.3f + %11.3f = %11.3f - char_rep_as_marker+rle\n", e0, e1, e0+e1);
}

static uint8_t* mtz(const uint8_t* src, int srclen)
{
  uint8_t t[256];
  for (int i = 0; i < 256; ++i)
    t[i] = i;

  uint8_t* dst = new uint8_t[srclen+1];
  for (int i = 0; i < srclen; ++i) {
    int c = src[i];
    int v0 = c, v1, k;
    for (k = 0; c != (v1=t[k]); ++k) {
      t[k] = v0;
      v0 = v1;
    }
    t[k] = v0;
    dst[i] = k;
  }
  dst[srclen] = dst[srclen] == 1 ? 2 : 1;
  return dst;
}

static void tst6(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtz(src0, srclen);
  int h[256] = {0};
  for (int i = 0; i < srclen; ++i)
    ++h[src[i]];
  delete [] src;
  double e0 = calc_entr(h, 256)/8;
  double e1 = 0/8;
  printf("%11.3f + %11.3f = %11.3f - mtz plain\n", e0, e1, e0+e1);
}

static void tst7(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtz(src0, srclen);
  int h[2][256] = {{0}};
  int x1 = 0, x2 = 0, x3 = 0;
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++h[0][c];
    int i0 = i;
    do ++i; while (src[i] == c);
    int rl = i - i0 - 1;
    ++h[1][rl<256 ? rl : 255];
    if (rl > 0 && c != 0) {
      ++x1;
      if (c == 1)
        ++x2;
      if (c == 2)
        ++x3;
    }
  }
  delete [] src;
  double e0 = calc_entr(h[0], 256)/8;
  double e1 = calc_entr(h[1], 256)/8;
  printf("%11.3f + %11.3f = %11.3f - mtz simple rle {%d %d %d}\n", e0, e1, e0+e1, x1, x2, x3);
}

static void tst8(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtz(src0, srclen);
  int h[2][256] = {{0}};
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++h[0][c];
    ++i;
    if (c == 0) {
      int i0 = i;
      while (src[i] == 0) ++i;
      int rl = i - i0;
      ++h[1][rl<256 ? rl : 255];
    }
  }
  delete [] src;

  double mne = srclen * 2;
  for (int rleThr = 0; rleThr < 8; ++rleThr) {
    if (rleThr > 0) {
      h[1][rleThr] = 0;
      h[0][0] += h[1][rleThr+1]*(rleThr+1);
    }
    double e0 = calc_entr(h[0], 256)/8;
    double e1 = calc_entr(h[1], 256)/8;
    double ee = e0 + e1;
    if (ee < mne) {
      printf("%11.3f + %11.3f = %11.3f - mtz zeros rle. Thr=%d\n", e0, e1, e0+e1, rleThr);
      mne = ee;
    }
  }
}

static void tst9(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtz(src0, srclen);
  src[srclen] = 255;

  double mne = srclen * 2;
  for (int xrleThr = 0; xrleThr < 8; ++xrleThr) {
    int h[2][512] = {{0}};
    for (int i = 0; i < srclen; ) {
      int c = src[i];
      int k = i + 1;
      while (src[k] == c) ++k;
      int rl = k - i - 1;
      i = k;

      ++h[0][c];
      if (c == 0) {
        ++h[1][rl<256 ? rl : 255];
      } else if (rl > xrleThr) {
        ++h[0][rl<256 ? rl+255 : 256+255];
      } else {
        h[0][c] += rl;
      }
    }

    for (int zrleThr = 0; zrleThr < 8; ++zrleThr) {
      if (zrleThr > 0) {
        h[1][zrleThr] = 0;
        h[0][0] += h[1][zrleThr+1]*(zrleThr+1);
      }

      double e1 = calc_entr(h[1], 256)/8;
      double e0 = calc_entr(h[0], 512)/8;
      double ee = e0 + e1;
      if (ee < mne) {
        printf("%11.3f + %11.3f = %11.3f - mtz zeros rle, ones intermixed. Zero RLE Thr=%d %d\n", e0, e1, e0+e1, xrleThr, zrleThr);
        mne = ee;
      }
    }
  }

  delete [] src;
}

static void tst10(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtz(src0, srclen);

  double mne = srclen * 2;
  // int mnx = 0, mny = 0;
  for (int xrleThr = 0; xrleThr < 16; ++xrleThr) {
    int h[3][256] = {{0}};
    for (int i = 0; i < srclen; ) {
      int c = src[i];
      int k = i + 1;
      while (src[k] == c) ++k;
      int rl = k - i - 1;
      i = k;

      ++h[0][c];
      if (c == 0) {
        ++h[1][rl<256 ? rl : 255];
      } else if (rl >= xrleThr) {
        ++h[2][rl<256 ? rl : 255];
        h[0][c] += xrleThr;
      } else {
        h[0][c] += rl;
      }
    }

    double e2 = calc_entr(h[2], 256)/8;
    for (int zrleThr = 0; zrleThr < 8; ++zrleThr) {
      if (zrleThr > 0) {
        h[1][zrleThr] = 0;
        h[0][0] += h[1][zrleThr+1]*(zrleThr+1);
      }
      double e0 = calc_entr(h[0], 256)/8;
      double e1 = calc_entr(h[1], 256)/8;
      double ee = e0 + e1 + e2;
      if (ee < mne) {
        printf("%11.3f + %11.3f  + %11.3f = %11.3f - mtz zeros/nzeros rle. Thrs %d %d\n", e0, e1, e2, e0+e1+e2, xrleThr, zrleThr);
        mne = ee;
      }
    }
  }

  delete [] src;
}


// static void tst9(const uint8_t* src0, int srclen)
// {
  // uint8_t* src = mtz(src0, srclen);
  // int h[3][256] = {{0}};
  // for (int i = 0; i < srclen; ) {
    // int c = src[i];
    // ++h[0][c];
    // ++i;
    // if (c < 8) {
      // int i0 = i;
      // while (src[i] == c) ++i;
      // int rl = i - i0;
      // ++h[c==0 ? 1 : 2][rl<256 ? rl : 255];
    // }
  // }
  // delete [] src;
  // double e0 = calc_entr(h[0], 256)/8;
  // double e1 = calc_entr(h[1], 256)/8;
  // double e2 = calc_entr(h[2], 256)/8;
  // printf("%11.3f + %11.3f + %11.3f = %11.3f - mtz zeros & one rle\n", e0, e1, e2, e0+e1+e2);
// }

