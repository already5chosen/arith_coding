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
static void tst11(const uint8_t* src, int srclen);
static void tst12(const uint8_t* src, int srclen);
static void tst13(const uint8_t* src, int srclen);
static void tst14(const uint8_t* src, int srclen);
static void tst15(const uint8_t* src, int srclen);
static void tst16(const uint8_t* src, int srclen);
static void tst17(const uint8_t* src, int srclen);
static void tst18(const uint8_t* src, int srclen);
static void tst19(const uint8_t* src, int srclen);
static void tst20(const uint8_t* src, int srclen);
static void tst21(const uint8_t* src, int srclen);
static void tst22(const uint8_t* src, int srclen);
static void tst23(const uint8_t* src, int srclen);
static void tst24(const uint8_t* src, int srclen);
static void tst25(const uint8_t* src, int srclen);
static void tst26(const uint8_t* src, int srclen);
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
      const size_t TILE_SIZE = 1024*1024*4;
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
            tst11(&bwtOut.at(0), tilelen);
            tst2(&bwtOut.at(0), tilelen);
            tst21(&bwtOut.at(0), tilelen);
            tst22(&bwtOut.at(0), tilelen);
            tst3(&bwtOut.at(0), tilelen);
            tst4(&bwtOut.at(0), tilelen);
            tst5(&bwtOut.at(0), tilelen);
            tst12(&bwtOut.at(0), tilelen);
            tst13(&bwtOut.at(0), tilelen);
            tst14(&bwtOut.at(0), tilelen);
            tst15(&bwtOut.at(0), tilelen);
            tst16(&bwtOut.at(0), tilelen);
            tst17(&bwtOut.at(0), tilelen);
            tst6(&bwtOut.at(0), tilelen);
            tst7(&bwtOut.at(0), tilelen);
            tst8(&bwtOut.at(0), tilelen);
            tst18(&bwtOut.at(0), tilelen);
            tst25(&bwtOut.at(0), tilelen);
            tst26(&bwtOut.at(0), tilelen);
            tst23(&bwtOut.at(0), tilelen);
            tst19(&bwtOut.at(0), tilelen);
            tst24(&bwtOut.at(0), tilelen);
            tst20(&bwtOut.at(0), tilelen);
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

static double calc_entr(int* h, int hlen, int* pTot = 0)
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
  if (tot > 0)
    res += log2(tot)*tot;
  if (pTot)
    *pTot = tot;
  return res;
}

static int nbits_base3(unsigned x)
{
  int n = 0;
  do {
    ++n;
    x /= 3;
  } while (x > 0);
  ++n;
  return n*2;
}

static void tst1(const uint8_t* src, int srclen)
{
  int h[256] = {0};
  for (int i = 0; i < srclen; ++i)
    ++h[src[i]];
  double e0 = calc_entr(h, 256)/8;
  double e1 = 0/8;
  double e2 = 0/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - plain\n", e0, e1, e2, e0+e1+e2);
}

static void tst2(const uint8_t* src, int srclen)
{
  int ha[256][16] = {{0}};
  int hb[256] = {0};
  int long_cnt_bits = 0;
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    int i0 = i;
    do ++i; while (src[i] == c);
    int rl = i - i0 - 1;
    if (rl < 15) {
      ++ha[c][rl];
    } else {
      ++ha[c][15];
      if (rl >= 255) {
        long_cnt_bits += nbits_base3(rl - 255);
        rl = 255;
      }
      ++hb[rl];
    }
  }
  double mne = srclen * 2;
  for (int thr = 0; thr < 15; ++thr) {
    int hc[256];
    int hr[256];
    for (int r = 0; r < 256; ++r)
      hr[r] = hb[r];
    for (int c = 0; c < 256; ++c) {
      int cnt = ha[c][15]; // long runs are counted as one character
      for (int r = 0; r < thr; ++r) {
        cnt += ha[c][r]*(r+1); // below thr: count each individual character
      }

      cnt     += ha[c][thr]*(thr+1); // at thr: count each individual character
      hr[thr] += ha[c][thr];         // at thr: account for run in hr

      for (int r = thr + 1; r < 15; ++r) {
        cnt   += ha[c][r]*(thr+1); // above thr: count first thr+1 charcters
        hr[r] += ha[c][r];         // above thr: account for run in hr
      }
      hc[c] = cnt;
    }

    double e0 = calc_entr(hc, 256)/8;
    double e1 = calc_entr(hr, 256)/8;
    double e2 = double(long_cnt_bits)/8;
    double ee = e0 + e1 + e2;
    if (ee < mne) {
      printf("%11.3f + %11.3f + %8.3f = %11.3f - simple rle. Thr %d\n", e0, e1, e2, e0+e1+e2, thr);
      mne = ee;
    }
  }
}

static void tst21(const uint8_t* src, int srclen)
{
  int h[257] = {0};
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    int i0 = i;
    do ++i; while (src[i] == c);
    int rl = i - i0;
    int cnt0 = h[0];
    int cnt1 = h[c+1];
    do {
      if (rl & 1)
        ++cnt1;
      else
        ++cnt0;
      rl >>= 1;
    } while (rl != 0);
    h[0]   = cnt0;
    h[c+1] = cnt1;
  }

  double e0 = calc_entr(h, 257)/8;
  double e1 = 0/8;
  double e2 = 0/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - rle. Similar to RUNA/RUNB\n", e0, e1, e2, e0+e1+e2);
}

static void tst22(const uint8_t* src, int srclen)
{
  int h[258] = {0};
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++h[c+2];
    int i0 = i;
    do ++i; while (src[i] == c);
    int rl = i - i0 - 1;
    while (rl != 0) {
      ++h[rl & 1];
      rl >>= 1;
    }
  }

  double e0 = calc_entr(h, 258)/8;
  double e1 = 0/8;
  double e2 = 0/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - rle. Intermixed RUNA/RUNB. %d %d\n", e0, e1, e2, e0+e1+e2, h[0], h[1]);
}

static void tst11(const uint8_t* src, int srclen)
{
  int hc[256] = {0};
  int hr[2]   = {0};
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    int i0 = i;
    do ++i; while (src[i] == c);
    int rl = i - i0 - 1;
    ++hc[c];
    ++hr[0];
    hr[1] += rl;
  }
  double e0 = calc_entr(hc, 256)/8;
  double e1 = calc_entr(hr, 2)/8;
  double e2 = 0 / 8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - 1-bit rle.\n", e0, e1, e2, e0+e1+e2);
}

static void tst3(const uint8_t* src, int srclen)
{
  int h[512] = {0};
  int long_cnt_bits = 0;
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++h[c];
    int i0 = i;
    do ++i; while (src[i] == c);
    int rl = i - i0;
    if (rl > 1) {
      int rl_code = 254 + rl;
      if (rl_code >= 511) {
        long_cnt_bits += nbits_base3(rl_code - 511);
        rl_code = 511;
      }
      ++h[rl_code];
    }
  }
  double e0 = calc_entr(h, 512)/8;
  double e1 = 0;
  double e2 = double(long_cnt_bits)/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - intermixed rle\n", e0, e1, e2, e0+e1+e2);
}


static void tst4(const uint8_t* src, int srclen)
{
  int h[2][257] = {{0}};
  int long_cnt_bits = 0;
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++h[0][c];
    int i0 = i;
    do ++i; while (src[i] == c);
    int rl = i - i0;
    if (rl > 1) {
      ++h[0][256];
      int rl_code = rl-2;
      if (rl_code >= 256) {
        long_cnt_bits += nbits_base3(rl_code - 256);
        rl_code = 256;
      }
      ++h[1][rl_code];
    }
  }
  double e0 = calc_entr(h[0], 257)/8;
  double e1 = calc_entr(h[1], 257)/8;
  double e2 = double(long_cnt_bits)/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - char&rle_marker+rle\n", e0, e1, e2, e0+e1+e2);
}

static void tst5(const uint8_t* src, int srclen)
{
  int h[2][256] = {{0}};
  int long_cnt_bits = 0;
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++h[0][c];
    int i0 = i;
    do ++i; while (src[i] == c);
    int rl = i - i0;
    if (rl > 1) {
      ++h[0][c];
      int rl_code = rl-2;
      if (rl_code >= 255) {
        long_cnt_bits += nbits_base3(rl_code - 255);
        rl_code = 255;
      }
      ++h[1][rl_code];
    }
  }
  double e0 = calc_entr(h[0], 256)/8;
  double e1 = calc_entr(h[1], 256)/8;
  double e2 = double(long_cnt_bits)/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - char_rep_as_marker+rle\n", e0, e1, e2, e0+e1+e2);
}

static void tst12(const uint8_t* src, int srclen)
{
  int h[258] = {0};
  int z0 = 0, z1 = 1;
  for (int i = 0; i < srclen;  ++i) {
    int c0 = src[i];
    int c = 256;
    if (c0 != z0) {
      c  = (c0 == z1) ? 257 : c0;
      z1 = z0;
    }
    z0 = c0;

    ++h[c];
  }
  double e0 = calc_entr(h, 258)/8;
  double e1 = 0/8;
  double e2 = 0/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf2 plain\n", e0, e1, e2, e0+e1+e2);
}

static void tst13(const uint8_t* src, int srclen)
{
  int hc[258] = {0};
  int hr[256] = {0};
  int long_cnt_bits = 0;
  int z0 = 0, z1 = 1;
  int rl = 0;
  for (int i = 0; i < srclen;  ++i) {
    int c0 = src[i];
    int c = 256;
    if (c0 != z0) {
      c  = (c0 == z1) ? 257 : c0;
      z1 = z0;
    }
    z0 = c0;

    if (c != 256) {
      if (rl != 0) {
        --rl;
        if (rl >= 255) {
          long_cnt_bits += nbits_base3(rl - 255);
          rl = 255;
        }
        ++hr[rl];
      }
      rl = 0;
    }
    if (rl == 0)
      ++hc[c];
    if (c == 256)
      ++rl;
  }
  if (rl != 0) {
    --rl;
    if (rl >= 255) {
      long_cnt_bits += nbits_base3(rl - 255);
      rl = 255;
    }
    ++hr[rl];
    rl = 0;
  }

  double e0 = calc_entr(hc, 258)/8;
  double e1 = calc_entr(hr, 256)/8;
  double e2 = double(long_cnt_bits)/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf2 zero rle %d\n", e0, e1, e2, e0+e1+e2, hc[256]);
}

static void tst14(const uint8_t* src, int srclen)
{
  int h[259] = {0};
  int z0 = 0, z1 = 1, z2 = 2;
  for (int i = 0; i < srclen;  ++i) {
    int c0 = src[i];
    int c = 256;
    if (c0 != z0) {
      c = 257;
      if (c0 != z1) {
        c = 258;
        if (c0 != z2) {
          c = c0;
        }
        z2 = z1;
      }
      z1 = z0;
    }
    z0 = c0;

    ++h[c];
  }
  double e0 = calc_entr(h, 259)/8;
  double e1 = 0/8;
  double e2 = 0/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf3 plain\n", e0, e1, e2, e0+e1+e2);
}

static void tst15(const uint8_t* src, int srclen)
{
  int hc[259] = {0};
  int hr[256] = {0};
  int long_cnt_bits = 0;
  int z0 = 0, z1 = 1, z2 = 2;
  int rl = 0;
  for (int i = 0; i < srclen;  ++i) {
    int c0 = src[i];
    int c = 256;
    if (c0 != z0) {
      c = 257;
      if (c0 != z1) {
        c = 258;
        if (c0 != z2) {
          c = c0;
        }
        z2 = z1;
      }
      z1 = z0;
    }
    z0 = c0;

    if (c != 256) {
      if (rl != 0) {
        --rl;
        if (rl >= 255) {
          long_cnt_bits += nbits_base3(rl - 255);
          rl = 255;
        }
        ++hr[rl];
      }
      rl = 0;
    }
    if (rl == 0)
      ++hc[c];
    if (c == 256)
      ++rl;
  }
  if (rl != 0) {
    --rl;
    if (rl >= 255) {
      long_cnt_bits += nbits_base3(rl - 255);
      rl = 255;
    }
    ++hr[rl];
    rl = 0;
  }

  double e0 = calc_entr(hc, 259)/8;
  double e1 = calc_entr(hr, 256)/8;
  double e2 = double(long_cnt_bits)/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf3 zero rle\n", e0, e1, e2, e0+e1+e2);
}

static void tst16(const uint8_t* src, int srclen)
{
  int h[264] = {0};
  uint64_t z = 0x0706050403020100ull;
  for (int i = 0; i < srclen;  ++i) {
    int c0 = src[i];
    uint64_t s = z;
    uint64_t msk = uint64_t(-1);
    int c = c0;
    for (int k = 256; msk != 0; ++k) {
      msk <<= 8;
      if ((s & 255)==c0) {
        c = k;
        break;
      }
      s >>= 8;
    }
    z = (z & msk) | ((z << 8) & ~msk) | c0;

    ++h[c];
  }
  double e0 = calc_entr(h, 264)/8;
  double e1 = 0/8;
  double e2 = 0/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf8 plain\n", e0, e1, e2, e0+e1+e2);
}

static void tst17(const uint8_t* src, int srclen)
{
  int hc[264] = {0};
  int hr[256] = {0};
  int long_cnt_bits = 0;
  uint64_t z = 0x0706050403020100ull;
  int rl = 0;
  for (int i = 0; i < srclen;  ++i) {
    int c0 = src[i];
    uint64_t s = z;
    uint64_t msk = uint64_t(-1);
    int c = c0;
    for (int k = 256; msk != 0; ++k) {
      msk <<= 8;
      if ((s & 255)==c0) {
        c = k;
        break;
      }
      s >>= 8;
    }
    z = (z & msk) | ((z << 8) & ~msk) | c0;

    if (c != 256) {
      if (rl != 0) {
        --rl;
        if (rl >= 255) {
          long_cnt_bits += nbits_base3(rl - 255);
          rl = 255;
        }
        ++hr[rl];
      }
      rl = 0;
    }
    if (rl == 0)
      ++hc[c];
    if (c == 256)
      ++rl;
  }
  if (rl != 0) {
    --rl;
    if (rl >= 255) {
      long_cnt_bits += nbits_base3(rl - 255);
      rl = 255;
    }
    ++hr[rl];
    rl = 0;
  }

  double e0 = calc_entr(hc, 264)/8;
  double e1 = calc_entr(hr, 256)/8;
  double e2 = double(long_cnt_bits)/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf8 zero rle\n", e0, e1, e2, e0+e1+e2);
}

static uint8_t* mtf(const uint8_t* src, int srclen)
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
  uint8_t* src = mtf(src0, srclen);
  int h[256] = {0};
  for (int i = 0; i < srclen; ++i)
    ++h[src[i]];
  delete [] src;
  double e0 = calc_entr(h, 256)/8;
  double e1 = 0/8;
  double e2 = 0/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf plain\n", e0, e1, e2, e0+e1+e2);
}

static void tst7(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtf(src0, srclen);
  int h[2][256] = {{0}};
  int long_cnt_bits = 0;
  int x1 = 0, x2 = 0, x3 = 0;
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++h[0][c];
    int i0 = i;
    do ++i; while (src[i] == c);
    int rl = i - i0 - 1;
    if (rl >= 255) {
      long_cnt_bits += nbits_base3(rl - 255);
      rl = 255;
    }
    ++h[1][rl];
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
  double e2 = double(long_cnt_bits)/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf simple rle {%d %d %d}\n", e0, e1, e2, e0+e1+e2, x1, x2, x3);
}

static void tst8(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtf(src0, srclen);

  double mne = srclen * 2;
  for (int rlMax = 255; rlMax >= 15; rlMax -= 16) {
    int h[2][256] = {{0}};
    int long_cnt_bits = 0;
    for (int i = 0; i < srclen; ) {
      int c = src[i];
      ++h[0][c];
      ++i;
      if (c == 0) {
        int i0 = i;
        while (src[i] == 0) ++i;
        int rl = i - i0;
        if (rl >= rlMax) {
          long_cnt_bits += nbits_base3(rl - rlMax);
          rl = rlMax;
        }
        ++h[1][rl];
      }
    }

    for (int rleThr = 0; rleThr < 8; ++rleThr) {
      if (rleThr > 0) {
        h[1][rleThr-1] = 0;
        for (int r = rleThr; r < 256; ++r)
          h[0][0] += h[1][r];
      }
      double e0 = calc_entr(h[0], 256)/8;
      double e1 = calc_entr(h[1], 256)/8;
      double e2 = double(long_cnt_bits)/8;
      double ee = e0 + e1 + e2;
      if (ee < mne) {
        printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf zeros rle. rlMax = %d. Thr=%d\n", e0, e1, e2, e0+e1+e2, rlMax, rleThr);
        mne = ee;
      }
    }
  }

  delete [] src;
}

static void tst18(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtf(src0, srclen);

  int h[257] = {0};
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++i;
    if (c == 0) {
      int i0 = i;
      while (src[i] == 0) ++i;
      unsigned rl = i - i0 + 1;
      // encode to RUNA=1/RUNB=0
      do {
        ++h[rl & 1];
        rl = (rl - 1) / 2;
      } while (rl != 0);
    } else {
      ++h[c+1];
    }
  }

  int tot;
  double e0 = calc_entr(h, 257, &tot)/8;
  double e1 = 0/8;
  double e2 = 0/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf zeros rle encoded with RUNA/RUNB %.3f %.3f %.3f %.3f  %.2f\n", e0, e1, e2, e0+e1+e2
    , double(h[0])/tot
    , double(h[1])/tot
    , double(h[3])/tot
    , double(h[4])/tot
    , double(srclen)/tot
    );

  delete [] src;
}

static void tst25(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtf(src0, srclen);

  int h[2][257] = {{0}};
  int tIdx = 0;
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++i;
    if (c == 0) {
      int i0 = i;
      while (src[i] == 0) ++i;
      unsigned rl = i - i0 + 1;
      // encode to RUNA=1/RUNB=0
      do {
        ++h[tIdx][rl & 1];
        rl = (rl - 1) / 2;
        tIdx = 1;
      } while (rl != 0);
    } else {
      ++h[tIdx][c+1];
      tIdx = 0;
    }
  }

  double e0 = calc_entr(h[0], 257)/8;
  double e1 = calc_entr(h[1], 257)/8;
  double e2 = 0/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf zeros rle encoded with RUNA/RUNB. 2 tables.\n", e0, e1, e2, e0+e1+e2
    );

  delete [] src;
}

static void tst26(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtf(src0, srclen);

  int h[3][257] = {{0}};
  int tIdx = 2;
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++i;
    if (c == 0) {
      int i0 = i;
      while (src[i] == 0) ++i;
      unsigned rl = i - i0 + 1;
      // encode to RUNA=1/RUNB=0
      do {
        int idx = rl & 1;
        ++h[tIdx][idx];
        tIdx = idx;
        rl = (rl - 1) / 2;
      } while (rl != 0);
    } else {
      ++h[tIdx][c+1];
      tIdx = 2;
    }
  }

  double e0 = calc_entr(h[2], 257)/8;
  double e1 = calc_entr(h[1], 257)/8;
  double e2 = calc_entr(h[0], 257)/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf zeros rle encoded with RUNA/RUNB. 3 tables.\n", e0, e1, e2, e0+e1+e2
    );

  delete [] src;
}

static const int VAL_RANGE = 1 << 14;
static void quantize_histogram(int* q, const int* h, int len)
{
  // find sum and maximum of h
  int maxH   = 0;
  int totCnt = 0;
  for (int i = 0; i < len; ++i) {
    totCnt += h[i];
    if (maxH < h[i])
      maxH = h[i];
  }

  // 1st pass - translate characters with counts that are significantly lower than maxCnt
  unsigned thr = maxH - maxH/8;
  unsigned remCnt   = totCnt;
  unsigned remRange = VAL_RANGE;
  for (int c = 0; c < len; ++c) {
    unsigned cnt = h[c];
    if (cnt < thr) {
      unsigned range = 0;
      if (cnt != 0) {
        // calculate range from full histogram
        range = (uint64_t(cnt)*(VAL_RANGE*2) + totCnt)/(totCnt*2);
        if (range == 0)
          range = 1;
        remCnt   -= cnt;
        remRange -= range;
      }
      q[c] = range;
    }
  }
  // 2nd pass - translate characters with higher counts
  for (int c = 0; remCnt != 0; ++c) {
    unsigned cnt = h[c];
    if (cnt >= thr) {
      // Calculate range from the remaining range and count
      // It is non-ideal, but this way we distribute the worst rounding errors
      // relatively evenly among higher ranges, where it has the smallest impact
      unsigned range = (uint64_t(cnt)*(remRange*2) + remCnt)/(remCnt*2);
      // (range < VAL_RANGE) is guaranteed , because we already handled the case of repetition of the same character
      remCnt   -= cnt;
      remRange -= range;
      q[c] = range;
    }
  }
}

static double calc_quantize_entropy(const int* q, const int* h, int len)
{
  // calculate entropy after quantization
  double e = 0;
  for (int c = 0; c < 257; ++c) {
    unsigned cnt = h[c];
    if (cnt)
      e -= log2(q[c]/double(VAL_RANGE))*cnt;
  }
  return e;
}

static void tst23(const uint8_t* src0, int srclen)
{
  const int CHUNK_SZ = 8192;
  uint8_t* src = mtf(src0, srclen);

  int nChunks = 0;
  double e00 = 0;
  double e0 = 0;
  double e1 = 0;
  int q0[257] = {0};
  for (int i = 0; i < srclen; ) {
    int h[257] = {0};
    for (int nsym = 0; i < srclen && nsym < CHUNK_SZ; ++nsym) {
      int c = src[i];
      ++i;
      if (c == 0) {
        int i0 = i;
        while (src[i] == 0) ++i;
        unsigned rl = i - i0 + 1;
        // encode to RUNA=1/RUNB=0
        do {
          ++h[rl & 1];
          rl = (rl - 1) / 2;
        } while (rl != 0);
      } else {
        ++h[c+1];
      }
    }

    e00 += calc_entr(h, 257)/8;

    int q[257];
    quantize_histogram(q, h, 257);
    e0 += calc_quantize_entropy(q, h, 257)/8;

    // encode q[] - q0[] with code, similar to RUNA/RUNB
    int hd[4] = {0};
    for (int c = 0; c < 257; ++c) {
      int dq = q[c];
      if (q0[c] > 31 || q0[c] < -31)
        dq -= q0[c];
      q0[c] = q[c];
      ++hd[dq >= 0];
      dq = dq < 0 ? -dq : dq;
      dq += 1;
      while (dq != 1) {
        ++hd[(dq&1) + 2];
        dq >>= 1;
      }
    }
    e1 += calc_entr(hd, 4)/8;

    ++nChunks;
  }

  double e2 = 0/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf zeros rle encoded with RUNA/RUNB. %d * %dB chunks  %11.3f\n", e0, e1, e2, e0+e1+e2, nChunks, CHUNK_SZ, e00);

  delete [] src;
}

static void tst19(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtf(src0, srclen);

  int h[259] = {0};
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++i;
    if (c == 0) {
      int i0 = i;
      while (src[i] == 0) ++i;
      unsigned rl = i - i0 + 1;
      // encode to RUNA=1/RUNB=0
      do {
        ++h[rl & 1];
        rl = (rl - 1) / 2;
      } while (rl != 0);
    } else {
      ++h[c+3];
      int i0 = i;
      while (src[i] == c) ++i;
      unsigned rl = i - i0;
      if (rl > 0) {
        // encode to RUNC=3/RUND=2
        do {
          ++h[(rl & 1)+2];
          rl = (rl - 1) / 2;
        } while (rl != 0);
      }
    }
  }

  double e0 = calc_entr(h, 259)/8;
  double e1 = 0/8;
  double e2 = 0/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf zeros rle encoded with RUNA/RUNB. Non-zero repetitions encoded with  RUNC/RUND\n", e0, e1, e2, e0+e1+e2);

  delete [] src;
}

static void tst24(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtf(src0, srclen);

  int h[258] = {0};
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    int i0 = i;
    ++i;
    while (src[i] == c) ++i;
    unsigned rl = i - i0;
    if (c == 0) {
      // encode to RUNA=1/RUNB=0
      do {
        ++h[rl & 1];
        rl = (rl - 1) / 2;
      } while (rl != 0);
    } else {
      // encode to ci=c+2/RUNC=2
      do {
        int idx = (rl & 1) != 0 ? c + 2: 2;
        ++h[idx];
        rl /= 2;
      } while (rl != 0);
    }
  }

  double e0 = calc_entr(h, 258)/8;
  double e1 = 0/8;
  double e2 = 0/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf zeros rle encoded with RUNA/RUNB. Non-zero repetitions encoded with  c/RUNC %d/%d\n", e0, e1, e2, e0+e1+e2, h[3], h[2]);

  delete [] src;
}

static void tst20(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtf(src0, srclen);

  int h[259] = {0};
  for (int i = 0; i < srclen; ) {
    int c = src[i];
    ++i;
    if (c == 0) {
      int i0 = i;
      while (src[i] == 0) ++i;
      unsigned rl = i - i0 + 1;
      // encode to RUNA=1/RUNB=0
      do {
        ++h[rl & 1];
        rl = (rl - 1) / 2;
      } while (rl != 0);
    } else {
      ++h[c+3];
      int i0 = i;
      while (src[i] == c) ++i;
      unsigned rl = i - i0;
      if (rl > 0) {
        // encode first bit to RUNC=3/RUND=2
        ++h[(rl & 1)+2];
        rl = (rl - 1) / 2;
        // encode remaining bits to RUNA=1/RUNB=0
        while (rl != 0) {
          ++h[rl & 1];
          rl = (rl - 1) / 2;
        }
      }
    }
  }

  double e0 = calc_entr(h, 259)/8;
  double e1 = 0/8;
  double e2 = 0/8;
  printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf zeros rle encoded with RUNA/RUNB. Non-zero repetitions encoded with  RUNC/RUND/RUNA/RUNB\n", e0, e1, e2, e0+e1+e2);

  delete [] src;
}

static void tst9(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtf(src0, srclen);
  src[srclen] = 255;

  double mne = srclen * 2;
  for (int xrleThr_i = 0; xrleThr_i < 16; ++xrleThr_i) {
    int xrleThr = xrleThr_i==0 ? 15 : xrleThr_i - 1;
    int h[2][512] = {{0}};
    int long_cnt_bits = 0;
    for (int i = 0; i < srclen; ) {
      int c = src[i];
      int k = i + 1;
      while (src[k] == c) ++k;
      int rl = k - i - 1;
      i = k;

      ++h[0][c];
      if (c == 0) {
        if (rl >= 255) {
          long_cnt_bits += nbits_base3(rl - 255);
          rl = 255;
        }
        ++h[1][rl];
      } else if (rl > xrleThr) {
        if (rl >= 255) {
          long_cnt_bits += nbits_base3(rl - 255);
          rl = 255;
        }
        ++h[0][rl+256];
      } else {
        h[0][c] += rl;
      }
    }

    for (int zrleThr = 0; zrleThr < 8; ++zrleThr) {
      if (zrleThr > 0) {
        h[1][zrleThr-1] = 0;
        for (int r = zrleThr; r < 256; ++r)
          h[0][0] += h[1][r];
      }

      double e1 = calc_entr(h[1], 256)/8;
      double e0 = calc_entr(h[0], 512)/8;
      double e2 = double(long_cnt_bits)/8;
      double ee = e0 + e1 + e2;
      if (ee < mne) {
        printf("%11.3f + %11.3f + %8.3f = %11.3f - mtf zeros rle, others intermixed. Zero RLE Thr=%d %d\n", e0, e1, e2, e0+e1+e2, xrleThr, zrleThr);
        mne = ee;
      }
    }
  }

  delete [] src;
}

static void tst10(const uint8_t* src0, int srclen)
{
  uint8_t* src = mtf(src0, srclen);

  double mne = srclen * 2;
  // int mnx = 0, mny = 0;
  for (int xrleThr_i = 0; xrleThr_i < 16; ++xrleThr_i) {
    int xrleThr = xrleThr_i==0 ? 15 : xrleThr_i - 1;
    int h[3][256] = {{0}};
    int long_cnt_bits = 0;
    for (int i = 0; i < srclen; ) {
      int c = src[i];
      int k = i + 1;
      while (src[k] == c) ++k;
      int rl = k - i - 1;
      i = k;

      ++h[0][c];
      if (c == 0) {
        if (rl >= 255) {
          long_cnt_bits += nbits_base3(rl - 255);
          rl = 255;
        }
        ++h[1][rl];
      } else if (rl >= xrleThr) {
        if (rl >= 255) {
          long_cnt_bits += nbits_base3(rl - 255);
          rl = 255;
        }
        ++h[2][rl];
        h[0][c] += xrleThr;
      } else {
        h[0][c] += rl;
      }
    }

    double e2 = double(long_cnt_bits)/8;
    double e3 = calc_entr(h[2], 256)/8;
    for (int zrleThr = 0; zrleThr < 8; ++zrleThr) {
      if (zrleThr > 0) {
        h[1][zrleThr-1] = 0;
        for (int r = zrleThr; r < 256; ++r)
          h[0][0] += h[1][r];
      }
      double e0 = calc_entr(h[0], 256)/8;
      double e1 = calc_entr(h[1], 256)/8;
      double ee = e0 + e1 + e2 + e3;
     if (ee < mne) {
        printf("%11.3f + %11.3f + %8.3f + %11.3f = %11.3f - mtf zeros/nzeros rle. Thrs %d %d\n", e0, e1, e2, e3, e0+e1+e2+e3, xrleThr, zrleThr);
       mne = ee;
     }
    }
  }

  delete [] src;
}


// static void tst9(const uint8_t* src0, int srclen)
// {
  // uint8_t* src = mtf(src0, srclen);
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
  // printf("%11.3f + %11.3f + %11.3f = %11.3f - mtf zeros & one rle\n", e0, e1, e2, e0+e1+e2);
// }
