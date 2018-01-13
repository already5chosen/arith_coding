#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>
#include <x86intrin.h>

#include "bwt_e.h"
#include "arithmetic_encode.h"

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
            bwtOut.resize(tilelen);
            int primary_i = bwt(&bwtOut.at(0), inptile, tilelen);
            hdr[6+0] = (primary_i >> 0) % 256;
            hdr[6+1] = (primary_i >> 8) % 256;
            hdr[6+2] = (primary_i >>16) % 256;
            t1 = __rdtsc();
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


