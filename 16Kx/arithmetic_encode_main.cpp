#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>
#include <x86intrin.h>

#include "arithmetic_encode.h"

int main(int argz, char** argv)
{
  if (argz < 3) {
    printf(
      "Usage:\n"
      "ari_encode input-file output-file [-v]\n"
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
          dst.clear();
          double info[8];
          uint64_t t0 = __rdtsc();
          int ressz = arithmetic_encode(&dst, inptile, tilelen, vFlag ? info : 0);
          uint64_t t1 = __rdtsc();
          if (vFlag)
            printf("%7u -> %7u. Model %7.3f. Coded %10.0f. Entropy %11.3f (%11.3f). %8.0f clocks. %4.1f clocks/char\n"
              ,unsigned(tilelen)
              ,ressz < 0 ? 0 : (ressz == 0 ? unsigned(tilelen) : unsigned(ressz))
              ,info[1]/8
              ,info[2]/8
              ,info[0]/8
              ,info[3]/8
              ,double(t1-t0)
              ,double(t1-t0)/tilelen
              );
          uint8_t hdr[6];
          hdr[0] = uint8_t(tilelen >> 0);
          hdr[1] = uint8_t(tilelen >> 8);
          hdr[2] = uint8_t(tilelen >> 16);
          uint8_t* pRes = 0;
          if (ressz > 0) {
            // normal compression
            hdr[3] = uint8_t(ressz >> 0);
            hdr[4] = uint8_t(ressz >> 8);
            hdr[5] = uint8_t(ressz >> 16);
            pRes = &dst.at(0);
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
          size_t wrlen = fwrite(hdr, 1, 6, fpout);
          if (wrlen != 6) {
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


