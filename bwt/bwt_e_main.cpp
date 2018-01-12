#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>
#include <x86intrin.h>

#include "bwt_e.h"

int main(int argz, char** argv)
{
  if (argz < 3) {
    printf(
      "Usage:\n"
      "bwt_e input-file output-file [-v]\n"
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
          int ressz = tilelen+3;
          dst.resize(ressz);
          uint64_t t0 = __rdtsc();
          int primary_i = bwt(&dst.at(0), inptile, tilelen);
          dst[tilelen+0] = (primary_i >> 0) % 256;
          dst[tilelen+1] = (primary_i >> 8) % 256;
          dst[tilelen+2] = (primary_i >>16) % 256;
          uint64_t t1 = __rdtsc();
          if (vFlag)
            printf("%7d chars. %8.0f clocks. %4.1f clocks/char\n"
              ,int(tilelen)
              ,double(t1-t0)
              ,double(t1-t0)/tilelen
              );
          size_t wrlen = fwrite(&dst.at(0), 1, ressz, fpout);
          if (wrlen != size_t(ressz)) {
            perror(outfilename);
            ret = 1;
            break;
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


