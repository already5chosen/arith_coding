#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>
#include <x86intrin.h>

#include "bwt_mtf_rle_e.h"

static void storeAs3octets(uint8_t* dst, unsigned val);
int main(int argz, char** argv)
{
  if (argz < 3) {
    printf(
      "Usage:\n"
      "bwt_mtf_rle_e input-file output-file [-v]\n"
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
      std::vector<int32_t> tmpDst;
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
          tmpDst.resize(tilelen);
          uint64_t t0 = __rdtsc();
          int bwtPrimaryIndex, nRuns;
          int ressz = bwt_mtf_rle(  // return length of destination array in octets
            &tmpDst.at(0), // srclen*uint32_t
            &bwtPrimaryIndex,
            &nRuns,     // number of runs in the destination array
            inptile, tilelen);
          uint8_t hdr[9];
          storeAs3octets(&hdr[0], tilelen);
          storeAs3octets(&hdr[3], ressz);
          storeAs3octets(&hdr[6], bwtPrimaryIndex);
          uint64_t t1 = __rdtsc();
          if (vFlag)
            printf("%7d->%7d chars. %10.0f clocks. %6.1f clocks/char\n"
              ,int(tilelen)
              ,ressz
              ,double(t1-t0)
              ,double(t1-t0)/tilelen
              );
          size_t wrlen = fwrite(hdr, 1, 9, fpout);
          if (wrlen != 9) {
            perror(outfilename);
            ret = 1;
            break;
          }
          wrlen = fwrite(&tmpDst.at(0), 1, ressz, fpout);
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

static void storeAs3octets(uint8_t* dst, unsigned val)
{
  dst[0] = (val >> 0) % 256;
  dst[1] = (val >> 8) % 256;
  dst[2] = (val >>16) % 256;
}

