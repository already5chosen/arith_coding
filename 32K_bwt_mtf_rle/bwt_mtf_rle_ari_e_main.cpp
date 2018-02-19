#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>
#include <x86intrin.h>

#include "bwt_mtf_rle_e.h"
#include "arithmetic_encode.h"

static void storeAs3octets(uint8_t* dst, unsigned val);
static bool isSingleCharacter(const uint8_t* src, int len);
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
          uint8_t hdr[9];
          storeAs3octets(&hdr[0], tilelen);
          hdr[5] = 255; // default to special case
          void* pRes = 0;
          size_t hdrlen = 6;
          int ressz  = 0;
          if (!isSingleCharacter(inptile, tilelen)) {
            tmpDst.resize(tilelen+256);
            uint64_t t0 = __rdtsc();
            bwt_sort(&tmpDst.at(0), inptile, tilelen);
            bwt_mtf_rle_meta_t meta;
            ressz = bwt_reorder_mtf_rle(  // return length of destination array in octets
              &tmpDst.at(0), // both input and output
              inptile, tilelen, &meta);
            storeAs3octets(&hdr[3], ressz);
            storeAs3octets(&hdr[6], meta.bwtPrimaryIndex);
            uint64_t t1 = __rdtsc();
            if (vFlag)
              printf("%7d->%7d chars. %10.0f clocks. %6.1f clocks/char\n"
                ,int(tilelen)
                ,ressz
                ,double(t1-t0)
                ,double(t1-t0)/tilelen
                );
            hdrlen = 9;
            pRes = &tmpDst.at(0);
          } else {
            // input consists of repetition of the same character
            hdr[3] = 1;
            hdr[4] = inptile[0];
          }
          // write result to file
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

static void storeAs3octets(uint8_t* dst, unsigned val)
{
  dst[0] = (val >> 0) % 256;
  dst[1] = (val >> 8) % 256;
  dst[2] = (val >>16) % 256;
}

static bool isSingleCharacter(const uint8_t* src, int len) // len  > 0
{
  uint8_t c0 = src[0];
  for (int i = 1; i < len; ++i) {
    if (src[i] != c0)
      return false;
  }
  return true; // repetition of single character
}