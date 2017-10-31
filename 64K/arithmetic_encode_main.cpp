#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

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
          double entropy;
          arithmetic_encode(&dst, inptile, tilelen, vFlag ? &entropy : 0);
          size_t ressz = dst.size();
          if (vFlag)
            printf("%7u -> %7u. Entropy %.4f\n", unsigned(tilelen), unsigned(ressz), entropy/8);
          uint8_t hdr[6];
          hdr[0] = uint8_t(tilelen >> 0);
          hdr[1] = uint8_t(tilelen >> 8);
          hdr[2] = uint8_t(tilelen >> 16);
          hdr[3] = uint8_t(ressz >> 0);
          hdr[4] = uint8_t(ressz >> 8);
          hdr[5] = uint8_t(ressz >> 16);
          size_t wrlen = fwrite(hdr, 1, 6, fpout);
          if (wrlen != 6) {
            perror(outfilename);
            ret = 1;
            break;
          }
          wrlen = fwrite(&dst.at(0), 1, ressz, fpout);
          if (wrlen != ressz) {
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


