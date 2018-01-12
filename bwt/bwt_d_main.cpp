#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>
#include <x86intrin.h>

#include "bwt_d.h"

int main(int argz, char** argv)
{
  if (argz < 3) {
    printf(
      "Usage:\n"
      "bwt_d input-file output-file [-v]\n"
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
      const size_t MAX_DST_TILE_SIZE = 1024*1024;
      const size_t MAX_SRC_TILE_SIZE = MAX_DST_TILE_SIZE + 3;
      std::vector<uint8_t> src(MAX_SRC_TILE_SIZE);
      std::vector<uint8_t> dst;
      ret = 1;

      while (true) {
        size_t srctilelen = fread(&src.at(0), 1, MAX_SRC_TILE_SIZE, fpinp);
        if (srctilelen <= 3) {
          if (ferror(fpinp))
            perror(inpfilename);
          else if (srctilelen == 0)
            ret = 0; // success
          else
            fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Short tile.": "");
          break;
        }
        int32_t dsttilelen = srctilelen - 3;
        int32_t first_i = (src[dsttilelen+2]*256+src[dsttilelen+1])*256+src[dsttilelen+0];
        if (first_i >= dsttilelen) {
          fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " first_i out of range.": "");
          break;
        }

        dst.resize(dsttilelen);
        uint64_t t0 = __rdtsc();
        ibwt(&dst.at(0), &src.at(0), dsttilelen, first_i);
        uint64_t t1 = __rdtsc();

        if (vFlag)
          printf(
            "%7d chars. %8.0f clocks. %4.1f clocks/char\n"
            ,dsttilelen
            ,double(t1-t0)
            ,double(t1-t0)/dsttilelen
            );

        size_t wrlen = fwrite(&dst.at(0), 1, dsttilelen, fpout);
        if (wrlen != size_t(dsttilelen)) {
          perror(outfilename);
          break;
        }
      }

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


