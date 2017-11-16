#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>
#include <x86intrin.h>

#include "arithmetic_decode.h"

int main(int argz, char** argv)
{
  if (argz < 3) {
    printf(
      "Usage:\n"
      "ari_decode input-file output-file [-v]\n"
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
      const size_t MAX_TILE_SIZE = 1024*1024;
      const size_t MAX_CODE_SIZE = MAX_TILE_SIZE + 1024;
      std::vector<uint8_t> src;
      std::vector<uint8_t> dst;
      ret = 1;

      while (true) {
        uint8_t hdr[6];
        size_t hdrlen = fread(hdr, 1, sizeof(hdr), fpinp);
        if (hdrlen != sizeof(hdr)) {
          if (ferror(fpinp))
            perror(inpfilename);
          else if (hdrlen == 0)
            ret = 0; // success
          else
            fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Short tile header.": "");
          break;
        }

        size_t tilelen = (hdr[2]*256u+hdr[1])*256u+hdr[0];
        size_t codelen = (hdr[5]*256u+hdr[4])*256u+hdr[3];

        if (tilelen < 1 || tilelen > MAX_TILE_SIZE) {
          fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Illegal tile length.": "");
          break;
        }

        if (hdr[5] == 255) {
          // special cases
          switch (hdr[3]) {
            case 0:
              // not compressible
              codelen = tilelen;
              break;
            case 1:
              // input consists of repetition of the same character
              codelen = 0;
              break;
            default:
              fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Illegal escape sequence in section header.": "");
              goto end_loop;
          }
        } else {
          if (codelen < 1 || codelen > MAX_CODE_SIZE) {
            fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Illegal code section length.": "");
            break;
          }
        }

        if (codelen != 0) {
          src.resize(codelen);
          size_t inplen = fread(&src.at(0), 1, codelen, fpinp);
          if (inplen < codelen) {
            if (feof(fpinp))
              fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Code section is shorter than specified in tile header.": "");
            else
              perror(inpfilename);
            break;
          }
        }

        int info[8]={0};
        uint8_t *pDst = 0;
        uint64_t t0 = __rdtsc();
        if (hdr[5] == 255) {
          // special cases
          switch (hdr[3]) {
            case 0:
              // not compressible
              pDst    = &src.at(0);
              info[0] = 0;
              info[1] = codelen;
              break;
            case 1:
              // input consists of repetition of the same character
              dst.resize(tilelen);
              pDst = &dst.at(0);
              memset(pDst, hdr[4], tilelen);
              info[0] = 8;
              info[1] = 1;
              break;
            default:
              break;
          }
        } else {
          dst.resize(tilelen);
          pDst = &dst.at(0);
          int declen = arithmetic_decode(pDst, tilelen, &src.at(0), codelen, vFlag ? info : 0);
          if (size_t(declen) != tilelen) {
            fprintf(stderr, "%s: %s invalid.\n", argv[0], inpfilename);
            if (vFlag) {
              if (declen < 0)
                fprintf(stderr, "Decoder parsing failure. Error %d.\n", declen);
              else
                fprintf(stderr, "Uncompressed section is shorter than specified in tile header. %d < %d.\n", declen, int(tilelen));
            }
            break;
          }
        }
        uint64_t t1 = __rdtsc();

        if (vFlag)
          printf(
            "%7u -> %7u. Model %7.3f. Coded %7d. %8.0f clocks. %4.1f clocks/char"
            #ifdef ENABLE_PERF_COUNT
            " %6d %6d %6d"
            #endif
            "\n"
            ,unsigned(codelen)
            ,unsigned(tilelen)
            ,info[0]/8.0
            ,info[1]
            ,double(t1-t0)
            ,double(t1-t0)/tilelen
            #ifdef ENABLE_PERF_COUNT
            ,info[2]
            ,info[3]
            ,info[4]
            #endif
            );

        if (pDst) {
          size_t wrlen = fwrite(pDst, 1, tilelen, fpout);
          if (wrlen != tilelen) {
            perror(outfilename);
            break;
          }
        }
      }
      end_loop:
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


