#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>
#include <x86intrin.h>

#include "arithmetic_decode.h"
#include "bwt_d.h"

int main(int argz, char** argv)
{
  if (argz < 3) {
    printf(
      "Usage:\n"
      "bwt_ari_decode input-file output-file [-v]\n"
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
      std::vector<uint8_t> ibwtInp;
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
        size_t codelen = ((hdr[5] & 63)*256u+hdr[4])*256u+hdr[3];

        if (tilelen < 1 || tilelen > MAX_TILE_SIZE) {
          fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Illegal tile length.": "");
          break;
        }

        int encType = 0;
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
          encType = (hdr[5] >> 6) & 3;
          if (encType > 2) {
            fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Illegal encode type.": "");
            break;
          }
        }

        if (codelen != 0) {
          size_t codelenEx = (encType == 2) ? codelen + 3 : codelen;
          src.resize(codelenEx);
          size_t inplen = fread(&src.at(0), 1, codelenEx, fpinp);
          if (inplen < codelenEx) {
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
        uint64_t t1 = t0;
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
          int declen = 0;
          int32_t primary_i = 0;
          if (encType == 2) {
            // BWT+arithmetic code with repetitions
            primary_i = (src[2]*256+src[1])*256+src[0];
            if (primary_i >= tilelen) {
              fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " primary_i out of range.": "");
              break;
            }
            ibwtInp.resize(tilelen);
            declen = arithmetic_decode(&ibwtInp.at(0), tilelen, &src.at(3), codelen, true, vFlag ? info : 0);
          } else {
            dst.resize(tilelen);
            declen = arithmetic_decode(&dst.at(0), tilelen, &src.at(0), codelen, encType != 0, vFlag ? info : 0);
          }
          t1  = __rdtsc();
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
          if (encType == 2) {
            dst.resize(tilelen);
            ibwt(&dst.at(0), &ibwtInp.at(0), tilelen, primary_i);
          }
          pDst = &dst.at(0);
        }
        uint64_t t2 = __rdtsc();

        if (vFlag)
          printf(
            "%7u -> %7u. Model %7.3f. Coded %7d. %9.0f clocks. %5.1f + %5.1f = %5.1f clocks/char"
            #ifdef ENABLE_PERF_COUNT
            " %6d %6d %6d"
            #endif
            "\n"
            ,unsigned(codelen)
            ,unsigned(tilelen)
            ,info[0]/8.0
            ,info[1]
            ,double(t2-t0)
            ,double(t1-t0)/tilelen
            ,double(t2-t1)/tilelen
            ,double(t2-t0)/tilelen
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


