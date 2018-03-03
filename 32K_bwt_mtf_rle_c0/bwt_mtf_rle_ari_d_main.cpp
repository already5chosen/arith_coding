#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>
#include <x86intrin.h>

#include "arithmetic_decode.h"
#include "bwt_d.h"

static uint32_t loadFrom3octets(const uint8_t* src);
int main(int argz, char** argv)
{
  if (argz < 3) {
    printf(
      "Usage:\n"
      "bwt_mtf_rle_d input-file output-file [-v]\n"
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
      std::vector<uint8_t> src;
      std::vector<uint8_t> ibwtInp;
      std::vector<int32_t> ibwtIdx;
      std::vector<uint8_t> dst;
      ret = 1;

      while (true) {
        uint8_t hdr[6];
        size_t hdrlen = fread(hdr, 1, 6, fpinp);
        if (hdrlen != sizeof(hdr)) {
          if (ferror(fpinp))
            perror(inpfilename);
          else if (hdrlen == 0)
            ret = 0; // success
          else
            fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Short tile header.": "");
          break;
        }
        size_t tilelen = loadFrom3octets(&hdr[0]);
        size_t codelen = loadFrom3octets(&hdr[3]);

        if (tilelen < 1 || tilelen > MAX_TILE_SIZE) {
          fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Illegal tile length.": "");
          break;
        }
        
        size_t srclen = codelen + 3;
        if (hdr[5] == 255) {
          // special cases
          switch (hdr[3]) {
            case 0:
              // not compressible
              srclen = codelen = tilelen;
              break;
            case 1:
              // input consists of repetition of the same character
              srclen = 0;
              break;
            default:
              fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Illegal escape sequence in section header.": "");
              goto end_loop;
          }
        } else {
          if (codelen < 1 || codelen > tilelen) {
            fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Illegal code section length.": "");
            break;
          }
        }
        
        if (srclen != 0) {
          if (srclen > src.size())
            src.resize(srclen);
          size_t inplen = fread(&src.at(0), 1, srclen, fpinp);
          if (inplen < srclen) {
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
              if (tilelen > dst.size())
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
          unsigned bwtPrimaryIndex = loadFrom3octets(&src.at(0));
          if (bwtPrimaryIndex >= tilelen) {
            fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " primary_i out of range.": "");
            break;
          }
          
          if (tilelen > ibwtInp.size())
            ibwtInp.resize(tilelen);
          int32_t histogram[256];
          int declen = arithmetic_decode(&ibwtInp.at(0), tilelen, histogram, &src.at(3), codelen, vFlag ? info : 0);
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
          // arithmetic decode/rle+mtf decode succeed
          if (tilelen > ibwtIdx.size())
            ibwtIdx.resize(tilelen);
          if (tilelen > dst.size())
            dst.resize(tilelen);
          
          t1 = __rdtsc();
          pDst = &dst.at(0);
          ibwt(
            pDst,
            &ibwtIdx.at(0),
            &ibwtInp.at(0),
            tilelen,
            bwtPrimaryIndex,
            histogram);
        }
        uint64_t t2 = __rdtsc();

        if (vFlag)
          printf(
            "%7u -> %7u. Model %7.3f. Coded %7d. %9.0f clocks. %5.1f+%4.1f=%5.1f clocks/char"
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

static uint32_t loadFrom3octets(const uint8_t* src)
{
  return (src[2]*256 + src[1])*256 + src[0];
}
