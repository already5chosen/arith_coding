#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>
#include <x86intrin.h>

#include "bwt_mtf_rle_d.h"

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
        uint8_t hdr[9];
        size_t hdrlen = fread(hdr, 1, 9, fpinp);
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
        unsigned bwtPrimaryIndex = loadFrom3octets(&hdr[6]);

        if (tilelen < 1 || tilelen > MAX_TILE_SIZE) {
          fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Illegal tile length.": "");
          break;
        }
        if (codelen < 1 || codelen > tilelen*2) {
          fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Illegal code length.": "");
          break;
        }
        if (bwtPrimaryIndex >= tilelen) {
          fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " primary_i out of range.": "");
          break;
        }

        src.resize(codelen);
        size_t inplen = fread(&src.at(0), 1, codelen, fpinp);
        if (inplen < codelen) {
          if (feof(fpinp))
            fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Code section is shorter than specified in tile header.": "");
          else
            perror(inpfilename);
          break;
        }

        ibwtIdx.resize(tilelen);
        ibwtInp.resize(tilelen);
        dst.resize(tilelen);

        uint64_t t0 = __rdtsc();
        int err = irle_imtf_ibwt(
          &dst.at(0),
          &ibwtInp.at(0),
          &ibwtIdx.at(0),
          tilelen,
          bwtPrimaryIndex,
          &src.at(0),
          codelen);
        if (err != 0) {
          fprintf(stderr, "%s: %s invalid.\n", argv[0], inpfilename);
          if (vFlag)
            fprintf(stderr, "Decoder parsing failure. Error %d.\n", err);
          break;
        }
        uint64_t t1 = __rdtsc();

        if (vFlag)
          printf(
            "*%7d->%7d chars. %10.0f clocks. %6.1f clocks/char\n"
            ,int(codelen)
            ,int(tilelen)
            ,double(t1-t0)
            ,double(t1-t0)/tilelen
            );

        size_t wrlen = fwrite(&dst.at(0), 1, tilelen, fpout);
        if (wrlen != size_t(tilelen)) {
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

static uint32_t loadFrom3octets(const uint8_t* src)
{
  return (src[2]*256 + src[1])*256 + src[0];
}
