#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

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

        if (codelen < 1 || codelen > MAX_CODE_SIZE) {
          fprintf(stderr, "%s: %s invalid.%s\n", argv[0], inpfilename, vFlag ? " Illegal code section length.": "");
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

        dst.resize(tilelen);
        int info[8];
        int declen = arithmetic_decode(&dst.at(0), tilelen, &src.at(0), codelen, vFlag ? info : 0);
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

        if (vFlag)
          printf("%7u -> %7u. Model %.3f. Coded %d.\n", unsigned(codelen), unsigned(tilelen), info[0]/8.0, info[1]);
        size_t wrlen = fwrite(&dst.at(0), 1, tilelen, fpout);
        if (wrlen != tilelen) {
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


