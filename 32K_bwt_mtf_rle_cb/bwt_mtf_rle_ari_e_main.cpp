#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>
#include <x86intrin.h>

#include "bwt_sort.h"
#include "bwt_mtf_rle_e.h"
#include "arithmetic_encode.h"

static void storeAs3octets(uint8_t* dst, unsigned val);
static bool isSingleCharacter(const uint8_t* src, int len);
int main(int argz, char** argv)
{
  if (argz < 3) {
    printf(
      "Usage:\n"
      "bwt_mtf_rle_e input-file output-file [-v[x]]\n"
      );
    return 1;
  }
  char* inpfilename = argv[1];
  char* outfilename = argv[2];
  bool vFlag = (argz > 3) && (strncmp("-v", argv[3], 2)==0);
  bool xFlag = vFlag ? argv[3][2]=='x' : false;

  int ret = 1;
  FILE* fpinp = fopen(inpfilename, "rb");
  if (fpinp) {
    size_t inpfilename_len = strlen(inpfilename);
    char* nametag = inpfilename_len < 4 ? inpfilename : &inpfilename[inpfilename_len-4];
    FILE* fpout = fopen(outfilename, "wb");
    if (fpout) {
      const size_t TILE_SIZE = 1024*1024;
      uint8_t* inptile = new uint8_t[TILE_SIZE+8];
      std::vector<uint64_t> tmpDst;
      ret = 0;
      int tile_i = 0;
      bool done = false;
      while (!done) {
        ++tile_i;
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
            size_t tmpDstLen = ((tilelen*3+256*3)*sizeof(int32_t))/sizeof(uint64_t);
            if (tmpDstLen > tmpDst.size())
              tmpDst.resize(tmpDstLen);

            uint64_t t0 = __rdtsc();
            bwt_sort(&tmpDst.at(0), inptile, tilelen);
            uint64_t t1 = __rdtsc();

            int32_t* bwtIdx = reinterpret_cast<int32_t*>(&tmpDst.at(0));
            uint32_t* encContext = reinterpret_cast<uint32_t*>(&bwtIdx[tilelen+256]);
            arithmetic_encode_init_context(encContext, tilelen);
            int bwtPrimaryIndex;
            bwt_reorder_mtf_rle(
              bwtIdx,
              inptile,
              tilelen,
              &bwtPrimaryIndex,
              arithmetic_encode_chunk_callback,
              encContext);
            uint64_t t2 = __rdtsc();

            uint8_t* ariEncDst = reinterpret_cast<uint8_t*>(bwtIdx);
            double info[64];
            ressz = arithmetic_encode(encContext, ariEncDst, tilelen, vFlag ? info : 0);
            uint64_t t3 = __rdtsc();
            if (vFlag) {
              printf(
               "%4s:%d "
               "%7u->%7u. Model %9.3f. Coded %9.0f. Entropy %11.3f (%11.3f). %10.0f clks. %6.1f+%5.1f+%5.1f=%6.1f clk/char"
               ,nametag
               ,tile_i
               ,unsigned(tilelen)
               ,ressz < 0 ? 0 : (ressz == 0 ? unsigned(tilelen) : unsigned(ressz))
               ,info[1]/8
               ,info[2]/8
               ,info[0]/8
               ,info[3]/8
               ,double(t3-t0)
               ,double(t1-t0)/tilelen
               ,double(t2-t1)/tilelen
               ,double(t3-t2)/tilelen
               ,double(t3-t0)/tilelen
             );
             if (xFlag) {
                printf("\n");
                for (int i = 0; i < 9; ++i) {
                  printf(" %4.0f %6.0f %5.0f %4.2f;%s"
                    ,info[4+9*0+i]
                    ,info[4+9*1+i]
                    ,info[4+9*2+i]/8
                    ,info[4+9*1+i] != 0 ? info[4+9*2+i]/info[4+9*1+i] : 0.0
                    ,i==3 ? "\n" : ""
                    );
                }
                printf("\n");
             } else {
              printf(" (%.0f)\n"
               ,info[4]
               );
             }
            }
            if (ressz > 0) {
              // normal compression
              storeAs3octets(&hdr[3], ressz);
              storeAs3octets(&hdr[6], bwtPrimaryIndex);
              hdrlen = 9;
              pRes = ariEncDst;
            } else {
              // not compressible
              hdr[3] = 0;
              hdr[4] = 0;
              pRes   = inptile;
              ressz  = tilelen;
            }
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