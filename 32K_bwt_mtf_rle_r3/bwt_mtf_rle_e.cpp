#include <cstring>

#include "bwt_mtf_rle_e.h"

static uint8_t* insertZeroRun(uint8_t* dst, unsigned zRunLen)
{ // encode length of zero run with ternary code
  do {
    zRunLen -= 1;
    int c = zRunLen % 3;
    *dst++ = c;
    zRunLen /= 3;
  } while (zRunLen > 0);
  return dst;
}

// return value - the length of destination array in octets
int bwt_reorder_mtf_rle(
 int32_t*            idx_dst, // srclen*uint32_t,
                              // on input  - result of bwt_sort
                              // on output - result of BWT followed by move-to-front and by zero-run-len encoding
 const uint8_t*      src,
 int                 srclen,
 bwt_mtf_rle_meta_t* pMeta,
 void                (*chunkCallback)(void* context, const uint8_t* chunk, int nSymbols),
 void*               chunkCallbackContext)
{
  // initialize calback machinery
  const int charPerChnck = 4096;

  // initialize move-to-front encoder table
  uint8_t t[256];
  for (int i = 0; i < 256; ++i)
    t[i] = i;

  uint8_t* dst = reinterpret_cast<uint8_t*>(idx_dst);
  uint8_t* chunk = dst;
  unsigned zRunLen = 0;
  // int dbg_i = 0;
  for (int i = 0; i < srclen; ++i) {
    int k = idx_dst[i];
    if (k == 0)
      k = srclen;
    if (k == 1)
      pMeta->bwtPrimaryIndex = i;
    uint8_t c = src[k-1];

    // move-to-front encoder
    // printf("[%3d]=%3d\n", dbg_i, c); ++dbg_i;
    int v1 = t[0];
    if (c == v1) {
      // c already at front - count length of zero run
      ++zRunLen;
      // printf("[%3d]=%3d\n", dbg_i, 0); ++dbg_i;
    } else {
      if (zRunLen != 0) {
        dst = insertZeroRun(dst, zRunLen);
        zRunLen = 0;
      }
      t[0] = c;
      int v0 = v1, k;
      for (k = 1; c != (v1=t[k]); ++k) {
        t[k] = v0;
        v0 = v1;
      }
      t[k] = v0;
      int mtfC = k;
      // printf("[%3d]=%3d\n", dbg_i, mtfC); ++dbg_i;
      int outC = mtfC + 2;
      *dst = outC;       // range [1..252] encoded as x+2
      if (mtfC >= 253) { // range [253..255] encoded as a pair {255,x}
        *dst++ = 255;
        *dst = mtfC;
      }
      ++dst;
      if (dst-chunk >= charPerChnck) {
        chunkCallback(chunkCallbackContext, chunk, dst - chunk);
        chunk = dst;
      }
    }
  }
  if (zRunLen != 0)
    dst = insertZeroRun(dst, zRunLen);
  if (dst - chunk != 0)
    chunkCallback(chunkCallbackContext, chunk, dst - chunk);

  return dst - reinterpret_cast<uint8_t*>(idx_dst); // return the length of destination array in octets
}
