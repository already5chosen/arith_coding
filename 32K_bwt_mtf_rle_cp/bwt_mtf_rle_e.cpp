#include <cstring>

#include "bwt_mtf_rle_e.h"

static int insertZeroRun(uint8_t* dst, unsigned zRunLen)
{ // encode length of zero run by method similar to bzip2' RUNA/RUNB
  uint8_t* dst0 = dst;
  zRunLen += 1;
  do {
    int c = zRunLen % 2;
    *dst++ = c;
    zRunLen /= 2;
  } while (zRunLen > 1);
  return dst-dst0;
}

// return value - the length of destination array in octets
int bwt_reorder_mtf_rle(
 int32_t*            idx_dst, // srclen*uint32_t,
                              // on input  - result of bwt_sort
                              // on output - result of BWT followed by move-to-front and by zero-run-len encoding
 const uint8_t*      src,
 int                 srclen,
 int*                pBwtPrimaryIndex,
 int                 (*chunkCallback)(void* context, const uint8_t* chunk, int nSymbols, int nRuns),
 void*               chunkCallbackContext)
{
  // initialize callback machinery
  const int symbolsPerChunk = chunkCallback(chunkCallbackContext, 0, 0, 0);

  // initialize move-to-front encoder table
  uint8_t t[256];
  for (int i = 0; i < 256; ++i)
    t[i] = i;

  int symbolsInChunk = 0;
  uint8_t* dst = reinterpret_cast<uint8_t*>(idx_dst);
  uint8_t* chunk = dst;
  unsigned zRunLen = 0;
  for (int i = 0; i < srclen; ++i) {
    int k = idx_dst[i];
    if (k == 0)
      k = srclen;
    if (k == 1)
      *pBwtPrimaryIndex = i;
    uint8_t c = src[k-1];

    // move-to-front encoder
    int v1 = t[0];
    if (c == v1) {
      // c already at front - count length of zero run
      ++zRunLen;
    } else {
      if (zRunLen != 0) {
        int nSym = insertZeroRun(dst, zRunLen);
        dst += nSym;
        zRunLen = 0;
        symbolsInChunk += nSym;
        if (symbolsInChunk >= symbolsPerChunk) {
          symbolsInChunk -= symbolsPerChunk;
          int chunklen = dst - chunk - symbolsInChunk;
          chunkCallback(chunkCallbackContext, chunk, chunklen, symbolsPerChunk);
          chunk += chunklen;
        }
      }
      t[0] = c;
      int v0 = v1, k;
      for (k = 1; c != (v1=t[k]); ++k) {
        t[k] = v0;
        v0 = v1;
      }
      t[k] = v0;
      int mtfC = k;
      int outC = mtfC + 1;
      *dst = outC;       // range [1..253] encoded as x+1
      if (mtfC >= 254) { // range [254..255] encode as a pair {255,x}
        *dst++ = 255;
        *dst = mtfC;
      }
      ++dst;
      ++symbolsInChunk;
      if (symbolsInChunk == symbolsPerChunk) {
        chunkCallback(chunkCallbackContext, chunk, dst - chunk, symbolsPerChunk);
        chunk = dst;
        symbolsInChunk = 0;
      }
    }
  }
  if (zRunLen != 0) {
    int nSym = insertZeroRun(dst, zRunLen);
    dst += nSym;
    symbolsInChunk += nSym;
  }
  if (symbolsInChunk != 0) {
    chunkCallback(chunkCallbackContext, chunk, dst - chunk, symbolsInChunk);
  }

  return dst - reinterpret_cast<uint8_t*>(idx_dst); // return the length of destination array in octets
}
