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

void bwt_reorder_mtf_rle(
 const int32_t*  idx, // result of bwt_sort
 const uint8_t*  src,
 int             srclen,
 int*            pBwtPrimaryIndex,
 void           (*chunkCallback)(void* context, const uint8_t* chunk, int chunklen),
 void*           chunkCallbackContext)
{
  // initialize move-to-front encoder table
  uint8_t t[256];
  for (int i = 0; i < 256; ++i)
    t[i] = i;
  unsigned zRunLen = 0;

  const int NOMINAL_OUT_BUF_LEN = 4000;
  uint8_t  outBuf[NOMINAL_OUT_BUF_LEN+32];
  uint8_t* pOut = outBuf;
  for (int i = 0; i < srclen; ++i) {
    int k = idx[i];
    if (k == 0)
      k = srclen;
    if (k == 1)
      *pBwtPrimaryIndex = i;
    uint8_t c = src[k-1];

    // move-to-front encoder
    int v1 = t[0];
    if (c != v1) {
      // c is not at front
      if (zRunLen != 0)
        pOut += insertZeroRun(pOut, zRunLen);

      t[0] = c;
      int v0 = v1, k;
      for (k = 1; c != (v1=t[k]); ++k) {
        t[k] = v0;
        v0 = v1;
      }
      t[k] = v0;
      int mtfC = k;
      int outC = mtfC + 1;
      *pOut = outC;       // range [1..253] encoded as x+1
      if (mtfC >= 254) { // range [254..255] encode as a pair {255,x}
        *pOut++ = 255;
        *pOut = mtfC;
      }
      ++pOut;
      if (pOut-outBuf >= NOMINAL_OUT_BUF_LEN) {
        chunkCallback(chunkCallbackContext, outBuf, pOut - outBuf);
        pOut = outBuf;
      }
      zRunLen = -1;
    }
    ++zRunLen;
  }

  if (zRunLen != 0)
    pOut += insertZeroRun(pOut, zRunLen);

  if (pOut - outBuf != 0)
    chunkCallback(chunkCallbackContext, outBuf, pOut - outBuf);
}
