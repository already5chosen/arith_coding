#include <cstring>

#include "bwt_mtf_rle_e.h"

// return the length of destination array in octets
static int insertZeroRun(uint8_t* dst, unsigned zRunLen)
{ // encode length of zero run by method similar to bzip2' RUNA/RUNB
  uint8_t* dst0 = dst;
  zRunLen += 1;
  do {
    int c = zRunLen % 2;
    *dst++ = c;
    zRunLen /= 2;
  } while (zRunLen > 1);

  // reverse order
  uint8_t* p0 = dst0;
  uint8_t* p1 = dst - 1;
  while (p0 < p1) {
    uint8_t v0 = *p0;
    uint8_t v1 = *p1;
    *p0++ = v1;
    *p1-- = v0;
  }

  return dst-dst0;
}

// return the length of destination array in octets
int bwt_reorder_mtf_rle(
 int32_t*            idx_dst, // srclen*uint32_t,
                              // on input  - result of bwt_sort
                              // on output - result of BWT followed by move-to-front and by zero-run-len encoding
 const uint8_t*      src,
 int                 srclen,
 bwt_mtf_rle_meta_t* pMeta)
{
  // initialize histogram
  memset(pMeta->histogram, 0, sizeof(pMeta->histogram));

  // initialize move-to-front encoder table
  uint8_t t[256];
  for (int i = 0; i < 256; ++i)
    t[i] = i;
  unsigned zRunLen = 0;

  uint8_t* dst = reinterpret_cast<uint8_t*>(idx_dst);
  for (int i = 0; i < srclen; ++i) {
    int k = idx_dst[i];
    if (k == 0)
      k = srclen;
    if (k == 1)
      pMeta->bwtPrimaryIndex = i;
    uint8_t c = src[k-1];

    // move-to-front encoder
    int v1 = t[0];
    if (c != v1) {
      // c is not at front
      if (zRunLen != 0) {
        int nz = insertZeroRun(dst, zRunLen);
        pMeta->histogram[257] += 1;  // count # of zero runs
        pMeta->histogram[0]   += nz; // count # zero characters
        dst += nz;
      }

      t[0] = c;
      int v0 = v1, k;
      for (k = 1; c != (v1=t[k]); ++k) {
        t[k] = v0;
        v0 = v1;
      }
      t[k] = v0;
      int mtfC = k;
      pMeta->histogram[mtfC] += 1;
      int outC = mtfC + 1;
      *dst = outC;       // range [1..253] encoded as x+1
      if (mtfC >= 254) { // range [254..255] encoded as a pair {255,x}
        *dst++ = 255;
        *dst = mtfC;
      }
      ++dst;
      zRunLen = -1;
    }
    ++zRunLen;
  }

  if (zRunLen != 0) {
    int nz = insertZeroRun(dst, zRunLen);
    pMeta->histogram[257] += 1;  // count # of zero runs
    pMeta->histogram[0]   += nz; // count # zero characters
    dst += nz;
  }
  pMeta->histogram[256] = (zRunLen != 0) ? // # of characters preceded by zero characters
    pMeta->histogram[0] - 1:
    pMeta->histogram[0];

  return dst - reinterpret_cast<uint8_t*>(idx_dst);
}
