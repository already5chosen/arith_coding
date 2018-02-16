#include <cstdint>
#include <cstring>
// #include <cstdio>

#include "bwt_mtf_rle_d.h"

// return 0 for success, negative number on failure
static int irle_imtf(
 uint8_t*       dst,     // [dstlen]
 int            dstlen,
 int32_t        histogram[256],
 const uint8_t* src,
 int            srclen)
{
  // initialize move-to-front decoder table
  uint8_t t[256];
  for (int i = 0; i < 256; ++i)
    t[i] = i;

  uint32_t rlAcc = 0;
  uint32_t rlMsb = 1;
  for (int i = 0; i < srclen; ++i) {
    int srcC = src[i];
    if (srcC < 2) {
      // zero run
      rlAcc |= (-srcC) & rlMsb;
      rlMsb += rlMsb;
      if (rlMsb == 0)
        return -1; // zero run too long (A)
    } else {
      if (rlMsb > 1) {
        // insert zero run
        uint32_t rl = rlAcc + rlMsb - 1;
        rlMsb = 1;
        rlAcc = 0;
        if (rl >= uint32_t(dstlen))
          return -2; // zero run too long (B)
        int c0 = t[0];
        histogram[c0] += rl;
        memset(dst, c0, rl);
        dst += rl;
        dstlen -= rl;
      }
      
      if (dstlen == 0) 
        return -3; // decoded section is longer than expected
      dstlen -= 1;
      
      int mtfC = srcC - 1;
      if (srcC == 255) {
        ++i;
        if (i == srclen)
          return -4; // end of input in the middle of escape sequence
        mtfC = src[i];
        if (mtfC < 254)
          return -5; // bad escape sequence
      }
      int c = t[mtfC];
      histogram[c] += 1;
      *dst++  = c;
      
      // update move-to-front decoder table
      memmove(&t[1], &t[0], mtfC);
      t[0] = c;
    }
  }
  if (rlMsb > 1) {
    // insert last zero run
    uint32_t rl = rlAcc + rlMsb - 1;
    if (rl > uint32_t(dstlen))
      return -6; // zero run too long (C)
    int c0 = t[0];
    histogram[c0] += rl;
    memset(dst, c0, rl);
    dstlen -= rl;
  }
  
  if (dstlen > 0)
      return -7; // decoded section is shorter than expected

  return 0;
}

static void ibwt(uint8_t* dst, const uint8_t* src, int32_t* ibwtIdx, int len, int bwtPrimaryIndex, int32_t histogram[256])
{
  // integrate histogram
  int32_t acc = 0;
  for (int i = 0; i < 256; ++i) {
    int32_t h = histogram[i];
    histogram[i] = acc;
    acc += h;
  }

  for (int i = 0; i < len; ++i) {
    int c = src[i];
    ibwtIdx[histogram[c]] = i;
    ++histogram[c];
  }

  int index = bwtPrimaryIndex;
  for (int i = 0; i < len; ++i) {
    dst[i] = src[index];
    index  = ibwtIdx[index];
  }
}

// return 0 for success, negative number on failure
int irle_imtf_ibwt(
 uint8_t*       dst,     // [dstlen]
 uint8_t*       ibwtInp, // [dstlen]
 int32_t*       ibwtIdx, // [dstlen]
 int            dstlen,
 int            bwtPrimaryIndex,
 const uint8_t* src,
 int            srclen)
{
  int32_t histogram[256]={0};
  int err = irle_imtf(ibwtInp, dstlen, histogram, src, srclen);
  if (err == 0)
    ibwt(dst, ibwtInp, ibwtIdx, dstlen, bwtPrimaryIndex, histogram);
  return err;
}

