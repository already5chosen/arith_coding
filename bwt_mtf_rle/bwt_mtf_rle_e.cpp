// #include <cstdio>
// #include <cmath>
#include <cstring>
#include <algorithm>

#include "bwt_mtf_rle_e.h"

class bwt_compare {
public:
  const uint8_t* m_src;
  int   m_srclen;
  bool operator() (int32_t a, int32_t b)  const
  {
    uint8_t ca = m_src[a];
    uint8_t cb = m_src[b];
    if (ca != cb)
      return ca < cb;
    // ca==cb;
    if (a < b) {
      if (m_srclen > b + 1) {
        int cmp = memcmp(&m_src[a+1], &m_src[b+1], m_srclen - b - 1);
        if (cmp != 0)
          return cmp < 0;
      }
      int cmp = memcmp(&m_src[a+m_srclen-b], &m_src[0], b-a);
      if (cmp != 0)
        return cmp < 0;
      cmp = memcmp(&m_src[0], &m_src[b-a], a);
      if (cmp != 0)
        return cmp < 0;
      return true;
    } else if (b < a) {
      if (m_srclen > a + 1) {
        int cmp = memcmp(&m_src[a+1], &m_src[b+1], m_srclen - a - 1);
        if (cmp != 0)
          return cmp < 0;
      }
      int cmp = memcmp(&m_src[0], &m_src[b+m_srclen-a], a-b);
      if (cmp != 0)
        return cmp < 0;
      cmp = memcmp(&m_src[a-b], &m_src[0], b);
      if (cmp != 0)
        return cmp < 0;
      return false;
    }
    return false;
  }
};

static uint8_t* insertZeroRun(uint8_t* dst, unsigned zRunLen)
{ // encode length of zero run by method similar to bzip2' RUNA/RUNB
  zRunLen += 1;
  do {
    *dst++ = zRunLen % 2;
    zRunLen /= 2;
  } while (zRunLen > 1);
  return dst;
}

int bwt_mtf_rle(  // return the length of destination array in octets
 int32_t*  tmp_dst, // srclen*uint32_t
 int*      pBwtPrimaryIndex,
 int*      pNruns,  // number of runs in the destination array
 const     uint8_t* src, 
 int       srclen)
{
  // BWT sort
  for (int i = 0; i < srclen; ++i)
    tmp_dst[i] = i;
  bwt_compare cmp;
  cmp.m_src = src;
  cmp.m_srclen = srclen;
  std::sort(&tmp_dst[0], &tmp_dst[srclen], cmp);

  // initialize move-to-front encoder table  
  uint8_t t[256];
  for (int i = 0; i < 256; ++i)
    t[i] = i;
  
  int nRuns = 0;
  uint8_t* dst = reinterpret_cast<uint8_t*>(tmp_dst);
  unsigned zRunLen = 0;
  for (int i = 0; i < srclen; ++i) {
    int k = tmp_dst[i];
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
      ++nRuns;
      if (zRunLen != 0) {
        ++nRuns;
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
      *dst = mtfC + 1; // range [1..253] encoded as x+1
      if (mtfC >= 254) { // range [254..255] encode as a pair {255,x}
        *dst++ = 255;
        *dst = mtfC;
      }
      ++dst;
    }
  }  
  if (zRunLen != 0) {
    ++nRuns;
    dst = insertZeroRun(dst, zRunLen);
  }
  *pNruns = nRuns;
  return dst - reinterpret_cast<uint8_t*>(tmp_dst); // return the length of destination array in octets
}
