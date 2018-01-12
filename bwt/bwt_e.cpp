// #include <cstdio>
// #include <cmath>
#include <cstring>
#include <algorithm>

#include "bwt_e.h"

class bwt_compare {
public:
  const uint8_t* m_src;
  int   m_srclen;
  bool operator() (int a, int b)  const
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

void bwt(uint8_t* dst, const uint8_t* src, int srclen)
{
  int* idx = new int[srclen];
  for (int i = 0; i < srclen; ++i)
    idx[i] = i;
  bwt_compare cmp;
  cmp.m_src = src;
  cmp.m_srclen = srclen;
  std::sort(&idx[0], &idx[srclen], cmp);
  int first_i = 0;
  for (int i = 0; i < srclen; ++i) {
    int k = idx[i];
    if (k == 0)
      k = srclen;
    dst[i] = src[k-1];
    if (k == 1)
      first_i = i;
  }
  delete [] idx;
  dst[srclen+0] = (first_i >> 0) % 256;
  dst[srclen+1] = (first_i >> 8) % 256;
  dst[srclen+2] = (first_i >>16) % 256;

}