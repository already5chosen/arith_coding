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

void bwt_sort(int32_t* dst, const uint8_t* src, int srclen)
{ // BWT sort
  for (int i = 0; i < srclen; ++i)
    dst[i] = i;
  bwt_compare cmp;
  cmp.m_src = src;
  cmp.m_srclen = srclen;
  std::sort(&dst[0], &dst[srclen], cmp);
}
