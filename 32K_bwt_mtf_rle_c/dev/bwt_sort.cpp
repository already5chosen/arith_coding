#include <cstring>
#include <algorithm>

#include "bwt_mtf_rle_e.h"

class bwt_compare {
public:
  const int* m_src;
  int   m_srclen;
  int   m_dist;
  bool operator() (int32_t a, int32_t b)  const
  {
    a += m_dist;
    b += m_dist;
    a = a >= m_srclen ? a - m_srclen : a;
    b = b >= m_srclen ? b - m_srclen : b;
    return m_src[a] < m_src[b];
  }
};

void bwt_sort(int32_t* dst, const uint8_t* src, int srclen)
{ // BWT sort
  int h[256]={0};
  for (int i = 0; i < srclen; ++i)
    ++h[src[i]];
  int h0[256];
  int acc=0;
  for (int i = 0; i < 256; ++i) {
    h0[i] = acc;
    acc += h[i];
  }
  memcpy(h, h0, sizeof(h));

  int* ord = new int[srclen*2];
  for (int i = 0; i < srclen; ++i) {
    unsigned c = src[i];
    int hv  = h[c];
    ord[i]  = h0[c];
    dst[hv] = i;
    h[c]    = hv + 1;
  }

  int* ord1 = &ord[srclen];
  bwt_compare cmp;
  cmp.m_src = ord;
  cmp.m_srclen = srclen;
  for (int dist = 1; dist < srclen; dist += dist) {
    cmp.m_dist = dist;
    for (int beg = 0; beg < srclen-1; ) {
      int v0 = beg;
      int end;
      for (end = beg + 1; end < srclen; ++end) {
        if (ord[dst[end]] != v0)
          break;
      }
      if (end - beg > 1) {
        std::sort(&dst[beg], &dst[end], cmp);
        int i0 = beg;
        int v0 = -1;
        for (int i = beg; i < end; ++i) {
          int indx0 = dst[i];
          int indx = indx0 + dist;
          if (indx >= srclen)
            indx -= srclen;
          int v = ord[indx];
          if (v != v0)
            i0 = i;
          v0 = v;
          ord1[i-beg] = i0;
        }
        for (int i = beg; i < end; ++i) 
          ord[dst[i]] = ord1[i-beg];
      }
      beg = end;
    }
  }
  delete [] ord;
}
