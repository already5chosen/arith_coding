#include <cstring>
#include <algorithm>

#include "bwt_mtf_rle_e.h"

class bwt_compare {
public:
  const int* m_src;
  int   m_srclen;
  int   m_dist;
  int get(int32_t indx) const {
    int32_t a0 = indx + m_dist;
    int32_t a1 = indx + (m_dist - m_srclen);
    return  m_src[a1 < 0 ? a0 : a1];
  }
  bool operator() (int32_t a, int32_t b) const {
    return get(a) < get(b);
  }
};

void bwt_sort(
  int32_t*       dst_tmp,  // length = srclen*3 
  const uint8_t* src, 
  int            srclen)
{ // BWT sort
  int h[256]={0};
  for (int i = 0; i < srclen; ++i)
    ++h[src[i]];

  int32_t* dst = &dst_tmp[0];
  int32_t* ord = &dst_tmp[srclen*1];
  int32_t* wrk = &dst_tmp[srclen*2];
  int h0[256];
  int acc=0;
  int wn = 0;
  for (int i = 0; i < 256; ++i) {
    h0[i] = acc;
    int hv = h[i];
    if (hv > 1) {
      // record segment that requires further sorting
      wrk[wn+0] = acc;
      wrk[wn+1] = hv;
      wn += hv;
    }
    acc += hv;
  }
  memcpy(h, h0, sizeof(h));

  for (int i = 0; i < srclen; ++i) {
    unsigned c = src[i];
    int hv  = h[c];
    ord[i]  = h0[c];
    dst[hv] = i;
    h[c]    = hv + 1;
  }

  bwt_compare cmp;
  cmp.m_src = ord;
  cmp.m_srclen = srclen;
  for (int dist = 1; dist < srclen && wn > 0; dist += dist) {
    cmp.m_dist = dist;
    int wi = 0;
    for (int ri = 0; ri < wn; ) {
      int beg = wrk[ri+0];
      int len = wrk[ri+1];
      std::sort(&dst[beg], &dst[beg+len], cmp);
      int i0 = 0;
      int v0 = cmp.get(dst[beg]);
      wrk[ri] = 0;
      for (int i = 1; i < len; ++i) {
        int v = cmp.get(dst[beg+i]);
        i0 = (v != v0) ? i : i0;
        v0 = v;
        wrk[ri+i] = i0;
      }
      ord[dst[beg]] = beg;
      for (int i = 1; i < len; ++i)
        ord[dst[beg+i]] = wrk[ri+i]+beg;
      i0 = 0;
      v0 = 0;
      for (int i = 1; i < len; ++i) {
        int v = wrk[ri+i];
        if (v != v0) {
          int seglen = i - i0;
          if (seglen > 1) {
            // record segment that requires further sorting
            wrk[wi+0] = beg+i0;
            wrk[wi+1] = seglen;
            wi += seglen;
          }
          i0 = i;
          v0 = v;
        }
      }
      if (len-i0 > 1) {
        // record segment that requires further sorting
        int seglen = len-i0;
        wrk[wi+0] = beg+i0;
        wrk[wi+1] = seglen;
        wi += seglen;
      }
      ri += len;
    }
    wn = wi;
  }
}
