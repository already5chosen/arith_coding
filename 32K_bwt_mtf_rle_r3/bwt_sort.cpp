#include <cstring>
#include <algorithm>

#include "bwt_mtf_rle_e.h"

// prepare_bwt_sort append first 8 characters to the end of the source
static void prepare_bwt_sort(uint8_t* src, int srclen)
{
  for (int i = 0; i < 8; ++i)
    src[srclen+i] = src[i];
}

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

static inline uint64_t load_8c(const uint8_t* src) {
  uint64_t r;
  memcpy(&r, src, sizeof(r));
  return r;
}

void bwt_sort(
  int32_t* dst_tmp, // length = srclen*3
  uint8_t* src,     // length = srclen+8, characters [srclen..srclen+7] modified, others preserved
  int      srclen)
{ // BWT sort
  int h[256]={0};
  for (int i = 0; i < srclen; ++i)
    ++h[src[i]];

  prepare_bwt_sort(src, srclen);

  int32_t* dst = &dst_tmp[0];
  int32_t* ord = &dst_tmp[srclen*1];
  int32_t* wrk = &dst_tmp[srclen*2];

  int h0[256];
  int acc=0;
  for (int i = 0; i < 256; ++i) {
    h0[i] = acc;
    int hv = h[i];
    acc += hv;
  }

  // sort by 8th character, result in ord[]
  memcpy(h, h0, sizeof(h));
  for (int i = 0; i < srclen; ++i) {
    unsigned c = src[i+7];
    int hv  = h[c];
    ord[hv] = i;
    h[c]    = hv + 1;
  }

  // sort by 7th, 6th, ... and 1st character, result in dst[]
  for (int offs = 6; offs >= 0; --offs) {
    int32_t* inp = offs % 2 == 0 ? ord : dst;
    int32_t* out = offs % 2 == 0 ? dst : ord;
    memcpy(h, h0, sizeof(h));
    for (int i = 0; i < srclen; ++i) {
      int ii = inp[i];
      unsigned c = src[ii+offs];
      int hv  = h[c];
      out[hv] = ii;
      h[c]    = hv + 1;
    }
  }

  // build ord[] and wrk[]
  int wn = 0;
  int i0  = 0;
  int ii0 = dst[0];
  uint64_t cc0 = load_8c(&src[ii0]);
  ord[ii0] = 0;
  for (int i = 1; i < srclen; ++i) {
    int ii = dst[i];
    uint64_t cc = load_8c(&src[ii]);
    if (cc != cc0) {
      int seglen = i - i0;
      if (seglen > 1) {
        // record segment that requires further sorting
        wrk[wn+0] = i0;
        wrk[wn+1] = seglen;
        wn += seglen;
      }
      i0 = i;
    }
    ord[ii] = i0;
    cc0 = cc;
  }
  if (srclen - i0 > 1) {
    // record last segment that requires further sorting
    int seglen = srclen - i0;
    wrk[wn+0] = i0;
    wrk[wn+1] = seglen;
    wn += seglen;
  }

  bwt_compare cmp;
  cmp.m_src = ord;
  cmp.m_srclen = srclen;
  for (int dist = 8; dist < srclen && wn > 0; dist += dist) {
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
