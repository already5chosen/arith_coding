#include <cstring>
#include <algorithm>

#include "bwt_mtf_rle_e.h"

#define PRF_BWT_SORT 0
#if PRF_BWT_SORT
 #include <cstdio>
 #include <x86intrin.h>
 #define PRF_Trace  ts[tsi++] = __rdtsc()
#else
 #define PRF_Trace
#endif

// prepare_bwt_sort append first 8 characters to the end of the source
static void prepare_bwt_sort(uint8_t* src, int srclen)
{
  for (int i = 0; i < 8; ++i)
    src[srclen+i] = src[i];
}

static int32_t bwt_qsort_median_of_3(const int32_t* x, int len, const int32_t* ord)
{
  int32_t y0 = x[0];
  int32_t y1 = x[len/2];
  int32_t val0 = ord[y0];
  int32_t val1 = ord[y1];
  if (val1 < val0) {
    val0 = val1;
    val1 = ord[y0];
  }
  int32_t val2 = ord[x[len-1]];
  if (val2 <= val0)
    return val0;
  if (val2 <= val1)
    return val2;
  return val1;
}

static void bwt_qsort(int32_t* x, int len, const int32_t* ord, int32_t* wrk)
{
  #if 0
  if (len > 16) {
    bwt_compare cmp;
    cmp.m_src = ord;
    std::sort(x, x + len, cmp);
    return;
  }
  #endif

  while (len > 16) {
    int32_t valm = bwt_qsort_median_of_3(x, len, ord);
    #if 0
    int32_t* wrk_l = &wrk[0];
    int32_t* wrk_h = &wrk[len-1];
    for (int i = 0; i < len; ++i) {
      int32_t y = x[i];
      int32_t val = ord[y];
      *wrk_l = y;
      *wrk_h = y;
      wrk_l += (val < valm);
      wrk_h -= (val > valm);
    }
    int32_t* wrk_m = wrk_l;
    int32_t* end_m = wrk_h + 1;
    for (int i = 0; wrk_m != end_m; ++i) {
      int32_t y = x[i];
      int32_t val = ord[y];
      *wrk_m = y;
      wrk_m += (val == valm);
    }
    memcpy(x, wrk, len*sizeof(*x));

    int l = wrk_l - wrk;
    int m = end_m - wrk_l;
    int h = len - m - l;
    if (l < h) {
      bwt_qsort(x, l, ord, wrk);
      x += l + m;
      len = h;
    } else {
      bwt_qsort(x + l + m, h, ord, wrk);
      len = l;
    }
    #else
    int32_t* pl = x;
    int32_t* ph = wrk;
    for (int i = 0; i < len; ++i) {
      int32_t y = x[i];
      int32_t val = ord[y];
      *pl = y;
      *ph = y;
      pl += (val <= valm);
      ph += (val >  valm);
    }
    int h = ph - wrk;
    if (h > 0) {
      memcpy(pl, wrk, h*sizeof(*x));
      int l = len - h;
      if (l < h) {
        bwt_qsort(x, l, ord, wrk);
        x += l;
        len = h;
      } else {
        bwt_qsort(x + l, h, ord, wrk);
        len = l;
      }
    } else {
      // ord[x[]] <= valm
      int32_t* pl = x;
      int32_t* pm = wrk;
      for (int i = 0; i < len; ++i) {
        int32_t y = x[i];
        int32_t val = ord[y];
        *pl = y;
        *pm = y;
        pl += (val != valm);
        pm += (val == valm);
      }
      int m = pm - wrk;
      memcpy(pl, wrk, m*sizeof(*x));
      len -= m;
    }
    #endif
  }

  // sort by straight insertion
  int32_t val0 = ord[x[0]];
  for (int beg = 1; beg < len; ++beg) {
    int32_t y   = x[beg];
    int32_t val = ord[y];
    if (val < val0) {
      val0 = val;
      for (int i = 0; i < beg; ++i) {
        int32_t y1 = x[i];
        x[i] = y;
        y = y1;
      }
      x[beg] = y;
    } else {
      int i;
      for (i = beg-1; ; --i) {
        int32_t yi  = x[i];
        if (ord[yi] <= val)
          break;
        x[i+1] = yi;
      }
      x[i+1] = y;
    }
  }
}

class bwt_compare {
public:
  const int* m_src;
  bool operator() (int32_t a, int32_t b) const {
    return m_src[a] < m_src[b];
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
#if PRF_BWT_SORT
  uint64_t ts[64];
  int tsi=0;
#endif

  PRF_Trace;

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
  PRF_Trace;

  bwt_compare cmp;
  cmp.m_src = ord;
  for (int dist = 8; dist < srclen && wn > 0; dist += dist) {
#if PRF_BWT_SORT
    uint64_t dt = ts[tsi-1];
#endif
    int wi = 0;
    for (int ri = 0; ri < wn; ) {
      int beg = wrk[ri+0];
      int len = wrk[ri+1];

      for (int i = beg; i < beg+len; ++i) {
        int32_t y = dst[i] + dist;
        if (y-srclen >= 0)
          y = y - srclen;
        dst[i] = y;
      }
#if PRF_BWT_SORT
      uint64_t t0 =  __rdtsc();
#endif
      // std::sort(&dst[beg], &dst[beg+len], cmp);
      bwt_qsort(&dst[beg], len, ord, &wrk[ri]);
#if PRF_BWT_SORT
      uint64_t t1 =  __rdtsc();
      dt += t1 - t0;
#endif
      int i0 = 0;
      int v0 = ord[dst[beg]];
      wrk[ri] = 0;
      for (int i = 1; i < len; ++i) {
        int v = ord[dst[beg+i]];
        i0 = (v != v0) ? i : i0;
        v0 = v;
        wrk[ri+i] = i0;
      }
      for (int i = beg; i < beg+len; ++i) {
        int32_t y = dst[i] - dist;
        if (y < 0)
          y = y + srclen;
        dst[i] = y;
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
#if PRF_BWT_SORT
    ts[tsi++] = dt;
#endif
    PRF_Trace;
  }

#if PRF_BWT_SORT
  double tot = ts[tsi-1]-ts[0];
  printf("%9.0f:", tot);
  for (int i = 0; i < tsi-1; ++i)
    printf(" %I64u", (ts[i+1]-ts[i])/100);
  printf(" : ");
  for (int i = 0; i < tsi-1; ++i)
    printf(" %.2f", (ts[i+1]-ts[i])*(100./tot));
  printf("\n");
#endif
}
