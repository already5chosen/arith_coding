#include <cstring>
#include "bwt_mtf_rle_e.h"

// prepare_bwt_sort append first 8 characters to the end of the source
static void prepare_bwt_sort(uint8_t* src, int srclen)
{
  for (int i = 0; i < 8; ++i)
    src[srclen+i] = src[i];
}

static void bwt_mergesort_core(
  int32_t*       dst,
  unsigned       len,
  const int32_t* src1, // [len/2]
  const int32_t* src2, // [len - len/2]
  const int32_t* ord)
{
  unsigned hlen1 = len / 2;
  unsigned hlen2 = len - hlen1;
  const int32_t* end1 = src1 + hlen1;
  const int32_t* end2 = src2 + hlen2;
  #if 1
  int32_t y1 = *src1;
  int32_t y2 = *src2;
  int32_t val1 = ord[y1];
  int32_t val2 = ord[y2];
  --dst;
  for (;;) {
    ++dst;
    if (val1 <= val2) {
      *dst = y1;
      ++src1;
      if (src1 == end1)
        break;
      y1   = *src1;
      val1 = ord[y1];
    } else {
      *dst = y2;
      ++src2;
      if (src2 == end2)
        break;
      y2   = *src2;
      val2 = ord[y2];
    }
  }
  ++dst;
  for (; src1 != end1; ++src1, ++dst)
    *dst = *src1;
  for (; src2 != end2; ++src2, ++dst)
    *dst = *src2;
  #else
  int32_t val1_max = ord[end1[-1]];
  int32_t val2_max = ord[end2[-1]];
  int32_t y1 = *src1++;
  int32_t y2 = *src2++;
  int32_t val1 = ord[y1];
  int32_t val2 = ord[y2];
  for (;;) {
    if (val1_max > val2) {
      while (val1 <= val2) {
        *dst++ = y1;
        y1 = *src1++;
        val1 = ord[y1];
      }
    } else {
      *dst++ = y1;
      while (src1 != end1)
        *dst++ = *src1++;
      *dst++ = y2;
      while (src2 != end2)
        *dst++ = *src2++;
      break;
    }
    if (val2_max > val1) {
      while (val2 <= val1) {
        *dst++ = y2;
        y2 = *src2++;
        val2 = ord[y2];
      }
    } else {
      *dst++ = y2;
      while (src2 != end2)
        *dst++ = *src2++;
      *dst++ = y1;
      while (src1 != end1)
        *dst++ = *src1++;
      break;
    }
  }
  #endif
}

static void bwt_mergesort(int32_t* x, unsigned len, const int32_t* ord, int32_t* wrk)
{
  if (len > 16) {
    unsigned hlen1 = len / 2;
    unsigned qlen1 = hlen1 / 2;
    bwt_mergesort(&x[0],     qlen1,       ord, wrk);
    bwt_mergesort(&x[qlen1], hlen1-qlen1, ord, wrk);
    bwt_mergesort_core(&wrk[0], hlen1, &x[0], &x[qlen1], ord);

    unsigned hlen2 = len - hlen1;
    unsigned qlen3 = hlen2 / 2;
    bwt_mergesort(&x[hlen1],       qlen3,       ord, &wrk[hlen1]);
    bwt_mergesort(&x[hlen1+qlen3], hlen2-qlen3, ord, &wrk[hlen1]);
    bwt_mergesort_core(&wrk[hlen1], hlen2, &x[hlen1], &x[hlen1+qlen3], ord);

    bwt_mergesort_core(&x[0], len, &wrk[0], &wrk[hlen1], ord);
    return;
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

  for (int dist = 8; dist < srclen && wn > 0; dist += dist) {
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
      bwt_mergesort(&dst[beg], len, ord, &wrk[ri]);
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
  }
}
