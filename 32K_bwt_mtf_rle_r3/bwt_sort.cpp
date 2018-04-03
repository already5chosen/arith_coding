#include <cstring>
#include "bwt_mtf_rle_e.h"

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
  while (len > 16) {
    int32_t valm = bwt_qsort_median_of_3(x, len, ord);
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
      int l = len - h;
      if (l < h) {
        bwt_qsort(x, l, ord, wrk+h);
        x += l;
        len = h;
      } else {
        bwt_qsort(wrk, h, ord, x+l);
        len = l;
      }
    } else {
      // ord[x[]] <= valm
      pl = x;
      ph = wrk;
      for (int i = 0; i < len; ++i) {
        int32_t y = x[i];
        int32_t val = ord[y];
        *pl = y;
        *ph = y;
        pl += (val != valm);
        ph += (val == valm);
      }
      h = ph - wrk;
      len -= h;
    }
    memcpy(pl, wrk, h*sizeof(*x));
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

struct flat_rec_t {
  int32_t key, dat;
};

static int32_t bwt_flat_qsort_median_of_3(const flat_rec_t* x, unsigned len)
{
  int32_t key0 = x[0].key;
  int32_t key1 = x[len/2].key;
  if (key1 < key0) {
    key0 = key1;
    key1 = x[0].key;
  }
  int32_t key2 = x[len-1].key;
  if (key2 <= key0)
    return key0;
  if (key2 <= key1)
    return key2;
  return key1;
}
static void bwt_flat_qsort(flat_rec_t* x, int len, flat_rec_t* wrk)
{
  while (len > 16) {
    int32_t keym = bwt_flat_qsort_median_of_3(x, len);
    flat_rec_t* pl = x;
    flat_rec_t* ph = wrk;
    for (int i = 0; i < len; ++i) {
      flat_rec_t r = x[i];
      *pl = r;
      *ph = r;
      int32_t key = r.key;
      pl += (key <= keym);
      ph += (key >  keym);
    }
    int h = ph - wrk;
    if (h > 0) {
      int l = len - h;
      if (l < h) {
        bwt_flat_qsort(x, l, wrk+h);
        x += l;
        len = h;
      } else {
        bwt_flat_qsort(wrk, h, x+l);
        len = l;
      }
    } else {
      // ord[x[]] <= keym
      pl = x;
      ph = wrk;
      for (int i = 0; i < len; ++i) {
        flat_rec_t r = x[i];
        *pl = r;
        *ph = r;
        int32_t key = r.key;
        pl += (key != keym);
        ph += (key == keym);
      }
      h = ph - wrk;
      len -= h;
    }
    memcpy(pl, wrk, h*sizeof(*x));
  }

  // sort by straight insertion
  int32_t key0 = x[0].key;
  for (int beg = 1; beg < len; ++beg) {
    flat_rec_t r = x[beg];
    int32_t key  = r.key;
    if (key < key0) {
      key0 = r.key;
      for (int i = 0; i < beg; ++i) {
        flat_rec_t r1 = x[i];
        x[i] = r;
        r = r1;
      }
      x[beg] = r;
    } else {
      int i;
      for (i = beg-1; ; --i) {
        flat_rec_t ri  = x[i];
        if (ri.key <= key)
          break;
        x[i+1] = ri;
      }
      x[i+1] = r;
    }
  }
}

static inline uint64_t load_8c(const uint8_t* src) {
  uint64_t r;
  memcpy(&r, src, sizeof(r));
  return r;
}

static void bwt_sort_flatten(flat_rec_t* dst, const int32_t* src, int len, const int32_t* ord, int srclen, int dist)
{
  for (int i = 0; i < len; ++i) {
    int32_t y = src[i];
    int32_t idx = y + dist;
    if (idx-srclen >= 0)
      idx = idx - srclen;
    dst[i].key = ord[idx];
    dst[i].dat = y;
  }
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
        wn += 2;
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
    wn += 2;
  }

  for (int dist = 8; dist < srclen && wn > 0; dist += dist) {
    // copy list of segments that require sorting to the upper part of wrk[]
    int ri = srclen - wn;
    if (ri != 0)
      memmove(&wrk[ri], &wrk[0], wn*sizeof(wrk[0]));

    wn = 0;
    while (ri != srclen) {
      int beg = wrk[ri+0];
      int len = wrk[ri+1];
      ri += 2;

      if (len == 2) {
        // very common simple case
        int32_t y0 = dst[beg+0];
        int32_t y1 = dst[beg+1];
        int32_t i0 = y0 + dist; if (i0-srclen >= 0) i0 = i0-srclen;
        int32_t i1 = y1 + dist; if (i1-srclen >= 0) i1 = i1-srclen;
        int32_t v0 = ord[i0];
        int32_t v1 = ord[i1];
        if (v1 != v0) {
          if (v1 < v0) {
            dst[beg+0] = y1;
            dst[beg+1] = y0;
            y1 = y0;
          }
          ++ord[y1];
        } else {
          // order still unresolved
          wrk[wn+0] = beg;
          wrk[wn+1] = len;
          wn += 2;
        }
      } else if (len*4 <= ri-wn ) {
        // there is sufficient space in the wrk[] to use flat structures
        flat_rec_t* pWrk = reinterpret_cast<flat_rec_t*>(&wrk[wn]);
        bwt_sort_flatten(pWrk, &dst[beg], len, ord, srclen, dist);
        bwt_flat_qsort(pWrk, len, &pWrk[len]);

        int v0 = pWrk[0].key;
        int y0 = pWrk[0].dat;
        dst[beg] = y0;
        ord[y0] = beg;
        int i0 = 0;
        for (int i = 1; i < len; ++i) {
          int v = pWrk[i].key;
          int y = pWrk[i].dat;
          if (v != v0) {
            v0 = v;
            int seglen = i - i0;
            if (seglen > 1) {
              // record segment that requires further sorting
              wrk[wn+0] = beg+i0;
              wrk[wn+1] = seglen;
              wn += 2;
            }
            i0 = i;
          }
          dst[beg+i] = y;
          ord[y] = beg+i0;
        }
        if (len-i0 > 1) {
          // record segment that requires further sorting
          int seglen = len-i0;
          wrk[wn+0] = beg+i0;
          wrk[wn+1] = seglen;
          wn += 2;
        }
      } else {
        int32_t* pWrk = &wrk[wn];
        for (int i = beg; i < beg+len; ++i) {
          int32_t y = dst[i] + dist;
          if (y-srclen >= 0)
            y = y - srclen;
          dst[i] = y;
        }
        bwt_qsort(&dst[beg], len, ord, pWrk);
        int i0 = 0;
        int v0 = ord[dst[beg]];
        pWrk[0] = 0;
        for (int i = 1; i < len; ++i) {
          int v = ord[dst[beg+i]];
          i0 = (v != v0) ? i : i0;
          v0 = v;
          pWrk[i] = i0;
        }
        for (int i = beg; i < beg+len; ++i) {
          int32_t y = dst[i] - dist;
          if (y < 0)
            y = y + srclen;
          dst[i] = y;
        }

        ord[dst[beg]] = beg;
        for (int i = 1; i < len; ++i)
          ord[dst[beg+i]] = pWrk[i]+beg;
        i0 = 0;
        v0 = 0;
        for (int i = 1; i < len; ++i) {
          int v = pWrk[i];
          if (v != v0) {
            int seglen = i - i0;
            if (seglen > 1) {
              // record segment that requires further sorting
              wrk[wn+0] = beg+i0;
              wrk[wn+1] = seglen;
              wn += 2;
            }
            i0 = i;
            v0 = v;
          }
        }
        if (len-i0 > 1) {
          // record segment that requires further sorting
          int seglen = len-i0;
          wrk[wn+0] = beg+i0;
          wrk[wn+1] = seglen;
          wn += 2;
        }
      }
    }
  }
}
