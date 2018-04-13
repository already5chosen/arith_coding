#include <cstring>
#include "bwt_sort.h"

// prepare_bwt_sort append first 8 characters to the end of the source
static void prepare_bwt_sort(uint8_t* src, int srclen)
{
  for (int i = 0; i < 8; ++i)
    src[srclen+i] = src[i];
}

static void integrate_histograms(unsigned h[3][256])
{ // integrate histograms
  unsigned acc0=0, acc1=0, acc2=0;
  for (int i = 0; i < 256; ++i) {
    unsigned val0 = h[0][i];
    unsigned val1 = h[1][i];
    unsigned val2 = h[2][i];
    h[0][i] = acc0;
    h[1][i] = acc1;
    h[2][i] = acc2;
    acc0 += val0;
    acc1 += val1;
    acc2 += val2;
  }
}

static void bwt_sort_indirect_sort_and_postprocess(int32_t* x, int len, const int32_t* ord, int srclen, int dist, int32_t* wrk)
{
  // sort by radix sort.
  // Since ord[] < 2**24 it is sufficient to sort by 3 8-bit digits

  // offset x[] by dist, calculate histograms
  unsigned h[3][256] = {{0}};
  for (int i = 0; i < len; ++i) {
    int32_t y = x[i] + dist;
    if (y-srclen >= 0)
      y = y - srclen;
    uint32_t val = ord[y];
    x[i] = y;
    ++h[0][val & 0xFF];
    ++h[1][(val >> 8) & 0xFF];
    ++h[2][(val >> 16)& 0xFF];
  }

  integrate_histograms(h);

  // sorting - 1st pass
  for (int i = 0; i < len; ++i) {
    int32_t y = x[i];
    unsigned c = uint32_t(ord[y]) & 0xFF;
    unsigned idx = h[0][c];
    wrk[idx] = y;
    h[0][c] = idx + 1;
  }

  // sorting - 2nd pass
  for (int i = 0; i < len; ++i) {
    int32_t y = wrk[i];
    unsigned c = (uint32_t(ord[y]) >> 8) & 0xFF;
    unsigned idx = h[1][c];
    x[idx] = y;
    h[1][c] = idx + 1;
  }

  // sorting - 3rd pass
  for (int i = 0; i < len; ++i) {
    int32_t y = x[i];
    unsigned c = (uint32_t(ord[y]) >> 16) & 0xFF;
    unsigned idx = h[2][c];
    wrk[idx] = y;
    h[2][c] = idx + 1;
  }

  // undo offset, fill wrk[] with next ord[]
  int i0 = 0;
  int32_t v0 = -1;
  for (int i = 0; i < len; ++i) {
    int32_t y = wrk[i];
    int32_t v = ord[y];
    i0 = (v != v0) ? i : i0;
    v0 = v;
    wrk[i] = i0;
    y -= dist;
    if (y < 0)
      y = y + srclen;
    x[i] = y;
  }
}


#if 0
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
#endif

#if 0
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
#endif

#if 0
static uint64_t bwt_flat_qsort_median_of_3(const flat_rec_t* x, unsigned len)
{
  uint64_t key0 = x[0].get();
  uint64_t key1 = x[len/2].get();
  if (key1 < key0) {
    uint64_t tmp = key0;
    key0 = key1;
    key1 = tmp;
  }
  uint64_t key2 = x[len-1].get();
  if (key2 <= key0)
    return key0;
  if (key2 <= key1)
    return key2;
  return key1;
}
static void bwt_flat_qsort(flat_rec_t* x, int len, flat_rec_t* wrk)
{
  while (len > 16) {
    uint64_t valm = bwt_flat_qsort_median_of_3(x, len);
    flat_rec_t* pl = x;
    flat_rec_t* ph = wrk;
    for (int i = 0; i < len; ++i) {
      uint64_t val = x[i].get();
      pl->set(val);
      ph->set(val);
      pl += (val <= valm);
      ph += (val >  valm);
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
      // ord[x[]] <= valm
      pl = x;
      ph = wrk;
      for (int i = 0; i < len; ++i) {
        uint64_t val = x[i].get();
        pl->set(val);
        ph->set(val);
        pl += (val != valm);
        ph += (val == valm);
      }
      h = ph - wrk;
      len -= h;
    }
    memcpy(pl, wrk, h*sizeof(*x));
  }

  // sort by straight insertion
  uint64_t val0 = x[0].get();
  for (int beg = 1; beg < len; ++beg) {
    uint64_t val = x[beg].get();
    if (val < val0) {
      val0 = val;
      for (int i = 0; i < beg; ++i) {
        uint64_t vali = x[i].get();
        x[i].set(val);
        val = vali;
      }
      x[beg].set(val);
    } else {
      int i;
      for (i = beg-1; ; --i) {
        uint64_t vali = x[i].get();
        if (vali <= val)
          break;
        x[i+1].set(vali);
      }
      x[i+1].set(val);
    }
  }
}
#endif

static void bwt_flat_mergesort_core(
  uint64_t*       dst,
  unsigned        len,
  const uint64_t* src1, // [len/2]
  const uint64_t* src2) // [len - len/2]
{
  unsigned hlen1 = len / 2;
  unsigned hlen2 = len - hlen1;
  const uint64_t* end1 = src1 + hlen1;
  const uint64_t* end2 = src2 + hlen2;
  #if 0
  uint64_t val1 = src1->get();
  uint64_t val2 = src2->get();
  --dst;
  for (;;) {
    ++dst;
    if (val1 <= val2) {
      dst->set(val1);
      ++src1;
      if (src1 == end1)
        break;
      val1 = src1->get();
    } else {
      dst->set(val2);
      ++src2;
      if (src2 == end2)
        break;
      val2 = src2->get();
    }
  }
  ++dst;
  for (; src1 != end1; ++src1, ++dst)
    *dst = *src1;
  for (; src2 != end2; ++src2, ++dst)
    *dst = *src2;
  #elif 0
  uint64_t val1_max = end1[-1];
  uint64_t val2_max = end2[-1];
  uint64_t val1 = *src1++;
  uint64_t val2 = *src2++;
  for (;;) {
    if (val1_max > val2) {
      while (val1 <= val2) {
        *dst++ = val1;
        val1 = *src1++;
      }
    } else {
      --src1;
      while (src1 != end1)
        *dst++ = *src1++;
      --src2;
      while (src2 != end2)
        *dst++ = *src2++;
      break;
    }
    if (val2_max > val1) {
      while (val2 <= val1) {
        *dst++ = val2;
        val2 = *src2++;
      }
    } else {
      --src2;
      while (src2 != end2)
        *dst++ = *src2++;
      --src1;
      while (src1 != end1)
        *dst++ = *src1++;
      break;
    }
  }
  #elif 1
  uint64_t val1_max = end1[-1];
  uint64_t val2_max = end2[-1];
  if (val2_max < val1_max) {
    val1_max = val2_max;
    const uint64_t* tmp = src1;
    src1 = src2;
    src2 = tmp;
    tmp = end1;
    end1 = end2;
    end2 = tmp;
  }
  // src1 will reach the end first
  uint64_t val1 = *src1++;
  uint64_t val2 = *src2++;
  for (;;) {
    if (val1_max > val2) {
      while (val1 <= val2) {
        *dst++ = val1;
        val1 = *src1++;
      }
    } else {
      --src1;
      while (src1 != end1)
        *dst++ = *src1++;
      --src2;
      while (src2 != end2)
        *dst++ = *src2++;
      break;
    }
    while (val2 <= val1) {
      *dst++ = val2;
      val2 = *src2++;
    }
  }
  #else
  uint64_t val1_max = end1[-1].get();
  uint64_t val2_max = end2[-1].get();
  if (val2_max < val1_max) {
    const flat_rec_t* tmp = src1;
    src1 = src2;
    src2 = tmp;
    tmp = end1;
    end1 = end2;
    end2 = tmp;
  }
  // src1 will reach the end first
  do {
    uint64_t val1 = src1->get();
    uint64_t val2 = src2->get();
    uint64_t val = val1 <= val2 ? val1 : val2;
    dst->set(val);
    int sel1 = (val1 <= val2);
    src1 += sel1;
    src2 += sel1 ^ 1;
    ++dst;
  } while (src1 != end1);
  while (src2 != end2)
    *dst++ = *src2++;
  #endif
}

static void bwt_flat_straigh_insertion_sort(
  uint64_t* x,
  unsigned  len)
{
  // sort by straight insertion
  uint64_t val0 = x[0];
  for (int beg = 1; beg < len; ++beg) {
    uint64_t val = x[beg];
    if (val < val0) {
      val0 = val;
      for (int i = 0; i < beg; ++i) {
        uint64_t vali = x[i];
        x[i] = val;
        val = vali;
      }
      x[beg] = val;
    } else {
      int i;
      for (i = beg-1; ; --i) {
        uint64_t vali = x[i];
        if (vali <= val)
          break;
        x[i+1] = vali;
      }
      x[i+1] = val;
    }
  }
}

static void bwt_flat_mergesort(
  uint64_t* x,
  unsigned  len,
  uint64_t* wrk) // wrk[] length = (len+1)/2
{
  if (len > 16) {
    unsigned hlen1 = len / 2;
    unsigned hlen2 = len - hlen1;

    unsigned qlen3 = hlen2 / 2;
    bwt_flat_mergesort(&x[hlen1],       qlen3,       &wrk[0]);
    bwt_flat_mergesort(&x[hlen1+qlen3], hlen2-qlen3, &wrk[0]);
    bwt_flat_mergesort_core(&wrk[0], hlen2, &x[hlen1], &x[hlen1+qlen3]);

    unsigned qlen1 = hlen1 / 2;
    bwt_flat_mergesort(&x[0],     qlen1,       &x[hlen1]);
    bwt_flat_mergesort(&x[qlen1], hlen1-qlen1, &x[hlen1]);
    bwt_flat_mergesort_core(&x[hlen2], hlen1, &x[0], &x[qlen1]);

    bwt_flat_mergesort_core(&x[0], len, &x[hlen2], &wrk[0]);
    return;
  }

  // sort by straight insertion
  bwt_flat_straigh_insertion_sort(x, len);
}

static void bwt_sort_flatten(uint64_t* dst, const int32_t* src, int len, const int32_t* ord, int srclen, int dist)
{
  for (int i = 0; i < len; ++i) {
    int32_t y = src[i];
    int32_t idx = y + dist;
    if (idx-srclen >= 0)
      idx = idx - srclen;
    dst[i] = (uint64_t(uint32_t(ord[idx])) << 32) | uint32_t(y);
  }
}

static void bwt_sort_flatten_and_count(
  uint64_t*      dst,
  const int32_t* src,
  int            len,
  const int32_t* ord,
  int            srclen,
  int            dist,
  unsigned       h[3][256])
{
  for (int i = 0; i < len; ++i) {
    int32_t y = src[i];
    int32_t idx = y + dist;
    if (idx-srclen >= 0)
      idx = idx - srclen;
    uint32_t key = ord[idx];
    dst[i] = (uint64_t(key) << 32) | uint32_t(y);
    ++h[0][key & 0xFF];
    ++h[1][(key >> 8) & 0xFF];
    ++h[2][(key >> 16) & 0xFF];
  }
}

#if 0
static void bwt_flatten_and_radixsort(
  const int32_t* src,
  int            len,
  const int32_t* ord,
  int            srclen,
  int            dist,
  uint64_t*      dst,
  uint64_t*      wrk)
{
  unsigned h[3][256] = {{0}};
  bwt_sort_flatten_and_count(wrk, src, len, ord, srclen, dist, h);

  integrate_histograms(h);

  // sorting - 1st pass. wrk->dst
  for (int i = 0; i < len; ++i) {
    uint64_t val = wrk[i];
    unsigned c = (val >> (8*4)) & 0xFF;
    unsigned idx = h[0][c];
    dst[idx] = val;
    h[0][c] = idx + 1;
  }

  // sorting - 2nd pass. dst->wrk
  for (int i = 0; i < len; ++i) {
    uint64_t val = dst[i];
    unsigned c = (val >> (8*5)) & 0xFF;
    unsigned idx = h[1][c];
    wrk[idx] = val;
    h[1][c] = idx + 1;
  }

  // sorting - 3rd pass. wrk->dst
  for (int i = 0; i < len; ++i) {
    uint64_t val = wrk[i];
    unsigned c = (val >> (8*6)) & 0xFF;
    unsigned idx = h[2][c];
    dst[idx] = val;
    h[2][c] = idx + 1;
  }
}
#endif

static void bwt_flat_radixsort_core(uint64_t* dst, const uint64_t* src, unsigned len, unsigned h[256], unsigned shift)
{
#if 0
  for (unsigned i = 0; i < len; ++i) {
    uint64_t val = src[i];
    unsigned c = (val >> shift) & 0xFF;
    unsigned idx = h[c];
    dst[idx] = val;
    h[c] = idx + 1;
  }
#else
  if (len % 2) {
    uint64_t val = src[0];
    unsigned c = (val >> shift) & 0xFF;
    unsigned idx = h[c];
    dst[idx] = val;
    h[c] = idx + 1;
    ++src;
  }
  unsigned hlen = len / 2;
  for (unsigned i = 0; i < hlen; ++i) {
    uint64_t val0 = src[i*2+0];
    uint64_t val1 = src[i*2+1];
    unsigned c0 = (val0 >> shift) & 0xFF;
    unsigned c1 = (val1 >> shift) & 0xFF;
    unsigned idx0 = h[c0];
    unsigned idx1 = h[c1];

    h[c0] = idx0 + 1;
    idx1 = (c0 == c1) ? idx0 + 1 : idx1;
    h[c1] = idx1 + 1;

    dst[idx0] = val0;
    dst[idx1] = val1;
  }
#endif
}

static void bwt_flatten_and_radixsort(
  const int32_t* src,
  int            len,
  const int32_t* ord,
  int            srclen,
  int            dist,
  uint64_t*      dst,
  uint64_t*      wrk)
{
  unsigned h[3][256] = {{0}};
  bwt_sort_flatten_and_count(wrk, src, len, ord, srclen, dist, h);
  integrate_histograms(h);
  // sorting - 1st pass. wrk->dst
  bwt_flat_radixsort_core(dst, wrk, len, h[0], 8*4);
  // sorting - 2nd pass. dst->wrk
  bwt_flat_radixsort_core(wrk, dst, len, h[1], 8*5);
  // sorting - 3rd pass. wrk->dst
  bwt_flat_radixsort_core(dst, wrk, len, h[2], 8*6);
}

static void bwt_flatten_and_sort(
  int32_t*       src,
  unsigned       len,
  const int32_t* ord,
  int            srclen,
  int            dist,
  uint64_t*      wrk,
  int            wrklen) // wrklen >= (len*3+1)/2
{
  if (len <= 512) {
    bwt_sort_flatten(wrk, src, len, ord, srclen, dist);
    if (len <= 32)
      bwt_flat_straigh_insertion_sort(wrk, len);
    else
      bwt_flat_mergesort(wrk, len, &wrk[len]);
  } else {
    if (len*2 <= wrklen) {
      bwt_flatten_and_radixsort(src, len, ord, srclen, dist, &wrk[0], &wrk[len]);
    } else {
      unsigned hlen1 = len/2;
      unsigned hlen2 = len-hlen1;
      bwt_flatten_and_radixsort(&src[0],     hlen1, ord, srclen, dist, &wrk[hlen2], &wrk[len]);
      bwt_flatten_and_radixsort(&src[hlen1], hlen2, ord, srclen, dist, &wrk[len],   &wrk[0]);
      bwt_flat_mergesort_core(&wrk[0], len, &wrk[hlen2], &wrk[len]);
    }
  }
}

static inline uint64_t load_8c(const uint8_t* src) {
  uint64_t r;
  memcpy(&r, src, sizeof(r));
  return r;
}

void bwt_sort(
  uint64_t* dst_tmp, // length = ((srclen*3+256)*sizeof(int32_t))/sizeof(uint64_t)
  uint8_t* src,     // length = srclen+8, characters [srclen..srclen+7] modified, others preserved
  int      srclen)
{ // BWT sort
  int h[256]={0};
  for (int i = 0; i < srclen; ++i)
    ++h[src[i]];

  prepare_bwt_sort(src, srclen);

  int32_t* dst = reinterpret_cast<int32_t*>(&dst_tmp[0]);
  int32_t* ord = &dst[srclen];
  uint64_t* wrk = &dst_tmp[(srclen*2*sizeof(int32_t))/sizeof(uint64_t)];

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
        wrk[wn] = (uint64_t(i0) << 32) + uint32_t(seglen);
        wn += 1;
      }
      i0 = i;
    }
    ord[ii] = i0;
    cc0 = cc;
  }
  if (srclen - i0 > 1) {
    // record last segment that requires further sorting
    int seglen = srclen - i0;
    wrk[wn] = (uint64_t(i0) << 32) + uint32_t(seglen);
    wn += 1;
  }

  int wrklen = ((srclen*3+256)*sizeof(int32_t))/sizeof(uint64_t) - (srclen*2*sizeof(int32_t))/sizeof(uint64_t);
  for (int dist = 8; dist < srclen && wn > 0; dist += dist) {
    // copy list of segments that require sorting to the upper part of wrk[]
    int ri = wrklen - wn;
    if (ri != 0)
      memmove(&wrk[ri], &wrk[0], wn*sizeof(wrk[0]));

    wn = 0;
    while (ri != wrklen) {
      int beg = wrk[ri] >> 32;
      int len = uint32_t(wrk[ri]);
      ri += 1;

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
          wrk[wn] = (uint64_t(beg) << 32) + uint32_t(len);
          wn += 1;
        }
      } else if (len*3+2 <= (ri-wn)*2) {
        // there is sufficient space in the wrk[] to use flat structures
        uint64_t* pWrk = &wrk[wn];
        bwt_flatten_and_sort(&dst[beg], len, ord, srclen, dist, pWrk, ri-wn);

        int v0 = uint32_t(pWrk[0] >> 32);
        int y0 = uint32_t(pWrk[0]);
        dst[beg] = y0;
        ord[y0] = beg;
        int i0 = 0;
        for (int i = 1; i < len; ++i) {
          int v = uint32_t(pWrk[i] >> 32);
          int y = uint32_t(pWrk[i]);
          if (v != v0) {
            v0 = v;
            int seglen = i - i0;
            if (seglen > 1) {
              // record segment that requires further sorting
              wrk[wn] = (uint64_t(beg+i0) << 32) + uint32_t(seglen);
              wn += 1;
            }
            i0 = i;
          }
          dst[beg+i] = y;
          ord[y] = beg+i0;
        }
        if (len-i0 > 1) {
          // record segment that requires further sorting
          int seglen = len-i0;
          wrk[wn] = (uint64_t(beg+i0) << 32) + uint32_t(seglen);
          wn += 1;
        }
      } else {
        int32_t* pWrk = reinterpret_cast<int32_t*>(&wrk[wn]);
        bwt_sort_indirect_sort_and_postprocess(&dst[beg], len, ord, srclen, dist, pWrk);

        ord[dst[beg]] = beg;
        int i0 = 0;
        int32_t v0 = 0;
        for (int i = 1; i < len; ++i) {
          int32_t v = pWrk[i];
          ord[dst[beg+i]] = v+beg;
          if (v != v0) {
            int seglen = i - i0;
            if (seglen > 1) {
              // record segment that requires further sorting
              wrk[wn] = (uint64_t(beg+i0) << 32) + uint32_t(seglen);
              wn += 1;
            }
            i0 = i;
            v0 = v;
          }
        }
        if (len-i0 > 1) {
          // record segment that requires further sorting
          int seglen = len-i0;
          wrk[wn] = (uint64_t(beg+i0) << 32) + uint32_t(seglen);
          wn += 1;
        }
      }
    }
  }
}
