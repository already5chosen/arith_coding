#include <cstdint>
#include <cstring>
#include <algorithm>
// #include <cstdio>
// #include <cmath>
// #include <cctype>

#include "bwt_mtf_rle_d.h"

void ibwt(uint8_t* dst, const uint8_t* src, int len, int first_i)
{
  int32_t histogram[256]={0};
  for (int i = 0; i < len; ++i)
    ++histogram[src[i]];
  // integrate histogram
  int32_t acc = 0;
  for (int i = 0; i < 256; ++i) {
    int32_t h = histogram[i];
    histogram[i] = acc;
    acc += h;
  }

  int* T = new int[len];
  for (int i = 0; i < len; ++i) {
    int c = src[i];
    T[histogram[c]] = i;
    ++histogram[c];
  }

  int index = first_i;
  for (int i = 0; i < len; ++i) {
    dst[i] = src[index];
    index  = T[index];
  }

  delete [] T;
}