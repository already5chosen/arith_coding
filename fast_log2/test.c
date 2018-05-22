#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <x86intrin.h>

double fast_log2(uint32_t x);
int main(void)
{
  const int ARR_SZ = 1 << 24;
  double* res = malloc(ARR_SZ*sizeof(double));
  if (!res)
    return 1;
  
  memset(res, 0, ARR_SZ*sizeof(double));
  uint64_t t0 = __rdtsc();
  for (int i = 0; i < ARR_SZ; ++i)
    res[i] = fast_log2(i+1);
  uint64_t t1 = __rdtsc();
  
  double err_sum1=0;
  double err_sum2=0;
  double err_min=1e100;
  double err_max=-1e100;
  int min_i = 0, max_i = 0;
  for (int i = 0; i < ARR_SZ; ++i) {
    double diff = res[i] - log2(i+1);
    err_sum1 += diff;
    err_sum2 += diff*diff;
    if (diff < err_min) { err_min = diff; min_i = i;}
    if (diff > err_max) { err_max = diff; max_i = i;}
  }
  printf("min %e at %d, max %e at %d, mean %e, std %e. %.2f clck/it\n"
    ,err_min, min_i + 1
    ,err_max, max_i + 1
    , err_sum1 / ARR_SZ
    , sqrt(err_sum2*ARR_SZ - err_sum1*err_sum1)/ARR_SZ
    , (double)(t1-t0)/ARR_SZ
    );
  
  free(res);
  return 0;
}
