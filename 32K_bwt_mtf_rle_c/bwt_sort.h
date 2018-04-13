#include <stdint.h>

void bwt_sort(
  uint64_t* dst_tmp, // length = ((srclen*3+256)*sizeof(int32_t))/sizeof(uint64_t)
  uint8_t*  src,     // length = srclen+8, characters [srclen..srclen+7] modified, others preserved
  int       srclen);
