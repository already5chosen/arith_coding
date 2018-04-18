#include <stdint.h>

void bwt_reorder_mtf_rle(
 const int32_t*  idx, // result of bwt_sort
 const uint8_t*  src,
 int             srclen,
 int*            pBwtPrimaryIndex,
 void           (*chunkCallback)(void* context, const uint8_t* chunk, int chunklen),
 void*           chunkCallbackContext);
