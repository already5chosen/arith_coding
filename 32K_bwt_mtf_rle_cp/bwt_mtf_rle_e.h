#include <stdint.h>

int bwt_reorder_mtf_rle(  // return the length of destination array in octets
 int32_t*            idx_dst, // srclen*uint32_t,
                              // on input  - result of bwt_sort
                              // on output - result of BWT followed by move-to-front and by zero-run-len encoding
 const uint8_t*      src,
 int                 srclen,
 int*                pBwtPrimaryIndex,
 int                 (*chunkCallback)(void* context, const uint8_t* chunk, int chunklen, int nSymbols),
 void*               chunkCallbackContext);
