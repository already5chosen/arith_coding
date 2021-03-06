#include <stdint.h>

void bwt_sort(
  uint64_t* dst_tmp, // length = ((srclen*3+256)*sizeof(int32_t))/sizeof(uint64_t)
  uint8_t* src,     // length = srclen+8, characters [srclen..srclen+7] modified, others preserved
  int      srclen);

typedef struct {
 int32_t bwtPrimaryIndex;
} bwt_mtf_rle_meta_t;

int bwt_reorder_mtf_rle(  // return the length of destination array in octets
 int32_t*            idx_dst, // srclen*uint32_t,
                              // on input  - result of bwt_sort
                              // on output - result of BWT followed by move-to-front and by zero-run-len encoding
 const uint8_t*      src,
 int                 srclen,
 bwt_mtf_rle_meta_t* pMeta,
 void                (*chunkCallback)(void* context, const uint8_t* chunk, int nSymbols),
 void*               chunkCallbackContext);

