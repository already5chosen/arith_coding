#include <stdint.h>

typedef struct {
 int32_t  bwtPrimaryIndex;
 uint32_t histogram[260]; 
 // histogram[0]      - non-zero characters that come after zero run
 // histogram[1..255] - non-zero characters
 // histogram[256]    - non-first RUNA of zero run
 // histogram[257]    - non-first RUNB of zero run
 // histogram[258]    - first RUNA of zero run
 // histogram[259]    - first RUNB of zero run
} bwt_mtf_rle_meta_t;

// return the length of destination array in octets
int bwt_reorder_mtf_rle(
 int32_t*            idx_dst, // srclen*uint32_t,
                              // on input  - result of bwt_sort
                              // on output - result of BWT followed by move-to-front and by zero-run-len encoding
 const uint8_t*      src,
 int                 srclen,
 bwt_mtf_rle_meta_t* pMeta);
