#include <stdint.h>

int bwt_mtf_rle(  // return the length of destination array in octets
 int32_t*  tmp_dst, // srclen*uint32_t
 int*      pBwtPrimaryIndex,
 int*      pNruns,  // number of runs in the destination array
 const     uint8_t* src, 
 int       srclen);
