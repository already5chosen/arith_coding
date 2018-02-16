#include <cstdint>

// return 0 for success, negative number on failure
int irle_imtf_ibwt(
 uint8_t*       dst,     // [dstlen]
 uint8_t*       ibwtInp, // [dstlen]
 int32_t*       ibwtIdx, // [dstlen]
 int            dstlen,
 int            bwtPrimaryIndex,
 const uint8_t* src, 
 int            srclen);
