#include <cstdint>

int arithmetic_decode(uint8_t* dst, int dstlen, int32_t histogram[256], const uint8_t* src, int srclen, int* pInfo=0);
