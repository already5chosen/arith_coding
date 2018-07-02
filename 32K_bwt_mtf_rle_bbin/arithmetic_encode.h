#include <stdint.h>

void arithmetic_encode_init_tables();

// return value:
//  0 - not compressible, because all input characters have approximately equal probability or because input is too short
// >0 - the length of compressed buffer
int arithmetic_encode(uint8_t* dst, const uint8_t* src, int srclen, uint32_t srcHistogram[259], int origlen, uint8_t* tmp, double* pInfo=0);
