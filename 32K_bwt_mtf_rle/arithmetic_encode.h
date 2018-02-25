#include <stdint.h>
#include <vector>

void arithmetic_encode_init_context(uint32_t* context);
void arithmetic_encode_chunk_callback(void* context, const uint8_t* chunk, int chunklen);

// return value:
//  0 - not compressible, because all input characters have approximately equal probability or because input is too short
// >0 - the length of compressed buffer
int arithmetic_encode(uint32_t* context, uint8_t* dst, const uint8_t* src, int srclen, int origlen, double* pInfo=0);
