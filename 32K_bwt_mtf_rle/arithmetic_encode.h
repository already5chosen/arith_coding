#include <stdint.h>
#include <vector>

int arithmetic_encode_get_context_length(int tilelen); // return # of 32-bit words
int arithmetic_encode_chunk_callback(void* context, const uint8_t* chunk, int nSymbols);

// return value:
// -1 - input consists of repetition of the same character
//  0 - not compressible, because all input characters have approximately equal probability or because input is too short
// >0 - the length of compressed buffer
int arithmetic_encode(std::vector<uint8_t>* dst, const uint8_t* src, int srclen, double* pInfo=0);
