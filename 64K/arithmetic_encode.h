#include <stdint.h>
#include <vector>

void arithmetic_encode(std::vector<uint8_t>* dst, const uint8_t* src, int srclen, double* pEntropy=0);
