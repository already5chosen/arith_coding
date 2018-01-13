#include <stdint.h>
#include <vector>

void   arithmetic_encode_calc_histograms(unsigned histograms[2][260], const uint8_t* src, int srclen);
double arithmetic_encode_calc_entropy(const unsigned histogram[260], int srclen);

// return value:
//  0 - not compressible, because all input characters have approximately equal probability or because input is too short
// >0 - the length of compressed buffer
int arithmetic_encode(
 std::vector<uint8_t>* dst, 
 const uint8_t* src, int srclen, 
 const unsigned histogram[260], bool rep, double* pInfo=0);
