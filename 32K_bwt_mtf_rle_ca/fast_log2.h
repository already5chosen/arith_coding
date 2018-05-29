#include <stdint.h>

#ifdef __cplusplus  
extern "C"
#endif
double fast_log2(uint32_t x); // x > 0, precision: ~10.5 digits (max error  -4.3e-011)
