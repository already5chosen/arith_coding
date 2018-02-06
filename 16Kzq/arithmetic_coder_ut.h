#include <cstdint>

void quantized_histogram_to_range(uint16_t* c2range, unsigned maxC, const uint8_t* qh, unsigned range_scale);
