#include <cstdint>

void quantized_histogram_to_range(uint16_t* c2range, unsigned len, const uint8_t* qh, unsigned range_scale);
