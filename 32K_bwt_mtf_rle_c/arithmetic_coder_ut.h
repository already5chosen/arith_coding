#include <cstdint>

int quantized_histogram_to_range(uint16_t* c2range, unsigned len, const uint8_t* qh, unsigned range_scale);
unsigned quantized_histogram_pair_to_range_qh_scale9(unsigned qh, unsigned range_scale);
unsigned quantized_histogram_pair_to_range_qh_scale255(unsigned qh, unsigned range_scale);
