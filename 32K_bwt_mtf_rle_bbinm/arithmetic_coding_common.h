#include <cstdint>
#include "arithmetic_coder_cfg.h"

void arithmetic_coding_common_init_tables();
void build_qh_encode_table(uint16_t c2lo[ARITH_CODER_QH_SCALE+2], const uint8_t modelQh[MODEL_QH_LEN]);
void initialize_qh_mtf_table(uint8_t dst[][ARITH_CODER_QH_SCALE+1], int nRows);


