#include <cstring>

#include "arithmetic_coding_common.h"
#include "arithmetic_coder_cfg.h"

static const unsigned VAL_RANGE = 1u << RANGE_BITS;
static uint32_t model_qh2h_tab[QQH_SCALE] = {
  0,
    46927670,   185659716,   410132882,   710536612,
  1073741824,  1483874706,  1923010482,  2371956814,
  2811092590,  3221225472,  3584430684,  3884834414,
  4109307580,  4248039626,
};
static uint32_t model_qh_tail_tab[ARITH_CODER_QH_SCALE+1];
static uint8_t  qhMtf0[ARITH_CODER_QH_SCALE+1];

void arithmetic_coding_common_init_tables()
{
  uint64_t val = uint64_t(1) << 32;
  uint64_t sum = val;
  for (int i = MODEL_QH_LEN+1; i <= ARITH_CODER_QH_SCALE; ++i) {
    val -= val / ARITH_CODER_QH_DECEY;
    sum += val;
    model_qh_tail_tab[i] = (val << 32)/sum;
  }

  for (int i = 0; i <= ARITH_CODER_QH_SCALE/2; ++i)
    qhMtf0[i] = ARITH_CODER_QH_SCALE/2 - i;
  for (int i = ARITH_CODER_QH_SCALE/2+1; i <= ARITH_CODER_QH_SCALE; ++i)
    qhMtf0[i] = i;
}

void build_qh_encode_table(uint16_t c2lo[ARITH_CODER_QH_SCALE+2], const uint8_t modelQh[MODEL_QH_LEN])
{
  unsigned ranges[ARITH_CODER_QH_SCALE+1]={0};
  unsigned remVal = VAL_RANGE;
  for (int i = 0; i < MODEL_QH_LEN; ++i) {
    uint32_t qhVal = modelQh[i];
    if (qhVal != QQH_SCALE) {
      unsigned ra = (uint64_t(model_qh2h_tab[qhVal]) * remVal + (uint32_t(1)<<31)) >> 32;
      if (remVal-ra < ARITH_CODER_QH_SCALE-i)
        ra = remVal - (ARITH_CODER_QH_SCALE-i);
      ranges[i] = ra;
      remVal   -= ra;
    } else {
      ranges[i] = remVal;
      remVal = 0;
      break;
    }
  }
  if (remVal > 0) {
    for (int i = ARITH_CODER_QH_SCALE; i > MODEL_QH_LEN; --i) {
      unsigned ra = (uint64_t(model_qh_tail_tab[i]) * remVal + (uint32_t(1)<<31)) >> 32;
      if (ra == 0) ra = 1;
      ranges[i] = ra;
      remVal -= ra;
    }
    ranges[MODEL_QH_LEN] = remVal;
  }

  unsigned acc = 0;
  for (int i = 0; i < ARITH_CODER_QH_SCALE+1; ++i) {
    c2lo[i] = acc;
    acc += ranges[i];
  }
  c2lo[ARITH_CODER_QH_SCALE+1] = VAL_RANGE;
}

void initialize_qh_mtf_table(uint8_t dst[][ARITH_CODER_QH_SCALE+1], int nRows)
{
  for (int i = 0; i < nRows; ++i)
    memcpy(dst[i], qhMtf0, sizeof(qhMtf0));
}
