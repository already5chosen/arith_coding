#include <cstdint>
#include <cstring>
#include <cstdio>
#include <cmath>

#include "arithmetic_decode.h"
#include "arithmetic_coder_cfg.h"

static const unsigned VAL_RANGE = 1u << RANGE_BITS;

static uint16_t qh2h_tab[ARITH_CODER_QH_SCALE+1];

void arithmetic_decode_init_tables()
{
  qh2h_tab[ARITH_CODER_QH_SCALE] = VAL_RANGE;
  double M_PI = 3.1415926535897932384626433832795;
  for (int i = 1; i < ARITH_CODER_QH_SCALE; ++i)
    qh2h_tab[i] = int(round((sin((i*2-ARITH_CODER_QH_SCALE)*(M_PI/2/ARITH_CODER_QH_SCALE))+1.0)*(VAL_RANGE/2)));
}

static int decode(
  uint8_t*       dst0,
  int            dstlen,
  int32_t        histogram[256],
  const uint8_t* src,
  int            srclen)
{
  uint8_t  cntrs[258] = {0};
  uint8_t  prevQh[258];
  uint16_t currH[258];
  memset(prevQh, ARITH_CODER_QH_SCALE/2, sizeof(prevQh));

  // initialize move-to-front decoder table
  uint8_t mtf_t[256];
  for (int i = 0; i < 256; ++i)
    mtf_t[i] = i;

  // initialize RLE decoder
  uint32_t rlInc = 1;

  // initialize arithmetic decoder
  const uint64_t MIN_RANGE = uint64_t(1) << (33-RANGE_BITS);

  uint8_t tmpbuf[16] = {0};
  bool useTmpbuf = false;
  uint64_t value = 0; // scaled by 2**64
  for (int k = 0; k < 8; ++k) {
    value <<= 8;
    if (srclen > 0)
      value |= *src++;
    --srclen;
  }
  uint64_t range = uint64_t(1) << (64-RANGE_BITS);   // scaled by 2**(64-rb).  Maintained in [2**(33-rb)+1..2*(64-rb)]

  unsigned prevC0 = 0;
  uint8_t* dst = dst0;
  while (dstlen > 0) {
    unsigned tLo = 0, tHi = 256;
    do {
      unsigned tMid = (tLo*3 + tHi)/4;
      unsigned idx = (tMid < 2)  & prevC0 ? tMid : tMid + 2; // separate statistics for non-MS characters of zero run
      int hVal = VAL_RANGE/2;
      enum { QH_DECODE_NONE, QH_DECODE_ZERO, QH_DECODE_SIGN, QH_DECODE_DIFF_1, QH_DECODE_DIFF_N };
      int qhDecodeState = QH_DECODE_NONE;
      unsigned qhDecodeDiff, qhDecodeRa, qhDecodeUp;
      int cntr = cntrs[idx];
      if (cntr != 0) {
        ++cntr;
        if (cntr == ARITH_CODER_CNT_MAX)
          cntr = 0;
        cntrs[idx] = cntr;
        hVal = currH[idx];
      } else {
        qhDecodeState = QH_DECODE_ZERO;
      }

      unsigned b;
      for (;;) {
        b = 0;
        if (hVal != 0) {
          b = 1;
          if (hVal != VAL_RANGE) {
            b = (value >= range*(VAL_RANGE-hVal));
            // keep decoder in sync with encoder
            uint64_t cLo = b == 0 ? 0                : VAL_RANGE - hVal;
            uint64_t cHi = b == 0 ? VAL_RANGE - hVal : VAL_RANGE;
            value -= range * cLo;
            range *= (cHi-cLo);  // at this point range is scaled by 2**64 - the same scale as value
            uint64_t nxtRange = range >> RANGE_BITS;

            if (nxtRange <= MIN_RANGE) {
              #ifdef ENABLE_PERF_COUNT
              ++pModel->m_renormalizationCount;
              #endif
              if (srclen < 8) {
                if (!useTmpbuf && srclen > 0) {
                  memcpy(tmpbuf, src, srclen);
                  src = tmpbuf;
                  useTmpbuf = true;
                }
              }

              uint32_t threeOctets =
                (uint32_t(src[2])       ) |
                (uint32_t(src[1]) << 1*8) |
                (uint32_t(src[0]) << 2*8);
              value = (value << 24) + threeOctets;
              src    += 3;
              srclen -= 3;
              nxtRange = range << (24 - RANGE_BITS);
              if (srclen < -7)
                return dst-dst0;
              if (value > ((nxtRange<<RANGE_BITS)-1)) {
                return -103; // should not happen
              }
            }
            range = nxtRange;
          }
        }
        if (qhDecodeState == QH_DECODE_NONE)
          break;

        unsigned prevQhVal = prevQh[idx];
        unsigned qhVal = prevQhVal;
        switch (qhDecodeState) {
          case QH_DECODE_ZERO:
          {
            if (b == 0) {
              qhDecodeState = QH_DECODE_NONE;
            } else {
              qhDecodeState = QH_DECODE_SIGN;
              if (prevQhVal == 0) {
                qhDecodeUp = 1;
                qhDecodeState = QH_DECODE_DIFF_1;
              } else if (prevQhVal == ARITH_CODER_QH_SCALE) {
                qhDecodeUp = 0;
                qhDecodeState = QH_DECODE_DIFF_1;
              }
            }
          } break;

          case QH_DECODE_SIGN:
          {
            qhDecodeUp = b;
            unsigned ra = qhDecodeUp ? ARITH_CODER_QH_SCALE - prevQhVal : prevQhVal;
            qhDecodeState = QH_DECODE_DIFF_1;
            if (ra == 1) {
              qhVal = qhDecodeUp ? prevQhVal+1 : prevQhVal-1;
              qhDecodeState = QH_DECODE_NONE;
            }
          } break;

          case QH_DECODE_DIFF_1:
          {
            qhDecodeDiff = 2;
            qhDecodeState = QH_DECODE_DIFF_N;
            if (b == 0) {
              qhVal = qhDecodeUp ? prevQhVal+1 : prevQhVal-1;
              qhDecodeState = QH_DECODE_NONE;
            } else {
              unsigned ra = qhDecodeUp ? ARITH_CODER_QH_SCALE - prevQhVal : prevQhVal;
              qhDecodeRa = ra - 1;
              if (qhDecodeRa == 1) {
                qhVal = qhDecodeUp ? prevQhVal+2 : prevQhVal-2;
                qhDecodeState = QH_DECODE_NONE;
              }
            }
          } break;

          default:
          {
            unsigned halfRa = qhDecodeRa / 2;
            if (b == 0) {
              qhDecodeRa = halfRa;
            } else {
              qhDecodeDiff += halfRa;
              qhDecodeRa -= halfRa;
            }
            if (qhDecodeRa==1) {
              qhVal = qhDecodeUp ? prevQhVal+qhDecodeDiff : prevQhVal-qhDecodeDiff;
              qhDecodeState = QH_DECODE_NONE;
            }
          } break;
        }

        if (qhDecodeState == QH_DECODE_NONE) {
          // printf("%3d %2d qh\n", idx, qhVal);
          prevQh[idx] = qhVal;
          currH[idx] = hVal = qh2h_tab[qhVal];
          cntrs[idx] = 1;
        }
      }
      tLo = (b == 0) ? tLo  : tMid + 1;
      tHi = (b == 0) ? tMid : tHi;
    } while (tLo != tHi);
    unsigned c = tLo;
    // arithmetic decode of one character done
    prevC0 = (c < 2);

    // RLE decode
    if (c < 2) {
      // zero run
      uint32_t delta_rl = c+rlInc;
      rlInc += delta_rl;

      if (delta_rl > uint32_t(dstlen))
        return -26; // zero run too long (A)

      // insert zero run
      int cm0 = mtf_t[0];
      histogram[cm0] += delta_rl;
      memset(dst, cm0, delta_rl);
      dst    += delta_rl;
      dstlen -= delta_rl;
    } else {
      rlInc = 1;
      // MTF decode
      int mtfC = c - 1;
      int dstC = mtf_t[mtfC];
      histogram[dstC] += 1;
      *dst++ = dstC;
      dstlen -= 1;
      // update move-to-front decoder table
      memmove(&mtf_t[1], &mtf_t[0], mtfC);
      mtf_t[0] = dstC;
    }
  }

  return dst - dst0;
}

// Arithmetic decode followed by RLE and MTF decode
int arithmetic_decode(
  uint8_t*       dst,
  int            dstlen,
  int32_t        histogram[256],
  const uint8_t* src,
  int            srclen,
  int*           pInfo)
{
  if (pInfo) {
    pInfo[0] = 0;
    pInfo[1] = 0;
    pInfo[2] = 0;
    pInfo[3] = 0;
    pInfo[4] = 0;
  }
  if (srclen < 1)
    return -1;

  memset(histogram, 0, sizeof(histogram[0])*256);
  int textlen = decode(dst, dstlen, histogram, src, srclen);
  return textlen;
}
