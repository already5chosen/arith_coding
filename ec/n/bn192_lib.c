/*
 Replacement for following BN_xxx functions

BN_add
BN_copy
BN_is_odd
BN_is_one
BN_is_zero
BN_mod_add_quick
BN_mod_exp
BN_mod_inverse
BN_mod_lshift1_quick
BN_mod_lshift_quick
BN_mod_mul
BN_mod_sqr
BN_mod_sub_quick
BN_mul
BN_nist_mod_192
BN_nnmod
BN_rshift1
BN_sqr
BN_ucmp
BN_zero
*/

#include <string.h>
#include "bn192_lib.h"

enum {
  BNLIB_NBYTES = 24,
};


void bn192_cleanup(void)
{
}

void bn192_copy(bn_word_t* result, const bn_word_t* src)
{
  memcpy(result, src, sizeof(bn_t));
}

void bn192_zero(bn_t result)
{
  memset(result, 0, sizeof(bn_t));
}

void bn192_one(bn_t result)
{
  memset(result, 0, sizeof(bn_t));
  result[0] = 1;
}

int bn192_is_odd(const bn_t a)
{
  return a[0] & 1;
}

int bn192_is_one(const bn_t a)
{
  if (a[0] != 1)
    return 0; // is not one
  for (int i = 1; i < ECDSA_NWORDS; ++i)
    if (a[i] != 0)
      return 0; // is not one
  return 1; // is one
}

int bn192_is_zero(const bn_t a)
{
  for (int i = 0; i < ECDSA_NWORDS; ++i)
    if (a[i] != 0)
      return 0; // is not zero
  return 1; // is zero
}

int bn192_ucmp(const bn_t a, const bn_t b)
{
  for (int i = ECDSA_NWORDS-1; i >= 0; --i) {
    bn_word_t av = a[i];
    bn_word_t bv = b[i];
    if (av != bv)
      return av < bv ? -1 : 1; // is non equal
  }
  return 0; // is equal
}

// return 1 when a+b_n >= 2^192 otherwise return 0
int bn192_is_ge_n(const bn_t a, const bn_ofn_t b_n)
{
  for (int i = ECDSA_NWORDS-1; i >= ECDSA_OFn_NWORDS; --i)
    if (a[i] != (bn_word_t)-1)
      return 0;
  bn_word_t carry = 0;
  for (int i = 0; i < ECDSA_OFn_NWORDS; ++i) {
    bn_word_t sum1 = a[i]   + carry;
    bn_word_t sum2 = b_n[i] + sum1;
    carry = (sum1 < carry) | (sum2 < sum1);
  }
  return carry;
}

void bn192_add_rshift1_n(bn_t result, const bn_t asrc, const bn_ofn_t b_n)
{
  bn_word_t a = asrc[0];
  bn_word_t b = b_n[0];
  bn_word_t r = a - b;
  bn_word_t borrow = a < b;
  r ^= (r & 1);
  for (int i = 1; i < ECDSA_OFn_NWORDS; ++i) {
    a = asrc[i];
    b = b_n[i];
    bn_word_t dif1 = a - borrow;
    bn_word_t dif = dif1 - b;
    borrow = (a < borrow) | (dif1 < b);
    bn_word_t lsb = dif & 1;
    r |= lsb;
    result[i-1] = (r >> 1) | (r << (BN192LIB_BITS_PER_WORD-1));
    r = dif ^ lsb;
  }
  for (int i = ECDSA_OFn_NWORDS; i < ECDSA_NWORDS; ++i) {
    a = asrc[i];
    bn_word_t dif = a - borrow;
    borrow = (a < borrow);
    bn_word_t lsb = dif & 1;
    r |= lsb;
    result[i-1] = (r >> 1) | (r << (BN192LIB_BITS_PER_WORD-1));
    r = dif ^ lsb;
  }
  r |= borrow ^ 1;
  result[ECDSA_NWORDS-1] = (r >> 1) | (r << (BN192LIB_BITS_PER_WORD-1));
}

// bn192_sub_n_core
// acc = acc - (2^192-asrc)
// return - borrow out
static bn_word_t bn192_sub_n_core(bn_t acc, const bn_ofn_t asrc)
{
  bn_word_t a = asrc[0];
  bn_word_t r = acc[0] + a;
  acc[0] = r;
  bn_word_t carry = r < a;
  for (int i = 1; i < ECDSA_OFn_NWORDS; ++i) {
    a = asrc[i];
    bn_word_t sum1 = acc[i] + carry;
    bn_word_t sum2 = sum1 + a;
    acc[i] = sum2;
    carry  = (sum1 < carry) | (sum2 < a);
  }
  for (int i = ECDSA_OFn_NWORDS; i < ECDSA_NWORDS; ++i) {
    uint32_t r = acc[i] + carry;
    acc[i] = r;
    carry = (r < carry);
  }
  return carry ^ 1;
}

// bn192_mod_add_quick_n
// result = (a + b) mod m where m=2^192-m_n, a < m, b < m, 0 < m_n < 2^96
void bn192_mod_add_quick_n(bn_t result, const bn_t a, const bn_t b, const bn_ofn_t m_n)
{
  bn_word_t carry = 0; // carry of a[]+b[]
  for (int i = 0; i < ECDSA_NWORDS; ++i) {
    bn_word_t sum1 = a[i] + carry;
    bn_word_t sum2 = b[i] + sum1;
    result[i] = sum2;
    carry = (sum1 < carry) | (sum2 < sum1);
  }
  if (!carry)
    carry = bn192_is_ge_n(result, m_n);
  if (carry)
    bn192_sub_n_core(result, m_n);
}

// bn192_mod_sub_quick_n
// result = (a - b) mod m where m=2^192-m_n, a < m, b < m, 0 < m_n < 2^96
void bn192_mod_sub_quick_n(bn_t result, const bn_t a, const bn_t b, const bn_ofn_t m_n)
{
  bn_word_t borrow = 0; // borrow of a[]-b[]
  for (int i = 0; i < ECDSA_NWORDS; ++i) {
    bn_word_t av = a[i];
    bn_word_t bv = b[i];
    bn_word_t dif = av - borrow;
    result[i] = dif - bv;
    borrow = (av < borrow) | (dif < bv);
  }
  if (borrow) {
    borrow = 0;
    for (int i = 0; i < ECDSA_OFn_NWORDS; ++i) {
      bn_word_t rv = result[i];
      bn_word_t mv = m_n[i];
      bn_word_t dif = rv - borrow;
      result[i] = dif - mv;
      borrow = (rv < borrow) | (dif < mv);
    }
    if (borrow) {
      for (int i = ECDSA_OFn_NWORDS; result[i] == 0; ++i) {
        bn_word_t rv = result[i];
        result[i] = rv -1;
        if (rv)
          break;
      }
    }
  }
}

// bn192_mod_inverse_n - solves (a*x) mod n == 1, where n is a prime number
//
// a   - number in range [2:n-1]
// n_n - inverse prime. n = 2^192-n_n is a prime number in range [2^192-2^96+1:2^192-1].
//
void bn192_mod_inverse_n(bn_t result, const bn_t a, const bn_ofn_t n_n)
{
  //
  // An inverse is computed as pow(a, n-2) mod n
  // They say, it works due to Fermat's Little Theorem.
  // I don't know if this method is the fastest, but it
  // certainly is one of the simplest.
  //

  bn_ofn_t not_exp; // exp[] = n_n+1 = ~(-n-2)
  bn_word_t add_w = 1;
  for (int i = 0; i < ECDSA_OFn_NWORDS; ++i) {
    bn_word_t nw = n_n[i] + add_w;
    not_exp[i] = nw;
    add_w  = nw < add_w;
  }

  int nz = 0;
  bn_t prod;
  for (int wi = ECDSA_NWORDS-1; wi >= 0; --wi)  {
    bn_word_t exp_w = wi < ECDSA_OFn_NWORDS ? not_exp[wi] : 0;
    for (int bi = 0; bi < BN192LIB_BITS_PER_WORD; ++bi) {
      if (nz)
        bn192_nist_mod_192_sqr_n(prod, prod, n_n);
      if (((exp_w >> (BN192LIB_BITS_PER_WORD-1)) & 1)==0) {
        if (nz)
          bn192_nist_mod_192_mul_n(prod, prod, a, n_n);
        else
          memcpy(prod, a, sizeof(bn_t));
        nz = 1;
      }
      exp_w += exp_w;
    }
  }
  memcpy(result, prod, sizeof(bn_t));
}



// bn192_mod_lshift1_quick_n
// result = (a + a) mod m, where a < m, m = 2^192-m_n, m_n < 2^96
void bn192_mod_lshift1_quick_n(bn_t result, const bn_t a, const bn_ofn_t m_n)
{
  bn192_mod_add_quick_n(result, a, a, m_n);
}

void bn192_mod_lshift_quick_n(bn_t result, const bn_t a, int n, const bn_ofn_t m_n)
{
  while (n > 0) {
    bn192_mod_lshift1_quick_n(result, a, m_n);
    a = result;
    --n;
  }
}

#ifndef MULTIPLICATION_ALGO
  #define MULTIPLICATION_ALGO 0
#endif

// Multiplication Algorithm variants
// MULTIPLICATION_ALGO  a*b           modulo reduction  Use mulx  Use mul
// 0   N/A
// 1                    convolution   word-by-word      Yes       Yes
// 2   N/A
// 3                    word-by-word  word-by-word      Yes       Yes
// 4   N/A
// 5                    convolution   word-by-word      No        Yes
// 6   N/A
// 7                    word-by-word  word-by-word      No        Yes
// 8   N/A
// 9   N/A
// 10  N/A
// 11                   bit-by-bit    bit-by-bit        No        No
// 12  N/A
// 13                   convolution   word-by-word      No        No
// 14  N/A
// 15                   word-by-word  word-by-word      No        No


#if (MULTIPLICATION_ALGO & 2)==0
// mulx_core - result = a * b
static void bn192_mulx_core(bn_word_t result[ECDSA_NWORDS*2], const bn_t a, const bn_t b)
{
#if (MULTIPLICATION_ALGO & 12)==12
  uint32_t acc0 = 0, acc1 = 0, acc2 = 0, acc3 = 0, acc4 = 0, acc5 = 0;
  for (int ri = 0;;) {
    int aiBeg = ri < ECDSA_NWORDS ? 0 : ri - ECDSA_NWORDS + 2;
    do {
      uint32_t acc6, acc7, acc8, acc9;
      if (ri != ECDSA_NWORDS*2-2) {
        int ai = aiBeg;
        int bi = ri - ai;
        int ni = bi - ai + 1;
        acc6 = acc7 = acc8 = acc9 = 0;
        do {
          uint32_t aw0 = a[ai+0];
          uint32_t aw1 = a[ai+1];
          uint32_t bw  = b[bi];
          uint32_t a0 = aw0 & 0xFFFF;
          uint32_t a2 = aw0 >> 16;
          uint32_t a4 = aw1 & 0xFFFF;
          uint32_t a6 = aw1 >> 16;
          int bit_i = 8;
          goto loop_entry;
          do {
            a0 += a0;
            a2 += a2;
            a4 += a4;
            a6 += a6;
            loop_entry:
            if (bw & 1) {
              acc0 += a0;
              acc2 += a2;
              acc4 += a4;
              acc6 += a6;
            }
            if (bw & (1u << 8)) {
              acc1 += a0;
              acc3 += a2;
              acc5 += a4;
              acc7 += a6;
            }
            if (bw & (1u << 16)) {
              acc2 += a0;
              acc4 += a2;
              acc6 += a4;
              acc8 += a6;
            }
            if (bw & (1u << 24)) {
              acc3 += a0;
              acc5 += a2;
              acc7 += a4;
              acc9 += a6;
            }
            bw >>= 1;
            --bit_i;
          } while (bit_i != 0);

          ai += 2;
          bi -= 2;
          ni -= 2;
        } while (ni > 0);
      }

      acc4 += acc1 >> 24; acc1 <<=  8; acc0 += acc1; acc4 += acc0 < acc1;
      acc4 += acc2 >> 16; acc2 <<= 16; acc0 += acc2; acc4 += acc0 < acc2;
      acc4 += acc3 >>  8; acc3 <<= 24; acc0 += acc3; acc4 += acc0 < acc3;

      result[ri] = acc0;
      if (ri == ECDSA_NWORDS*2-2)
        goto done;

      acc0 = acc4;
      acc1 = acc5;
      acc2 = acc6;
      acc3 = acc7;
      acc4 = acc8;
      acc5 = acc9;

      ++ri;
    } while (ri & 1);
  }
  done:
  result[ECDSA_NWORDS*2-1] = acc4 + (acc5 << 8);
#else
  uint64_t mxc = 0;
  for (int ai_beg_inc = 0; ai_beg_inc < 2; ++ai_beg_inc) {
    int ri     = ai_beg_inc==0 ? 0              : ECDSA_NWORDS-1;
    int ri_end = ai_beg_inc==0 ? ECDSA_NWORDS-1 : ECDSA_NWORDS*2-1;
    for (int ai_beg = 0; ri < ri_end; ai_beg += ai_beg_inc, ++ri) {
#if (MULTIPLICATION_ALGO & 4)==0
      uint64_t acc_l = mxc;
      uint64_t acc_h = 0;
      for (int ai = ai_beg; ai <= ri - ai_beg; ++ai) {
        uint64_t mx = (uint64_t)a[ai] * b[ri-ai];
        acc_l += (uint32_t)mx;
        acc_h += (uint32_t)(mx>>32);
      }
      result[ri] = (uint32_t)acc_l;
      mxc = acc_h + (uint32_t)(acc_l>>32);
#else
      uint32_t acc_l = (uint32_t)mxc;
      uint64_t acc_h = (uint32_t)(mxc >> 32);
      uint64_t acc_m = 0;
      for (int ai = ai_beg; ai <= ri - ai_beg; ++ai) {
        uint32_t ax = a[ai];
        uint32_t bx = b[ri-ai];
        uint32_t ah = ax >> 16, al = ax & 0xFFFF;
        uint32_t bh = bx >> 16, bl = bx & 0xFFFF;
        uint32_t ll = al * bl;
        uint32_t hh = ah * bh;
        uint32_t lh = al * bh;
        uint32_t hl = ah * bl;
        acc_l += ll;
        hh += (acc_l < ll); // no carry out from hh, because hh can't be UINT32_MAX
        acc_h += hh;
        acc_m += lh;
        acc_m += hl;
      }
      acc_m = (acc_m << 16) + acc_l;
      result[ri] = (uint32_t)acc_m;
      mxc = acc_h + (uint32_t)(acc_m>>32);
#endif
    }
  }
  result[ECDSA_NWORDS*2-1] = (uint32_t)mxc;
#endif
}
#endif

#if (MULTIPLICATION_ALGO & 3)==3 && (MULTIPLICATION_ALGO & 12)!=8
// mulw_core - multiply bn_t number by word
// result - a * b
static
#if (MULTIPLICATION_ALGO & 12) != 12
inline
#endif
void bn192_mulw_core(bn_word_t result[ECDSA_NWORDS+1], const bn_t a, bn_word_t b)
{
#if (MULTIPLICATION_ALGO & 12)==0
  // mul and mulx
  uint64_t mx = (uint64_t)a[0] * b;
  result[0] = (bn_word_t)mx;
  bn_word_t mh =  (uint32_t)(mx>>32);
  for (int i = 1; i < ECDSA_NWORDS; ++i) {
    mx = (uint64_t)a[i] * b + mh;
    result[i] = (bn_word_t)mx;
    mh = (uint32_t)(mx>>32);
  }
  result[ECDSA_NWORDS] = mh;
#elif  (MULTIPLICATION_ALGO & 12)==4
  // mul, no mulx
  uint32_t ax = a[0];
  uint32_t bh = b >> 16;
  uint32_t bl = b & 0xFFFF;
  uint32_t ml = ax * b;
  uint32_t ah = ax >> 16;
  uint32_t al = ax & 0xFFFF;
  uint32_t mh = ah * bh;
  result[0] = ml;
  uint32_t ahXbl = ah * bl;
  uint32_t alXbh = al * bh;
  mh += ahXbl >> 16;
  mh += alXbh >> 16;
  uint32_t hl = (ahXbl & 0xFFFF) + (alXbh & 0xFFFF);
  mh += hl >> 16;
  mh += (hl << 16) > ml;
  for (int i = 1; i < ECDSA_NWORDS; ++i) {
    ax = a[i];
    ml = ax * b;
    ah = ax >> 16;
    al = ax & 0xFFFF;
    uint32_t mhi = ah * bh;
    ahXbl = ah * bl;
    alXbh = al * bh;
    mhi += ahXbl >> 16;
    mhi += alXbh >> 16;
    hl = (ahXbl & 0xFFFF) + (alXbh & 0xFFFF);
    mhi += hl >> 16;
    mhi += (hl << 16) > ml;
    ml += mh;
    result[i] = ml;
    mh = mhi + (mh > ml);
  }
  result[ECDSA_NWORDS] = mh;
#else
  // no mul, no mulx
  uint32_t mh = 0;
  for (int i = 0; i < ECDSA_NWORDS; i += 2) {
    // process 2 words per iteration
    const uint32_t MSK32  = (uint32_t)-1;
    const uint32_t MSK8   = MSK32 >> (32-8);
    const uint32_t MSK24  = MSK32 >> (32-24);
    uint32_t a0 = a[i+0];
    uint32_t a1 = a[i+1];
    // split a1:a0 into 24-16-24-bit sub-words
    uint32_t a00 = a0 & MSK24;
    uint32_t a24 = (a0 >> 24) | ((a1  & MSK8) << 8);
    uint32_t a40 = a1 >> 8;
    uint32_t bb = b;
    uint32_t m00=0, m08=0, m16=0;
    uint32_t m24=0, m32=0, m40=0;
    uint32_t m48=0, m56=0, m64=0;
    int bi = 8;
    goto entry;
    do {
      a00 += a00;
      a24 += a24;
      a40 += a40;
      entry:
      if (bb & 1) {
        m00 += a00;
        m24 += a24;
        m40 += a40;
      }
      if (bb & ((uint32_t)1<<8)) {
        m08 += a00;
        m32 += a24;
        m48 += a40;
      }
      if (bb & ((uint32_t)1<<16)) {
        m16 += a00;
        m40 += a24;
        m56 += a40;
      }
      if (bb & ((uint32_t)1<<24)) {
        m24 += a00;
        m48 += a24;
        m64 += a40;
      }
      bb >>= 1;
    } while (--bi);

    m00 += mh;
    m32 += m00 < mh;

    m00 += (m08 <<  8); m32 += m00 < (m08 <<  8); m32 += (m08 >> 24);
    m00 += (m16 << 16); m32 += m00 < (m16 << 16); m32 += (m16 >> 16);
    m00 += (m24 << 24); m32 += m00 < (m24 << 24); m32 += (m24 >>  8);
    m32 += (m40 <<  8); m64 += m32 < (m40 <<  8); m64 += (m40 >> 24);
    m32 += (m48 << 16); m64 += m32 < (m48 << 16); m64 += (m48 >> 16);
    m32 += (m56 << 24); m64 += m32 < (m56 << 24); m64 += (m56 >>  8);

    result[i+0] = m00;
    result[i+1] = m32;
    mh = m64;
  }
  result[ECDSA_NWORDS] = mh;
#endif
}
#endif

#if (MULTIPLICATION_ALGO & 1)==1 && (MULTIPLICATION_ALGO & 12)!=8
// bn192_mulsubw_n_core - multiply bn_t number = (2^192-a_n) by word and subtract product
//                        from accumulator
// acc = acc - (2^192-a_n) * b
// return MS word of accumulator
static inline bn_word_t bn192_mulsubw_n_core(bn_word_t acc[ECDSA_NWORDS+1], const bn_ofn_t a_n, bn_word_t b)
{
  uint32_t mh = 0;
#if (MULTIPLICATION_ALGO & 4)==0
  for (int i = 0; i < ECDSA_OFn_NWORDS; ++i) {
    uint64_t mx = (uint64_t)a_n[i] * b + mh; // at most UINT32_MAX << 32
    uint32_t accw = acc[i] + (uint32_t)mx;
    acc[i] = accw;
    mh = (uint32_t)(mx>>32) + (accw < (uint32_t)mx); // no carry out because when MS word of mx=UINT32_MAX then LS word of mx = 0
  }
#elif (MULTIPLICATION_ALGO & 12)==12
  for (int i = 0; i < ECDSA_OFn_NWORDS; ++i) {
    uint32_t aw = a_n[i];
    uint64_t mx;
    if (aw < 3) {
      if (aw == 0) {
        mx = 0;
      } else if (aw == 1) {
        mx = b;
      } else { // aw == 2
        mx  = b;
        mx += b;
      }
      // uint32_t mxh, mxl;
    } else {
      // Optimize for size, because this case is less common
      // It only happens during calculations by modulo group->order
      uint32_t bl = b & 0xFFFF;
      uint32_t bh = b >> 16;
      uint32_t ll = 0, lh = 0, hl = 0, hh = 0;
      int bi = 16;
      goto loop_entry;
      do {
        bl += bl;
        bh += bh;
        loop_entry:
        if (aw & 1) {
          ll += bl;
          lh += bh;
        }
        if (aw & ((uint32_t)1 << 16)) {
          hl += bl;
          hh += bh;
        }
        aw >>= 1;
      } while (--bi);
      lh += hl;
      hh += (lh < hl) << 16;
      hh += lh >> 16;
      lh <<= 16;
      ll += lh;
      hh += ll < lh;
      mx = ((uint64_t)hh << 32) | ll;
    }
    mx += mh;
    uint32_t accw = acc[i] + (uint32_t)mx;
    acc[i] = accw;
    mh = (uint32_t)(mx>>32) + (accw < (uint32_t)mx); // no carry out because when MS word of mx=UINT32_MAX then LS word of mx = 0
  }
#else
  uint32_t bh = b >> 16;
  uint32_t bl = b & 0xFFFF;
  for (int i = 0; i < ECDSA_OFn_NWORDS; ++i) {
    uint32_t ax = a_n[i];
    uint32_t ll = ax * b;
    uint32_t ah = ax >> 16;
    uint32_t al = ax & 0xFFFF;
    uint32_t hh = ah * bh;
    uint32_t hl = ah * bl;
    uint32_t lh = al * bh;
    uint32_t ml = mh + ll;
    hh += (ml < ll);
    hl += lh;
    hh += (hl < lh) << 16;
    hh += hl >> 16;
    hh += (hl << 16) > ll;
    uint32_t accw = acc[i] + ml;
    acc[i] = accw;
    mh = hh + (accw < ml); // no carry out because when hh=UINT32_MAX then ml = 0
  }
#endif
  for (int i = ECDSA_OFn_NWORDS; i < ECDSA_NWORDS; ++i) {
    uint32_t accw = acc[i] + mh;
    acc[i] = accw;
    mh = (accw < mh); // carry out
  }
  mh = b - mh;
  return acc[ECDSA_NWORDS] - mh;
}
#endif

#if (MULTIPLICATION_ALGO & 3)==3
// bn192_add_core
// result = result + a
// return - carry out
static inline bn_word_t bn192_add_core(bn_t result, const bn_t a)
{
  uint64_t sum = result[0];
  sum += a[0];
  result[0] = (bn_word_t)sum;
  bn_word_t carry = (uint32_t)(sum>>32);
  for (int i = 1; i < ECDSA_NWORDS; ++i) {
    sum = result[i];
    sum += a[i];
    sum += carry;
    result[i] = (bn_word_t)sum;
    carry = (uint32_t)(sum>>32);
  }
  return carry;
}
#endif

#if MULTIPLICATION_ALGO==11
// bn192_dbl_core
// acc += acc
// return - carry out
static inline bn_word_t bn192_dbl_core(bn_t acc)
{
  bn_word_t carry = 0;
  for (int i = 0; i < ECDSA_NWORDS; ++i) {
    bn_word_t w = acc[i];
    bn_word_t w2 = w + w;
    acc[i] = w2 + carry;
    carry = w2 < w;
  }
  return carry;
}
#endif

// bn192_nist_mod_192_mul
// result = (a * b) mod (2^192-m) where a < (2^192-m), b < (2^192-m), m < 2^96
void bn192_nist_mod_192_mul_n(bn_t result, const bn_t a, const bn_t b, const bn_ofn_t m_n)
{
#if MULTIPLICATION_ALGO==11
  bn_t acc = {0};
  for (int i = ECDSA_NWORDS-1; i >= 0; --i) {
    bn_word_t bw = b[i];
    for (int bi = 0; bi < BN192LIB_BITS_PER_WORD; ++bi) {
      bn_word_t carry = bn192_dbl_core(acc); // acc *= 2;
      while (carry)
        carry -= bn192_sub_n_core(acc, m_n);
      if (bw & (1u << 31)) {
        carry = bn192_add_core(acc, a);
        while (carry)
          carry -= bn192_sub_n_core(acc, m_n);
      }
      bw += bw;
    }
  }
  // copy and conditionally last subtract
  bn192_nist_mod_192_n(result, acc, m_n);

#else

#if (MULTIPLICATION_ALGO & 3)==0
 #error "MULTIPLICATION_ALGO % 4 == 0 is not supported"
#endif

#if (MULTIPLICATION_ALGO & 3)==2
 #error "MULTIPLICATION_ALGO % 4 == 2 is not supported"
#endif

#if (MULTIPLICATION_ALGO & 3)==1
  bn_word_t acc[ECDSA_NWORDS*2]; // buffer for a[]*b[];
  bn192_mulx_core(acc, a, b);    // multiply
  for (int i = ECDSA_NWORDS-1; i >= 0; --i) {
    bn_word_t msw = acc[ECDSA_NWORDS+i];
    if (msw > 1)
      msw = bn192_mulsubw_n_core(&acc[i], m_n, msw);
    while (msw != 0)
      msw -= bn192_sub_n_core(&acc[i], m_n);
  }
  // copy and conditionally last subtract
  bn192_nist_mod_192_n(result, acc, m_n);
#endif

#if (MULTIPLICATION_ALGO & 3)==3
  bn_word_t acc[ECDSA_NWORDS*2] = {0};
  for (int i = ECDSA_NWORDS-1; i >= 0; --i) {
    bn_word_t mx[ECDSA_NWORDS+1];
    bn192_mulw_core(mx, a, b[i]);
    acc[i] = mx[0];
    bn_word_t carry = bn192_add_core(&acc[i+1], &mx[1]);
    while (carry)
      carry -= bn192_sub_n_core(&acc[i+1], m_n);
    bn_word_t msw = acc[ECDSA_NWORDS+i];
    if (msw > 1)
      msw = bn192_mulsubw_n_core(&acc[i], m_n, msw);
    while (msw != 0)
      msw -= bn192_sub_n_core(&acc[i], m_n);
  }
  // copy and conditionally last subtract
  bn192_nist_mod_192_n(result, acc, m_n);
#endif

#endif
}

void bn192_nist_mod_192_sqr_n(bn_t result, const bn_t a, const bn_ofn_t m_n)
{
  bn192_nist_mod_192_mul_n(result, a, a, m_n);
}

void bn192_rshift1(bn_t result, const bn_t a)
{
  bn_word_t msb = 0;
  for (int i = ECDSA_NWORDS-1; i >= 0; --i) {
    bn_word_t av = a[i];
    result[i] = (av >> 1) | msb;
    msb = av << (BN192LIB_BITS_PER_WORD-1);
  }
}

void bn192_nist_mod_192_n(bn_t result, const bn_t a, const bn_ofn_t m_n)
{
  if (a[ECDSA_NWORDS-1] == (bn_word_t)-1) {
    bn_t tmp;
    memcpy(tmp, a, sizeof(bn_t));
    if (!bn192_sub_n_core(tmp, m_n)) {
      memcpy(result, tmp, sizeof(bn_t));
      return;
    }
  }
  memcpy(result, a, sizeof(bn_t));
}
