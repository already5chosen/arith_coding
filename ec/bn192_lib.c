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

void bn192_add_rshift1(bn_t result, const bn_t a, const bn_t b)
{
  bn_word_t av = a[0];
  bn_word_t r = av + b[0];
  bn_word_t carry = r < av;
  for (int i = 1; i < ECDSA_NWORDS; ++i) {
    av = a[i];
    bn_word_t sum1 = av + b[i];
    bn_word_t sum2 = sum1 + carry;
    result[i-1] = (r >> 1) | (sum2 << (BN192LIB_BITS_PER_WORD-1));
    carry = (sum1 < av) | (sum2 < sum1);
    r = sum2;
  }
  result[ECDSA_NWORDS-1] = (r >> 1) | (carry << (BN192LIB_BITS_PER_WORD-1));
}

// bn192_mod_add_quick
// result = (a + b) mod m where a < m, b < m, m > 2^191
void bn192_mod_add_quick(bn_t result, const bn_t a, const bn_t b, const bn_t m)
{
  bn_word_t carry  = 0; // carry of a[]+b[]
  bn_word_t borrow = 0; // borrow of a[]+b[]-m[]
  bn_t result2; // a[]+b[]-m[]
  for (int i = 0; i < ECDSA_NWORDS; ++i) {
    bn_word_t av = a[i];
    bn_word_t sum1 = av + b[i];
    bn_word_t sum2 = sum1 + carry;
    bn_word_t mv = m[i];
    result[i] = sum2;
    carry = (sum1 < av) | (sum2 < sum1);

    bn_word_t diff = sum2 - borrow;
    borrow = (sum2 < borrow) | (diff < mv);
    result2[i] = diff - mv;
  }
  if (carry || !borrow)
    memcpy(result, result2, sizeof(bn_t));
}

// bn192_mod_sub_quick
// result = (a - b) mod m where a < m, b < m, m > 2^191
void bn192_mod_sub_quick(bn_t result, const bn_t a, const bn_t b, const bn_t m)
{
  bn_word_t borrow = 0; // borrow of a[]-b[]
  bn_word_t carry  = 0; // carry  of a[]-b[]+m[]
  bn_t result2; // a[]-b[]+m[]
  for (int i = 0; i < ECDSA_NWORDS; ++i) {
    bn_word_t av = a[i];
    bn_word_t bv = b[i];
    bn_word_t mv = m[i];
    bn_word_t dif1 = av - borrow;
    bn_word_t dif2 = dif1 - bv;
    borrow = (av < borrow) | (dif1 < bv);

    bn_word_t sum1 = dif2 + mv;
    bn_word_t sum2 = sum1 + carry;
    carry = (sum1 < mv) | (sum2 < sum1);

    result[i]  = dif2;
    result2[i] = sum2;
  }
  if (borrow)
    memcpy(result, result2, sizeof(bn_t));
}

// bn192_mod_inverse - solves (a*x) mod n == 1, where n is a prime number
//
// a - number in range [2:n-1]
// n - prime number in range [2^191:2^192-1].
// An implementation is efficient when m is close to 2^192
//
void bn192_mod_inverse(bn_t result, const bn_t a, const bn_t n)
{
  //
  // An inverse is computed as pow(a, n-2) mod n
  // They say, it works due to Fermat's Little Theorem.
  // I don't know if this method is the fastest, but it
  // certainly is one of the simplest.
  //

  bn_t exp; // exp[] = n[2]-2
  bn_word_t sub_w = 2;
  for (int i = 0; i < ECDSA_NWORDS; ++i) {
    bn_word_t nw = n[i];
    exp[i] = nw - sub_w;
    sub_w  = nw < sub_w;
  }

  int nz = 0;
  bn_t prod;
  for (int wi = ECDSA_NWORDS-1; wi >= 0; --wi)  {
    bn_word_t exp_w = exp[wi];
    for (int bi = 0; bi < BN192LIB_BITS_PER_WORD; ++bi) {
      if (nz)
        bn192_nist_mod_192_sqr(prod, prod, n);
      if ((exp_w >> (BN192LIB_BITS_PER_WORD-1)) & 1) {
        if (nz)
          bn192_nist_mod_192_mul(prod, prod, a, n);
        else
          memcpy(prod, a, sizeof(bn_t));
        nz = 1;
      }
      exp_w += exp_w;
    }
  }
  memcpy(result, prod, sizeof(bn_t));
}

// bn192_mod_lshift1_quick
// result = (a + a) mod m where a < m, m > 2^191
void bn192_mod_lshift1_quick(bn_t result, const bn_t a, const bn_t m)
{
  bn_word_t carry  = 0; // carry of a[]+a[]
  bn_word_t borrow = 0; // borrow of a[]+a[]-m[]
  bn_t result2; // a[]+a[]-m[]
  for (int i = 0; i < ECDSA_NWORDS; ++i) {
    bn_word_t av = a[i];
    bn_word_t mv = m[i];
    bn_word_t sum = av + av + carry;
    carry = av >> (BN192LIB_BITS_PER_WORD-1);

    bn_word_t diff = sum - borrow;
    borrow = (sum < borrow) | (diff < mv);
    result[i] = sum;
    result2[i] = diff - mv;
  }
  if (carry || !borrow)
    memcpy(result, result2, sizeof(bn_t));
}

void bn192_mod_lshift_quick(bn_t result, const bn_t a, int n, const bn_t m)
{
  while (n > 0) {
    bn192_mod_lshift1_quick(result, a, m);
    a = result;
    --n;
  }
}

void bn192_nist_mod_192(bn_t result, const bn_t a, const bn_t m)
{
  memcpy(result, a, sizeof(bn_t));
  if (bn192_ucmp(a, m) >= 0) {
    // subtract
    bn_word_t borrow = 0;
    for (int i = 0; i < ECDSA_NWORDS; ++i) {
      bn_word_t av = a[i];
      bn_word_t bv = m[i];
      bn_word_t dif1 = av - borrow;
      bn_word_t dif2 = dif1 - bv;
      borrow = (av < borrow) | (dif1 < bv);
      result[i] = dif2;
    }
  }
}

#ifdef NONCONV_MULTIPLY

// mulw_core - multiply bn_t number by word
// result - a * b
static inline void bn192_mulw_core(bn_word_t result[ECDSA_NWORDS+1], const bn_t a, bn_word_t b)
{
#ifdef NO_MULX
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
  uint64_t mx = (uint64_t)a[0] * b;
  result[0] = (bn_word_t)mx;
  bn_word_t mh =  (uint32_t)(mx>>32);
  for (int i = 1; i < ECDSA_NWORDS; ++i) {
    mx = (uint64_t)a[i] * b + mh;
    result[i] = (bn_word_t)mx;
    mh = (uint32_t)(mx>>32);
  }
  result[ECDSA_NWORDS] = mh;
#endif
}

// add_core
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

// sub_core
// result = result - a
// return - borrow out
static inline bn_word_t bn192_sub_core(bn_t result, const bn_t a)
{
  uint64_t dif = result[0];
  dif -= a[0];
  result[0] = (bn_word_t)dif;
  bn_word_t borrow = (uint32_t)(dif>>32) & 1;
  for (int i = 1; i < ECDSA_NWORDS; ++i) {
    dif = result[i];
    dif -= a[i];
    dif -= borrow;
    result[i] = (bn_word_t)dif;
    borrow = (uint32_t)(dif>>32) & 1;
  }
  return borrow;
}

// bn192_nist_mod_192_mul
// result = (a * b) mod m where a < m, b < m, m > 2^191
// An implementation is efficient when m is close to 2^192
void bn192_nist_mod_192_mul(bn_t result, const bn_t a, const bn_t b, const bn_t m)
{
  bn_word_t acc[ECDSA_NWORDS*2] = {0};
  for (int i = ECDSA_NWORDS-1; i >= 0; --i) {
    bn_word_t mx[ECDSA_NWORDS+1];
    bn192_mulw_core(mx, a, b[i]);
    acc[i] = mx[0];
    bn_word_t carry = bn192_add_core(&acc[i+1], &mx[1]);
    while (carry)
      carry -= bn192_sub_core(&acc[i+1], m);
    bn_word_t msw = acc[i+ECDSA_NWORDS];
    if (msw > 1) {
      bn192_mulw_core(mx, m, msw);
      msw -= mx[ECDSA_NWORDS];
      msw -= bn192_sub_core(&acc[i], mx);
    }
    while (msw)
      msw -= bn192_sub_core(&acc[i], m);
  }
  // copy and conditionally last subtract
  bn192_nist_mod_192(result, acc, m);
}

#else

// mulx_core - result = a * b
static void bn192_mulx_core(bn_word_t result[ECDSA_NWORDS*2], const bn_t a, const bn_t b)
{
#ifdef NO_MULX
  uint32_t mxc = 0;
  for (int ri = 0; ri < ECDSA_NWORDS*4-1; ++ri) {
    int ai_beg = 0;
    int bi_beg = ri;
    if (bi_beg >= ECDSA_NWORDS*2) {
      bi_beg = ECDSA_NWORDS*2-1;
      ai_beg = ri - bi_beg;
    }
    uint32_t acc_l = mxc;
    uint32_t acc_h = 0;
    for (int ai = ai_beg, bi = bi_beg; ai <= bi_beg; ++ai, --bi) {
      uint32_t mx = (uint32_t)((const uint16_t*)a)[ai] * ((const uint16_t*)b)[bi];
      acc_l += mx & 0xFFFF;
      acc_h += mx >> 16;
    }
    ((uint16_t*)result)[ri] = (uint16_t)acc_l;
    mxc = acc_h + (acc_l>>16);
  }
  ((uint16_t*)result)[ECDSA_NWORDS*4-1] = (uint16_t)mxc;
#else
  uint64_t mxc = 0;
  for (int ri = 0; ri < ECDSA_NWORDS*2-1; ++ri) {
    int ai_beg = 0;
    int bi_beg = ri;
    if (bi_beg >= ECDSA_NWORDS) {
      bi_beg = ECDSA_NWORDS-1;
      ai_beg = ri - bi_beg;
    }
    uint64_t acc_l = mxc;
    uint64_t acc_h = 0;
    for (int ai = ai_beg, bi = bi_beg; ai <= bi_beg; ++ai, --bi) {
      uint64_t mx = (uint64_t)a[ai] * b[bi];
      acc_l += (uint32_t)mx;
      acc_h += (uint32_t)(mx>>32);
    }
    result[ri] = (uint32_t)acc_l;
    mxc = acc_h + (uint32_t)(acc_l>>32);
  }
  result[ECDSA_NWORDS*2-1] = (uint32_t)mxc;
#endif
}

// bn192_nist_mod_192_mul
// result = (a * b) mod m where a < m, b < m, m > 2^191
// An implementation is efficient when m is close to 2^192
void bn192_nist_mod_192_mul(bn_t result, const bn_t a, const bn_t b, const bn_t m)
{
  bn_word_t aXb[ECDSA_NWORDS*2]; // buffer for a[]*b[];
  bn192_mulx_core(aXb, a, b);    // multiply
  while (!bn192_is_zero(&aXb[ECDSA_NWORDS])) {
    bn_word_t dXm[ECDSA_NWORDS*2]; // buffer for aXb_h[]*m[];
    // multiply and subtract
    bn192_mulx_core(dXm, &aXb[ECDSA_NWORDS], m);
    bn_word_t borrow = 0;
    for (int i = 0; i < ECDSA_NWORDS*2; ++i) {
      bn_word_t av = aXb[i];
      bn_word_t bv = dXm[i];
      bn_word_t dif1 = av - borrow;
      bn_word_t dif2 = dif1 - bv;
      borrow = (av < borrow) | (dif1 < bv);
      aXb[i] = dif2;
    }
  }

  // copy and conditionally last subtract
  bn192_nist_mod_192(result, aXb, m);
}
#endif

void bn192_nist_mod_192_sqr(bn_t result, const bn_t a, const bn_t m)
{
  bn192_nist_mod_192_mul(result, a, a, m);
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
