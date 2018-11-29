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


#ifndef MULTIPLICATION_ALGO
  #define MULTIPLICATION_ALGO 1
#endif

// Multiplication Algorithm variants
// MULTIPLICATION_ALGO  a*b           modulo reduction  Use mulx  Use mul
// 0                    convolution   convolution       Yes       Yes
// 1                    convolution   word-by-word      Yes       Yes
// 2 - N/A
// 3                    word-by-word  word-by-word      Yes       Yes

#if (MULTIPLICATION_ALGO & 2) == 0
// mulx_core - result = a * b
static void bn192_mulx_core(bn_word_t result[ECDSA_NWORDS*2], const bn_t a, const bn_t b)
{
#if 1
  bn_doubleword_t mxc = 0;
  for (int ri = 0; ri < ECDSA_NWORDS*2-1; ++ri) {
    int ai_beg = 0;
    int bi_beg = ri;
    if (bi_beg >= ECDSA_NWORDS) {
      bi_beg = ECDSA_NWORDS-1;
      ai_beg = ri - bi_beg;
    }
    bn_doubleword_t acc_l = mxc;
    bn_doubleword_t acc_h = 0;
    for (int ai = ai_beg, bi = bi_beg; ai <= bi_beg; ++ai, --bi) {
      bn_doubleword_t mx = (bn_doubleword_t)a[ai] * b[bi];
      acc_l += (bn_word_t)mx;
      acc_h += (bn_word_t)(mx>>(BN192LIB_BYTES_PER_WORD*8));
    }
    result[ri] = (bn_word_t)acc_l;
    mxc = acc_h + (bn_word_t)(acc_l>>(BN192LIB_BYTES_PER_WORD*8));
  }
  result[ECDSA_NWORDS*2-1] = (bn_word_t)mxc;
#else
  unsigned __int128 acc = (unsigned __int128)a[0]*b[0];
  uint64_t w0 = (uint64_t)acc, w1x = (uint64_t)(acc >> 64), w1y, w1z;
  result[0] = w0;

  acc = (unsigned __int128)a[0]*b[1] + w1x; w0 = (uint64_t)acc; w1x = (uint64_t)(acc >> 64);
  acc = (unsigned __int128)a[1]*b[0] + w0;  w0 = (uint64_t)acc; w1y = (uint64_t)(acc >> 64);
  result[1] = w0;
  acc = (unsigned __int128)w1x + w1y;

  acc += (unsigned __int128)a[0]*b[2];      w0 = (uint64_t)acc; w1x = (uint64_t)(acc >> 64);
  acc =  (unsigned __int128)a[1]*b[1] + w0; w0 = (uint64_t)acc; w1y = (uint64_t)(acc >> 64);
  acc =  (unsigned __int128)a[2]*b[0] + w0; w0 = (uint64_t)acc; w1z = (uint64_t)(acc >> 64);
  result[2] = w0;
  acc = ((unsigned __int128)w1x + w1y) + w1z; w0 = (uint64_t)acc; w1x = (uint64_t)(acc >> 64);

  acc =  (unsigned __int128)a[1]*b[2] + w0; w0 = (uint64_t)acc; w1y = (uint64_t)(acc >> 64);
  acc =  (unsigned __int128)a[2]*b[1] + w0; w0 = (uint64_t)acc; w1z = (uint64_t)(acc >> 64);
  result[3] = w0;
  acc = ((unsigned __int128)w1x + w1y) + w1z;

  acc += (unsigned __int128)a[2]*b[2];      w0 = (uint64_t)acc; w1x = (uint64_t)(acc >> 64);
  result[4] = w0;
  result[5] = w1x;
#endif
}
#endif

#if 0 // MULTIPLICATION_ALGO == 1
// sqrx_core - result = a * a
static void bn192_sqrx_core(bn_word_t result[ECDSA_NWORDS*2], const bn_t a)
{
  bn_doubleword_t acc = (bn_doubleword_t)a[0]*a[0];
  bn_word_t w0 = (bn_word_t)acc, w1x = (bn_word_t)(acc >> (BN192LIB_BYTES_PER_WORD*8)), w1y, w1z;
  result[0] = w0;

  bn_doubleword_t mx = (bn_doubleword_t)a[0]*a[1];
  acc = mx + w1x; w0 = (bn_word_t)acc; w1x = (bn_word_t)(acc >> (BN192LIB_BYTES_PER_WORD*8));
  acc = mx + w0;  w0 = (bn_word_t)acc; w1y = (bn_word_t)(acc >> (BN192LIB_BYTES_PER_WORD*8));
  result[1] = w0;
  acc = (bn_doubleword_t)w1x + w1y;

  mx = (bn_doubleword_t)a[1]*a[1];
  acc += mx;      w0 = (bn_word_t)acc; w1x = (bn_word_t)(acc >> (BN192LIB_BYTES_PER_WORD*8));
  mx = (bn_doubleword_t)a[0]*a[2];
  acc =  mx + w0; w0 = (bn_word_t)acc; w1y = (bn_word_t)(acc >> (BN192LIB_BYTES_PER_WORD*8));
  acc =  mx + w0; w0 = (bn_word_t)acc; w1z = (bn_word_t)(acc >> (BN192LIB_BYTES_PER_WORD*8));
  result[2] = w0;
  acc = ((bn_doubleword_t)w1x + w1y) + w1z; w0 = (bn_word_t)acc; w1x = (bn_word_t)(acc >> (BN192LIB_BYTES_PER_WORD*8));

  mx = (bn_doubleword_t)a[1]*a[2];
  acc =  mx + w0; w0 = (bn_word_t)acc; w1y = (bn_word_t)(acc >> (BN192LIB_BYTES_PER_WORD*8));
  acc =  mx + w0; w0 = (bn_word_t)acc; w1z = (bn_word_t)(acc >> (BN192LIB_BYTES_PER_WORD*8));
  result[3] = w0;
  acc = ((bn_doubleword_t)w1x + w1y) + w1z;

  mx = (bn_doubleword_t)a[2]*a[2];
  acc += mx;      w0 = (bn_word_t)acc; w1x = (bn_word_t)(acc >> (BN192LIB_BYTES_PER_WORD*8));
  result[4] = w0;
  result[5] = w1x;
}
#endif

#if (MULTIPLICATION_ALGO & 1)==1
// sub_core
// acc = acc - a
// return - borrow out
static inline bn_word_t bn192_sub_core(bn_t acc, const bn_t a)
{
  bn_word_t rw = acc[0];
  bn_word_t aw = a[0];
  acc[0] = rw - aw;
  bn_word_t borrow = rw < aw;
  for (int i = 1; i < ECDSA_NWORDS; ++i) {
    rw = acc[i];
    aw = a[i];
    bn_word_t dif = rw - aw;
    acc[i] = dif - borrow;
    borrow = (rw < aw) | (dif < borrow);
  }
  return borrow;
}
#endif

#if MULTIPLICATION_ALGO == 3
// bn192_add_core
// result = result + a
// return - carry out
static inline bn_word_t bn192_add_core(bn_t result, const bn_t a)
{
  bn_doubleword_t sum = result[0];
  sum += a[0];
  result[0] = (bn_word_t)sum;
  bn_word_t carry = (bn_word_t)(sum>>(BN192LIB_BYTES_PER_WORD*8));
  for (int i = 1; i < ECDSA_NWORDS; ++i) {
    sum = result[i];
    sum += a[i];
    sum += carry;
    result[i] = (bn_word_t)sum;
    carry = (bn_word_t)(sum>>(BN192LIB_BYTES_PER_WORD*8));
  }
  return carry;
}
#endif

#if MULTIPLICATION_ALGO == 3
// mulw_core - multiply bn_t number by word
// result - a * b
static inline void bn192_mulw_core(bn_word_t result[ECDSA_NWORDS+1], const bn_t a, bn_word_t b)
{
  bn_doubleword_t mx = (bn_doubleword_t)a[0] * b;
  result[0] = (bn_word_t)mx;
  bn_word_t mh =  (bn_word_t)(mx>>(BN192LIB_BYTES_PER_WORD*8));
  for (int i = 1; i < ECDSA_NWORDS; ++i) {
    mx = (bn_doubleword_t)a[i] * b + mh;
    result[i] = (bn_word_t)mx;
    mh = (bn_word_t)(mx>>(BN192LIB_BYTES_PER_WORD*8));
  }
  result[ECDSA_NWORDS] = mh;
}
#endif

#if (MULTIPLICATION_ALGO & 1)==1
// bn192_mulsubw_core - multiply bn_t number by word and subtract product
//                      from accumulator
// acc = acc - a * b
// return MS word of accumulator
static inline bn_word_t bn192_mulsubw_core(bn_word_t acc[ECDSA_NWORDS+1], const bn_t a, bn_word_t b)
{
  bn_word_t mh = 0;
  for (int i = 0; i < ECDSA_NWORDS; ++i) {
    bn_doubleword_t mx = (bn_doubleword_t)a[i] * b + mh;
    bn_word_t accw = acc[i];
    acc[i] = accw - (bn_word_t)mx;
    mh = (bn_word_t)(mx>>(BN192LIB_BYTES_PER_WORD*8)) + (accw < (bn_word_t)mx); // no carry out because when MS word of mx=UINT64_MAX then LS word of mx = 0
  }
  return acc[ECDSA_NWORDS] - mh;
}
#endif

// bn192_nist_mod_192_mul
// result = (a * b) mod m where a < m, b < m, m > 2^191
// An implementation is efficient when m is close to 2^192
void bn192_nist_mod_192_mul(bn_t result, const bn_t a, const bn_t b, const bn_t m)
{
#if MULTIPLICATION_ALGO == 0
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
#endif

#if MULTIPLICATION_ALGO == 1
  bn_word_t acc[ECDSA_NWORDS*2]; // buffer for a[]*b[];
  bn192_mulx_core(acc, a, b);    // multiply
  for (int i = ECDSA_NWORDS-1; i >= 0; --i) {
    bn_word_t msw = acc[ECDSA_NWORDS+i];
    if (msw > 1)
      msw = bn192_mulsubw_core(&acc[i], m, msw);
    while (msw != 0)
      msw -= bn192_sub_core(&acc[i], m);
  }
  // copy and conditionally last subtract
  bn192_nist_mod_192(result, acc, m);
#endif

#if MULTIPLICATION_ALGO == 3
  bn_word_t acc[ECDSA_NWORDS*2] = {0};
  for (int i = ECDSA_NWORDS-1; i >= 0; --i) {
    bn_word_t mx[ECDSA_NWORDS+1];
    bn192_mulw_core(mx, a, b[i]);
    acc[i] = mx[0];
    bn_word_t carry = bn192_add_core(&acc[i+1], &mx[1]);
    while (carry)
      carry -= bn192_sub_core(&acc[i+1], m);
    bn_word_t msw = acc[ECDSA_NWORDS+i];
    if (msw > 1)
      msw = bn192_mulsubw_core(&acc[i], m, msw);
    while (msw != 0)
      msw -= bn192_sub_core(&acc[i], m);
  }
  // copy and conditionally last subtract
  bn192_nist_mod_192(result, acc, m);
#endif
}

void bn192_nist_mod_192_sqr(bn_t result, const bn_t a, const bn_t m)
{
#if 0 // MULTIPLICATION_ALGO == 1
  bn_word_t acc[ECDSA_NWORDS*2]; // buffer for a[]*b[];
  bn192_sqrx_core(acc, a);    // multiply
  for (int i = ECDSA_NWORDS-1; i >= 0; --i) {
    bn_word_t msw = acc[ECDSA_NWORDS+i];
    if (msw > 1)
      msw = bn192_mulsubw_core(&acc[i], m, msw);
    while (msw != 0)
      msw -= bn192_sub_core(&acc[i], m);
  }
  // copy and conditionally last subtract
  bn192_nist_mod_192(result, acc, m);
#else
  bn192_nist_mod_192_mul(result, a, a, m);
#endif
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
