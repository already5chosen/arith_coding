#ifndef INCLUDED_BN192_H
#define INCLUDED_BN192_H

#include <stdint.h>

/*
 Replacement for following BN_xxx functions

BN_add
BN_copy
BN_is_odd
BN_is_one
BN_is_zero
BN_mod_add_quick
BN_mod_inverse
BN_mod_lshift1_quick
BN_mod_lshift_quick
BN_mod_sub_quick
BN_mul
BN_nist_mod_192
BN_one
BN_rshift1
BN_sqr
BN_ucmp
BN_zero
*/

enum {
  ECDSA_NBYTES = 24,
  ECDSA_NBITS  = ECDSA_NBYTES*8,
  BN192LIB_BYTES_PER_WORD = 8,
  BN192LIB_BITS_PER_WORD  = BN192LIB_BYTES_PER_WORD*8,
  ECDSA_NWORDS = ECDSA_NBYTES/BN192LIB_BYTES_PER_WORD,
};

typedef uint64_t bn_word_t, bn_t[ECDSA_NWORDS];
typedef unsigned __int128 bn_doubleword_t;

int bn192_is_odd(const bn_t a);
int bn192_is_one(const bn_t a);
int bn192_is_zero(const bn_t a);
int bn192_ucmp(const bn_t a, const bn_t b);

void bn192_copy(bn_t result, const bn_t src);
void bn192_zero(bn_t result);
void bn192_one(bn_t result);

void bn192_add_rshift1(bn_t result, const bn_t a, const bn_t b);
void bn192_mod_add_quick(bn_t result, const bn_t a, const bn_t b, const bn_t m);
void bn192_mod_sub_quick(bn_t result, const bn_t a, const bn_t b, const bn_t m);
void bn192_mod_inverse(bn_t result, const bn_t a, const bn_t n);
void bn192_mod_lshift1_quick(bn_t result, const bn_t a, const bn_t m);
void bn192_mod_lshift_quick(bn_t result, const bn_t a, int n, const bn_t m);
void bn192_nist_mod_192    (bn_t result, const bn_t a, const bn_t m);
void bn192_nist_mod_192_mul(bn_t result, const bn_t a, const bn_t b, const bn_t m);
void bn192_nist_mod_192_sqr(bn_t result, const bn_t a, const bn_t m);
void bn192_rshift1(bn_t result, const bn_t a);

#endif
