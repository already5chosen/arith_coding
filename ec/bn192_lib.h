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
BN_one
BN_rshift1
BN_sqr
BN_ucmp
BN_zero
*/

int bn192_init(void);
void bn192_cleanup(void);

int bn192_add_rshift1(uint8_t* result, const uint8_t* a, const uint8_t* b);
int bn192_copy(uint8_t* result, const uint8_t* src);
int bn192_is_odd(const uint8_t* a);
int bn192_is_one(const uint8_t* a);
int bn192_is_zero(const uint8_t* a);
int bn192_mod_add_quick(uint8_t* result, const uint8_t* a, const uint8_t* b, const uint8_t* m);
int bn192_mod_sub_quick(uint8_t* result, const uint8_t* a, const uint8_t* b, const uint8_t* m);
int bn192_mod_exp(uint8_t* result, const uint8_t* a, const uint8_t* p, const uint8_t* m, BN_CTX *ctx);
int bn192_mod_inverse(uint8_t* result, const uint8_t* a, const uint8_t* n, BN_CTX *ctx);
int bn192_mod_lshift1_quick(uint8_t* result, const uint8_t* a, const uint8_t* m);
int bn192_mod_lshift_quick(uint8_t* result, const uint8_t* a, int n, const uint8_t* m);
int bn192_mod_mul(uint8_t* result, const uint8_t* a, const uint8_t* b, const uint8_t* m, BN_CTX *ctx);
int bn192_mod_sqr(uint8_t* result, const uint8_t* a, const uint8_t* m, BN_CTX *ctx);
int bn192_nist_mod_192_mul(uint8_t* result, const uint8_t* a, const uint8_t* b, const uint8_t* m, BN_CTX *ctx);
int bn192_nist_mod_192_sqr(uint8_t* result, const uint8_t* a, const uint8_t* m, BN_CTX *ctx);
int bn192_nnmod(uint8_t* result, const uint8_t* m, const uint8_t* d, BN_CTX *ctx);
int bn192_rshift1(uint8_t* result, const uint8_t* a);
int bn192_one(uint8_t* result);
int bn192_ucmp(const uint8_t* a, const uint8_t* b);
int bn192_zero(uint8_t* result);

#endif
