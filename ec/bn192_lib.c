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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <openssl/bn.h>
#include <openssl/err.h>

#include "ecerr.h"
#include "bn192_lib.h"

enum {
  BNLIB_NBYTES = 24,
};

static BIGNUM* st_tmp[4]={0};

int bn192_init(void)
{
  for (int i = 0; i < sizeof(st_tmp)/sizeof(st_tmp[0]); ++i) {
    if ((st_tmp[i] = BN_new())==NULL) {
      bn192_cleanup();
      return 0;
    }
  }
  return 1;
}

void bn192_cleanup(void)
{
  for (int i = 0; i < sizeof(st_tmp)/sizeof(st_tmp[0]); ++i) {
    if (st_tmp[i]) {
      BN_free(st_tmp[i]);
      st_tmp[i] = 0;
    }
  }
}

int bn192_copy(uint8_t* result, const uint8_t* src)
{
  memcpy(result, src, BNLIB_NBYTES*sizeof(*result));
  return 1;
}

int bn192_zero(uint8_t* result)
{
  memset(result, 0, BNLIB_NBYTES*sizeof(*result));
  return 1;
}

int bn192_one(uint8_t* result)
{
  memset(result, 0, BNLIB_NBYTES*sizeof(*result));
  result[0] = 1;
  return 1;
}

int bn192_is_odd(const uint8_t* a)
{
  return a[0] & 1;
}

int bn192_is_one(const uint8_t* a)
{
  if (a[0] != 1)
    return 0; // is not one
  for (int i = 1; i < BNLIB_NBYTES; ++i)
    if (a[i] != 0)
      return 0; // is not one
  return 1; // is one
}

int bn192_is_zero(const uint8_t* a)
{
  for (int i = 0; i < BNLIB_NBYTES; ++i)
    if (a[i] != 0)
      return 0; // is not zero
  return 1; // is zero
}

int bn192_ucmp(const uint8_t* a, const uint8_t* b)
{
  for (int i = BNLIB_NBYTES-1; i >= 0; --i) {
    int av = a[i];
    int bv = b[i];
    if (av != bv)
      return av < bv ? -1 : 1; // is non equal
  }
  return 0; // is equal
}

int bn192_add_rshift1(uint8_t* result, const uint8_t* a, const uint8_t* b)
{
  if (!BN_lebin2bn(a, BNLIB_NBYTES, st_tmp[0])) return 0;
  if (!BN_lebin2bn(b, BNLIB_NBYTES, st_tmp[1])) return 0;

  if (!BN_add    (st_tmp[0], st_tmp[0], st_tmp[1])) return 0;
  if (!BN_rshift1(st_tmp[0], st_tmp[0]))            return 0;

  return BN_bn2lebinpad(st_tmp[0], result, BNLIB_NBYTES);
}

int bn192_mod_add_quick(uint8_t* result, const uint8_t* a, const uint8_t* b, const uint8_t* m)
{
  if (!BN_lebin2bn(a, BNLIB_NBYTES, st_tmp[0])) return 0;
  if (!BN_lebin2bn(b, BNLIB_NBYTES, st_tmp[1])) return 0;
  if (!BN_lebin2bn(m, BNLIB_NBYTES, st_tmp[2])) return 0;

  if (!BN_mod_add_quick(st_tmp[0], st_tmp[0], st_tmp[1], st_tmp[2])) return 0;

  return BN_bn2lebinpad(st_tmp[0], result, BNLIB_NBYTES);
}

int bn192_mod_sub_quick(uint8_t* result, const uint8_t* a, const uint8_t* b, const uint8_t* m)
{
  if (!BN_lebin2bn(a, BNLIB_NBYTES, st_tmp[0])) return 0;
  if (!BN_lebin2bn(b, BNLIB_NBYTES, st_tmp[1])) return 0;
  if (!BN_lebin2bn(m, BNLIB_NBYTES, st_tmp[2])) return 0;

  if (!BN_mod_sub_quick(st_tmp[0], st_tmp[0], st_tmp[1], st_tmp[2])) return 0;

  return BN_bn2lebinpad(st_tmp[0], result, BNLIB_NBYTES);
}

int bn192_mod_exp(uint8_t* result, const uint8_t* a, const uint8_t* p, const uint8_t* m, BN_CTX *ctx)
{
  if (!BN_lebin2bn(a, BNLIB_NBYTES, st_tmp[0])) return 0;
  if (!BN_lebin2bn(p, BNLIB_NBYTES, st_tmp[1])) return 0;
  if (!BN_lebin2bn(m, BNLIB_NBYTES, st_tmp[2])) return 0;

  if (!BN_mod_exp(st_tmp[0], st_tmp[0], st_tmp[1], st_tmp[2], ctx)) return 0;

  return BN_bn2lebinpad(st_tmp[0], result, BNLIB_NBYTES);
}

int bn192_mod_inverse(uint8_t* result, const uint8_t* a, const uint8_t* n, BN_CTX *ctx)
{
  if (!BN_lebin2bn(a, BNLIB_NBYTES, st_tmp[0])) return 0;
  if (!BN_lebin2bn(n, BNLIB_NBYTES, st_tmp[1])) return 0;

  if (!BN_mod_inverse(st_tmp[0], st_tmp[0], st_tmp[1], ctx)) return 0;

  return BN_bn2lebinpad(st_tmp[0], result, BNLIB_NBYTES);
}

int bn192_mod_lshift1_quick(uint8_t* result, const uint8_t* a, const uint8_t* m)
{
  if (!BN_lebin2bn(a, BNLIB_NBYTES, st_tmp[0])) return 0;
  if (!BN_lebin2bn(m, BNLIB_NBYTES, st_tmp[1])) return 0;

  if (!BN_mod_lshift1_quick(st_tmp[0], st_tmp[0], st_tmp[1])) return 0;

  return BN_bn2lebinpad(st_tmp[0], result, BNLIB_NBYTES);
}

int bn192_mod_lshift_quick(uint8_t* result, const uint8_t* a, int n, const uint8_t* m)
{
  if (!BN_lebin2bn(a, BNLIB_NBYTES, st_tmp[0])) return 0;
  if (!BN_lebin2bn(m, BNLIB_NBYTES, st_tmp[1])) return 0;

  if (!BN_mod_lshift_quick(st_tmp[0], st_tmp[0], n, st_tmp[1])) return 0;

  return BN_bn2lebinpad(st_tmp[0], result, BNLIB_NBYTES);
}

int bn192_mod_mul(uint8_t* result, const uint8_t* a, const uint8_t* b, const uint8_t* m, BN_CTX *ctx)
{
  if (!BN_lebin2bn(a, BNLIB_NBYTES, st_tmp[0])) return 0;
  if (!BN_lebin2bn(b, BNLIB_NBYTES, st_tmp[1])) return 0;
  if (!BN_lebin2bn(m, BNLIB_NBYTES, st_tmp[2])) return 0;

  if (!BN_mod_mul(st_tmp[0], st_tmp[0], st_tmp[1], st_tmp[2], ctx)) return 0;

  return BN_bn2lebinpad(st_tmp[0], result, BNLIB_NBYTES);
}

int bn192_mod_sqr(uint8_t* result, const uint8_t* a, const uint8_t* m, BN_CTX *ctx)
{
  if (!BN_lebin2bn(a, BNLIB_NBYTES, st_tmp[0])) return 0;
  if (!BN_lebin2bn(m, BNLIB_NBYTES, st_tmp[1])) return 0;

  if (!BN_mod_sqr(st_tmp[0], st_tmp[0], st_tmp[1], ctx)) return 0;

  return BN_bn2lebinpad(st_tmp[0], result, BNLIB_NBYTES);
}

int bn192_nist_mod_192_mul(uint8_t* result, const uint8_t* a, const uint8_t* b, const uint8_t* m, BN_CTX *ctx)
{
  if (!BN_lebin2bn(a, BNLIB_NBYTES, st_tmp[0])) return 0;
  if (!BN_lebin2bn(b, BNLIB_NBYTES, st_tmp[1])) return 0;
  if (!BN_lebin2bn(m, BNLIB_NBYTES, st_tmp[2])) return 0;

  if (!BN_mul         (st_tmp[0], st_tmp[0], st_tmp[1], ctx)) return 0;
  if (!BN_nist_mod_192(st_tmp[0], st_tmp[0], st_tmp[2], ctx)) return 0;

  return BN_bn2lebinpad(st_tmp[0], result, BNLIB_NBYTES);
}

int bn192_nist_mod_192_sqr(uint8_t* result, const uint8_t* a, const uint8_t* m, BN_CTX *ctx)
{
  if (!BN_lebin2bn(a, BNLIB_NBYTES, st_tmp[0])) return 0;
  if (!BN_lebin2bn(m, BNLIB_NBYTES, st_tmp[1])) return 0;

  if (!BN_sqr         (st_tmp[0], st_tmp[0],            ctx)) return 0;
  if (!BN_nist_mod_192(st_tmp[0], st_tmp[0], st_tmp[1], ctx)) return 0;

  return BN_bn2lebinpad(st_tmp[0], result, BNLIB_NBYTES);
}

int bn192_nnmod(uint8_t* result, const uint8_t* m, const uint8_t* d, BN_CTX *ctx)
{
  if (!BN_lebin2bn(m, BNLIB_NBYTES, st_tmp[0])) return 0;
  if (!BN_lebin2bn(d, BNLIB_NBYTES, st_tmp[1])) return 0;

  if (!BN_nnmod(st_tmp[0], st_tmp[0], st_tmp[1], ctx)) return 0;

  return BN_bn2lebinpad(st_tmp[0], result, BNLIB_NBYTES);
}

int bn192_rshift1(uint8_t* result, const uint8_t* a)
{
  if (!BN_lebin2bn(a, BNLIB_NBYTES, st_tmp[0])) return 0;

  if (!BN_rshift1(st_tmp[0], st_tmp[0]))        return 0;

  return BN_bn2lebinpad(st_tmp[0], result, BNLIB_NBYTES);
}
