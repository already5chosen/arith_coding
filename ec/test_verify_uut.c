// #include "ecs_locl.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <openssl/ec.h>
#include <openssl/obj_mac.h>
#include <openssl/err.h>

enum {
  ECDSA_NBYTES = 24,
  ECDSA_NBITS  = ECDSA_NBYTES*8,
};
static EC_GROUP* st_group = NULL;
static BIGNUM* st_pub_x = NULL;
static BIGNUM* st_pub_y = NULL;
static EC_POINT* st_pub_key = NULL;
static BIGNUM*   st_group_p = NULL;
static BIGNUM*   st_group_a = NULL;
static BIGNUM*   st_group_b = NULL;


void uut_cleanup(void)
{
  if (st_group_p) {
    BN_free(st_group_p);
    st_group_p = 0;
  }
  if (st_group_a) {
    BN_free(st_group_a);
    st_group_a = 0;
  }
  if (st_group_b) {
    BN_free(st_group_b);
    st_group_b = 0;
  }
  if (st_group) {
    EC_GROUP_free(st_group);
    st_group = 0;
  }
  if (st_pub_x) {
    BN_free(st_pub_x);
    st_pub_x = 0;
  }
  if (st_pub_y) {
    BN_free(st_pub_y);
    st_pub_y = 0;
  }
  if (st_pub_key) {
    EC_POINT_free(st_pub_key);
    st_pub_key = 0;
  }
}

int uut_init(void) {
  st_group_p = BN_new();
  st_group_a = BN_new();
  st_group_b = BN_new();
  st_pub_x   = BN_new();
  st_pub_y   = BN_new();

  if (!st_group_p || !st_group_a
    || !st_group_b || !st_pub_x || !st_pub_x)
    goto err;
  st_group = EC_GROUP_new_by_curve_name(NID_X9_62_prime192v1);
  if (!st_group)
    goto err;
  if (!EC_GROUP_get_curve(st_group, st_group_p,
                      st_group_a, st_group_b, 0))
    goto err;
  st_pub_key = EC_POINT_new(st_group);
  if (!st_pub_key)
    goto err;

  return 1;

  err:
  uut_cleanup();
  return 0;
}

int uut_set_public_key(unsigned char key[2][24])
{
  if (!BN_lebin2bn(key[0], 24, st_pub_x))
    return 0;
  if (!BN_lebin2bn(key[1], 24, st_pub_y))
    return 0;
  if (!EC_POINT_set_affine_coordinates(
    st_group, st_pub_key, st_pub_x, st_pub_y, 0))
    return 0;

  return 1;
}

static int ec_group_do_inverse_ord(
 const EC_GROUP* group,
 BIGNUM*         res,
 const BIGNUM*   x,
 BN_CTX*         ctx)
{
    BIGNUM *e = NULL;
    BN_CTX *new_ctx = NULL;
    int ret = 0;

    if (ctx == NULL && (ctx = new_ctx = BN_CTX_secure_new()) == NULL)
        return 0;

    BN_CTX_start(ctx);
    if ((e = BN_CTX_get(ctx)) == NULL)
        goto err;

    const BIGNUM *order = EC_GROUP_get0_order(group);
    if (!order)
      goto err;

    if (!BN_set_word(e, 2))
        goto err;
    if (!BN_sub(e, order, e))
        goto err;

    if (!BN_mod_exp(res, x, e, order, ctx))
        goto err;

    ret = 1;

 err:
    if (ctx != NULL)
        BN_CTX_end(ctx);
    BN_CTX_free(new_ctx);
    return ret;

}

// derived from ec_GFp_nist_field_mul()
// in \openssl-1.1.1\crypto\ec\ecp_nist.c
static int ec_GFp_nist_field_mul(
  const EC_GROUP* group,
  BIGNUM*         r,
  const BIGNUM*   a,
  const BIGNUM*   b,
  BN_CTX*         ctx)
{
    int ret = 0;
    BN_CTX *ctx_new = NULL;

    if (!group || !r || !a || !b) {
        ECerr(EC_F_EC_GFP_NIST_FIELD_MUL, ERR_R_PASSED_NULL_PARAMETER);
        goto err;
    }
    if (!ctx)
        if ((ctx_new = ctx = BN_CTX_new()) == NULL)
            goto err;

    if (!BN_mul(r, a, b, ctx))
        goto err;
    if (!BN_nist_mod_192(r, r, st_group_p, ctx))
        goto err;

    ret = 1;
 err:
    BN_CTX_free(ctx_new);
    return ret;
}

// derived from ec_GFp_nist_field_sqr()
// in \openssl-1.1.1\crypto\ec\ecp_nist.c
static int ec_GFp_nist_field_sqr(
  const EC_GROUP* group,
  BIGNUM*         r,
  const BIGNUM*   a,
  BN_CTX*         ctx)
{
    int ret = 0;
    BN_CTX *ctx_new = NULL;

    if (!group || !r || !a) {
        ECerr(EC_F_EC_GFP_NIST_FIELD_SQR, EC_R_PASSED_NULL_PARAMETER);
        goto err;
    }
    if (!ctx)
        if ((ctx_new = ctx = BN_CTX_new()) == NULL)
            goto err;

    if (!BN_sqr(r, a, ctx))
        goto err;
    if (!BN_nist_mod_192(r, r, st_group_p, ctx))
        goto err;

    ret = 1;
 err:
    BN_CTX_free(ctx_new);
    return ret;
}


// derived from ec_GFp_simple_dbl()
// in \openssl-1.1.1\crypto\ec\ecp_smpl.c
static int ec_GFp_simple_dbl(const EC_GROUP *group, EC_POINT *r, const EC_POINT *a,
                      BN_CTX *ctx)
{
    int (*field_mul) (const EC_GROUP *, BIGNUM *, const BIGNUM *,
                      const BIGNUM *, BN_CTX *);
    int (*field_sqr) (const EC_GROUP *, BIGNUM *, const BIGNUM *, BN_CTX *);
    const BIGNUM *p;
    BN_CTX *new_ctx = NULL;
    BIGNUM *n0, *n1, *n2, *n3;
    BIGNUM *aX, *aY, *aZ;
    BIGNUM *rX, *rY, *rZ;
    int ret = 0;

    field_mul = ec_GFp_nist_field_mul;
    field_sqr = ec_GFp_nist_field_sqr;
    p = st_group_p;

    if (ctx == NULL) {
        ctx = new_ctx = BN_CTX_new();
        if (ctx == NULL)
            return 0;
    }

    BN_CTX_start(ctx);
    rX = BN_CTX_get(ctx);
    rY = BN_CTX_get(ctx);
    rZ = BN_CTX_get(ctx);
    aX = BN_CTX_get(ctx);
    aY = BN_CTX_get(ctx);
    aZ = BN_CTX_get(ctx);
    n0 = BN_CTX_get(ctx);
    n1 = BN_CTX_get(ctx);
    n2 = BN_CTX_get(ctx);
    n3 = BN_CTX_get(ctx);
    if (n3 == NULL)
        goto err;

    if (EC_POINT_is_at_infinity(group, a)) {
        BN_zero(rX);
        BN_zero(rY);
        BN_zero(rZ);
        goto done;
    }

   if (!EC_POINT_get_Jprojective_coordinates_GFp(
          group, a, aX, aY, aZ, ctx))
        goto err;

    /*
     * Note that in this function we must not read components of 'a' once we
     * have written the corresponding components of 'r'. ('r' might the same
     * as 'a'.)
     */

    /* n1 */
  //  if (a->Z_is_one) {
  //      if (!field_sqr(group, n0, a->X, ctx))
  //          goto err;
  //      if (!BN_mod_lshift1_quick(n1, n0, p))
  //          goto err;
  //      if (!BN_mod_add_quick(n0, n0, n1, p))
  //          goto err;
  //      if (!BN_mod_add_quick(n1, n0, group->a, p))
  //          goto err;
  //      /* n1 = 3 * X_a^2 + a_curve */
  //  } else if (group->a_is_minus3) {
  //      if (!field_sqr(group, n1, a->Z, ctx))
  //          goto err;
  //      if (!BN_mod_add_quick(n0, a->X, n1, p))
  //          goto err;
  //      if (!BN_mod_sub_quick(n2, a->X, n1, p))
  //          goto err;
  //      if (!field_mul(group, n1, n0, n2, ctx))
  //          goto err;
  //      if (!BN_mod_lshift1_quick(n0, n1, p))
  //          goto err;
  //      if (!BN_mod_add_quick(n1, n0, n1, p))
  //          goto err;
  //      /*-
  //       * n1 = 3 * (X_a + Z_a^2) * (X_a - Z_a^2)
  //       *    = 3 * X_a^2 - 3 * Z_a^4
  //       */
  //  } else {
        if (!field_sqr(group, n0, aX, ctx))
            goto err;
        if (!BN_mod_lshift1_quick(n1, n0, p))
            goto err;
        if (!BN_mod_add_quick(n0, n0, n1, p))
            goto err;
        if (!field_sqr(group, n1, aZ, ctx))
            goto err;
        if (!field_sqr(group, n1, n1, ctx))
            goto err;
        if (!field_mul(group, n1, n1, st_group_a, ctx))
            goto err;
        if (!BN_mod_add_quick(n1, n1, n0, p))
            goto err;
        /* n1 = 3 * X_a^2 + a_curve * Z_a^4 */
  //  }

    /* Z_r */
  // if (a->Z_is_one) {
  //     if (!BN_copy(n0, a->Y))
  //         goto err;
  // } else {
        if (!field_mul(group, n0, aY, aZ, ctx))
            goto err;
  //  }
    if (!BN_mod_lshift1_quick(rZ, n0, p))
        goto err;
  //  r->Z_is_one = 0;
    /* Z_r = 2 * Y_a * Z_a */

    /* n2 */
    if (!field_sqr(group, n3, aY, ctx))
        goto err;
    if (!field_mul(group, n2, aX, n3, ctx))
        goto err;
    if (!BN_mod_lshift_quick(n2, n2, 2, p))
        goto err;
    /* n2 = 4 * X_a * Y_a^2 */

    /* X_r */
    if (!BN_mod_lshift1_quick(n0, n2, p))
        goto err;
    if (!field_sqr(group, rX, n1, ctx))
        goto err;
    if (!BN_mod_sub_quick(rX, rX, n0, p))
        goto err;
    /* X_r = n1^2 - 2 * n2 */

    /* n3 */
    if (!field_sqr(group, n0, n3, ctx))
        goto err;
    if (!BN_mod_lshift_quick(n3, n0, 3, p))
        goto err;
    /* n3 = 8 * Y_a^4 */

    /* Y_r */
    if (!BN_mod_sub_quick(n0, n2, rX, p))
        goto err;
    if (!field_mul(group, n0, n1, n0, ctx))
        goto err;
    if (!BN_mod_sub_quick(rY, n0, n3, p))
        goto err;
    /* Y_r = n1 * (n2 - X_r) - n3 */

 done:
    if (!EC_POINT_set_Jprojective_coordinates_GFp(
      group, r, rX, rY, rZ, ctx))
      goto err;

    ret = 1;

 err:
    BN_CTX_end(ctx);
    BN_CTX_free(new_ctx);
    return ret;
}

// derived from ec_GFp_simple_add()
// in \openssl-1.1.1\crypto\ec\ecp_smpl.c
static int ec_GFp_simple_add(const EC_GROUP *group, EC_POINT *r, const EC_POINT *a,
                      const EC_POINT *b, BN_CTX *ctx)
{
    int (*field_mul) (const EC_GROUP *, BIGNUM *, const BIGNUM *,
                      const BIGNUM *, BN_CTX *);
    int (*field_sqr) (const EC_GROUP *, BIGNUM *, const BIGNUM *, BN_CTX *);
    const BIGNUM *p;
    BN_CTX *new_ctx = NULL;
    BIGNUM *n0, *n1, *n2, *n3, *n4, *n5, *n6;
    BIGNUM *aX, *aY, *aZ;
    BIGNUM *bX, *bY, *bZ;
    BIGNUM *rX, *rY, *rZ;
    int ret = 0;

    if (a == b)
        return ec_GFp_simple_dbl(group, r, a, ctx);
    if (EC_POINT_is_at_infinity(group, a))
        return EC_POINT_copy(r, b);
    if (EC_POINT_is_at_infinity(group, b))
        return EC_POINT_copy(r, a);

    field_mul = ec_GFp_nist_field_mul;
    field_sqr = ec_GFp_nist_field_sqr;
    p = st_group_p;

    if (ctx == NULL) {
        ctx = new_ctx = BN_CTX_new();
        if (ctx == NULL)
            return 0;
    }

    BN_CTX_start(ctx);
    aX = BN_CTX_get(ctx);
    aY = BN_CTX_get(ctx);
    aZ = BN_CTX_get(ctx);
    bX = BN_CTX_get(ctx);
    bY = BN_CTX_get(ctx);
    bZ = BN_CTX_get(ctx);
    rX = BN_CTX_get(ctx);
    rY = BN_CTX_get(ctx);
    rZ = BN_CTX_get(ctx);
    n0 = BN_CTX_get(ctx);
    n1 = BN_CTX_get(ctx);
    n2 = BN_CTX_get(ctx);
    n3 = BN_CTX_get(ctx);
    n4 = BN_CTX_get(ctx);
    n5 = BN_CTX_get(ctx);
    n6 = BN_CTX_get(ctx);
    if (n6 == NULL)
        goto end;

   if (!EC_POINT_get_Jprojective_coordinates_GFp(
    group, a, aX, aY, aZ, ctx))
        goto end;

   if (!EC_POINT_get_Jprojective_coordinates_GFp(
    group, b, bX, bY, bZ, ctx))
        goto end;

    /*
     * Note that in this function we must not read components of 'a' or 'b'
     * once we have written the corresponding components of 'r'. ('r' might
     * be one of 'a' or 'b'.)
     */

    /* n1, n2 */
  //  if (b->Z_is_one) {
  //      if (!BN_copy(n1, a->X))
  //          goto end;
  //      if (!BN_copy(n2, a->Y))
  //          goto end;
  //      /* n1 = X_a */
  //      /* n2 = Y_a */
  //  } else {
        if (!field_sqr(group, n0, bZ, ctx))
            goto end;
        if (!field_mul(group, n1, aX, n0, ctx))
            goto end;
        /* n1 = X_a * Z_b^2 */

        if (!field_mul(group, n0, n0, bZ, ctx))
            goto end;
        if (!field_mul(group, n2, aY, n0, ctx))
            goto end;
        /* n2 = Y_a * Z_b^3 */
  //  }

    /* n3, n4 */
  //  if (a->Z_is_one) {
  //      if (!BN_copy(n3, b->X))
  //          goto end;
  //      if (!BN_copy(n4, b->Y))
  //          goto end;
  //      /* n3 = X_b */
  //      /* n4 = Y_b */
  //  } else {
        if (!field_sqr(group, n0, aZ, ctx))
            goto end;
        if (!field_mul(group, n3, bX, n0, ctx))
            goto end;
        /* n3 = X_b * Z_a^2 */

        if (!field_mul(group, n0, n0, aZ, ctx))
            goto end;
        if (!field_mul(group, n4, bY, n0, ctx))
            goto end;
        /* n4 = Y_b * Z_a^3 */
  //  }

    /* n5, n6 */
    if (!BN_mod_sub_quick(n5, n1, n3, p))
        goto end;
    if (!BN_mod_sub_quick(n6, n2, n4, p))
        goto end;
    /* n5 = n1 - n3 */
    /* n6 = n2 - n4 */

    if (BN_is_zero(n5)) {
        if (BN_is_zero(n6)) {
            /* a is the same point as b */
            BN_CTX_end(ctx);
            ret = ec_GFp_simple_dbl(group, r, a, ctx);
            ctx = NULL;
            goto end;
        } else {
            /* a is the inverse of b */
            BN_zero(rZ);
  //          r->Z_is_one = 0;
            ret = 1;
            goto end;
        }
    }

    /* 'n7', 'n8' */
    if (!BN_mod_add_quick(n1, n1, n3, p))
        goto end;
    if (!BN_mod_add_quick(n2, n2, n4, p))
        goto end;
    /* 'n7' = n1 + n3 */
    /* 'n8' = n2 + n4 */

    /* Z_r */
  //  if (a->Z_is_one && b->Z_is_one) {
  //      if (!BN_copy(r->Z, n5))
  //          goto end;
  //  } else {
  //      if (a->Z_is_one) {
  //          if (!BN_copy(n0, b->Z))
  //              goto end;
  //      } else if (b->Z_is_one) {
  //          if (!BN_copy(n0, a->Z))
  //              goto end;
  //      } else {
            if (!field_mul(group, n0, aZ, bZ, ctx))
                goto end;
  //      }
        if (!field_mul(group, rZ, n0, n5, ctx))
            goto end;
  //  }
  //  r->Z_is_one = 0;
    /* Z_r = Z_a * Z_b * n5 */

    /* X_r */
    if (!field_sqr(group, n0, n6, ctx))
        goto end;
    if (!field_sqr(group, n4, n5, ctx))
        goto end;
    if (!field_mul(group, n3, n1, n4, ctx))
        goto end;
    if (!BN_mod_sub_quick(rX, n0, n3, p))
        goto end;
    /* X_r = n6^2 - n5^2 * 'n7' */

    /* 'n9' */
    if (!BN_mod_lshift1_quick(n0, rX, p))
        goto end;
    if (!BN_mod_sub_quick(n0, n3, n0, p))
        goto end;
    /* n9 = n5^2 * 'n7' - 2 * X_r */

    /* Y_r */
    if (!field_mul(group, n0, n0, n6, ctx))
        goto end;
    if (!field_mul(group, n5, n4, n5, ctx))
        goto end;               /* now n5 is n5^3 */
    if (!field_mul(group, n1, n2, n5, ctx))
        goto end;
    if (!BN_mod_sub_quick(n0, n0, n1, p))
        goto end;
    if (BN_is_odd(n0))
        if (!BN_add(n0, n0, p))
            goto end;
    /* now  0 <= n0 < 2*p,  and n0 is even */
    if (!BN_rshift1(rY, n0))
        goto end;
    /* Y_r = (n6 * 'n9' - 'n8' * 'n5^3') / 2 */

    if (!EC_POINT_set_Jprojective_coordinates_GFp(
      group, r, rX, rY, rZ, ctx))
      goto end;

    ret = 1;

 end:
    if (ctx)                    /* otherwise we already called BN_CTX_end */
        BN_CTX_end(ctx);
    BN_CTX_free(new_ctx);
    return ret;
}

#if 0
// derived from ec_GFp_simple_point_get_affine_coordinates()
// in \openssl-1.1.1\crypto\ec\ecp_smpl.c
int ec_GFp_simple_point_get_affine_coordinates(const EC_GROUP *group,
                                               const EC_POINT *point,
                                               BIGNUM *x, BIGNUM *y,
                                               BN_CTX *ctx)
{
    BN_CTX *new_ctx = NULL;
    BIGNUM *Z, *Z_1, *Z_2, *Z_3;
    const BIGNUM *Z_;
    int ret = 0;

    if (EC_POINT_is_at_infinity(group, point)) {
        ECerr(EC_F_EC_GFP_SIMPLE_POINT_GET_AFFINE_COORDINATES,
              EC_R_POINT_AT_INFINITY);
        return 0;
    }

    if (ctx == NULL) {
        ctx = new_ctx = BN_CTX_new();
        if (ctx == NULL)
            return 0;
    }

    BN_CTX_start(ctx);
    Z = BN_CTX_get(ctx);
    Z_1 = BN_CTX_get(ctx);
    Z_2 = BN_CTX_get(ctx);
    Z_3 = BN_CTX_get(ctx);
    if (Z_3 == NULL)
        goto err;

    /* transform  (X, Y, Z)  into  (x, y) := (X/Z^2, Y/Z^3) */

    if (group->meth->field_decode) {
        if (!group->meth->field_decode(group, Z, point->Z, ctx))
            goto err;
        Z_ = Z;
    } else {
        Z_ = point->Z;
    }

    if (BN_is_one(Z_)) {
        if (group->meth->field_decode) {
            if (x != NULL) {
                if (!group->meth->field_decode(group, x, point->X, ctx))
                    goto err;
            }
            if (y != NULL) {
                if (!group->meth->field_decode(group, y, point->Y, ctx))
                    goto err;
            }
        } else {
            if (x != NULL) {
                if (!BN_copy(x, point->X))
                    goto err;
            }
            if (y != NULL) {
                if (!BN_copy(y, point->Y))
                    goto err;
            }
        }
    } else {
        if (!BN_mod_inverse(Z_1, Z_, group->field, ctx)) {
            ECerr(EC_F_EC_GFP_SIMPLE_POINT_GET_AFFINE_COORDINATES,
                  ERR_R_BN_LIB);
            goto err;
        }

        if (group->meth->field_encode == 0) {
            /* field_sqr works on standard representation */
            if (!group->meth->field_sqr(group, Z_2, Z_1, ctx))
                goto err;
        } else {
            if (!BN_mod_sqr(Z_2, Z_1, group->field, ctx))
                goto err;
        }

        if (x != NULL) {
            /*
             * in the Montgomery case, field_mul will cancel out Montgomery
             * factor in X:
             */
            if (!group->meth->field_mul(group, x, point->X, Z_2, ctx))
                goto err;
        }

        if (y != NULL) {
            if (group->meth->field_encode == 0) {
                /*
                 * field_mul works on standard representation
                 */
                if (!group->meth->field_mul(group, Z_3, Z_2, Z_1, ctx))
                    goto err;
            } else {
                if (!BN_mod_mul(Z_3, Z_2, Z_1, group->field, ctx))
                    goto err;
            }

            /*
             * in the Montgomery case, field_mul will cancel out Montgomery
             * factor in Y:
             */
            if (!group->meth->field_mul(group, y, point->Y, Z_3, ctx))
                goto err;
        }
    }

    ret = 1;

 err:
    BN_CTX_end(ctx);
    BN_CTX_free(new_ctx);
    return ret;
}
#endif

/** Computes r = q1 * m1 + q2 * m2
 *  \param  group  underlying EC_GROUP object
 *  \param  r      EC_POINT object for result
 *                 (can't be the same as q1 or q2)
 *  \param  q1     EC_POINT object to be multiplied by m1
 *  \param  m1     BIGNUM positive multiplication factor for q1
 *  \param  q2     EC_POINT object to be multiplied by m2
 *  \param  m2     BIGNUM non-negative multiplication factor for q2
 *  \param  ctx    BN_CTX object (optional)
 *  \return 1 on success and 0 if an error occurred
 */
static int EC_POINT_muladd2(
  const EC_GROUP* group,
  EC_POINT*       r,
  const EC_POINT* q1,
  const BIGNUM*   m1,
  const EC_POINT* q2,
  const BIGNUM*   m2,
  BN_CTX*         ctx)
{
  unsigned char m_oct[2][ECDSA_NBYTES];
  if (BN_bn2lebinpad(m1, m_oct[0], ECDSA_NBYTES) != ECDSA_NBYTES) {
    fprintf(stderr, "BN_bn2lebinpad(s) fail.\n");
    return 0;
  }
  if (BN_bn2lebinpad(m2, m_oct[1], ECDSA_NBYTES) != ECDSA_NBYTES) {
    fprintf(stderr, "BN_bn2lebinpad(s) fail.\n");
    return 0;
  }

  const EC_POINT* q[2] = { q1, q2 };
  int nz = 0;
  for (int bit_i = ECDSA_NBITS-1; bit_i >= 0; --bit_i) {
    if (nz) {
      if (!ec_GFp_simple_dbl(group, r, r, ctx)) {
        ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_EC_LIB);
        return 0;
      }
    }
    for (int k = 0; k < 2; ++k) {
      if ((m_oct[k][bit_i/8] >> (bit_i%8)) & 1) {
        // bit is set
        if (nz) {
          if (!ec_GFp_simple_add(group, r, r, q[k], ctx)) {
            ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_EC_LIB);
            return 0;
          }
        } else {
          if (!EC_POINT_copy(r, q[k])) {
            ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_INTERNAL_ERROR);
            return 0;
          }
        }
        nz = 1;
      }
    }
  }
  return 1; // success
}


static int ossl_ecdsa_verify_sig(
  const unsigned char *dgst,
  const unsigned char signature_octets[2][ECDSA_NBYTES])
{
    int ret = -1;
    BN_CTX *ctx;
    BIGNUM *u1, *u2, *m, *X, *sig_r, *sig_s;
    EC_POINT *point = NULL;
    const EC_GROUP *group   = st_group;
    const EC_POINT *pub_key = st_pub_key;

    ctx = BN_CTX_new();
    if (ctx == NULL) {
        ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_MALLOC_FAILURE);
        return -1;
    }
    BN_CTX_start(ctx);
    sig_r = BN_CTX_get(ctx);
    sig_s = BN_CTX_get(ctx);
    u1 = BN_CTX_get(ctx);
    u2 = BN_CTX_get(ctx);
    m = BN_CTX_get(ctx);
    X = BN_CTX_get(ctx);
    if (X == NULL) {
        ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_BN_LIB);
        goto err;
    }

    if (!BN_lebin2bn(signature_octets[0], ECDSA_NBYTES, sig_r)) {
      ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_BN_LIB);
      goto err;
    }
    if (!BN_lebin2bn(signature_octets[1], ECDSA_NBYTES, sig_s)) {
      ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_BN_LIB);
      goto err;
    }

    const BIGNUM *order = EC_GROUP_get0_order(group);
    if (order == NULL) {
        ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_EC_LIB);
        goto err;
    }

    if (BN_is_zero(sig_r) || BN_is_negative(sig_r) ||
        BN_ucmp(sig_r, order) >= 0 || BN_is_zero(sig_s) ||
        BN_is_negative(sig_s) || BN_ucmp(sig_s, order) >= 0) {
        ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, EC_R_BAD_SIGNATURE);
        ret = 0;                /* signature is invalid */
        goto err;
    }
    // calculate tmp1 = inv(S) mod order
    if (!ec_group_do_inverse_ord(group, u2, sig_s, ctx)) {
        ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_BN_LIB);
        goto err;
    }

    //
    // There is no need to truncate digest,
    // because it is always supplied by caller
    // with correct length
    //

    // digest -> m
    if (!BN_bin2bn(dgst, ECDSA_NBYTES, m)) {
        ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_BN_LIB);
        goto err;
    }

    // u1 = m * tmp mod order
    if (!BN_mod_mul(u1, m, u2, order, ctx)) {
        ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_BN_LIB);
        goto err;
    }
    // u2 = r * w mod q
    if (!BN_mod_mul(u2, sig_r, u2, order, ctx)) {
        ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_BN_LIB);
        goto err;
    }

    if ((point = EC_POINT_new(group)) == NULL) {
        ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_MALLOC_FAILURE);
        goto err;
    }
    // point = generator*u1 + pub_key*u2
    if (!EC_POINT_muladd2(group, point,
      EC_GROUP_get0_generator(group), u1,
      pub_key,                        u2,
      ctx)) {
        ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_EC_LIB);
        goto err;
    }

    if (!EC_POINT_get_affine_coordinates(group, point, X, NULL, ctx)) {
        ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_EC_LIB);
        goto err;
    }

    if (!BN_nnmod(u1, X, order, ctx)) {
        ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_BN_LIB);
        goto err;
    }
    /*  if the signature is correct u1 is equal to sig_r */
    ret = (BN_ucmp(u1, sig_r) == 0);
 err:
    BN_CTX_end(ctx);
    BN_CTX_free(ctx);
    EC_POINT_free(point);
    return ret;
}

int uut_verify(
 const unsigned char digest[24],
 const unsigned char signature_octets[2][24],
 int* res)
{
  *res = 0;
  int ret = ossl_ecdsa_verify_sig(digest, signature_octets);
  if (ret >= 0) {
    *res = ret;
  } else {
    unsigned long err = ERR_get_error();
    char errstr[512];
    ERR_error_string_n(err, errstr, sizeof(errstr));
    fprintf(stderr, "ECDSA_do_verify: error %lu.\n%s\n", err, errstr);
  }
  return (ret >= 0);
}