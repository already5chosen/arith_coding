#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <openssl/err.h>

#include "ecerr.h"
#include "bn192_lib.h"

typedef struct {
  bn_t X;
  bn_t Y;
  bn_t Z;        // Jacobian projective coordinates: (X, Y, Z) represents (X/Z^2, Y/Z^3) if Z != 0
  int  Z_is_one; // enable optimized point arithmetics for special case
} ec_point_t;

typedef struct {
  bn_t field;
  bn_t a;
  // bn_t b; // it seems, as long as we trust public key, parameter b is not used in verification process
  bn_t order;
  ec_point_t generator;
} ec_group_t;

static ec_group_t st_group = {
 { // p
  0xFFFFFFFFFFFFFFFF,
  0xFFFFFFFFFFFFFFFE,
  0xFFFFFFFFFFFFFFFF,
 },
 { // a
  0xFFFFFFFFFFFFFFFC,
  0xFFFFFFFFFFFFFFFE,
  0xFFFFFFFFFFFFFFFF,
 },
 { // order
  0x146BC9B1B4D22831,
  0xFFFFFFFF99DEF836,
  0xFFFFFFFFFFFFFFFF,
 },
  // generator
 {
  { // .X
  0xF4FF0AFD82FF1012,
  0x7CBF20EB43A18800,
  0x188DA80EB03090F6,
  },
  { // .Y
  0x73f977a11e794811,
  0x631011ed6b24cdd5,
  0x07192b95ffc8da78,
  },
  { 1 }, // .Z
  1,     // .Z_is_one
 },
};

static ec_point_t st_pub_key;

void uut_cleanup(void)
{
}

int uut_init(void) {
  return 1;
}

int uut_set_public_key(unsigned char key[2][24])
{
  bn192_copy(st_pub_key.X, (bn_word_t*)key[0]); // non-portable, depends on LE byte order
  bn192_copy(st_pub_key.Y, (bn_word_t*)key[1]); // non-portable, depends on LE byte order
  bn192_one(st_pub_key.Z);
  st_pub_key.Z_is_one = 1;

  return 1;
}

static int ec_point_is_at_infinity(const ec_point_t* a) {
  return bn192_is_zero(a->Z);
}

// derived from ec_GFp_simple_dbl()
// in \openssl-1.1.1\crypto\ec\ecp_smpl.c
static void ec_point_dbl(
  const ec_group_t* group,
  ec_point_t*       r,
  const ec_point_t* a)
{
    if (ec_point_is_at_infinity(a)) {
      bn192_zero(r->Z);
      r->Z_is_one = 0;
      return;
    }

    bn_t n0, n1, n2, n3;
    bn_t rX;

    /*
     * Note that in this function we must not read components of 'a' once we
     * have written the corresponding components of 'r'. ('r' might the same
     * as 'a'.)
     */

    /* n1 */
    if (a->Z_is_one) {
        bn192_nist_mod_192_sqr(n0, a->X, group->field);
        bn192_mod_lshift1_quick(n1, n0, group->field);
        bn192_mod_add_quick(n0, n0, n1, group->field);
        bn192_mod_add_quick(n1, n0, group->a, group->field);
        /* n1 = 3 * X_a^2 + a_curve */
  //  } else if (group->a_is_minus3) {
    } else {
        bn192_nist_mod_192_sqr(n1, a->Z, group->field);
        bn192_mod_add_quick(n0, a->X, n1, group->field);
        bn192_mod_sub_quick(n2, a->X, n1, group->field);
        bn192_nist_mod_192_mul(n1, n0, n2, group->field);
        bn192_mod_lshift1_quick(n0, n1, group->field);
        bn192_mod_add_quick(n1, n0, n1, group->field);
        /*-
         * n1 = 3 * (X_a + Z_a^2) * (X_a - Z_a^2)
         *    = 3 * X_a^2 - 3 * Z_a^4
         */
  //  } else {
  //      bn192_nist_mod_192_sqr(n0, a->X, group->field);
  //      bn192_mod_lshift1_quick(n1, n0, group->field);
  //      bn192_mod_add_quick(n0, n0, n1, group->field);
  //      bn192_nist_mod_192_sqr(n1, a->Z, group->field);
  //      bn192_nist_mod_192_sqr(n1, n1, group->field);
  //      bn192_nist_mod_192_mul(n1, n1, group->a, group->field);
  //      bn192_mod_add_quick(n1, n1, n0, group->field);
  //      /* n1 = 3 * X_a^2 + a_curve * Z_a^4 */
    }

    /* Z_r */
    if (a->Z_is_one) {
        bn192_copy(n0, a->Y);
    } else {
        bn192_nist_mod_192_mul(n0, a->Y, a->Z, group->field);
    }
    bn192_mod_lshift1_quick(r->Z, n0, group->field);
    r->Z_is_one = 0;
    /* Z_r = 2 * Y_a * Z_a */

    /* n2 */
    bn192_nist_mod_192_sqr(n3, a->Y, group->field);
    bn192_nist_mod_192_mul(n2, a->X, n3, group->field);
    bn192_mod_lshift_quick(n2, n2, 2, group->field);
    /* n2 = 4 * X_a * Y_a^2 */

    /* X_r */
    bn192_mod_lshift1_quick(n0, n2, group->field);
    bn192_nist_mod_192_sqr(rX, n1, group->field);
    bn192_mod_sub_quick(r->X, rX, n0, group->field);
    /* X_r = n1^2 - 2 * n2 */

    /* n3 */
    bn192_nist_mod_192_sqr(n0, n3, group->field);
    bn192_mod_lshift_quick(n3, n0, 3, group->field);
    /* n3 = 8 * Y_a^4 */

    /* Y_r */
    bn192_mod_sub_quick(n0, n2, r->X, group->field);
    bn192_nist_mod_192_mul(n0, n1, n0, group->field);
    bn192_mod_sub_quick(r->Y, n0, n3, group->field);
    /* Y_r = n1 * (n2 - X_r) - n3 */
}

// derived from ec_GFp_simple_add()
// in \openssl-1.1.1\crypto\ec\ecp_smpl.c
static void ec_point_add(
  const ec_group_t* group,
  ec_point_t*       r,
  const ec_point_t* a,
  const ec_point_t* b)
{
    if (a == b) {
        ec_point_dbl(group, r, a);
        return;
    }
    if (ec_point_is_at_infinity(a)) {
        *r = *b;
        return;
    }
    if (ec_point_is_at_infinity(b)) {
        *r = *a;
        return;
    }

    bn_t n0, n1, n2, n3, n4, n5, n6;

    /*
     * Note that in this function we must not read components of 'a' or 'b'
     * once we have written the corresponding components of 'r'. ('r' might
     * be one of 'a' or 'b'.)
     */

    /* n1, n2 */
    if (b->Z_is_one) {
        bn192_copy(n1, a->X);
        bn192_copy(n2, a->Y);
        /* n1 = X_a */
        /* n2 = Y_a */
    } else {
        bn192_nist_mod_192_sqr(n0, b->Z, group->field);
        bn192_nist_mod_192_mul(n1, a->X, n0, group->field);
        /* n1 = X_a * Z_b^2 */

        bn192_nist_mod_192_mul(n0, n0, b->Z, group->field);
        bn192_nist_mod_192_mul(n2, a->Y, n0, group->field);
        /* n2 = Y_a * Z_b^3 */
    }

    /* n3, n4 */
    if (a->Z_is_one) {
        bn192_copy(n3, b->X);
        bn192_copy(n4, b->Y);
        /* n3 = X_b */
        /* n4 = Y_b */
    } else {
        bn192_nist_mod_192_sqr(n0, a->Z, group->field);
        bn192_nist_mod_192_mul(n3, b->X, n0, group->field);
        /* n3 = X_b * Z_a^2 */

        bn192_nist_mod_192_mul(n0, n0, a->Z, group->field);
        bn192_nist_mod_192_mul(n4, b->Y, n0, group->field);
        /* n4 = Y_b * Z_a^3 */
    }

    /* n5, n6 */
    bn192_mod_sub_quick(n5, n1, n3, group->field);
    bn192_mod_sub_quick(n6, n2, n4, group->field);
    /* n5 = n1 - n3 */
    /* n6 = n2 - n4 */

    if (bn192_is_zero(n5)) {
        if (bn192_is_zero(n6)) {
            /* a is the same point as b */
            ec_point_dbl(group, r, a);
        } else {
            /* a is the inverse of b */
            bn192_zero(r->Z);
            r->Z_is_one = 0;
        }
        return;
    }

    /* 'n7', 'n8' */
    bn192_mod_add_quick(n1, n1, n3, group->field);
    bn192_mod_add_quick(n2, n2, n4, group->field);
    /* 'n7' = n1 + n3 */
    /* 'n8' = n2 + n4 */

    /* Z_r */
    if (a->Z_is_one && b->Z_is_one) {
        bn192_copy(r->Z, n5);
    } else {
        if (a->Z_is_one) {
            bn192_copy(n0, b->Z);
        } else if (b->Z_is_one) {
            bn192_copy(n0, a->Z);
        } else {
            bn192_nist_mod_192_mul(n0, a->Z, b->Z, group->field);
        }
        bn192_nist_mod_192_mul(r->Z, n0, n5, group->field);
    }
    r->Z_is_one = 0;
    /* Z_r = Z_a * Z_b * n5 */

    /* X_r */
    bn192_nist_mod_192_sqr(n0, n6, group->field);
    bn192_nist_mod_192_sqr(n4, n5, group->field);
    bn192_nist_mod_192_mul(n3, n1, n4, group->field);
    bn192_mod_sub_quick(r->X, n0, n3, group->field);
    /* X_r = n6^2 - n5^2 * 'n7' */

    /* 'n9' */
    bn192_mod_lshift1_quick(n0, r->X, group->field);
    bn192_mod_sub_quick(n0, n3, n0, group->field);
    /* n9 = n5^2 * 'n7' - 2 * X_r */

    /* Y_r */
    bn192_nist_mod_192_mul(n0, n0, n6, group->field);
    bn192_nist_mod_192_mul(n5, n4, n5, group->field);
    bn192_nist_mod_192_mul(n1, n2, n5, group->field);
    bn192_mod_sub_quick(n0, n0, n1, group->field);
    if (bn192_is_odd(n0)) {
      bn192_add_rshift1(r->Y, n0, group->field);
    } else {
      /* now  0 <= n0 < 2*p,  and n0 is even */
      bn192_rshift1(r->Y, n0);
    }
    /* Y_r = (n6 * 'n9' - 'n8' * 'n5^3') / 2 */
    return;
}

// derived from ec_GFp_simple_point_get_affine_coordinates()
// in \openssl-1.1.1\crypto\ec\ecp_smpl.c
static int ec_point_get_affine_coordinates(
  const ec_group_t* group,
  const ec_point_t* point,
  bn_t              x,
  bn_t              y)
{
    bn_t Z_1, Z_2, Z_3;
    const bn_word_t *Z_;

    if (ec_point_is_at_infinity(point)) {
        ECerr(EC_F_EC_GFP_SIMPLE_POINT_GET_AFFINE_COORDINATES,
              EC_R_POINT_AT_INFINITY);
        return 0;
    }

    /* transform  (X, Y, Z)  into  (x, y) := (X/Z^2, Y/Z^3) */

  //  if (group->meth->field_decode) {
  //      if (!group->meth->field_decode(group, Z, point->Z))
  //          goto err;
  //      Z_ = Z;
  //  } else {
        Z_ = point->Z;
  //  }

    if (bn192_is_one(Z_)) {
  //      if (group->meth->field_decode) {
  //          if (x != NULL) {
  //              if (!group->meth->field_decode(group, x, point->X))
  //                  goto err;
  //          }
  //          if (y != NULL) {
  //              if (!group->meth->field_decode(group, y, point->Y))
  //                  goto err;
  //          }
  //      } else {
           if (x != NULL) {
               bn192_copy(x, point->X);
           }
           if (y != NULL) {
               bn192_copy(y, point->Y);
           }
  //      }
    } else {
        bn192_mod_inverse(Z_1, Z_, group->field);
  //      if (group->meth->field_encode == 0) {
  //          /* field_sqr works on standard representation */
  //          if (!group->meth->field_sqr(group, Z_2, Z_1))
  //              goto err;
  //      } else {
            bn192_nist_mod_192_sqr(Z_2, Z_1, group->field);
  //      }

        if (x != NULL) {
            /*
             * in the Montgomery case, field_mul will cancel out Montgomery
             * factor in X:
             */
            bn192_nist_mod_192_mul(x, point->X, Z_2, group->field);
        }

        if (y != NULL) {
  //          if (group->meth->field_encode == 0) {
  //              /*
  //               * field_mul works on standard representation
  //               */
  //              if (!group->meth->field_mul(group, Z_3, Z_2, Z_1))
  //                  goto err;
  //          } else {
                bn192_nist_mod_192_mul(Z_3, Z_2, Z_1, group->field);
  //          }

            /*
             * in the Montgomery case, field_mul will cancel out Montgomery
             * factor in Y:
             */
            bn192_nist_mod_192_mul(y, point->Y, Z_3, group->field);
        }
    }
    return 1;
}


/** Computes r = q1 * m1 + q2 * m2
 *  \param  group  underlying ec_group_t object
 *  \param  r      ec_point_t object for result
 *                 (can't be the same as q1 or q2)
 *  \param  q1     ec_point_t object to be multiplied by m1
 *  \param  m1     BIGNUM positive multiplication factor for q1
 *  \param  q2     ec_point_t object to be multiplied by m2
 *  \param  m2     BIGNUM non-negative multiplication factor for q2
 */
static void ec_point_muladd2(
  const ec_group_t* group,
  ec_point_t*       r,
  const ec_point_t* q1,
  const bn_t        m1,
  const ec_point_t* q2,
  const bn_t        m2)
{
  const bn_word_t* m_oct[2] = {m1, m2};
  const ec_point_t* q[2] = { q1, q2 };
  int nz = 0;
  for (int bit_i = ECDSA_NBITS-1; bit_i >= 0; --bit_i) {
    if (nz)
      ec_point_dbl(group, r, r);

    for (int k = 0; k < 2; ++k) {
      if ((m_oct[k][bit_i/BN192LIB_BITS_PER_WORD] >> (bit_i%BN192LIB_BITS_PER_WORD)) & 1) {
        // bit is set
        if (nz) {
          ec_point_add(group, r, r, q[k]);
        } else {
          *r = *q[k];
        }
        nz = 1;
      }
    }
  }
}


static int ossl_ecdsa_verify_sig(
  const bn_t       digest,
  const bn_word_t* sig_r,
  const bn_word_t* sig_s)
{
    const ec_group_t* group   = &st_group;
    const ec_point_t* pub_key = &st_pub_key;
    const bn_word_t* order = group->order;

    if (bn192_is_zero(sig_r) || bn192_ucmp(sig_r, order) >= 0 ||
        bn192_is_zero(sig_s) || bn192_ucmp(sig_s, order) >= 0)
    {
        return 0;                /* signature is invalid */
    }

    // calculate tmp1 = inv(S) mod order
    bn_t u2;
    bn192_mod_inverse(u2, sig_s, group->order);

    //
    // There is no need to truncate digest,
    // because it is always supplied by caller
    // with correct length
    //

    // u1 = digest * tmp mod order
    bn_t u1;
    bn192_nist_mod_192_mul(u1, digest, u2, order);
    // u2 = r * w mod q
    bn192_nist_mod_192_mul(u2, sig_r, u2, order);

    // point = generator*u1 + pub_key*u2
    ec_point_t point;
    ec_point_muladd2(group, &point,
      &group->generator, u1,
      pub_key,           u2);

    bn_t X;
    if (!ec_point_get_affine_coordinates(group, &point, X, NULL)) {
        ECerr(EC_F_OSSL_ECDSA_VERIFY_SIG, ERR_R_EC_LIB);
        return -1; // internal error (point==Inf, it shouldn't happen)
    }

    bn192_nist_mod_192(u1, X, order);
    /*  if the signature is correct u1 is equal to sig_r */
    return (bn192_ucmp(u1, sig_r) == 0);
}

int uut_verify(
 const unsigned char digest_be[24],
 const unsigned char signature_octets[2][24],
 int* res)
{
  *res = 0;

  // digest -> m, reverse bytes order
  bn_t digest;
  for (int i = 0; i < ECDSA_NBYTES; ++i)
    ((uint8_t*)digest)[ECDSA_NBYTES-1-i] = digest_be[i]; // non-portable, depends on LE byte order
  const bn_word_t* sig_r = (const bn_word_t*)signature_octets[0]; // non-portable, depends on LE byte order
  const bn_word_t* sig_s = (const bn_word_t*)signature_octets[1]; // non-portable, depends on LE byte order
  int ret = ossl_ecdsa_verify_sig(digest, sig_r, sig_s);
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
