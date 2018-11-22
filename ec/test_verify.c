// #include "ecs_locl.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <openssl/ec.h>
#include <openssl/obj_mac.h>
#include <openssl/err.h>

int uut_init(void);
int uut_set_public_key(unsigned char key[2][24]);
int uut_verify(
 const unsigned char digest[24],
 const unsigned char signature_octets[2][24],
 int* res);
void uut_cleanup(void);

static int test(
  BN_CTX* ctx,
  EC_KEY* eckey,
  unsigned char digest_octets[24])
{
  ECDSA_SIG* sig = ECDSA_do_sign(digest_octets, 24, eckey);
  if (!sig) {
    fprintf(stderr, "ECDSA_do_sign() fail.\n");
    return 0;
  }

  int ret = 0;
  const BIGNUM* r = ECDSA_SIG_get0_r(sig);
  const BIGNUM* s = ECDSA_SIG_get0_s(sig);
  if (r && s) {
    unsigned char signature_octets[2][24];
    if (BN_bn2lebinpad(r, signature_octets[0], 24) != 24) {
      fprintf(stderr, "BN_bn2lebinpad(r) fail.\n");
      goto done;
    }
    if (BN_bn2lebinpad(s, signature_octets[1], 24) != 24) {
      fprintf(stderr, "BN_bn2lebinpad(s) fail.\n");
      goto done;
    }

    int res_1;
    if (!uut_verify(digest_octets, signature_octets, &res_1)) {
      fprintf(stderr, "uut_verify() fail.\n");
      goto done;
    }

    int bit_i = rand() % (2*24*8);
    signature_octets[bit_i/(24*8)][(bit_i%(24*8))/8] ^= 1 << (bit_i %8);
    int res_2;
    if (!uut_verify(digest_octets, signature_octets, &res_2)) {
      fprintf(stderr, "uut_verify() fail.\n");
      goto done;
    }
    signature_octets[bit_i/(24*8)][(bit_i%(24*8))/8] ^= 1 << (bit_i %8);

    bit_i = rand() % (24*8);
    digest_octets[bit_i/8] ^= 1 << (bit_i %8);
    int res_3;
    if (!uut_verify(digest_octets, signature_octets, &res_3)) {
      fprintf(stderr, "uut_verify() fail.\n");
      goto done;
    }
    digest_octets[bit_i/8] ^= 1 << (bit_i %8);

    if (res_1 != 1) {
      fprintf(stderr, "uut_verify() failure instead of success.\n");
      goto done;
    }

    if (res_2 != 0) {
      fprintf(stderr, "uut_verify() success instead of failure.\n");
      goto done;
    }

    if (res_3 != 0) {
      fprintf(stderr, "uut_verify() success instead of failure.\n");
      goto done;
    }

    ret = 1;
    done:;
  } else {
    fprintf(stderr, "ECDSA_SIG_get0_r() fail.\n");
  }
  ECDSA_SIG_free(sig);
  return ret;
}

static void free_bn(BIGNUM **buf, int buflen) {
  for (int i = 0; i < buflen; ++i)
    BN_free(buf[i]);
}

static int allocate_bn(BIGNUM **buf, int buflen) {
  for (int i = 0; i < buflen; ++i) {
    BIGNUM *a = BN_new();
    if (!a) {
      free_bn(buf, i);
      return 0;
    }
    buf[i] = a;
  }
  return 1;
}

int main(int argz, char** argv)
{
  int nKeys    = 100;
  int nDigests = 100;
  if (argz > 1) {
    int val = strtol(argv[1], 0, 0);
    if (val > 0)
      nKeys = val;
    if (argz > 2) {
      val = strtol(argv[2], 0, 0);
      if (val > 0)
        nDigests = val;
    }
  }

  BN_CTX* ctx = BN_CTX_new();
  if (!ctx) {
    fprintf(stderr, "BN_CTX_new() fail.\n");
    return 1;
  }

  EC_KEY* eckey = EC_KEY_new_by_curve_name(NID_X9_62_prime192v1);
  if (eckey) {
    BIGNUM *bn[3];
    if (allocate_bn(bn, 3)) {
      if (uut_init()) {
        for (int key_i = 0; key_i < nKeys; ++key_i) {
          if (!EC_KEY_generate_key(eckey)) {
            fprintf(stderr, "EC_KEY_generate_key() fail at iteration %d.\n", key_i);
            goto done;
          }
          const EC_POINT* pub_key = EC_KEY_get0_public_key(eckey);
          if (!pub_key) {
            fprintf(stderr, "EC_KEY_get0_public_key() fail at iteration %d.\n", key_i);
            goto done;
          }
          const EC_GROUP* group = EC_KEY_get0_group(eckey);
          if (!group) {
            fprintf(stderr, "EC_KEY_get0_group() fail at iteration %d.\n", key_i);
            goto done;
          }
          BIGNUM *pub_x = bn[0];
          BIGNUM *pub_y = bn[1];
          if (!EC_POINT_get_affine_coordinates(
            group, pub_key, pub_x, pub_y, ctx)) {
            fprintf(stderr, "EC_POINT_get_affine_coordinates() fail.\n");
            goto done;
          }
          unsigned char pub_key_octets[2][24];
          if (BN_bn2lebinpad(pub_x, pub_key_octets[0], 24) != 24) {
            fprintf(stderr, "BN_bn2lebinpad(pub_x) fail.\n");
            goto done;
          }
          if (BN_bn2lebinpad(pub_y, pub_key_octets[1], 24) != 24) {
            fprintf(stderr, "BN_bn2lebinpad(pub_y) fail.\n");
            goto done;
          }
          if (!uut_set_public_key(pub_key_octets)) {
            fprintf(stderr, "uut_set_public_key() fail.\n");
            goto done;
          }

          for (int dgst_i = 0; dgst_i < nDigests; ++dgst_i) {
            BIGNUM *digest = bn[2];
            if (BN_rand(digest, 192, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY)) {
              unsigned char digest_octets[24];
              if (BN_bn2lebinpad(digest, digest_octets, 24) == 24) {
                if (!test(ctx, eckey, digest_octets)) {
                  fprintf(stderr, "Test fail at iteration %d.%d.\n", key_i, dgst_i);
                  goto done;
                }
              } else {
                fprintf(stderr, "BN_bn2lebinpad() fail at iteration %d.%d.\n", key_i, dgst_i);
                goto done;
              }
            } else {
              fprintf(stderr, "BN_rand() fail at iteration %d.%d.\n", key_i, dgst_i);
              goto done;
            }
          }
        }
        done:
        uut_cleanup();
      }
      free_bn(bn, 3);
    } else {
      fprintf(stderr, "BN_new() fail.\n");
    }
    EC_KEY_free(eckey);
  } else {
    fprintf(stderr, "EC_KEY_new_by_curve_name() fail.\n");
  }
  BN_CTX_free(ctx);
  return 0;
}

