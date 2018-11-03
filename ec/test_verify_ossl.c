#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <openssl/ec.h>
#include <openssl/err.h>
#include <openssl/obj_mac.h>

static EC_KEY* st_eckey = NULL;

void uut_cleanup(void)
{
  if (st_eckey) {
    EC_KEY_free(st_eckey);
    st_eckey = 0;
  }
}

int uut_init(void) {
  st_eckey = EC_KEY_new_by_curve_name(NID_X9_62_prime192v1);
  if (!st_eckey)
    goto err;

  return 1;

  err:
  uut_cleanup();
  return 0;
}

int uut_set_public_key(unsigned char key[2][24])
{
  int ret = 0;
  BIGNUM* x = BN_lebin2bn(key[0], 24, NULL);
  if (x) {
    BIGNUM* y = BN_lebin2bn(key[1], 24, NULL);
    if (y) {
      ret = EC_KEY_set_public_key_affine_coordinates(st_eckey, x, y);
      BN_free(y);
    }
    BN_free(x);
  }
  return ret;
}

int uut_verify(
 const unsigned char digest_be[24],
 const unsigned char signature_octets[2][24],
 int* res)
{
  *res = 0;
  BIGNUM* sig_r = BN_lebin2bn(signature_octets[0], 24, 0);
  BIGNUM* sig_s = BN_lebin2bn(signature_octets[1], 24, 0);
  if (!sig_r || !sig_s) {
    fprintf(stderr, "BN_lebin2bn() fail. Error %lu\n", ERR_get_error());
    return 0;
  }

  ECDSA_SIG* sig = ECDSA_SIG_new();
  if (!sig) {
    fprintf(stderr, "ECDSA_SIG_new() fail. Error %lu\n", ERR_get_error());
    BN_free(sig_r);
    BN_free(sig_s);
    return 0;
  }

  if (!ECDSA_SIG_set0(sig, sig_r, sig_s)) {
    fprintf(stderr, "ECDSA_SIG_set0() fail. Error %lu\n", ERR_get_error());
    ECDSA_SIG_free(sig);
    BN_free(sig_r);
    BN_free(sig_s);
    return 0;
  }
  // now sig_r and sig_s are owned by sig

  int ret = ECDSA_do_verify(digest_be, 24, sig, st_eckey);
  if (ret >= 0) {
    *res = ret;
  } else {
    unsigned long err = ERR_get_error();
    char errstr[512];
    ERR_error_string_n(err, errstr, sizeof(errstr));
    fprintf(stderr, "ECDSA_do_verify: error %lu.\n%s\n", err, errstr);
  }
  ECDSA_SIG_free(sig);
  return (ret >= 0);
}
