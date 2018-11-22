#include <string.h>

#include "secp192k1_ecdsa_verify.h"

int uut_init(void)
{
  return 1;
}

int uut_set_public_key(unsigned char key_octets[2][24])
{
  bn_t key[2];
  memcpy(key, key_octets, sizeof(key));
  secp192k1_ecdsa_set_public_key(key);
  return 1;
}

int uut_verify(
 const unsigned char digest_be[24],
 const unsigned char signature_octets[2][24],
 int* pRes)
{
  // digest -> m, reverse bytes order
  bn_t digest;
  for (int i = 0; i < ECDSA_NBYTES; ++i)
    ((uint8_t*)digest)[ECDSA_NBYTES-1-i] = digest_be[i]; // non-portable, depends on LE byte order
  bn_t signature[2];
  memcpy(signature, signature_octets, sizeof(signature));
  int res = secp192k1_ecdsa_verify(digest, signature);
  if (pRes)
    *pRes = res;
  return res >= 0;
}

void uut_cleanup(void)
{
}
