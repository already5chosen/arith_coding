#include "bn192_lib.h"

void secp192k1_ecdsa_set_public_key(const bn_t key[2]);
int secp192k1_ecdsa_verify(const bn_t digest, const bn_t signature[2]);
