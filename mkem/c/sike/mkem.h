#pragma once

#include <stdint.h>
#include "sike.h"

#define MKEM_SYMBYTES 16

// We don't have any shared group data except that we use it to diversify
// the hashes.
#define MKEM_GROUPBYTES MKEM_SYMBYTES

#define MKEM_SECRETKEYBYTES (CRYPTO_SECRETKEYBYTES + MKEM_SYMBYTES)
#define MKEM_PUBLICKEYBYTES CRYPTO_PUBLICKEYBYTES

#define MKEM_SSBYTES MKEM_SYMBYTES

#define MKEM_CTDBYTES (MKEM_SYMBYTES)
#define MKEM_CTIBYTES (CRYPTO_CIPHERTEXTBYTES - MKEM_CTDBYTES)

// Generates a new group.
void mkem_group(uint8_t group[MKEM_GROUPBYTES]);

// Generates a new keypair for the given group.
void mkem_keypair(
    uint8_t pk[MKEM_PUBLICKEYBYTES],
    uint8_t sk[MKEM_SECRETKEYBYTES],
    const uint8_t group[MKEM_GROUPBYTES]
);

// Generates a shared secret, writes it to ss and encapsulates it
// for the len public keys given by pks.  The common part of all cipher
// texts is written to cti.  The public key dependent part of the
// ciphertext for pks[i] is written to cti[i], which must point to
// a buffer with MKEM_CTDBYTES of space.
void mkem_enc(
    uint8_t **ctds,
    uint8_t cti[MKEM_CTIBYTES],
    uint8_t ss[MKEM_SSBYTES],
    const uint8_t group[MKEM_GROUPBYTES],
    const uint8_t **pks,
    const unsigned int len
);

// Decapsulates the shared secret in the ciphertext (cti, ctd) for sk
// and write it to ss.  cti is the public key independent part of the cipher
// text and ctd is the public key dependent part.  
void mkem_dec(
    uint8_t ss[MKEM_SSBYTES],
    const uint8_t cti[MKEM_CTIBYTES],
    const uint8_t ctd[MKEM_CTDBYTES],
    const uint8_t group[MKEM_GROUPBYTES],
    const uint8_t pk[MKEM_PUBLICKEYBYTES],
    const uint8_t sk[MKEM_SECRETKEYBYTES]
);
