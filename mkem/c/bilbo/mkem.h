#pragma once

#include <stdint.h>
#include "params.h"
#include "namespace.h"

#define MKEM_SYMBYTES 16

// Group is simply the seed of A.
#define MKEM_GROUPBYTES BYTES_SEED_A

// Matrix s and mkem seed
#define MKEM_SECRETKEYBYTES 10256

// Matrix B
#define MKEM_PUBLICKEYBYTES 10240

#define MKEM_SSBYTES (MKEM_SYMBYTES)

#define MKEM_CTIBYTES 10240
#define MKEM_CTDBYTES 24


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
