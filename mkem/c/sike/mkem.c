#include <string.h>

#include "mkem.h"
#include "internal.h"
#include "randombytes.h"
#include "fips202.h"
#include "verify.h"
#include "util.h"

void mkem_group(uint8_t group[MKEM_GROUPBYTES])
{
    uint8_t buf[MKEM_GROUPBYTES+1];
    buf[0] = 4;
    randombytes(buf+1, MKEM_GROUPBYTES);
    shake128(group, MKEM_GROUPBYTES, buf, sizeof(buf));
}

void mkem_keypair(
    uint8_t pk[MKEM_PUBLICKEYBYTES],
    uint8_t sk[MKEM_SECRETKEYBYTES],
    const uint8_t group[MKEM_GROUPBYTES]
) {
    (void)(group);
    crypto_kem_keypair(pk, sk);
    randombytes(sk + CRYPTO_SECRETKEYBYTES, MKEM_SYMBYTES); // mkem seed
}

static void _mkem_enc(
    uint8_t **ctds,
    uint8_t cti[MKEM_CTIBYTES],
    const uint8_t group[MKEM_GROUPBYTES],
    const uint8_t **pks,
    const unsigned int len,
    const uint8_t pre_ss[MKEM_SSBYTES]
) {
    unsigned int i, j;
    uint8_t esk_seed[1 + MKEM_GROUPBYTES + MKEM_SSBYTES];
    uint8_t esk[SECRETKEY_A_BYTES];

    // Generate ephemeral secret key
    esk_seed[0] = 0;
    memcpy(esk_seed + 1, pre_ss, MKEM_SSBYTES);
    memcpy(esk_seed + 1 + MKEM_SSBYTES, group, MKEM_GROUPBYTES);
    shake128(
        (uint8_t*)esk,
        SECRETKEY_A_BYTES,
        esk_seed,
        sizeof(esk_seed)
    );
    esk[SECRETKEY_A_BYTES - 1] &= MASK_ALICE;

    EphemeralKeyGeneration_A(esk, cti);

    for (j = 0; j < len; j++) {
        uint8_t jinv[FP2_ENCODED_BYTES];
        uint8_t hjinv[MKEM_SYMBYTES];

        EphemeralSecretAgreement_A(esk, pks[j], jinv);
        shake128(hjinv, sizeof(hjinv), jinv, sizeof(jinv));
        for (i = 0; i < MKEM_SYMBYTES; i++)
            ctds[j][i] = pre_ss[i] ^ hjinv[i];
    }
}


void mkem_enc(
    uint8_t **ctds,
    uint8_t cti[MKEM_CTIBYTES],
    uint8_t ss[MKEM_SSBYTES],
    const uint8_t group[MKEM_GROUPBYTES],
    const uint8_t **pks,
    const unsigned int len
) {
    uint8_t buf[MKEM_SSBYTES+1];
    uint8_t* pre_ss = buf + 1;

    randombytes(pre_ss, MKEM_SSBYTES);
    _mkem_enc(ctds, cti, group, pks, len, pre_ss);

    buf[0] = 2;
    shake128(ss, MKEM_SSBYTES, buf, sizeof(buf));
}

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
) {
    unsigned int i;
    int ok;

    uint8_t ss_buf[MKEM_SSBYTES+1] = {0};
    uint8_t *pre_ss2 = ss_buf+1;
    uint8_t jinv[FP2_ENCODED_BYTES];
    uint8_t hjinv[MKEM_SYMBYTES];

    uint8_t ctip[MKEM_CTIBYTES];
    uint8_t ctdp[MKEM_CTIBYTES];
    uint8_t *ctdp_ptr = ctdp;
    uint8_t ssp_buf[1+MKEM_SYMBYTES+MKEM_CTIBYTES+MKEM_CTDBYTES];
    const uint8_t *seed = sk + CRYPTO_SECRETKEYBYTES;

    EphemeralSecretAgreement_B(sk + MSG_BYTES, cti, jinv);
    shake128(hjinv, sizeof(hjinv), jinv, sizeof(jinv));
    for (i = 0; i < MKEM_SYMBYTES; i++)
        pre_ss2[i] = ctd[i] ^ hjinv[i];

    _mkem_enc(&ctdp_ptr, ctip, group, &pk, 1, pre_ss2);

    // Compute shared secret in case of decryption error
    ssp_buf[0] = 3;
    memcpy(ssp_buf+1, seed, MKEM_SYMBYTES);
    memcpy(ssp_buf+1+MKEM_SYMBYTES, cti, MKEM_CTIBYTES);
    memcpy(ssp_buf+1+MKEM_SYMBYTES+MKEM_CTIBYTES, ctd, MKEM_CTDBYTES);
    shake128(ss, MKEM_SYMBYTES, ssp_buf, sizeof(ssp_buf));

    // Compute actual shared secret
    ss_buf[0] = 2;
    shake128(ss_buf, MKEM_SYMBYTES, ss_buf, sizeof(ss_buf));
    ok = 1 - (verify(ctdp, ctd, MKEM_CTDBYTES)
                | verify(ctip, cti, MKEM_CTIBYTES));
    cmov(ss, ss_buf, MKEM_SSBYTES, ok);
}
