#include <string.h>
#include <stdio.h>

#include "mkem.h"
#include "randombytes.h"
#include "fips202.h"
#include "verify.h"
#include "util.h"
#include "ntrulpr.h"

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
    int16_t ag[p];
    int8_t a[p];

    generator(ag, group);
    short_random(a);
    mult_small_ntt(ag, ag, a);

    crypto_encode_px3(sk, a);

    crypto_round_and_encode(pk, ag);

    randombytes(sk + (MKEM_SECRETKEYBYTES - MKEM_SYMBYTES), MKEM_SYMBYTES);
}

static void _mkem_enc(
    uint8_t **ctds,
    uint8_t cti[MKEM_CTIBYTES],
    const uint8_t group[MKEM_GROUPBYTES],
    const uint8_t **pks,
    const unsigned int len,
    uint8_t pre_ss[MKEM_SSBYTES]
) {
    unsigned int j;
    uint8_t shared_rnd_buf[1 + MKEM_GROUPBYTES + MKEM_SSBYTES];
    uint32_t shared_rnd[p] = {0};
    int8_t b[p];
    int16_t bg[p];
    int8_t r[I];

    crypto_decode_Ix2(r, pre_ss);

    // Sample b
    shared_rnd_buf[0] = 0;
    memcpy(shared_rnd_buf + 1, pre_ss, MKEM_SSBYTES);
    memcpy(shared_rnd_buf + 1 + MKEM_SSBYTES, group, MKEM_GROUPBYTES);
    shake128(
        (uint8_t*)shared_rnd,
        p*sizeof(uint32_t),
        shared_rnd_buf,
        sizeof(shared_rnd_buf)
    );

    short_fromlist(b, shared_rnd);

    // Compute B := Round(b G)
    generator(bg, group);

    mult_small_ntt(bg, bg, b);

    crypto_round_and_encode(cti, bg);

    for (j = 0; j < len; j++) {
        int16_t ba[p];
        int8_t t[I];

        // Compute b A
        crypto_decode_rounded(ba, pks[j]);
        mult_small_ntt(ba, ba, b);

        // Compute T
        top_add_r_q12(t, ba, r);
        crypto_encode_Ixtau(ctds[j], t);
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
    uint8_t buf[MKEM_SSBYTES + 1];
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
    int ok;

    int16_t ab[p];
    int8_t a[p];
    int8_t t[p];
    int8_t r[I];

    uint8_t ss_buf[MKEM_SSBYTES+1] = {0};
    uint8_t *mup = ss_buf+1;

    uint8_t ctip[MKEM_CTIBYTES];
    uint8_t ctdp[MKEM_CTIBYTES];
    uint8_t *ctdp_ptr = ctdp;
    uint8_t ssp_buf[1+MKEM_SYMBYTES+MKEM_CTIBYTES+MKEM_CTDBYTES];
    const uint8_t *seed = sk + (MKEM_SECRETKEYBYTES - MKEM_SYMBYTES);

    // Compute a B
    crypto_decode_rounded(ab, cti);

    crypto_decode_px3(a, sk);
    mult_small_ntt(ab, ab, a);

    // Compute r
    crypto_decode_Ixtau(t, ctd);
    sign_right_sub_4w1(r, t, ab);
    crypto_encode_Ix2(mup, r);

    _mkem_enc(&ctdp_ptr, ctip, group, &pk, 1, mup);

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
