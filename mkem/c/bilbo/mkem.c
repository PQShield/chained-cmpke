#include <string.h>
#include <stdio.h>

#include "namespace.h"
#include "mkem.h"
#include "randombytes.h"
#include "fips202.h"
#include "verify.h"
#include "common.h"

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
    unsigned int i;

    uint8_t seed_SE[MKEM_SYMBYTES+1];

    uint32_t S_rnd[2*PARAMS_N*PARAMS_NBAR] = {0};
    uint32_t *E_rnd = &S_rnd[PARAMS_N*PARAMS_NBAR];

    uint16_t  B[PARAMS_NBAR*PARAMS_N] = {0};
    uint16_t  S[PARAMS_NBAR*PARAMS_N] = {0};
    uint16_t  E[PARAMS_NBAR*PARAMS_N] = {0};
    
    seed_SE[0] = 0x5f;
    randombytes(seed_SE+1, MKEM_SYMBYTES);
    shake128(
        (uint8_t*)S_rnd,
        2*PARAMS_NBAR*PARAMS_N*sizeof(uint32_t),
        seed_SE,
        sizeof(seed_SE)
    );

    for (i = 0; i < 2*PARAMS_N*PARAMS_NBAR; i++)
        S_rnd[i] = PQCLEAN_FRODOKEM640SHAKE_OPT_LE_TO_UINT32(S_rnd[i]);
    
    PQCLEAN_FRODOKEM640SHAKE_OPT_sample_n(S, S_rnd, PARAMS_N * PARAMS_NBAR);
    PQCLEAN_FRODOKEM640SHAKE_OPT_sample_n(E, E_rnd, PARAMS_N * PARAMS_NBAR);
    PQCLEAN_FRODOKEM640SHAKE_OPT_mul_add_as_plus_e(B, S, E, group);

    // Compute B = A S + E
    PQCLEAN_FRODOKEM640SHAKE_OPT_pack(pk, MKEM_PUBLICKEYBYTES, B,
            PARAMS_NBAR*PARAMS_N, PARAMS_LOGQ);
    
    for (i = 0; i < PARAMS_N*PARAMS_NBAR; i++)
        S[i] = PQCLEAN_FRODOKEM640SHAKE_OPT_UINT16_TO_LE(S[i]);

    memcpy(sk, S, PARAMS_NBAR*PARAMS_N*sizeof(uint16_t));
    randombytes(sk + PARAMS_NBAR*PARAMS_N*sizeof(uint16_t), MKEM_SYMBYTES);
}

static void _mkem_enc(
    uint8_t **ctds,
    uint8_t cti[MKEM_CTIBYTES],
    const uint8_t group[MKEM_GROUPBYTES],
    const uint8_t **pks,
    const unsigned int len,
    uint8_t pre_ss[MKEM_SSBYTES]
) {
    unsigned int i, j;
    uint8_t shared_rnd_buf[1 + MKEM_GROUPBYTES + MKEM_SSBYTES];
    uint32_t shared_rnd[PARAMS_MBAR*PARAMS_N*2] = {0};
    uint32_t* Sp_rnd = shared_rnd;
    uint32_t* Ep_rnd = shared_rnd + PARAMS_MBAR*PARAMS_N;
    uint16_t Sp[PARAMS_MBAR*PARAMS_N] = {0};
    uint16_t Ep[PARAMS_MBAR*PARAMS_N] = {0};
    uint16_t Bp[PARAMS_MBAR*PARAMS_N] = {0};

    // Sample S' and E'
    shared_rnd_buf[0] = 0;
    memcpy(shared_rnd_buf + 1, pre_ss, MKEM_SSBYTES);
    memcpy(shared_rnd_buf + 1 + MKEM_SSBYTES, group, MKEM_GROUPBYTES);
    shake128(
        (uint8_t*)shared_rnd,
        PARAMS_MBAR*PARAMS_N*2*sizeof(uint32_t), 
        shared_rnd_buf,
        sizeof(shared_rnd_buf)
    );

    for (i = 0; i < PARAMS_MBAR*PARAMS_N*2; i++)
        shared_rnd[i] = PQCLEAN_FRODOKEM640SHAKE_OPT_LE_TO_UINT32(shared_rnd[i]);

    PQCLEAN_FRODOKEM640SHAKE_OPT_sample_n(Sp, Sp_rnd, PARAMS_N * PARAMS_MBAR);
    PQCLEAN_FRODOKEM640SHAKE_OPT_sample_n(Ep, Ep_rnd, PARAMS_N * PARAMS_MBAR);

    // Compute B' = S' A + E'
    PQCLEAN_FRODOKEM640SHAKE_OPT_mul_add_sa_plus_e(Bp, Sp, Ep, group);

    PQCLEAN_FRODOKEM640SHAKE_OPT_pack(cti, MKEM_CTIBYTES, Bp,
        PARAMS_N * PARAMS_MBAR, PARAMS_LOGQ);

    for (j = 0; j < len; j++) {
        uint8_t rnd_buf[1 + MKEM_GROUPBYTES + MKEM_SSBYTES + MKEM_PUBLICKEYBYTES];
        uint32_t rnd[PARAMS_MBAR * PARAMS_NBAR] = {0};
        uint16_t Epp[PARAMS_NBAR*PARAMS_MBAR] = {0};
        uint16_t B[PARAMS_N*PARAMS_NBAR] = {0};
        uint16_t V[PARAMS_MBAR*PARAMS_NBAR] = {0};
        uint16_t C[PARAMS_MBAR*PARAMS_NBAR] = {0};

        // Sample E''
        rnd_buf[0] = 1;
        memcpy(rnd_buf + 1, pre_ss, MKEM_SSBYTES);
        memcpy(rnd_buf + 1 + MKEM_SSBYTES, pks[j], MKEM_PUBLICKEYBYTES);
        memcpy(rnd_buf + 1 + MKEM_SSBYTES + MKEM_PUBLICKEYBYTES,
                group, MKEM_GROUPBYTES);

        shake128(
            (uint8_t*)rnd,
            PARAMS_NBAR*PARAMS_MBAR*sizeof(uint32_t),
            rnd_buf,
            sizeof(rnd_buf)
        );

        for (i = 0; i < PARAMS_MBAR*PARAMS_NBAR; i++)
            rnd[i] = PQCLEAN_FRODOKEM640SHAKE_OPT_LE_TO_UINT32(rnd[i]);

        PQCLEAN_FRODOKEM640SHAKE_OPT_sample_n(Epp, rnd, PARAMS_NBAR*PARAMS_MBAR);

        // Compute  V = S' B + E''
        PQCLEAN_FRODOKEM640SHAKE_OPT_unpack(B, PARAMS_N*PARAMS_NBAR, pks[j],
                MKEM_PUBLICKEYBYTES, PARAMS_LOGQ);
        PQCLEAN_FRODOKEM640SHAKE_OPT_mul_add_sb_plus_e(V, B, Sp, Epp);
        
        // c₂ = compress(V + encode(μ))
        PQCLEAN_FRODOKEM640SHAKE_OPT_key_encode(C, (uint16_t*)pre_ss);
        PQCLEAN_FRODOKEM640SHAKE_OPT_add(C, V, C);

        // Kyber compression
        //
        //  ⌈(2ᵈ/q)x⌋ mod⁺ 2ᵈ = ⌊((x << d) + q/2) / q⌋ mod⁺ 2ᵈ
        //
        // where d = DC2.
        for (int i = 0; i < PARAMS_MBAR*PARAMS_NBAR; i++)
            C[i] = (((C[i] << PARAMS_DC2) + (1 << (PARAMS_LOGQ - 1)))
                    >> PARAMS_LOGQ) & ((1 << PARAMS_DC2) - 1);


        PQCLEAN_FRODOKEM640SHAKE_OPT_pack(ctds[j], MKEM_CTDBYTES, C,
                PARAMS_MBAR*PARAMS_NBAR, PARAMS_DC2);
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
    unsigned int i;
    int ok;

    uint16_t S[PARAMS_N*PARAMS_NBAR] = {0};
    uint16_t Bp[PARAMS_MBAR*PARAMS_N] = {0};
    uint16_t C[PARAMS_MBAR*PARAMS_NBAR] = {0};
    uint16_t W[PARAMS_MBAR*PARAMS_NBAR] = {0};
    uint8_t ss_buf[MKEM_SSBYTES+1] = {0};
    uint8_t *mup = ss_buf+1;

    uint8_t ctip[MKEM_CTIBYTES];
    uint8_t ctdp[MKEM_CTIBYTES];
    uint8_t *ctdp_ptr = ctdp;
    uint8_t ssp_buf[1+MKEM_SYMBYTES+MKEM_CTIBYTES+MKEM_CTDBYTES];
    const uint8_t *seed = sk + 2*PARAMS_NBAR*PARAMS_N*sizeof(uint16_t);

    for (i = 0; i < PARAMS_NBAR*PARAMS_N; i++)
        S[i] = sk[2*i] | (sk[2*i+1] << 8);

    PQCLEAN_FRODOKEM640SHAKE_OPT_unpack(Bp, PARAMS_MBAR*PARAMS_N, cti,
            MKEM_CTIBYTES, PARAMS_LOGQ);
    PQCLEAN_FRODOKEM640SHAKE_OPT_unpack(C, PARAMS_MBAR*PARAMS_NBAR, ctd,
            MKEM_CTDBYTES, PARAMS_DC2);

    // Kyber decompression for C
    // 
    //   ⌈(q/2ᵈ)x⌋ = (qx + (1<<(d-1))) >> d = x << (log₂(q) - d)
    for (i = 0; i < PARAMS_NBAR*PARAMS_MBAR; i++)
        C[i] <<= PARAMS_LOGQ - PARAMS_DC2;

    // Compute μ' = Decode(C - B' S)
    PQCLEAN_FRODOKEM640SHAKE_OPT_mul_bs(W, Bp, S);
    PQCLEAN_FRODOKEM640SHAKE_OPT_sub(W, C, W);
    PQCLEAN_FRODOKEM640SHAKE_OPT_key_decode((uint16_t*)mup, W);

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
