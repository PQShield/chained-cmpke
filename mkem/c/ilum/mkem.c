#include "mkem.h"
#include "verify.h"
#include "polyvec.h"
#include "symmetric.h"
#include "indcpa.h"
#include "randombytes.h"

// TODO - symmetric sizes; is 128 bit fine everywhere?
//      - only sample the first half of e₂

void mkem_group(uint8_t group[MKEM_GROUPBYTES])
{
    randombytes(group, KYBER_SYMBYTES);
    hash_h(group, group, KYBER_SYMBYTES); // don't expose system RNG
}

void mkem_keypair(
    uint8_t pk[MKEM_PUBLICKEYBYTES],
    uint8_t sk[MKEM_SECRETKEYBYTES],
    const uint8_t group[MKEM_GROUPBYTES]
) {
    unsigned int i;

    uint8_t sigma[KYBER_SYMBYTES];

    uint8_t nonce = 0;
    polyvec a[KYBER_K], e, pkpv, skpv;

    randombytes(sigma, KYBER_SYMBYTES);

    gen_matrix(a, group, 0);

    for(i=0;i<KYBER_K;i++)
        poly_getnoise_eta1(&skpv.vec[i], sigma, nonce++);
    for(i=0;i<KYBER_K;i++)
        poly_getnoise_eta1(&e.vec[i], sigma, nonce++);

    polyvec_ntt(&skpv);
    polyvec_ntt(&e);

    // matrix-vector multiplication
    for(i=0;i<KYBER_K;i++) {
        polyvec_basemul_acc_montgomery(&pkpv.vec[i], &a[i], &skpv);
        poly_tomont(&pkpv.vec[i]);
    }

    polyvec_add(&pkpv, &pkpv, &e);
    polyvec_reduce(&pkpv);

    polyvec_tobytes(sk, &skpv);
    polyvec_tobytes(pk, &pkpv);
    randombytes(sk + KYBER_POLYVECBYTES, KYBER_SYMBYTES); // mkem seed
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
    uint8_t nonce = 0;
    uint8_t shared_coins[KYBER_SYMBYTES];

    poly dss; // decompressed shared secret
    polyvec at[KYBER_K]; // Aᵀ
    polyvec rs, e1;
    polyvec u;

    // Compute shared randomness as G₁(m) = SHAKE-256(0‖m‖ρ)
    uint8_t sc_buf[KYBER_SYMBYTES + MKEM_GROUPBYTES + 1];
    sc_buf[0] = 0;
    for (i = 0; i < KYBER_SYMBYTES; i++)
        sc_buf[i+1] = pre_ss[i];
    for (i = 0; i < MKEM_GROUPBYTES; i++)
        sc_buf[i+1+KYBER_SYMBYTES] = group[i];
    shake256(shared_coins, KYBER_SYMBYTES, sc_buf, sizeof(sc_buf));

    // Generate A
    gen_matrix(at, group, 1);

    // Decompress shared secret
    poly_frommsg(&dss, pre_ss); // dss = Decompress(Decode_1(ss), 1) 

    // Generate vectors r and e.
    for (i = 0; i < KYBER_K; i++)
        poly_getnoise_eta1(rs.vec+i, shared_coins, nonce++);
    polyvec_ntt(&rs);
    for (i = 0; i < KYBER_K; i++)
        poly_getnoise_eta2(e1.vec+i, shared_coins, nonce++);

    // Compute u = Aᵀ r + e₁
    for (i = 0; i < KYBER_K; i++)
        polyvec_basemul_acc_montgomery(&u.vec[i], &at[i], &rs);
    polyvec_invntt_tomont(&u);
    polyvec_add(&u, &u, &e1);
    polyvec_reduce(&u);
    polyvec_compress(cti, &u);

    // Now for the per public key computation.
    for (i = 0; i < len; i++) {
        uint8_t coins[KYBER_SYMBYTES];
        uint8_t c_buf[KYBER_SYMBYTES + MKEM_PUBLICKEYBYTES
                + MKEM_GROUPBYTES + 1];
        poly e2, v;
        polyvec t;

        // Compute coins
        c_buf[0] = 1;
        for (j = 0; j < KYBER_SYMBYTES; j++)
            c_buf[1+j] = pre_ss[j];
        for (j = 0; j < MKEM_PUBLICKEYBYTES; j++)
            c_buf[1+j+KYBER_SYMBYTES] = pks[i][j];
        for (j = 0; j < MKEM_GROUPBYTES; j++)
            c_buf[1+j+KYBER_SYMBYTES+MKEM_PUBLICKEYBYTES] = group[j];
        shake256(coins, KYBER_SYMBYTES, c_buf, sizeof(c_buf));

        // Generate e₂
        poly_getnoise_eta2(&e2, coins, 0);

        // Unpack public key
        polyvec_frombytes(&t, pks[i]);

        // Compute v = <t, r> + e₂ + dss
        polyvec_basemul_acc_montgomery(&v, &t, &rs);
        poly_invntt_tomont(&v);
        poly_add(&v, &v, &e2);
        poly_add(&v, &v, &dss);
        poly_reduce(&v);
        poly_compress(ctds[i], &v);
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
    uint8_t pre_ss[MKEM_SSBYTES];

    // Generate shared secret
    randombytes(pre_ss, KYBER_SYMBYTES);
    hash_h(pre_ss, pre_ss, KYBER_SYMBYTES); // Don't expose system RNG
    hash_h(ss, pre_ss, KYBER_SYMBYTES);

    _mkem_enc(ctds, cti, group, pks, len, pre_ss);
}
//
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
    uint8_t ctd2[MKEM_CTDBYTES];
    uint8_t *ctd2_ptr = ctd2;
    uint8_t cti2[MKEM_CTIBYTES];
    uint8_t buf[MKEM_CTIBYTES + MKEM_CTDBYTES + KYBER_SYMBYTES];

    polyvec u, s;
    poly v, pre_ss2;
    uint8_t pre_ss[MKEM_SSBYTES];
    int ok;

    // Decompress ciphertext
    poly_decompress(&v, ctd);       // v is the pk-dependent part
    polyvec_decompress(&u, cti);    // bigger u is pk-independent

    // Unpack secret key
    polyvec_frombytes(&s, sk);

    // Compute pre_ss = Compress(v-<s,u>, 1)
    polyvec_ntt(&u);
    polyvec_basemul_acc_montgomery(&pre_ss2, &s, &u);
    poly_invntt_tomont(&pre_ss2);
    poly_sub(&pre_ss2, &v, &pre_ss2);
    poly_reduce(&pre_ss2);
    poly_tomsg(pre_ss, &pre_ss2);

    // Reencapsulate
    _mkem_enc(&ctd2_ptr, cti2, group, &pk, 1, pre_ss);

    // Set ss = H'(seed, ct)
    for (i = 0; i < MKEM_CTIBYTES; i++)
        buf[i] = cti[i];
    for (i = 0; i < MKEM_CTDBYTES; i++)
        buf[MKEM_CTIBYTES+i] = ctd[i];
    for (i = 0; i < KYBER_SYMBYTES; i++)
        buf[MKEM_CTIBYTES+MKEM_CTDBYTES+i] = sk[KYBER_POLYVECBYTES+i];

    ok = 1 - (verify(ctd2, ctd, MKEM_CTDBYTES)
                | verify(cti2, cti, MKEM_CTIBYTES));

    // Set ss = H(pre_ss) if ct2 = ct.
    hash_h(pre_ss, pre_ss, KYBER_SYMBYTES);
    cmov(ss, pre_ss, KYBER_SYMBYTES, ok);
}
