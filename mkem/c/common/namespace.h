#pragma once

#ifndef PREFIX
#error missing prefix
#endif

#define NS_CONCAT(A,B) A ## B
#define NS_CONCAT2(A,B) NS_CONCAT(A,B)
#define PREFIXED(s) NS_CONCAT2(PREFIX, s)

#define shake128 PREFIXED(shake128)
#define sha3_128 PREFIXED(sha3_128)
#define shake128_init PREFIXED(shake128_init)
#define shake128_absorb PREFIXED(shake128_absorb)
#define shake128_finalize PREFIXED(shake128_finalize)
#define shake128_squeeze PREFIXED(shake128_squeeze)
#define shake128_absorb_once PREFIXED(shake128_absorb_once)
#define shake128_squeezeblocks PREFIXED(shake128_squeezeblocksabsorb_once)

#define shake256 PREFIXED(shake256)
#define sha3_256 PREFIXED(sha3_256)
#define shake256_init PREFIXED(shake256_init)
#define shake256_absorb PREFIXED(shake256_absorb)
#define shake256_finalize PREFIXED(shake256_finalize)
#define shake256_squeeze PREFIXED(shake256_squeeze)
#define shake256_absorb_once PREFIXED(shake256_absorb_once)
#define shake256_squeezeblocks PREFIXED(shake256_squeezeblocksabsorb_once)

#define verify PREFIXED(verify)
#define cmov PREFIXED(cmov)
#define randombytes PREFIXED(randombytes)

#define mkem_keypair PREFIXED(mkem_keypair)
#define mkem_group PREFIXED(mkem_group)
#define mkem_enc PREFIXED(mkem_enc)
#define mkem_dec PREFIXED(mkem_dec)
