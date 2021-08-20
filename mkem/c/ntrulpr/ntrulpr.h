#pragma once

#include <stdint.h>

#define p 757
#define q 7879
#define w 242
#define tau 4
#define tau0 3011
#define tau1 33
#define tau2 1995
#define tau3 1978
#define I 128
#define q12 3939


void crypto_decode_pxint32(uint32_t *x, const uint8_t *s);
void crypto_decode_pxint16(uint16_t *x, const unsigned char *s);
void crypto_encode_pxint16(unsigned char *s, const uint16_t *x);
void crypto_decode_Ixtau(int8_t *x, const unsigned char *s);
void crypto_encode_Ixtau(unsigned char *s, const int8_t *x);
void crypto_sort_int32(int32_t *x, long long n);
void crypto_sort_uint32(uint32_t *x, long long n);
void crypto_decode_Ix2(int8_t *x, const unsigned char *s);
void crypto_encode_Ix2(unsigned char *s, const int8_t *x);
void crypto_decode_rounded(int16_t *R0,const unsigned char *s);
void crypto_encode_rounded(unsigned char *out,const int16_t *R0);
void crypto_round_and_encode(unsigned char *out, const int16_t *a);
void crypto_encode_px3(unsigned char *s, const int8_t *f);
void crypto_decode_px3(int8_t* v, const unsigned char *s);
void mult_small(int16_t* f, int8_t *g);
void mult_small_ntt(int16_t*r, int16_t* f, int8_t *g);
void mult_small_prefix(int16_t* r, int16_t* f, int8_t *g);
void generator(int16_t *g, const uint8_t* seed);
void short_random(int8_t *out);
void top_add_r_q12(int8_t *t, const int16_t *x, const int8_t *r);
void short_fromlist(int8_t *out, const uint32_t* in);
void sign_right_sub_4w1(int8_t *r, const int8_t *t, const int16_t *x);
