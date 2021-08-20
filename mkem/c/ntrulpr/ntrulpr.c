#include "namespace.h"
#include "ntt_1907713_1536.h"
#include "ntrulpr.h"
#include "randombytes.h"
#include "fips202.h"

// This is a heavily modified version of ntrulpr653 from PQClean, which in turn
// is based on the factorized version from supercop.
//
// Changes:
//
//  - The "normal" parameters changed: p, q, w,  τ₀, τ₁, τ₂, τ₃.
//  - Accordingly we changed the strategy for Barrett reduction and
//    polynomial multiplication.
//  - We compute polynomial multiplication in GF(956929)/<2^1536+1>
//  - We changed tau to 4 (from 16) and I to 128 (from 256).
//  - We changed the Top to
// 
//      Top(x) = floor((τ₁ (x + τ₀) + 2^15)/2^16)
//
//    Plain NTRU Lpr has exponents (14, 15) instead of (15, 16).
//  - We use Shake128 instead of AES

#if -1 != ~0
#error We assume twos complement
#endif

#if q != 7879
#error Unsupported q
#endif

/* return -1 if x<0; otherwise return 0 */
static int int16_negative_mask(int16_t x) {
    uint16_t u = (uint16_t) x;
    u >>= 15;
    return -(int) u;
    /* alternative with gcc -fwrapv: */
    /* x>>15 compiles to CPU's arithmetic right shift */
}

// Reduces x modulo q=7879.
static int16_t freeze_u32(uint32_t x) {
    // This is Barrett reduction.
    //
    // Note that x mod q = x - floor(x/q)q.  We approximate 1/q by b/2^44
    // where b=2,232,794,269.  We have
    //
    //      |xb/2^44 - x/q| ≤ |x| |b/2^44 - 1/q| ≤ 2^-14.
    //
    // If x is not a multiple of q, then x/q is at least 1/q≈2^-13  apart from
    // an integer.  Hence floor(xb/2^44) = floor(x/q).
    //
    // For the remaining case, assume x is a multiple of q.  Note 1/q ≤ b/2^44.
    // Thus x/q ≤ xb/2^44 hence floor(x/q) = floor(xb/2^44) as well.
    //
    // Finally note that |xb| ≤ 2^63.05 and so floor(xb/2^44) = (xb)>>44
    // when computing with uint64s.
    return x - (uint32_t)(((uint64_t)x*2232794269) >> 44)*q;
}

// Reduce x modulo q=7879 provided |x| ≤ 390q = 3072810
static int16_t freeze_32(int32_t x) {
    // XXX we might also use two step reduction using only 32-bit mults:
    //      - first x -= q(x >> 13)
    //      - and then x -= q((17035 x)  >> 27).

    // Add 390q to turn positive.
    return freeze_u32((uint32_t)(x + 3072810 + q12)) - q12;
}

// Returns Top(x[i] + r[i] * (q-1)/2).
void top_add_r_q12(int8_t *t, const int16_t *x, const int8_t *r) {
    int i;
    for (i = 0; i < I; i++) {
        // XXX we could do with a lighter reduction here
        t[i] = (int8_t) ((tau1 *
                    (int32_t)(freeze_32(x[i] + r[i]*q12)
                        + tau0) + 32768) >> 16);
    }
}

// Returns sign bit of Right(t[i]) - x[i] + 4w + 1
void sign_right_sub_4w1(int8_t *r, const int8_t *t, const int16_t *x) {
    int i;
    for (i = 0; i < I; i++) {
        // XXX lighter reduction
        r[i] = (int8_t)-int16_negative_mask(
                freeze_32(
                   freeze_32(tau3 * (int32_t)t[i] - tau2)
                    - x[i] + 4 * w + 1));
    }
}

void short_random(int8_t *out) {
    uint32_t l[p];

    randombytes((uint8_t*)l, sizeof l);
    crypto_decode_pxint32(l, (uint8_t*)l);
    short_fromlist(out, l);
}

// See https://ntruprime.cr.yp.to/divergence-20180430.pdf
void short_fromlist(int8_t *out, const uint32_t* in) {
    uint32_t l[p];
    int i;

    // For the first w, we zero the last bit, so that the last two bits
    // are randomly 0 are 2, which will be interpreted as -1, 1.
    for (i = 0; i < w; i++)
        l[i] = in[i] & (uint32_t)-2;

    // For the remainder we'll set the last two bits to 1, which will
    // be interpreted as 0.
    for (i = w; i < p; i++)
        l[i] = (in[i] & (uint32_t) - 3) | 1; // set last bit and zero pen.ult. bit

    crypto_sort_uint32(l, p);
    for (i = 0; i < p; i++)
        out[i] = (int8_t)((l[i] & 3)-1);
}

void generator(int16_t *g, const uint8_t* seed) {
    uint32_t L[p];
    int i;

    shake128((uint8_t*)L, 4*p, seed, 16);
    crypto_decode_pxint32(L, (uint8_t*)L);
    for (i = 0; i < p; i++)
        g[i] = freeze_u32(L[i]) - q12;
}

// Only computes the first I coefficients of f*g
void mult_small_prefix(int16_t* r, int16_t* f, int8_t *g) {
    int32_t fg[I];  // lower coefficients of fg
    int32_t fg2[I]; // higher coefficients of fg
    int32_t result;
    int i, j;

    for (i = 0; i < I; ++i) {
        result = 0;
        for (j = 0; j <= i; ++j) {
            result += f[j] * (int32_t)g[i - j];
        }
        fg[i] = result;
    }

    for (i = p; i < p + I; ++i) {
        result = 0;
        for (j = i - p + 1; j < p; ++j) {
            result += f[j] * (int32_t)g[i - j];
        }
        fg2[i-p] = result;
    }

    // |fg[i]| ≤ w (q-1)/2.
    fg[I-1] += fg2[I-1];
    for (i = I-2; i >= 0; --i) {
        fg[i] += fg2[i];
        fg[i+1] += fg2[i];
    }

    // |fg[i]| ≤ 3 w (q-1)/2 ≈ 363q ≤ 390q.
    for (i = 0; i < I; i++) {
        r[i] = freeze_32(fg[i]);
    }

    // TODO check effect of L1 cache
}

void mult_small_ntt(int16_t* r, int16_t* f, int8_t *g) {
    int32_t f2[1536];
    int32_t g2[1536];
    int i;

    for (i = 0; i < p; i++) {
        f2[i] = f[i];
        g2[i] = g[i];
    }

    for (i = p; i < 1536; i++) {
        f2[i] = 0;
        g2[i] = 0;
    }

    ntt_1907713_1536(f2);
    ntt_1907713_1536(g2);
    ntt_1907713_1536_mul(f2, g2);
    ntt_1907713_1536_inverse(f2);

    r[0] = freeze_32(f2[0] + f2[p]);
    for (i = 1; i < p; i++) {
        r[i] = freeze_32(f2[i] + f2[i+p] +f2[i+p-1]);
    }
}

// Assumes f is frozen.
void mult_small(int16_t* f, int8_t *g) {
    int32_t fg[p + p - 1];
    int32_t result;
    int i, j;

    for (i = 0; i < p; ++i) {
        result = 0;
        for (j = 0; j <= i; ++j) {
            result += f[j] * (int32_t)g[i - j];
        }
        fg[i] = result;
    }

    for (i = p; i < p + p - 1; ++i) {
        result = 0;
        for (j = i - p + 1; j < p; ++j) {
            result += f[j] * (int32_t)g[i - j];
        }
        fg[i] = result;
    }

    // |fg[i]| ≤ w (q-1)/2.
    for (i = p + p - 2; i >= p; --i) {
        fg[i - p] += fg[i];
        fg[i - p + 1] += fg[i];
    }

    // |fg[i]| ≤ 3 w (q-1)/2 ≈ 363q ≤ 390q.
    for (i = 0; i < p; i++) {
        f[i] = freeze_32(fg[i]);
    }

    // TODO check effect of L1 cache
}

#if q >> 13 != 0
#error unsupported q
#endif

void crypto_round_and_encode(unsigned char *out, const int16_t *a) {
    int16_t x[p];
    int i;

    for (i = 0; i < p; ++i) {
        // Set x to nearest multiple of 3.  Works for |x|<2^14.
        x[i] = (int16_t) (3 * (((int32_t)a[i] * 10923 + 16384) >> 15));
    }

    crypto_encode_rounded(out, x);
}

#if p != 757 || q != 7879
#error Unsupported p and q
#endif

// CPU division instruction typically takes time depending on x.
// This software is designed to take time independent of x.
// Time still varies depending on m; user must ensure that m is constant.
// Time also varies on CPUs where multiplication is variable-time.
// There could be more CPU issues.
// There could also be compiler issues.
static void uint32_divmod_uint14(uint32_t *y, uint16_t *r, uint32_t x,
        uint16_t m)
{
  uint32_t v = 0x80000000;
  uint32_t qpart;
  uint32_t mask;

  v /= m;

  /* caller guarantees m > 0 */
  /* caller guarantees m < 16384 */
  /* vm <= 2^31 <= vm+m-1 */
  /* xvm <= 2^31 x <= xvm+x(m-1) */

  *y = 0;

  qpart = (x*(uint64_t)v)>>31;
  /* 2^31 qpart <= xv <= 2^31 qpart + 2^31-1 */
  /* 2^31 qpart m <= xvm <= 2^31 qpart m + (2^31-1)m */
  /* 2^31 qpart m <= 2^31 x <= 2^31 qpart m + (2^31-1)m + x(m-1) */
  /* 0 <= 2^31 newx <= (2^31-1)m + x(m-1) */
  /* 0 <= newx <= (1-1/2^31)m + x(m-1)/2^31 */
  /* 0 <= newx <= (1-1/2^31)(2^14-1) + (2^32-1)((2^14-1)-1)/2^31 */

  x -= qpart*m; *y += qpart;
  /* x <= 49146 */

  qpart = (x*(uint64_t)v)>>31;
  /* 0 <= newx <= (1-1/2^31)m + x(m-1)/2^31 */
  /* 0 <= newx <= m + 49146(2^14-1)/2^31 */
  /* 0 <= newx <= m + 0.4 */
  /* 0 <= newx <= m */

  x -= qpart*m; *y += qpart;
  /* x <= m */

  x -= m; *y += 1;
  mask = -(x>>31);
  x += mask&(uint32_t)m; *y += mask;
  /* x < m */

  *r = x;
}

static uint16_t uint32_mod_uint14(uint32_t x, uint16_t m)
{
  uint32_t quotient;
  uint16_t r;
  uint32_divmod_uint14(&quotient,&r,x,m);
  return r;
}

// Autogenerated by supercop's crypto_decode/653x1541/portable/decodegen.py
// with parameters "757 2627 2627 3939 True".  That's p, (q-1)/3+1, (q-1)/3+1,
// (q-1)/2 and "yes, we want to divide by 3".
void crypto_decode_rounded(int16_t *R0,const unsigned char *s)
{
  uint16_t R1[379],R2[190],R3[95],R4[48],R5[24],R6[12],R7[6],R8[3],R9[2],R10[1];
  long long i;
  uint16_t r0;
  uint32_t r1,r2;

  s += 1076;
  r1 = 0;
  r1 = (r1<<8)|*--s;
  r1 = (r1<<8)|*--s;
  r1 = uint32_mod_uint14(r1,1634); /* needed only for invalid inputs */
  R10[0] = r1;

  r2 = R10[0];
  r2 = (r2<<8)|*--s;
  uint32_divmod_uint14(&r1,&r0,r2,1408);
  R9[0] = r0;
  r1 = uint32_mod_uint14(r1,297); /* needed only for invalid inputs */
  R9[1] = r1;

  R8[2] = R9[1];
  r2 = R9[0];
  r2 = (r2<<8)|*--s;
  r2 = (r2<<8)|*--s;
  uint32_divmod_uint14(&r1,&r0,r2,9604);
  R8[0] = r0;
  r1 = uint32_mod_uint14(r1,9604); /* needed only for invalid inputs */
  R8[1] = r1;

  r2 = R8[2];
  r2 = (r2<<8)|*--s;
  uint32_divmod_uint14(&r1,&r0,r2,98);
  R7[4] = r0;
  r1 = uint32_mod_uint14(r1,775); /* needed only for invalid inputs */
  R7[5] = r1;
  for (i = 1;i >= 0;--i) {
    r2 = R8[i];
    uint32_divmod_uint14(&r1,&r0,r2,98);
    R7[2*i] = r0;
    r1 = uint32_mod_uint14(r1,98); /* needed only for invalid inputs */
    R7[2*i+1] = r1;
  }

  r2 = R7[5];
  r2 = (r2<<8)|*--s;
  uint32_divmod_uint14(&r1,&r0,r2,158);
  R6[10] = r0;
  r1 = uint32_mod_uint14(r1,1255); /* needed only for invalid inputs */
  R6[11] = r1;
  for (i = 4;i >= 0;--i) {
    r2 = R7[i];
    r2 = (r2<<8)|*--s;
    uint32_divmod_uint14(&r1,&r0,r2,158);
    R6[2*i] = r0;
    r1 = uint32_mod_uint14(r1,158); /* needed only for invalid inputs */
    R6[2*i+1] = r1;
  }

  r2 = R6[11];
  r2 = (r2<<8)|*--s;
  uint32_divmod_uint14(&r1,&r0,r2,3211);
  R5[22] = r0;
  r1 = uint32_mod_uint14(r1,100); /* needed only for invalid inputs */
  R5[23] = r1;
  for (i = 10;i >= 0;--i) {
    r2 = R6[i];
    r2 = (r2<<8)|*--s;
    r2 = (r2<<8)|*--s;
    uint32_divmod_uint14(&r1,&r0,r2,3211);
    R5[2*i] = r0;
    r1 = uint32_mod_uint14(r1,3211); /* needed only for invalid inputs */
    R5[2*i+1] = r1;
  }

  r2 = R5[23];
  r2 = (r2<<8)|*--s;
  r2 = (r2<<8)|*--s;
  uint32_divmod_uint14(&r1,&r0,r2,14506);
  R4[46] = r0;
  r1 = uint32_mod_uint14(r1,451); /* needed only for invalid inputs */
  R4[47] = r1;
  for (i = 22;i >= 0;--i) {
    r2 = R5[i];
    r2 = (r2<<8)|*--s;
    r2 = (r2<<8)|*--s;
    uint32_divmod_uint14(&r1,&r0,r2,14506);
    R4[2*i] = r0;
    r1 = uint32_mod_uint14(r1,14506); /* needed only for invalid inputs */
    R4[2*i+1] = r1;
  }

  R3[94] = R4[47];
  for (i = 46;i >= 0;--i) {
    r2 = R4[i];
    r2 = (r2<<8)|*--s;
    uint32_divmod_uint14(&r1,&r0,r2,1927);
    R3[2*i] = r0;
    r1 = uint32_mod_uint14(r1,1927); /* needed only for invalid inputs */
    R3[2*i+1] = r1;
  }

  r2 = R3[94];
  r2 = (r2<<8)|*--s;
  r2 = (r2<<8)|*--s;
  uint32_divmod_uint14(&r1,&r0,r2,11236);
  R2[188] = r0;
  r1 = uint32_mod_uint14(r1,2627); /* needed only for invalid inputs */
  R2[189] = r1;
  for (i = 93;i >= 0;--i) {
    r2 = R3[i];
    r2 = (r2<<8)|*--s;
    r2 = (r2<<8)|*--s;
    uint32_divmod_uint14(&r1,&r0,r2,11236);
    R2[2*i] = r0;
    r1 = uint32_mod_uint14(r1,11236); /* needed only for invalid inputs */
    R2[2*i+1] = r1;
  }

  R1[378] = R2[189];
  for (i = 188;i >= 0;--i) {
    r2 = R2[i];
    uint32_divmod_uint14(&r1,&r0,r2,106);
    R1[2*i] = r0;
    r1 = uint32_mod_uint14(r1,106); /* needed only for invalid inputs */
    R1[2*i+1] = r1;
  }

  R0[756] = 3*R1[378]-3939;
  for (i = 377;i >= 0;--i) {
    r2 = R1[i];
    r2 = (r2<<8)|*--s;
    r2 = (r2<<8)|*--s;
    uint32_divmod_uint14(&r1,&r0,r2,2627);
    R0[2*i] = 3*r0-3939;
    r1 = uint32_mod_uint14(r1,2627); /* needed only for invalid inputs */
    R0[2*i+1] = 3*r1-3939;
  }
}


// Autogenerated by supercop's crypto_encode/653x1541/portable/encodegen.py
// with parameters "757 2627 2627 3939 True".  That's p, (q-1)/3+1, (q-1)/3+1,
// (q-1)/2 and "yes, we want to divide by 3".
void crypto_encode_rounded(unsigned char *out,const int16_t *R0)
{
  // XXX: caller could overlap R with input
  uint16_t R[379];
  long i;
  uint16_t r0,r1;
  uint32_t r2;

  for (i = 0;i < 378;++i) {
    r0 = (((R0[2*i]+3939)&16383)*10923)>>15;
    r1 = (((R0[2*i+1]+3939)&16383)*10923)>>15;
    r2 = r0+r1*(uint32_t)2627;
    *out++ = r2; r2 >>= 8;
    *out++ = r2; r2 >>= 8;
    R[i] = r2;
  }
  R[378] = (((R0[756]+3939)&16383)*10923)>>15;

  for (i = 0;i < 189;++i) {
    r0 = R[2*i];
    r1 = R[2*i+1];
    r2 = r0+r1*(uint32_t)106;
    R[i] = r2;
  }
  R[189] = R[378];

  for (i = 0;i < 95;++i) {
    r0 = R[2*i];
    r1 = R[2*i+1];
    r2 = r0+r1*(uint32_t)11236;
    *out++ = r2; r2 >>= 8;
    *out++ = r2; r2 >>= 8;
    R[i] = r2;
  }

  for (i = 0;i < 47;++i) {
    r0 = R[2*i];
    r1 = R[2*i+1];
    r2 = r0+r1*(uint32_t)1927;
    *out++ = r2; r2 >>= 8;
    R[i] = r2;
  }
  R[47] = R[94];

  for (i = 0;i < 24;++i) {
    r0 = R[2*i];
    r1 = R[2*i+1];
    r2 = r0+r1*(uint32_t)14506;
    *out++ = r2; r2 >>= 8;
    *out++ = r2; r2 >>= 8;
    R[i] = r2;
  }

  for (i = 0;i < 11;++i) {
    r0 = R[2*i];
    r1 = R[2*i+1];
    r2 = r0+r1*(uint32_t)3211;
    *out++ = r2; r2 >>= 8;
    *out++ = r2; r2 >>= 8;
    R[i] = r2;
  }
  r0 = R[22];
  r1 = R[23];
  r2 = r0+r1*(uint32_t)3211;
  *out++ = r2; r2 >>= 8;
  R[11] = r2;

  for (i = 0;i < 6;++i) {
    r0 = R[2*i];
    r1 = R[2*i+1];
    r2 = r0+r1*(uint32_t)158;
    *out++ = r2; r2 >>= 8;
    R[i] = r2;
  }

  for (i = 0;i < 2;++i) {
    r0 = R[2*i];
    r1 = R[2*i+1];
    r2 = r0+r1*(uint32_t)98;
    R[i] = r2;
  }
  r0 = R[4];
  r1 = R[5];
  r2 = r0+r1*(uint32_t)98;
  *out++ = r2; r2 >>= 8;
  R[2] = r2;

  r0 = R[0];
  r1 = R[1];
  r2 = r0+r1*(uint32_t)9604;
  *out++ = r2; r2 >>= 8;
  *out++ = r2; r2 >>= 8;
  R[0] = r2;
  R[1] = R[2];

  r0 = R[0];
  r1 = R[1];
  r2 = r0+r1*(uint32_t)1408;
  *out++ = r2; r2 >>= 8;
  R[0] = r2;

  r0 = R[0];
  *out++ = r0; r0 >>= 8;
  *out++ = r0; r0 >>= 8;
}

#if p%4 != 1
#error unsupported p
#endif

void crypto_encode_px3(unsigned char *s, const int8_t *v) {
    uint8_t x;
    const uint8_t* f = (const uint8_t*)v;
    int i;

    for (i = 0; i < p / 4; ++i) {
        x = *f++ + 1;
        x += (*f++ + 1) << 2;
        x += (*f++ + 1) << 4;
        x += (*f++ + 1) << 6;
        *s++ = x;
    }
    x = *f++ + 1;
    *s++ = x;
}

void crypto_decode_px3(int8_t* v, const unsigned char *s) {
    uint8_t x;
    uint8_t* f = (uint8_t*)v;
    int i;

    for (i = 0; i < p / 4; ++i) {
        x = *s++;
        *f++ = (uint8_t) ((x & 3) - 1);
        x >>= 2;
        *f++ = (uint8_t) ((x & 3) - 1);
        x >>= 2;
        *f++ = (uint8_t) ((x & 3) - 1);
        x >>= 2;
        *f++ = (uint8_t) ((x & 3) - 1);
    }
    x = *s++;
    *f++ = (uint8_t) ((x & 3) - 1);
}

#if I%8 != 0
#error unsupported I
#endif

void crypto_encode_Ix2(unsigned char *s, const int8_t *x) {
    int i;
    for (i = 0; i < I / 8; ++i) {
        s[i] = 0;
    }
    for (i = 0; i < I; ++i) {
        s[i >> 3] |= (unsigned char) (((uint8_t)x[i] & 1) << (i & 7));
    }
}

void crypto_decode_Ix2(int8_t *x, const unsigned char *s) {
    int i;
    for (i = 0; i < I; ++i) {
        x[i] = 1 & (s[i >> 3] >> (i & 7));
    }
}

#if tau != 4
#error Unsupported tau
#endif

void crypto_encode_Ixtau(unsigned char *s, const int8_t *x) {
    int i;
    for (i = 0; i < I/4; ++i) {
        s[i] = x[4*i] | (x[4*i+1] << 2) | (x[4*i+2] << 4)
            | ((uint8_t)x[4*i+3] << 6);
    }
}

void crypto_decode_Ixtau(int8_t *x, const unsigned char *s) {
    int i;
    for (i = 0; i < I/4; ++i) {
        x[4*i] = s[i] & 3;
        x[4*i+1] = (s[i] >> 2) & 3;
        x[4*i+2] = (s[i] >> 4) & 3;
        x[4*i+3] = (s[i] >> 6);
    }
}

void crypto_encode_pxint16(unsigned char *s, const uint16_t *x) {
    int i;

    for (i = 0; i < p; ++i) {
        uint16_t u = *x++;
        *s++ = (unsigned char) u;
        *s++ = (unsigned char) (u >> 8);
    }
}

void crypto_decode_pxint16(uint16_t *x, const unsigned char *s) {
    int i;

    for (i = 0; i < p; ++i) {
        *x = (uint16_t)s[0] | ((uint16_t)s[1] << 8);
        x += 1;
        s += 2;
    }
}

void crypto_decode_pxint32(uint32_t *x, const unsigned char *s) {
    int i;

    for (i = 0; i < p; ++i) {
        *x = (uint32_t)s[0] | ((uint32_t)s[1] << 8)
            | ((uint32_t)s[2] << 16) | ((uint32_t)s[3] << 24);
        x += 1;
        s += 4;
    }
}

void crypto_sort_uint32(uint32_t *x, long long n) {
    long long j;
    for (j = 0; j < n; ++j) {
        x[j] ^= 0x80000000;
    }
    crypto_sort_int32((int32_t *)x, n);
    for (j = 0; j < n; ++j) {
        x[j] ^= 0x80000000;
    }
}

#define int32_MINMAX(a,b) \
    do { \
        int32_t ab = (b) ^ (a); \
        int32_t c = (int32_t)((int64_t)(b) - (int64_t)(a)); \
        c ^= ab & (c ^ (b)); \
        c >>= 31; \
        c &= ab; \
        (a) ^= c; \
        (b) ^= c; \
    } while(0)

/* assume 2 <= n <= 0x40000000 */
void crypto_sort_int32(int32_t *x, long long n) {
    int32_t top, k, m, r, i;
    long long j;

    top = 1;
    while (top < n - top) {
        top += top;
    }

    for (k = top; k >= 1; k >>= 1) {
        i = 0;
        while (i + 2 * k <= n) {
            for (j = i; j < i + k; ++j) {
                int32_MINMAX(x[j], x[j + k]);
            }
            i += 2 * k;
        }
        for (j = i; j < n - k; ++j) {
            int32_MINMAX(x[j], x[j + k]);
        }

        i = 0;
        j = 0;
        for (m = top; m > k; m >>= 1) {
            if (j != i) {
                for (;;) {
                    if (j == n - m) {
                        goto done;
                    }
                    int32_t a = x[j + k];
                    for (r = m; r > k; r >>= 1) {
                        int32_MINMAX(a, x[j + r]);
                    }
                    x[j + k] = a;
                    ++j;
                    if (j == i + k) {
                        i += 2 * k;
                        break;
                    }
                }
            }
            while (i + k <= n - m) {
                for (j = i; j < i + k; ++j) {
                    int32_t a = x[j + k];
                    for (r = m; r > k; r >>= 1) {
                        int32_MINMAX(a, x[j + r]);
                    }
                    x[j + k] = a;
                }
                i += 2 * k;
            }
            /* now i + k > n - m */
            j = i;
            while (j < n - m) {
                int32_t a = x[j + k];
                for (r = m; r > k; r >>= 1) {
                    int32_MINMAX(a, x[j + r]);
                }
                x[j + k] = a;
                ++j;
            }

done:
            ;
        }
    }
}
