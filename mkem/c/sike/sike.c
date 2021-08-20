/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: supersingular isogeny parameters and generation of functions for P434
*********************************************************************************************/  

#include <string.h>

#include "fips202.h"
#include "randombytes.h"

#include "sike.h" 
#include "internal.h"


/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: portable modular arithmetic for P434
*********************************************************************************************/


// Global constants
extern const uint64_t p434[NWORDS64_FIELD];
extern const uint64_t p434p1[NWORDS64_FIELD]; 
extern const uint64_t p434x2[NWORDS64_FIELD];  
extern const uint64_t p434x4[NWORDS64_FIELD];


void mp_sub434_p2(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision subtraction with correction with 2*p, c = a-b+2p. 
    unsigned int i, borrow = 0;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, a[i], b[i], borrow, c[i]); 
    }

    borrow = 0;
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(borrow, c[i], ((digit_t*)p434x2)[i], borrow, c[i]); 
    }
} 


void mp_sub434_p4(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision subtraction with correction with 4*p, c = a-b+4p. 
    unsigned int i, borrow = 0;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, a[i], b[i], borrow, c[i]); 
    }

    borrow = 0;
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(borrow, c[i], ((digit_t*)p434x4)[i], borrow, c[i]); 
    }
} 


void fpadd434(const digit_t* a, const digit_t* b, digit_t* c)
{ // Modular addition, c = a+b mod p434.
  // Inputs: a, b in [0, 2*p434-1] 
  // Output: c in [0, 2*p434-1] 
    unsigned int i, carry = 0;
    digit_t mask;

    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(carry, a[i], b[i], carry, c[i]); 
    }

    carry = 0;
    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(carry, c[i], ((digit_t*)p434x2)[i], carry, c[i]); 
    }
    mask = 0 - (digit_t)carry;

    carry = 0;
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(carry, c[i], ((digit_t*)p434x2)[i] & mask, carry, c[i]); 
    }
} 


void fpsub434(const digit_t* a, const digit_t* b, digit_t* c)
{ // Modular subtraction, c = a-b mod p434.
  // Inputs: a, b in [0, 2*p434-1] 
  // Output: c in [0, 2*p434-1] 
    unsigned int i, borrow = 0;
    digit_t mask;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, a[i], b[i], borrow, c[i]); 
    }
    mask = 0 - (digit_t)borrow;

    borrow = 0;
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(borrow, c[i], ((digit_t*)p434x2)[i] & mask, borrow, c[i]); 
    }
}


void fpneg434(digit_t* a)
{ // Modular negation, a = -a mod p434.
  // Input/output: a in [0, 2*p434-1] 
    unsigned int i, borrow = 0;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, ((digit_t*)p434x2)[i], a[i], borrow, a[i]); 
    }
}


void fpdiv2_434(const digit_t* a, digit_t* c)
{ // Modular division by two, c = a/2 mod p434.
  // Input : a in [0, 2*p434-1] 
  // Output: c in [0, 2*p434-1] 
    unsigned int i, carry = 0;
    digit_t mask;
        
    mask = 0 - (digit_t)(a[0] & 1);    // If a is odd compute a+p434
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(carry, a[i], ((digit_t*)p434)[i] & mask, carry, c[i]); 
    }

    mp_shiftr1(c, NWORDS_FIELD);
} 


void fpcorrection434(digit_t* a)
{ // Modular correction to reduce field element a in [0, 2*p434-1] to [0, p434-1].
    unsigned int i, borrow = 0;
    digit_t mask;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, a[i], ((digit_t*)p434)[i], borrow, a[i]); 
    }
    mask = 0 - (digit_t)borrow;

    borrow = 0;
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(borrow, a[i], ((digit_t*)p434)[i] & mask, borrow, a[i]); 
    }
}


void digit_x_digit(const digit_t a, const digit_t b, digit_t* c)
{ // Digit multiplication, digit * digit -> 2-digit result    
    register digit_t al, ah, bl, bh, temp;
    digit_t albl, albh, ahbl, ahbh, res1, res2, res3, carry;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*4), mask_high = (digit_t)(-1) << (sizeof(digit_t)*4);

    al = a & mask_low;                        // Low part
    ah = a >> (sizeof(digit_t) * 4);          // High part
    bl = b & mask_low;
    bh = b >> (sizeof(digit_t) * 4);

    albl = al*bl;
    albh = al*bh;
    ahbl = ah*bl;
    ahbh = ah*bh;
    c[0] = albl & mask_low;                   // C00

    res1 = albl >> (sizeof(digit_t) * 4);
    res2 = ahbl & mask_low;
    res3 = albh & mask_low;  
    temp = res1 + res2 + res3;
    carry = temp >> (sizeof(digit_t) * 4);
    c[0] ^= temp << (sizeof(digit_t) * 4);    // C01   

    res1 = ahbl >> (sizeof(digit_t) * 4);
    res2 = albh >> (sizeof(digit_t) * 4);
    res3 = ahbh & mask_low;
    temp = res1 + res2 + res3 + carry;
    c[1] = temp & mask_low;                   // C10 
    carry = temp & mask_high; 
    c[1] ^= (ahbh & mask_high) + carry;       // C11
}


void mp_mul(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{ // Multiprecision comba multiply, c = a*b, where lng(a) = lng(b) = nwords.   
    unsigned int i, j;
    digit_t t = 0, u = 0, v = 0, UV[2];
    unsigned int carry = 0;
    
    for (i = 0; i < nwords; i++) {
        for (j = 0; j <= i; j++) {
            MUL(a[j], b[i-j], UV+1, UV[0]); 
            ADDC(0, UV[0], v, carry, v); 
            ADDC(carry, UV[1], u, carry, u); 
            t += carry;
        }
        c[i] = v;
        v = u; 
        u = t;
        t = 0;
    }

    for (i = nwords; i < 2*nwords-1; i++) {
        for (j = i-nwords+1; j < nwords; j++) {
            MUL(a[j], b[i-j], UV+1, UV[0]); 
            ADDC(0, UV[0], v, carry, v); 
            ADDC(carry, UV[1], u, carry, u); 
            t += carry;
        }
        c[i] = v;
        v = u; 
        u = t;
        t = 0;
    }
    c[2*nwords-1] = v; 
}


void rdc_mont(digit_t* ma, digit_t* mc)
{ // Efficient Montgomery reduction using comba and exploiting the special form of the prime p434.
  // mc = ma*R^-1 mod p434x2, where R = 2^448.
  // If ma < 2^448*p434, the output mc is in the range [0, 2*p434-1].
  // ma is assumed to be in Montgomery representation.
    unsigned int i, j, carry, count = p434_ZERO_WORDS;
    digit_t UV[2], t = 0, u = 0, v = 0;

    for (i = 0; i < NWORDS_FIELD; i++) {
        mc[i] = 0;
    }

    for (i = 0; i < NWORDS_FIELD; i++) {
        for (j = 0; j < i; j++) {
            if (j < (i-p434_ZERO_WORDS+1)) { 
                MUL(mc[j], ((digit_t*)p434p1)[i-j], UV+1, UV[0]);
                ADDC(0, UV[0], v, carry, v); 
                ADDC(carry, UV[1], u, carry, u); 
                t += carry; 
            }
        }
        ADDC(0, v, ma[i], carry, v); 
        ADDC(carry, u, 0, carry, u); 
        t += carry; 
        mc[i] = v;
        v = u;
        u = t;
        t = 0;
    }    

    for (i = NWORDS_FIELD; i < 2*NWORDS_FIELD-1; i++) {
        if (count > 0) {
            count -= 1;
        }
        for (j = i-NWORDS_FIELD+1; j < NWORDS_FIELD; j++) {
            if (j < (NWORDS_FIELD-count)) { 
                MUL(mc[j], ((digit_t*)p434p1)[i-j], UV+1, UV[0]);
                ADDC(0, UV[0], v, carry, v); 
                ADDC(carry, UV[1], u, carry, u); 
                t += carry;
            }
        }
        ADDC(0, v, ma[i], carry, v); 
        ADDC(carry, u, 0, carry, u); 
        t += carry; 
        mc[i-NWORDS_FIELD] = v;
        v = u;
        u = t;
        t = 0;
    }
    ADDC(0, v, ma[2*NWORDS_FIELD-1], carry, v); 
    mc[NWORDS_FIELD-1] = v;
}


// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 
// Internally, the number of digits used to represent all these elements is obtained by approximating the number of bits to the immediately greater multiple of 32.
// For example, a 434-bit field element is represented with Ceil(434 / 64) = 7 64-bit digits or Ceil(434 / 32) = 14 32-bit digits.

//
// Curve isogeny system "SIDHp434". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p434^2), where A=6, B=1, C=1 and p434 = 2^216*3^137-1
//
         
const uint64_t p434[NWORDS64_FIELD]              = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFDC1767AE2FFFFFF, 
                                                     0x7BC65C783158AEA3, 0x6CFC5FD681C52056, 0x0002341F27177344 };
const uint64_t p434x2[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFB82ECF5C5FFFFFF,
                                                     0xF78CB8F062B15D47, 0xD9F8BFAD038A40AC, 0x0004683E4E2EE688 }; 
const uint64_t p434x4[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFC, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xF705D9EB8BFFFFFF, 
                                                     0xEF1971E0C562BA8F, 0xB3F17F5A07148159, 0x0008D07C9C5DCD11 }; 
const uint64_t p434p1[NWORDS64_FIELD]            = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xFDC1767AE3000000,
                                                     0x7BC65C783158AEA3, 0x6CFC5FD681C52056, 0x0002341F27177344 };  
const uint64_t p434x16p[2*NWORDS64_FIELD]        = { 0x0000000000000010, 0x0000000000000000, 0x0000000000000000, 0x47D130A3A0000000, 
                                                     0x873470F9D4EA2B80, 0x6074052FC75BF530, 0x54497C1B1D119772, 0xC55F373D2CDCA412, 
                                                     0x732CA2221C664B96, 0x6445AB96AF6359A5, 0x221708AB42ABE1B4, 0xAE3D3D0063244F01, 
                                                     0x18B920F2ECF68816, 0x0000004DB194809D }; 
// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000001000000 }; 
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0x58AEA3FDC1767AE3, 0xC520567BC65C7831, 0x1773446CFC5FD681, 0x0000000002341F27 };
// Alice's generator values {XPA0 + XPA1*i, XQA0 + xQA1*i, XRA0 + XRA1*i} in GF(p434^2), expressed in Montgomery representation
const uint64_t A_gen[6*NWORDS64_FIELD]           = { 0x05ADF455C5C345BF, 0x91935C5CC767AC2B, 0xAFE4E879951F0257, 0x70E792DC89FA27B1, 
                                                     0xF797F526BB48C8CD, 0x2181DB6131AF621F, 0x00000A1C08B1ECC4,    // XPA0
                                                     0x74840EB87CDA7788, 0x2971AA0ECF9F9D0B, 0xCB5732BDF41715D5, 0x8CD8E51F7AACFFAA, 
                                                     0xA7F424730D7E419F, 0xD671EB919A179E8C, 0x0000FFA26C5A924A,    // XPA1
                                                     0xFEC6E64588B7273B, 0xD2A626D74CBBF1C6, 0xF8F58F07A78098C7, 0xE23941F470841B03, 
                                                     0x1B63EDA2045538DD, 0x735CFEB0FFD49215, 0x0001C4CB77542876,    // XQA0
                                                     0xADB0F733C17FFDD6, 0x6AFFBD037DA0A050, 0x680EC43DB144E02F, 0x1E2E5D5FF524E374,
                                                     0xE2DDA115260E2995, 0xA6E4B552E2EDE508, 0x00018ECCDDF4B53E,    // XQA1
                                                     0x01BA4DB518CD6C7D, 0x2CB0251FE3CC0611, 0x259B0C6949A9121B, 0x60E17AC16D2F82AD, 
                                                     0x3AA41F1CE175D92D, 0x413FBE6A9B9BC4F3, 0x00022A81D8D55643,    // XRA0
                                                     0xB8ADBC70FC82E54A, 0xEF9CDDB0D5FADDED, 0x5820C734C80096A0, 0x7799994BAA96E0E4, 
                                                     0x044961599E379AF8, 0xDB2B94FBF09F27E2, 0x0000B87FC716C0C6 };  // XRA1
// Bob's generator values {XPB0, XQB0, XRB0 + XRB1*i} in GF(p434^2), expressed in Montgomery representation
const uint64_t B_gen[6*NWORDS64_FIELD]           = { 0x6E5497556EDD48A3, 0x2A61B501546F1C05, 0xEB919446D049887D, 0x5864A4A69D450C4F, 
                                                     0xB883F276A6490D2B, 0x22CC287022D5F5B9, 0x0001BED4772E551F,    // XPB0 
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000,    // XPB1
                                                     0xFAE2A3F93D8B6B8E, 0x494871F51700FE1C, 0xEF1A94228413C27C, 0x498FF4A4AF60BD62, 
                                                     0xB00AD2A708267E8A, 0xF4328294E017837F, 0x000034080181D8AE,    // XQB0
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000,    // XQB1
                                                     0x283B34FAFEFDC8E4, 0x9208F44977C3E647, 0x7DEAE962816F4E9A, 0x68A2BA8AA262EC9D, 
                                                     0x8176F112EA43F45B, 0x02106D022634F504, 0x00007E8A50F02E37,    // XRB0
                                                     0xB378B7C1DA22CCB1, 0x6D089C99AD1D9230, 0xEBE15711813E2369, 0x2B35A68239D48A53, 
                                                     0x445F6FD138407C93, 0xBEF93B29A3F6B54B, 0x000173FA910377D3 };  // XRB1
// Montgomery constant Montgomery_R2 = (2^448)^2 mod p434
const uint64_t Montgomery_R2[NWORDS64_FIELD]     = { 0x28E55B65DCD69B30, 0xACEC7367768798C2, 0xAB27973F8311688D, 0x175CC6AF8D6C7C0B,
                                                     0xABCD92BF2DDE347E, 0x69E16A61C7686D9A, 0x000025A89BCDD12A };                                                   
// Value one in Montgomery representation 
const uint64_t Montgomery_one[NWORDS64_FIELD]    = { 0x000000000000742C, 0x0000000000000000, 0x0000000000000000, 0xB90FF404FC000000, 
                                                     0xD801A4FB559FACD4, 0xE93254545F77410C, 0x0000ECEEA7BD2EDA };


// Fixed parameters for isogeny tree computation
const unsigned int strat_Alice[MAX_Alice-1] = { 
48, 28, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 13, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 
1, 1, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 
1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1 };

const unsigned int strat_Bob[MAX_Bob-1] = { 
66, 33, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 
2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 3, 1, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 
1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1 };
           
// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions
#define fpcopy                        fpcopy434
#define fpzero                        fpzero434
#define fpadd                         fpadd434
#define fpsub                         fpsub434
#define fpneg                         fpneg434
#define fpdiv2                        fpdiv2_434
#define fpcorrection                  fpcorrection434
#define fpmul_mont                    fpmul434_mont
#define fpsqr_mont                    fpsqr434_mont
#define fpinv_mont                    fpinv434_mont
#define fpinv_chain_mont              fpinv434_chain_mont
#define fpinv_mont_bingcd             fpinv434_mont_bingcd
#define fp2copy                       fp2copy434
#define fp2zero                       fp2zero434
#define fp2add                        fp2add434
#define fp2sub                        fp2sub434
#define mp_sub_p2                     mp_sub434_p2
#define mp_sub_p4                     mp_sub434_p4
#define sub_p4                        mp_sub_p4
#define fp2neg                        fp2neg434
#define fp2div2                       fp2div2_434
#define fp2correction                 fp2correction434
#define fp2mul_mont                   fp2mul434_mont
#define fp2sqr_mont                   fp2sqr434_mont
#define fp2inv_mont                   fp2inv434_mont
#define fp2inv_mont_bingcd            fp2inv434_mont_bingcd
#define fpequal_non_constant_time     fpequal434_non_constant_time
#define mp_add_asm                    mp_add434_asm
#define mp_subaddx2_asm               mp_subadd434x2_asm
#define mp_dblsubx2_asm               mp_dblsub434x2_asm


/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: core functions over GF(p) and GF(p^2)
*********************************************************************************************/


void clear_words(void* mem, digit_t nwords)
{ // Clear digits from memory. "nwords" indicates the number of digits to be zeroed.
  // This function uses the volatile type qualifier to inform the compiler not to optimize out the memory clearing.
    volatile digit_t *v = mem; 

    for (unsigned int i = 0; i < nwords; i++)
        v[i] = 0;
}


int8_t ct_compare(const uint8_t *a, const uint8_t *b, unsigned int len) 
{ // Compare two byte arrays in constant time.
  // Returns 0 if the byte arrays are equal, -1 otherwise.
    uint8_t r = 0;

    for (unsigned int i = 0; i < len; i++)
        r |= a[i] ^ b[i];

    return (-(int8_t)r) >> (8*sizeof(uint8_t)-1);
}


void ct_cmov(uint8_t *r, const uint8_t *a, unsigned int len, int8_t selector) 
{ // Conditional move in constant time.
  // If selector = -1 then load r with a, else if selector = 0 then keep r.

    for (unsigned int i = 0; i < len; i++)
        r[i] ^= selector & (a[i] ^ r[i]);
}


static void encode_to_bytes(const digit_t* x, unsigned char* enc, int nbytes)
{ // Encoding digits to bytes according to endianness
#ifdef _BIG_ENDIAN_
    int ndigits = nbytes / sizeof(digit_t);
    int rem = nbytes % sizeof(digit_t);

    for (int i = 0; i < ndigits; i++)
        ((digit_t*)enc)[i] = BSWAP_DIGIT(x[i]);
    if (rem) {
        digit_t ld = BSWAP_DIGIT(x[ndigits]);
        memcpy(enc + ndigits*sizeof(digit_t), (unsigned char*)&ld, rem);
    }
#else    
    memcpy(enc, (const unsigned char*)x, nbytes);
#endif
}


static void decode_to_digits(const unsigned char* x, digit_t* dec, int nbytes, int ndigits)
{ // Decoding bytes to digits according to endianness

    dec[ndigits - 1] = 0;
    memcpy((unsigned char*)dec, x, nbytes);
#ifdef _BIG_ENDIAN_
    for (int i = 0; i < ndigits; i++)
        dec[i] = BSWAP_DIGIT(dec[i]);
#endif
}


static void fp2_encode(const f2elm_t x, unsigned char *enc)
{ // Conversion of GF(p^2) element from Montgomery to standard representation, and encoding by removing leading 0 bytes
    f2elm_t t;

    from_fp2mont(x, t);
    encode_to_bytes(t[0], enc, FP2_ENCODED_BYTES / 2);
    encode_to_bytes(t[1], enc + FP2_ENCODED_BYTES / 2, FP2_ENCODED_BYTES / 2);
}


static void fp2_decode(const unsigned char *x, f2elm_t dec)
{ // Parse byte sequence back into GF(p^2) element, and conversion to Montgomery representation

    decode_to_digits(x, dec[0], FP2_ENCODED_BYTES / 2, NWORDS_FIELD);
    decode_to_digits(x + FP2_ENCODED_BYTES / 2, dec[1], FP2_ENCODED_BYTES / 2, NWORDS_FIELD);
    to_fp2mont(dec, dec);
}


void fpcopy(const felm_t a, felm_t c)
{ // Copy a field element, c = a.
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++)
        c[i] = a[i];
}


void fpzero(felm_t a)
{ // Zero a field element, a = 0.
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++)
        a[i] = 0;
}


void to_mont(const felm_t a, felm_t mc)
{ // Conversion to Montgomery representation,
  // mc = a*R^2*R^(-1) mod p = a*R mod p, where a in [0, p-1].
  // The Montgomery constant R^2 mod p is the global value "Montgomery_R2". 

    fpmul_mont(a, (digit_t*)&Montgomery_R2, mc);
}


void from_mont(const felm_t ma, felm_t c)
{ // Conversion from Montgomery representation to standard representation,
  // c = ma*R^(-1) mod p = a mod p, where ma in [0, p-1].
    digit_t one[NWORDS_FIELD] = {0};
    
    one[0] = 1;
    fpmul_mont(ma, one, c);
    fpcorrection(c);
}


void copy_words(const digit_t* a, digit_t* c, const unsigned int nwords)
{ // Copy wordsize digits, c = a, where lng(a) = nwords.
    unsigned int i;
        
    for (i = 0; i < nwords; i++)                      
        c[i] = a[i];
}

void fpmul_mont(const felm_t ma, const felm_t mb, felm_t mc)
{ // Multiprecision multiplication, c = a*b mod p.
    dfelm_t temp = {0};

    mp_mul(ma, mb, temp, NWORDS_FIELD);
    rdc_mont(temp, mc);
}


void fpsqr_mont(const felm_t ma, felm_t mc)
{ // Multiprecision squaring, c = a^2 mod p.
    dfelm_t temp = {0};

    mp_mul(ma, ma, temp, NWORDS_FIELD);
    rdc_mont(temp, mc);
}


void fpinv_mont(felm_t a)
{ // Field inversion using Montgomery arithmetic, a = a^(-1)*R mod p.
    felm_t tt;

    fpcopy(a, tt);
    fpinv_chain_mont(tt);
    fpsqr_mont(tt, tt);
    fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, a);
}


void fp2copy(const f2elm_t a, f2elm_t c)
{ // Copy a GF(p^2) element, c = a.
    fpcopy(a[0], c[0]);
    fpcopy(a[1], c[1]);
}


void fp2zero(f2elm_t a)
{ // Zero a GF(p^2) element, a = 0.
    fpzero(a[0]);
    fpzero(a[1]);
}


void fp2neg(f2elm_t a)
{ // GF(p^2) negation, a = -a in GF(p^2).
    fpneg(a[0]);
    fpneg(a[1]);
}


void fp2add(const f2elm_t a, const f2elm_t b, f2elm_t c)           
{ // GF(p^2) addition, c = a+b in GF(p^2).
    fpadd(a[0], b[0], c[0]);
    fpadd(a[1], b[1], c[1]);
}


void fp2sub(const f2elm_t a, const f2elm_t b, f2elm_t c)          
{ // GF(p^2) subtraction, c = a-b in GF(p^2).
    fpsub(a[0], b[0], c[0]);
    fpsub(a[1], b[1], c[1]);
}


void fp2div2(const f2elm_t a, f2elm_t c)          
{ // GF(p^2) division by two, c = a/2  in GF(p^2).
    fpdiv2(a[0], c[0]);
    fpdiv2(a[1], c[1]);
}


void fp2correction(f2elm_t a)
{ // Modular correction, a = a in GF(p^2).
    fpcorrection(a[0]);
    fpcorrection(a[1]);
}


static void mp_addfast(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision addition, c = a+b.    
#if (OS_TARGET == OS_WIN) || defined(GENERIC_IMPLEMENTATION) || (TARGET == TARGET_ARM)

    mp_add(a, b, c, NWORDS_FIELD);
    
#elif (OS_TARGET == OS_NIX)                 
    
    mp_add_asm(a, b, c);    

#endif
}


static void mp2_add(const f2elm_t a, const f2elm_t b, f2elm_t c)       
{ // GF(p^2) addition without correction, c = a+b in GF(p^2). 
    mp_addfast(a[0], b[0], c[0]);
    mp_addfast(a[1], b[1], c[1]);
}


static void mp2_sub_p2(const f2elm_t a, const f2elm_t b, f2elm_t c)       
{ // GF(p^2) subtraction with correction with 2*p, c = a-b+2p in GF(p^2).    
    mp_sub_p2(a[0], b[0], c[0]);  
    mp_sub_p2(a[1], b[1], c[1]);
}


unsigned int mp_add(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{ // Multiprecision addition, c = a+b, where lng(a) = lng(b) = nwords. Returns the carry bit.
    unsigned int i, carry = 0;
        
    for (i = 0; i < nwords; i++) {                      
        ADDC(carry, a[i], b[i], carry, c[i]);
    }

    return carry;
}


void fp2sqr_mont(const f2elm_t a, f2elm_t c)
{ // GF(p^2) squaring using Montgomery arithmetic, c = a^2 in GF(p^2).
  // Inputs: a = a0+a1*i, where a0, a1 are in [0, 2*p-1] 
  // Output: c = c0+c1*i, where c0, c1 are in [0, 2*p-1] 
    felm_t t1, t2, t3;
    
    mp_addfast(a[0], a[1], t1);                      // t1 = a0+a1 
    sub_p4(a[0], a[1], t2);                          // t2 = a0-a1
    mp_addfast(a[0], a[0], t3);                      // t3 = 2a0
    fpmul_mont(t1, t2, c[0]);                        // c0 = (a0+a1)(a0-a1)
    fpmul_mont(t3, a[1], c[1]);                      // c1 = 2a0*a1
}


unsigned int mp_sub(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{ // Multiprecision subtraction, c = a-b, where lng(a) = lng(b) = nwords. Returns the borrow bit.
    unsigned int i, borrow = 0;

    for (i = 0; i < nwords; i++)
        SUBC(borrow, a[i], b[i], borrow, c[i]);

    return borrow;
}


static void mp_subaddfast(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision subtraction followed by addition with p*2^MAXBITS_FIELD, c = a-b+(p*2^MAXBITS_FIELD) if a-b < 0, otherwise c=a-b. 
#if (OS_TARGET == OS_WIN) || defined(GENERIC_IMPLEMENTATION) || (TARGET == TARGET_ARM)
    felm_t t1;

    digit_t mask = 0 - (digit_t)mp_sub(a, b, c, 2*NWORDS_FIELD);
    for (int i = 0; i < NWORDS_FIELD; i++)
        t1[i] = ((digit_t*)PRIME)[i] & mask;
    mp_addfast((digit_t*)&c[NWORDS_FIELD], t1, (digit_t*)&c[NWORDS_FIELD]);

#elif (OS_TARGET == OS_NIX)               

    mp_subaddx2_asm(a, b, c);     

#endif
}


static void mp_dblsubfast(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision subtraction, c = c-a-b, where lng(a) = lng(b) = 2*NWORDS_FIELD.
#if (OS_TARGET == OS_WIN) || defined(GENERIC_IMPLEMENTATION) || (TARGET == TARGET_ARM)

    mp_sub(c, a, c, 2*NWORDS_FIELD);
    mp_sub(c, b, c, 2*NWORDS_FIELD);

#elif (OS_TARGET == OS_NIX)                 

    mp_dblsubx2_asm(a, b, c);

#endif
}


void fp2mul_mont(const f2elm_t a, const f2elm_t b, f2elm_t c)
{ // GF(p^2) multiplication using Montgomery arithmetic, c = a*b in GF(p^2).
  // Inputs: a = a0+a1*i and b = b0+b1*i, where a0, a1, b0, b1 are in [0, 2*p-1] 
  // Output: c = c0+c1*i, where c0, c1 are in [0, 2*p-1] 
    felm_t t1, t2;
    dfelm_t tt1, tt2, tt3; 
    
    mp_addfast(a[0], a[1], t1);                      // t1 = a0+a1
    mp_addfast(b[0], b[1], t2);                      // t2 = b0+b1
    mp_mul(a[0], b[0], tt1, NWORDS_FIELD);           // tt1 = a0*b0
    mp_mul(a[1], b[1], tt2, NWORDS_FIELD);           // tt2 = a1*b1
    mp_mul(t1, t2, tt3, NWORDS_FIELD);               // tt3 = (a0+a1)*(b0+b1)
    mp_dblsubfast(tt1, tt2, tt3);                    // tt3 = (a0+a1)*(b0+b1) - a0*b0 - a1*b1
    mp_subaddfast(tt1, tt2, tt1);                    // tt1 = a0*b0 - a1*b1 + p*2^MAXBITS_FIELD if a0*b0 - a1*b1 < 0, else tt1 = a0*b0 - a1*b1
    rdc_mont(tt3, c[1]);                             // c[1] = (a0+a1)*(b0+b1) - a0*b0 - a1*b1 
    rdc_mont(tt1, c[0]);                             // c[0] = a0*b0 - a1*b1
}


void fpinv_chain_mont(felm_t a)
{ // Chain to compute a^(p-3)/4 using Montgomery arithmetic.
    unsigned int i, j;
    
#if (NBITS_FIELD == 434)
    felm_t t[31], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(a, tt, t[0]);
    for (i = 0; i <= 29; i++) fpmul_mont(t[i], tt, t[i+1]);

    fpcopy(a, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[30], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[30], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[21], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[26], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[25], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[30], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[26], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[28], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (j = 0; j < 35; j++) {
        for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[30], tt, tt);
    }
    fpcopy(tt, a);   
    
#elif (NBITS_FIELD == 503)
    felm_t t[15], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(a, tt, t[0]);
    for (i = 0; i <= 13; i++) fpmul_mont(t[i], tt, t[i+1]);

    fpcopy(a, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 12; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (j = 0; j < 49; j++) {
        for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[14], tt, tt);
    }
    fpcopy(tt, a);

#elif (NBITS_FIELD == 610)
    felm_t t[31], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(a, tt, t[0]);
    for (i = 0; i <= 29; i++) fpmul_mont(t[i], tt, t[i+1]);

    fpcopy(a, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[30], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[25], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[28], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i = 0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[28], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[15], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[15], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[27], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[28], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[29], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[30], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[25], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[28], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[26], tt, tt);
    for (i = 0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (j = 0; j < 50; j++) {
        for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[30], tt, tt);
    }
    fpcopy(tt, a);    

#elif (NBITS_FIELD == 751)
    felm_t t[27], tt;
    
    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(a, tt, t[0]);
    fpmul_mont(t[0], tt, t[1]);
    fpmul_mont(t[1], tt, t[2]);
    fpmul_mont(t[2], tt, t[3]); 
    fpmul_mont(t[3], tt, t[3]);
    for (i = 3; i <= 8; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[9], tt, t[9]);
    for (i = 9; i <= 20; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[21], tt, t[21]); 
    for (i = 21; i <= 24; i++) fpmul_mont(t[i], tt, t[i+1]); 
    fpmul_mont(t[25], tt, t[25]);
    fpmul_mont(t[25], tt, t[26]);

    fpcopy(a, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[15], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[26], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[18], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[18], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[17], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[25], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[18], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[21], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[17], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[26], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (j = 0; j < 61; j++) {
        for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[26], tt, tt);
    }
    fpcopy(tt, a);  
#endif
}


void fp2inv_mont(f2elm_t a)
{// GF(p^2) inversion using Montgomery arithmetic, a = (a0-i*a1)/(a0^2+a1^2).
    f2elm_t t1;

    fpsqr_mont(a[0], t1[0]);                         // t10 = a0^2
    fpsqr_mont(a[1], t1[1]);                         // t11 = a1^2
    fpadd(t1[0], t1[1], t1[0]);                      // t10 = a0^2+a1^2
    fpinv_mont(t1[0]);                               // t10 = (a0^2+a1^2)^-1
    fpneg(a[1]);                                     // a = a0-i*a1
    fpmul_mont(a[0], t1[0], a[0]);
    fpmul_mont(a[1], t1[0], a[1]);                   // a = (a0-i*a1)*(a0^2+a1^2)^-1
}


void to_fp2mont(const f2elm_t a, f2elm_t mc)
{ // Conversion of a GF(p^2) element to Montgomery representation,
  // mc_i = a_i*R^2*R^(-1) = a_i*R in GF(p^2). 

    to_mont(a[0], mc[0]);
    to_mont(a[1], mc[1]);
}


void from_fp2mont(const f2elm_t ma, f2elm_t c)
{ // Conversion of a GF(p^2) element from Montgomery representation to standard representation,
  // c_i = ma_i*R^(-1) = a_i in GF(p^2).

    from_mont(ma[0], c[0]);
    from_mont(ma[1], c[1]);
}


void mp_shiftleft(digit_t* x, unsigned int shift, const unsigned int nwords)
{
    unsigned int i, j = 0;

    while (shift > RADIX) {
        j += 1;
        shift -= RADIX;
    }

    for (i = 0; i < nwords-j; i++) 
        x[nwords-1-i] = x[nwords-1-i-j];
    for (i = nwords-j; i < nwords; i++) 
        x[nwords-1-i] = 0;
    if (shift != 0) {
        for (j = nwords-1; j > 0; j--) 
            SHIFTL(x[j], x[j-1], shift, x[j], RADIX);
        x[0] <<= shift;
    }
}


void mp_shiftr1(digit_t* x, const unsigned int nwords)
{ // Multiprecision right shift by one.
    unsigned int i;

    for (i = 0; i < nwords-1; i++) {
        SHIFTR(x[i+1], x[i], 1, x[i], RADIX);
    }
    x[nwords-1] >>= 1;
}


void mp_shiftl1(digit_t* x, const unsigned int nwords)
{ // Multiprecision left shift by one.
    int i;

    for (i = nwords-1; i > 0; i--) {
        SHIFTL(x[i], x[i-1], 1, x[i], RADIX);
    }
    x[0] <<= 1;
}

#ifdef COMPRESS

static unsigned int is_felm_zero(const felm_t x)
{ // Is x = 0? return 1 (TRUE) if condition is true, 0 (FALSE) otherwise.
  // SECURITY NOTE: This function does not run in constant-time.
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++) {
        if (x[i] != 0) return 0;
    }
    return 1;
}

static unsigned int is_felm_one(const felm_t x)
{ // Is x = 0? return 1 (TRUE) if condition is true, 0 (FALSE) otherwise.
  // SECURITY NOTE: This function does not run in constant-time.
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++) {
        if (x[i] != 0) return 0;
    }
    return 1;
}

void mul3(unsigned char *a) 
{ // Computes a = 3*a
  // The input is assumed to be OBOB_BITS-2 bits long and stored in SECRETKEY_B_BYTES
    digit_t temp1[NWORDS_ORDER] = {0}, temp2[NWORDS_ORDER] = {0};
        
    decode_to_digits(a, temp1, SECRETKEY_B_BYTES, NWORDS_ORDER);
    mp_add(temp1, temp1, temp2, NWORDS_ORDER);               // temp2 = 2*a
    mp_add(temp1, temp2, temp1, NWORDS_ORDER);               // temp1 = 3*a
    encode_to_bytes(temp1, a, SECRETKEY_B_BYTES);
    
    clear_words((void*)temp1, NWORDS_ORDER);
    clear_words((void*)temp2, NWORDS_ORDER);
}


unsigned int mod3(digit_t* a) 
{ // Computes the input modulo 3
  // The input is assumed to be NWORDS_ORDER long 
    digit_t temp;
    hdigit_t *val = (hdigit_t*)a, r = 0;

    for (int i = (2*NWORDS_ORDER-1); i >= 0; i--) {
        temp = ((digit_t)r << (sizeof(hdigit_t)*8)) | (digit_t)val[i];
        r = temp % 3;
    }

    return r;
}


void fp2shl(const f2elm_t a, const int k, f2elm_t c) 
{  // c = (2^k)*a
   fp2copy(a, c);
   for (int j = 0; j < k; j++) {
      fp2add(c, c, c);
   }
}


void fp2_conj(const f2elm_t v, f2elm_t r)
{ // r = a - b*i where v = a + b*i
    fpcopy(v[0],r[0]);
    fpcopy(v[1],r[1]);
    
    if(!is_felm_zero(r[1])) {
        fpneg(r[1]);
    }
}


void sqr_Fp2_cycl(f2elm_t a, const felm_t one)
{ // Cyclotomic squaring on elements of norm 1, using a^(p+1) = 1.
     felm_t t0;
 
     fpadd(a[0], a[1], t0);              // t0 = a0 + a1
     fpsqr_mont(t0, t0);                 // t0 = t0^2
     fpsub(t0, one, a[1]);               // a1 = t0 - 1   
     fpsqr_mont(a[0], t0);               // t0 = a0^2
     fpadd(t0, t0, t0);                  // t0 = t0 + t0
     fpsub(t0, one, a[0]);               // a0 = t0 - 1
}


void cube_Fp2_cycl(f2elm_t a, const felm_t one)
{ // Cyclotomic cubing on elements of norm 1, using a^(p+1) = 1.
     felm_t t0;
   
     fpadd(a[0], a[0], t0);              // t0 = a0 + a0
     fpsqr_mont(t0, t0);                 // t0 = t0^2
     fpsub(t0, one, t0);                 // t0 = t0 - 1
     fpmul_mont(a[1], t0, a[1]);         // a1 = t0*a1
     fpsub(t0, one, t0);
     fpsub(t0, one, t0);                 // t0 = t0 - 2
     fpmul_mont(a[0], t0, a[0]);         // a0 = t0*a0
}






static bool is_zero(digit_t* a, unsigned int nwords)
{ // Check if multiprecision element is zero.
  // SECURITY NOTE: This function does not run in constant time.

    for (unsigned int i = 0; i < nwords; i++) {
        if (a[i] != 0) {
            return false;
        } 
    }

    return true;
}


unsigned char is_sqr_fp2(const f2elm_t a, felm_t s) 
{ // Test if a is a square in GF(p^2) and return 1 if true, 0 otherwise
  // If a is a quadratic residue, s will be assigned with a partially computed square root of a
    int i;
    felm_t a0,a1,z,temp;
    
    fpsqr_mont(a[0],a0);
    fpsqr_mont(a[1],a1);
    fpadd(a0,a1,z);
    
    fpcopy(z,s);
    for (i = 0; i < OALICE_BITS - 2; i++) {             
        fpsqr_mont(s, s);
    }
    for (i = 0; i < OBOB_EXPON; i++) {
        fpsqr_mont(s, temp);
        fpmul_mont(s, temp, s);
    }  
    fpsqr_mont(s,temp);          // s = z^((p+1)/4)
    fpcorrection(temp);
    fpcorrection(z);
    if (memcmp((unsigned char*)temp, (unsigned char*)z, NBITS_TO_NBYTES(NBITS_FIELD)) != 0)  // s^2 !=? z
        return 0;
    
    return 1;
}


void sqrt_Fp2(const f2elm_t u, f2elm_t y)
{ // Computes square roots of elements in (Fp2)^2 using Hamburg's trick. 
    felm_t t0, t1, t2, t3;
    digit_t *a  = (digit_t*)u[0], *b  = (digit_t*)u[1];
    unsigned int i;

    fpsqr_mont(a, t0);                   // t0 = a^2
    fpsqr_mont(b, t1);                   // t1 = b^2
    fpadd(t0, t1, t0);                   // t0 = t0+t1 
    fpcopy(t0, t1);
    for (i = 0; i < OALICE_BITS - 2; i++) {   // t = t3^((p+1)/4)
        fpsqr_mont(t1, t1);
    }
    for (i = 0; i < OBOB_EXPON; i++) {
        fpsqr_mont(t1, t0);
        fpmul_mont(t1, t0, t1);
    }  
    fpadd(a, t1, t0);                    // t0 = a+t1      
    fpdiv2(t0, t0);                      // t0 = t0/2 
    fpcopy(t0, t2);
    fpinv_chain_mont(t2);                // t2 = t0^((p-3)/4)      
    fpmul_mont(t0, t2, t1);              // t1 = t2*t0             
    fpmul_mont(t2, b, t2);               // t2 = t2*b       
    fpdiv2(t2, t2);                      // t2 = t2/2 
    fpsqr_mont(t1, t3);                  // t3 = t1^2              
    fpcorrection(t0);
    fpcorrection(t3);
           
    if (memcmp(t0, t3, NBITS_TO_NBYTES(NBITS_FIELD)) == 0) {
        fpcopy(t1, y[0]);
        fpcopy(t2, y[1]);
    } else {
        fpneg(t1);
        fpcopy(t2, y[0]);
        fpcopy(t1, y[1]);
    }
}


static void power2_setup(digit_t* x, int mark, const unsigned int nwords)
{ // Set up the value 2^mark.
    unsigned int i;

    for (i = 0; i < nwords; i++) x[i] = 0;

    i = 0;
    while (mark >= 0) {
        if (mark < RADIX) {
            x[i] = (digit_t)1 << mark;
        }
        mark -= RADIX;
        i += 1;
    }    
}


int8_t cmp_f2elm(const f2elm_t x, const f2elm_t y)
{ // Comparison of two GF(p^2) elements in constant time. 
  // Is x != y? return -1 if condition is true, 0 otherwise.
    f2elm_t a, b;      
    uint8_t r = 0;
    
    fp2copy(x, a);
    fp2copy(y, b);
    fp2correction(a);
    fp2correction(b);
    
    for (int i = NWORDS_FIELD-1; i >= 0; i--)
        r |= (a[0][i] ^ b[0][i]) | (a[1][i] ^ b[1][i]);

    return (-(int8_t)r) >> (8*sizeof(uint8_t)-1);
}


static unsigned int is_felm_even(const felm_t x)
{ // Is x even? return 1 (TRUE) if condition is true, 0 (FALSE) otherwise.
    return (unsigned int)((x[0] & 1) ^ 1);
}


static unsigned int is_felm_lt(const felm_t x, const felm_t y)
{ // Is x < y? return 1 (TRUE) if condition is true, 0 (FALSE) otherwise.
  // SECURITY NOTE: This function does not run in constant-time.

    for (int i = NWORDS_FIELD-1; i >= 0; i--) {
        if (x[i] < y[i]) { 
            return true;
        } else if (x[i] > y[i]) {
            return false;
        }
    }
    return false;
}


static unsigned int is_orderelm_lt(const digit_t *x, const digit_t *y)
{ // Is x < y? return 1 (TRUE) if condition is true, 0 (FALSE) otherwise.
  // SECURITY NOTE: This function does not run in constant-time.

    for (int i = NWORDS_ORDER-1; i >= 0; i--) {
        if (x[i] < y[i]) { 
            return true;
        } else if (x[i] > y[i]) {
            return false;
        }
    }
    return false;
}


static void fpinv_mont_bingcd_partial(const felm_t a, felm_t x1, unsigned int* k)
{ // Partial Montgomery inversion via the binary GCD algorithm.
    felm_t u, v, x2;
    unsigned int cwords;  // Number of words necessary for x1, x2

    fpcopy(a, u);
    fpcopy((digit_t*)PRIME, v);
    fpzero(x1); x1[0] = 1;
    fpzero(x2);
    *k = 0;

    while (!is_felm_zero(v)) {
        cwords = ((*k + 1) / RADIX) + 1;
        if ((cwords < NWORDS_FIELD)) {
            if (is_felm_even(v)) {
                mp_shiftr1(v, NWORDS_FIELD);
                mp_shiftl1(x1, cwords);
            } else if (is_felm_even(u)) {
                mp_shiftr1(u, NWORDS_FIELD);
                mp_shiftl1(x2, cwords);
            } else if (!is_felm_lt(v, u)) {
                mp_sub(v, u, v, NWORDS_FIELD);
                mp_shiftr1(v, NWORDS_FIELD);
                mp_add(x1, x2, x2, cwords);
                mp_shiftl1(x1, cwords);
            } else {
                mp_sub(u, v, u, NWORDS_FIELD);
                mp_shiftr1(u, NWORDS_FIELD);
                mp_add(x1, x2, x1, cwords);
                mp_shiftl1(x2, cwords);
            }
        } else {
            if (is_felm_even(v)) {
                mp_shiftr1(v, NWORDS_FIELD);
                mp_shiftl1(x1, NWORDS_FIELD);
            } else if (is_felm_even(u)) {
                mp_shiftr1(u, NWORDS_FIELD);
                mp_shiftl1(x2, NWORDS_FIELD);
            } else if (!is_felm_lt(v, u)) {
                mp_sub(v, u, v, NWORDS_FIELD);
                mp_shiftr1(v, NWORDS_FIELD);
                mp_add(x1, x2, x2, NWORDS_FIELD);
                mp_shiftl1(x1, NWORDS_FIELD);
            } else {
                mp_sub(u, v, u, NWORDS_FIELD);
                mp_shiftr1(u, NWORDS_FIELD);
                mp_add(x1, x2, x1, NWORDS_FIELD);
                mp_shiftl1(x2, NWORDS_FIELD);
            }
        }
        *k += 1;
    }

    if (is_felm_lt((digit_t*)PRIME, x1)) {
        mp_sub(x1, (digit_t*)PRIME, x1, NWORDS_FIELD);
    }
}


void fpinv_mont_bingcd(felm_t a)
{ // Field inversion via the binary GCD using Montgomery arithmetic, a = a^-1*r' mod p.
  // SECURITY NOTE: This function does not run in constant-time and is therefore only suitable for 
  //                operations not involving any secret data.
    felm_t x, t;
    unsigned int k;

    if (is_felm_zero(a) == true)
        return;

    fpinv_mont_bingcd_partial(a, x, &k);
    if (k <= MAXBITS_FIELD) { 
        fpmul_mont(x, (digit_t*)&Montgomery_R2, x);
        k += MAXBITS_FIELD;
    }
    fpmul_mont(x, (digit_t*)&Montgomery_R2, x);
    power2_setup(t, 2*MAXBITS_FIELD - k, NWORDS_FIELD);
    fpmul_mont(x, t, a);
}


void fp2inv_mont_bingcd(f2elm_t a)
{// GF(p^2) inversion using Montgomery arithmetic, a = (a0-i*a1)/(a0^2+a1^2)
 // This uses the binary GCD for inversion in fp and is NOT constant time!!!
    f2elm_t t1;

    fpsqr_mont(a[0], t1[0]);             // t10 = a0^2
    fpsqr_mont(a[1], t1[1]);             // t11 = a1^2
    fpadd(t1[0], t1[1], t1[0]);          // t10 = a0^2+a1^2
    fpinv_mont_bingcd(t1[0]);            // t10 = (a0^2+a1^2)^-1
    fpneg(a[1]);                         // a = a0-i*a1
    fpmul_mont(a[0], t1[0], a[0]);
    fpmul_mont(a[1], t1[0], a[1]);       // a = (a0-i*a1)*(a0^2+a1^2)^-1
}


void mont_n_way_inv(const f2elm_t* vec, const int n, f2elm_t* out)
{ // n-way simultaneous inversion using Montgomery's trick.
  // SECURITY NOTE: This function does not run in constant time.
  // Also, vec and out CANNOT be the same variable!
    f2elm_t t1;
    int i;

    fp2copy(vec[0], out[0]);                      // out[0] = vec[0]
    for (i = 1; i < n; i++) {
        fp2mul_mont(out[i-1], vec[i], out[i]);    // out[i] = out[i-1]*vec[i]
    }

    fp2copy(out[n-1], t1);                        // t1 = 1/out[n-1]
    fp2inv_mont_bingcd(t1);
    
    for (i = n-1; i >= 1; i--) {
        fp2mul_mont(out[i-1], t1, out[i]);        // out[i] = t1*out[i-1]
        fp2mul_mont(t1, vec[i], t1);              // t1 = t1*vec[i]
    }
    fp2copy(t1, out[0]);                          // out[0] = t1
}


void multiply(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{ // Multiprecision comba multiply, c = a*b, where lng(a) = lng(b) = nwords.
  // NOTE: a and c CANNOT be the same variable!
    unsigned int i, j, carry = 0;
    digit_t t = 0, u = 0, v = 0, UV[2];
    
    for (i = 0; i < nwords; i++) {
        for (j = 0; j <= i; j++) {
            MUL(a[j], b[i-j], UV+1, UV[0]);
            ADDC(0, UV[0], v, carry, v);
            ADDC(carry, UV[1], u, carry, u);
            t += carry;
        }
        c[i] = v;
        v = u;
        u = t;
        t = 0;
    }
    for (i = nwords; i < 2*nwords-1; i++) {
        for (j = i-nwords+1; j < nwords; j++) {
            MUL(a[j], b[i-j], UV+1, UV[0]);
            ADDC(0, UV[0], v, carry, v);
            ADDC(carry, UV[1], u, carry, u);
            t += carry;
        }
        c[i] = v;
        v = u;
        u = t;
        t = 0;
    }
    c[2*nwords-1] = v;
}


void Montgomery_neg(digit_t* a, digit_t* order)
{ // Modular negation, a = -a mod p.
  // Input/output: a in [0, 2*p-1] 
    unsigned int i, borrow = 0;
    
    for (i = 0; i < NWORDS_ORDER; i++) {
        SUBC(borrow, order[i], a[i], borrow, a[i]);
    }
}


void Montgomery_multiply_mod_order(const digit_t* ma, const digit_t* mb, digit_t* mc, const digit_t* order, const digit_t* Montgomery_rprime)
{ // Montgomery multiplication modulo the group order, mc = ma*mb*r' mod order, where ma,mb,mc in [0, order-1].
  // ma, mb and mc are assumed to be in Montgomery representation.
  // The Montgomery constant r' = -r^(-1) mod 2^(log_2(r)) is the value "Montgomery_rprime", where r is the order.  
  // Assume log_2(r) is a multiple of RADIX bits
    unsigned int i, cout = 0, bout = 0;
    digit_t mask, P[2*NWORDS_ORDER] = {0}, Q[2*NWORDS_ORDER] = {0}, temp[2*NWORDS_ORDER] = {0};

    multiply(ma, mb, P, NWORDS_ORDER);                 // P = ma * mb
    multiply(P, Montgomery_rprime, Q, NWORDS_ORDER);   // Q = P * r' mod 2^(log_2(r))
    multiply(Q, order, temp, NWORDS_ORDER);            // temp = Q * r
    cout = mp_add(P, temp, temp, 2*NWORDS_ORDER);      // (cout, temp) = P + Q * r     

    for (i = 0; i < NWORDS_ORDER; i++) {               // (cout, mc) = (P + Q * r)/2^(log_2(r))
        mc[i] = temp[NWORDS_ORDER+i];
    }

    // Final, constant-time subtraction     
    bout = mp_sub(mc, order, mc, NWORDS_ORDER);        // (cout, mc) = (cout, mc) - r
    mask = (digit_t)cout - (digit_t)bout;              // if (cout, mc) >= 0 then mask = 0x00..0, else if (cout, mc) < 0 then mask = 0xFF..F
    
    for (i = 0; i < NWORDS_ORDER; i++) {               // temp = mask & r
        temp[i] = (order[i] & mask);
    }
    
    mp_add(mc, temp, mc, NWORDS_ORDER);                //  mc = mc + (mask & r)
}


void to_Montgomery_mod_order(const digit_t* a, digit_t* mc, const digit_t* order, const digit_t* Montgomery_rprime, const digit_t* Montgomery_Rprime)
{ // Conversion of elements in Z_r to Montgomery representation, where the order r is up to NBITS_ORDER bits.
    Montgomery_multiply_mod_order(a, Montgomery_Rprime, mc, order, Montgomery_rprime);
}


void from_Montgomery_mod_order(const digit_t* ma, digit_t* c, const digit_t* order, const digit_t* Montgomery_rprime)
{ // Conversion of elements in Z_r from Montgomery to standard representation, where the order is up to NBITS_ORDER bits.
    digit_t one[NWORDS_ORDER] = {0};
    one[0] = 1;

    Montgomery_multiply_mod_order(ma, one, c, order, Montgomery_rprime);
}


static unsigned int is_zero_mod_order(const digit_t* x)
{ // Is x = 0? return 1 (TRUE) if condition is true, 0 (FALSE) otherwise
  // SECURITY NOTE: This function does not run in constant time.
    unsigned int i;

    for (i = 0; i < NWORDS_ORDER; i++) {
        if (x[i] != 0) return false;
    }
    return true;
}


static unsigned int is_even_mod_order(const digit_t* x)
{ // Is x even? return 1 (TRUE) if condition is true, 0 (FALSE) otherwise.
    return (unsigned int)((x[0] & 1) ^ 1);
}


static unsigned int is_lt_mod_order(const digit_t* x, const digit_t* y)
{ // Is x < y? return 1 (TRUE) if condition is true, 0 (FALSE) otherwise.
  // SECURITY NOTE: This function does not run in constant time.
    int i;

    for (i = NWORDS_ORDER-1; i >= 0; i--) {
        if (x[i] < y[i]) { 
            return true;
        } else if (x[i] > y[i]) {
            return false;
        }
    }
    return false;
}


static void Montgomery_inversion_mod_order_bingcd_partial(const digit_t* a, digit_t* x1, unsigned int* k, const digit_t* order)
{ // Partial Montgomery inversion modulo order.
    digit_t u[NWORDS_ORDER], v[NWORDS_ORDER], x2[NWORDS_ORDER] = {0};
    unsigned int cwords;  // number of words necessary for x1, x2

    copy_words(a, u, NWORDS_ORDER);
    copy_words(order, v, NWORDS_ORDER);
    copy_words(x2, x1, NWORDS_ORDER);
    x1[0] = 1;
    *k = 0;

    while (!is_zero_mod_order(v)) {
        cwords = ((*k + 1) / RADIX) + 1;
        if ((cwords < NWORDS_ORDER)) {
            if (is_even_mod_order(v)) {
                mp_shiftr1(v, NWORDS_ORDER);
                mp_shiftl1(x1, cwords);
            } else if (is_even_mod_order(u)) {
                mp_shiftr1(u, NWORDS_ORDER);
                mp_shiftl1(x2, cwords);
            } else if (!is_lt_mod_order(v, u)) {
                mp_sub(v, u, v, NWORDS_ORDER);
                mp_shiftr1(v, NWORDS_ORDER);
                mp_add(x1, x2, x2, cwords);
                mp_shiftl1(x1, cwords);
            } else {
                mp_sub(u, v, u, NWORDS_ORDER);
                mp_shiftr1(u, NWORDS_ORDER);
                mp_add(x1, x2, x1, cwords);
                mp_shiftl1(x2, cwords);
            }
        } else {
            if (is_even_mod_order(v)) {
                mp_shiftr1(v, NWORDS_ORDER);
                mp_shiftl1(x1, NWORDS_ORDER);
            } else if (is_even_mod_order(u)) {
                mp_shiftr1(u, NWORDS_ORDER);
                mp_shiftl1(x2, NWORDS_ORDER);
            } else if (!is_lt_mod_order(v, u)) {
                mp_sub(v, u, v, NWORDS_ORDER);
                mp_shiftr1(v, NWORDS_ORDER);
                mp_add(x1, x2, x2, NWORDS_ORDER);
                mp_shiftl1(x1, NWORDS_ORDER);
            } else {
                mp_sub(u, v, u, NWORDS_ORDER);
                mp_shiftr1(u, NWORDS_ORDER);
                mp_add(x1, x2, x1, NWORDS_ORDER);
                mp_shiftl1(x2, NWORDS_ORDER);
            }
        }
        *k += 1;
    }

    if (is_lt_mod_order(order, x1)) {
        mp_sub(x1, order, x1, NWORDS_ORDER);
    }
}


void Montgomery_inversion_mod_order_bingcd(const digit_t* a, digit_t* c, const digit_t* order, const digit_t* Montgomery_rprime, const digit_t* Montgomery_Rprime)
{// Montgomery inversion modulo order, c = a^(-1)*R mod order.
    digit_t x[NWORDS_ORDER], t[NWORDS_ORDER] = {0};
    unsigned int k;

    if (is_zero((digit_t*)a, NWORDS_ORDER) == true) {
        copy_words(t, c, NWORDS_ORDER);
        return;
    }

    Montgomery_inversion_mod_order_bingcd_partial(a, x, &k, order);
    if (k <= NBITS_ORDER) {
        Montgomery_multiply_mod_order(x, Montgomery_Rprime, x, order, Montgomery_rprime);
        k += NBITS_ORDER;
    }

    Montgomery_multiply_mod_order(x, Montgomery_Rprime, x, order, Montgomery_rprime);
    power2_setup(t, 2*NBITS_ORDER - k, NWORDS_ORDER);
    Montgomery_multiply_mod_order(x, t, c, order, Montgomery_rprime);
}


void inv_mod_orderA(const digit_t* a, digit_t* c)
{ // Inversion of an odd integer modulo an even integer of the form 2^m.
  // Algorithm 3: Explicit Quadratic Modular inverse modulo 2^m from Dumas'12: http://arxiv.org/pdf/1209.6626.pdf
  // If the input is invalid (even), the function outputs c = a.
    unsigned int i, f, s = 0;
    digit_t am1[NWORDS_ORDER] = {0};
    digit_t tmp1[NWORDS_ORDER] = {0};
    digit_t tmp2[2*NWORDS_ORDER] = {0};
    digit_t one[NWORDS_ORDER] = {0};
    digit_t order[NWORDS_ORDER] = {0};
    digit_t mask = (digit_t)((uint64_t)(-1) >> (NBITS_ORDER - OALICE_BITS));

    order[NWORDS_ORDER-1] = (digit_t)((uint64_t)1 << (64 - (NBITS_ORDER - OALICE_BITS)));  // Load most significant digit of Alice's order
    one[0] = 1;
        
    mp_sub(a, one, am1, NWORDS_ORDER);                   // am1 = a-1

    if (((a[0] & (digit_t)1) == 0) || (is_zero(am1, NWORDS_ORDER) == true)) {  // Check if the input is even or one 
        copy_words(a, c, NWORDS_ORDER);
        c[NWORDS_ORDER-1] &= mask;                       // mod 2^m
    } else { 
        mp_sub(order, am1, c, NWORDS_ORDER);
        mp_add(c, one, c, NWORDS_ORDER);                 // c = 2^m - a + 2

        copy_words(am1, tmp1, NWORDS_ORDER);
        while ((tmp1[0] & (digit_t)1) == 0) {
            s += 1;
            mp_shiftr1(tmp1, NWORDS_ORDER);
        }

        f = OALICE_BITS / s;
        for (i = 1; i < f; i <<= 1) {
            multiply(am1, am1, tmp2, NWORDS_ORDER);            // tmp2 = am1^2  
            copy_words(tmp2, am1, NWORDS_ORDER);
            am1[NWORDS_ORDER-1] &= mask;                       // am1 = tmp2 mod 2^m
            mp_add(am1, one, tmp1, NWORDS_ORDER);              // tmp1 = am1 + 1
            tmp1[NWORDS_ORDER-1] &= mask;                      // mod 2^m
            multiply(c, tmp1, tmp2, NWORDS_ORDER);             // c = c*tmp1
            copy_words(tmp2, c, NWORDS_ORDER);
            c[NWORDS_ORDER-1] &= mask;                         // mod 2^m
        }
    }
}


void recover_os(const f2elm_t X1, const f2elm_t Z1, const f2elm_t X2, const f2elm_t Z2, const f2elm_t x, const f2elm_t y, const f2elm_t A, f2elm_t X3, f2elm_t Y3, f2elm_t Z3)
{
    f2elm_t t0, t1, t2, t3;
    
    // X3 := 2*y*Z1*Z2*X1;
    // Y3 := Z2*((X1+x*Z1+2*A*Z1)*(X1*x+Z1)-2*A*Z1^2)-(X1-x*Z1)^2*X2;
    // Z3 := 2*y*Z1*Z2*Z1;
    
    fp2add(y, y, t0);
    fp2mul_mont(t0, Z1, t0);
    fp2mul_mont(t0, Z2, t0);       // t0 = 2*y*Z1*Z2
    fp2mul_mont(t0, Z1, Z3);       // Z3 = 2*y*Z1*Z2*Z1       
    fp2mul_mont(t0, X1, X3);       // X3 = 2*y*Z1*Z2*X1
    fp2add(A, A, t0);
    fp2mul_mont(t0, Z1, t0);       // t0 = 2*A*Z1  
    fp2mul_mont(x, Z1, t1);        // t1 = x*Z1  
    fp2add(X1, t1, t2);            // t2 = X1+x*Z1
    fp2sub(X1, t1, t1);            // t1 = X1-x*Z1
    fp2add(t0, t2, t3);            // t3 = X1+x*Z1+2*A*Z1
    fp2mul_mont(t0, Z1, t0);       // t0 = 2*A*Z1^2 
    fp2sqr_mont(t1, t1);           // t1 = (X1-x*Z1)^2
    fp2mul_mont(x, X1, t2);        // t2 = x*X1
    fp2add(t2, Z1, t2);            // t2 = X1*x+Z1
    fp2mul_mont(t2, t3, t2);       // t2 = (X1+x*Z1+2*A*Z1)*(X1*x+Z1)
    fp2sub(t2, t0, t0);            // t0 = (X1+x*Z1+2*A*Z1)*(X1*x+Z1)-2*A*Z1^2
    fp2mul_mont(t1, X2, t1);       // t1 = (X1-x*Z1)^2*X2
    fp2mul_mont(t0, Z2, t0);       // t0 = Z2*[(X1+x*Z1+2*A*Z1)*(X1*x+Z1)-2*A*Z1^2]
    fp2sub(t0, t1, Y3);            // Y3 = Z2*[(X1+x*Z1+2*A*Z1)*(X1*x+Z1)-2*A*Z1^2] - (X1-x*Z1)^2*X2
}
// Closing COMPRESSED
#endif

/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: elliptic curve and isogeny functions
*********************************************************************************************/


void xDBL(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24)
{ // Doubling of a Montgomery point in projective coordinates (X:Z).
  // Input: projective Montgomery x-coordinates P = (X1:Z1), where x1=X1/Z1 and Montgomery curve constants A+2C and 4C.
  // Output: projective Montgomery x-coordinates Q = 2*P = (X2:Z2).
    f2elm_t t0, t1;
    
    mp2_sub_p2(P->X, P->Z, t0);                     // t0 = X1-Z1
    mp2_add(P->X, P->Z, t1);                        // t1 = X1+Z1
    fp2sqr_mont(t0, t0);                            // t0 = (X1-Z1)^2 
    fp2sqr_mont(t1, t1);                            // t1 = (X1+Z1)^2 
    fp2mul_mont(C24, t0, Q->Z);                     // Z2 = C24*(X1-Z1)^2   
    fp2mul_mont(t1, Q->Z, Q->X);                    // X2 = C24*(X1-Z1)^2*(X1+Z1)^2
    mp2_sub_p2(t1, t0, t1);                         // t1 = (X1+Z1)^2-(X1-Z1)^2 
    fp2mul_mont(A24plus, t1, t0);                   // t0 = A24plus*[(X1+Z1)^2-(X1-Z1)^2]
    mp2_add(Q->Z, t0, Q->Z);                        // Z2 = A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2
    fp2mul_mont(Q->Z, t1, Q->Z);                    // Z2 = [A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2]*[(X1+Z1)^2-(X1-Z1)^2]
}


void xDBLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24, const int e)
{ // Computes [2^e](X:Z) on Montgomery curve with projective constant via e repeated doublings.
  // Input: projective Montgomery x-coordinates P = (XP:ZP), such that xP=XP/ZP and Montgomery curve constants A+2C and 4C.
  // Output: projective Montgomery x-coordinates Q <- (2^e)*P.
    int i;
    
    copy_words((digit_t*)P, (digit_t*)Q, 2*2*NWORDS_FIELD);

    for (i = 0; i < e; i++) {
        xDBL(Q, Q, A24plus, C24);
    }
}

#if (OALICE_BITS % 2 == 1)

void get_2_isog(const point_proj_t P, f2elm_t A, f2elm_t C)
{ // Computes the corresponding 2-isogeny of a projective Montgomery point (X2:Z2) of order 2.
  // Input:  projective point of order two P = (X2:Z2).
  // Output: the 2-isogenous Montgomery curve with projective coefficients A/C.
    
    fp2sqr_mont(P->X, A);                           // A = X2^2
    fp2sqr_mont(P->Z, C);                           // C = Z2^2
    mp2_sub_p2(C, A, A);                            // A = Z2^2 - X2^2
}


void eval_2_isog(point_proj_t P, point_proj_t Q)
{ // Evaluates the isogeny at the point (X:Z) in the domain of the isogeny, given a 2-isogeny phi.
  // Inputs: the projective point P = (X:Z) and the 2-isogeny kernel projetive point Q = (X2:Z2).
  // Output: the projective point P = phi(P) = (X:Z) in the codomain. 
    f2elm_t t0, t1, t2, t3;
    
    mp2_add(Q->X, Q->Z, t0);                        // t0 = X2+Z2
    mp2_sub_p2(Q->X, Q->Z, t1);                     // t1 = X2-Z2
    mp2_add(P->X, P->Z, t2);                        // t2 = X+Z
    mp2_sub_p2(P->X, P->Z, t3);                     // t3 = X-Z
    fp2mul_mont(t0, t3, t0);                        // t0 = (X2+Z2)*(X-Z)
    fp2mul_mont(t1, t2, t1);                        // t1 = (X2-Z2)*(X+Z)
    mp2_add(t0, t1, t2);                            // t2 = (X2+Z2)*(X-Z) + (X2-Z2)*(X+Z)
    mp2_sub_p2(t0, t1, t3);                         // t3 = (X2+Z2)*(X-Z) - (X2-Z2)*(X+Z)
    fp2mul_mont(P->X, t2, P->X);                    // Xfinal
    fp2mul_mont(P->Z, t3, P->Z);                    // Zfinal
}

#endif

void get_4_isog(const point_proj_t P, f2elm_t A24plus, f2elm_t C24, f2elm_t* coeff)
{ // Computes the corresponding 4-isogeny of a projective Montgomery point (X4:Z4) of order 4.
  // Input:  projective point of order four P = (X4:Z4).
  // Output: the 4-isogenous Montgomery curve with projective coefficients A+2C/4C and the 3 coefficients 
  //         that are used to evaluate the isogeny at a point in eval_4_isog().
    
    mp2_sub_p2(P->X, P->Z, coeff[1]);               // coeff[1] = X4-Z4
    mp2_add(P->X, P->Z, coeff[2]);                  // coeff[2] = X4+Z4
    fp2sqr_mont(P->Z, coeff[0]);                    // coeff[0] = Z4^2
    mp2_add(coeff[0], coeff[0], coeff[0]);          // coeff[0] = 2*Z4^2
    fp2sqr_mont(coeff[0], C24);                     // C24 = 4*Z4^4
    mp2_add(coeff[0], coeff[0], coeff[0]);          // coeff[0] = 4*Z4^2
    fp2sqr_mont(P->X, A24plus);                     // A24plus = X4^2
    mp2_add(A24plus, A24plus, A24plus);             // A24plus = 2*X4^2
    fp2sqr_mont(A24plus, A24plus);                  // A24plus = 4*X4^4
}


void eval_4_isog(point_proj_t P, f2elm_t* coeff)
{ // Evaluates the isogeny at the point (X:Z) in the domain of the isogeny, given a 4-isogeny phi defined 
  // by the 3 coefficients in coeff (computed in the function get_4_isog()).
  // Inputs: the coefficients defining the isogeny, and the projective point P = (X:Z).
  // Output: the projective point P = phi(P) = (X:Z) in the codomain. 
    f2elm_t t0, t1;
    
    mp2_add(P->X, P->Z, t0);                        // t0 = X+Z
    mp2_sub_p2(P->X, P->Z, t1);                     // t1 = X-Z
    fp2mul_mont(t0, coeff[1], P->X);                // X = (X+Z)*coeff[1]
    fp2mul_mont(t1, coeff[2], P->Z);                // Z = (X-Z)*coeff[2]
    fp2mul_mont(t0, t1, t0);                        // t0 = (X+Z)*(X-Z)
    fp2mul_mont(coeff[0], t0, t0);                  // t0 = coeff[0]*(X+Z)*(X-Z)
    mp2_add(P->X, P->Z, t1);                        // t1 = (X-Z)*coeff[2] + (X+Z)*coeff[1]
    mp2_sub_p2(P->X, P->Z, P->Z);                   // Z = (X-Z)*coeff[2] - (X+Z)*coeff[1]
    fp2sqr_mont(t1, t1);                            // t1 = [(X-Z)*coeff[2] + (X+Z)*coeff[1]]^2
    fp2sqr_mont(P->Z, P->Z);                        // Z = [(X-Z)*coeff[2] - (X+Z)*coeff[1]]^2
    mp2_add(t1, t0, P->X);                          // X = coeff[0]*(X+Z)*(X-Z) + [(X-Z)*coeff[2] + (X+Z)*coeff[1]]^2
    mp2_sub_p2(P->Z, t0, t0);                       // t0 = [(X-Z)*coeff[2] - (X+Z)*coeff[1]]^2 - coeff[0]*(X+Z)*(X-Z)
    fp2mul_mont(P->X, t1, P->X);                    // Xfinal
    fp2mul_mont(P->Z, t0, P->Z);                    // Zfinal
}


void xTPL(const point_proj_t P, point_proj_t Q, const f2elm_t A24minus, const f2elm_t A24plus)              
{ // Tripling of a Montgomery point in projective coordinates (X:Z).
  // Input: projective Montgomery x-coordinates P = (X:Z), where x=X/Z and Montgomery curve constants A24plus = A+2C and A24minus = A-2C.
  // Output: projective Montgomery x-coordinates Q = 3*P = (X3:Z3).
    f2elm_t t0, t1, t2, t3, t4, t5, t6;
                                    
    mp2_sub_p2(P->X, P->Z, t0);                     // t0 = X-Z 
    fp2sqr_mont(t0, t2);                            // t2 = (X-Z)^2           
    mp2_add(P->X, P->Z, t1);                        // t1 = X+Z 
    fp2sqr_mont(t1, t3);                            // t3 = (X+Z)^2
    mp2_add(P->X, P->X, t4);                        // t4 = 2*X
    mp2_add(P->Z, P->Z, t0);                        // t0 = 2*Z 
    fp2sqr_mont(t4, t1);                            // t1 = 4*X^2
    mp2_sub_p2(t1, t3, t1);                         // t1 = 4*X^2 - (X+Z)^2 
    mp2_sub_p2(t1, t2, t1);                         // t1 = 4*X^2 - (X+Z)^2 - (X-Z)^2
    fp2mul_mont(A24plus, t3, t5);                   // t5 = A24plus*(X+Z)^2 
    fp2mul_mont(t3, t5, t3);                        // t3 = A24plus*(X+Z)^4
    fp2mul_mont(A24minus, t2, t6);                  // t6 = A24minus*(X-Z)^2
    fp2mul_mont(t2, t6, t2);                        // t2 = A24minus*(X-Z)^4
    mp2_sub_p2(t2, t3, t3);                         // t3 = A24minus*(X-Z)^4 - A24plus*(X+Z)^4
    mp2_sub_p2(t5, t6, t2);                         // t2 = A24plus*(X+Z)^2 - A24minus*(X-Z)^2
    fp2mul_mont(t1, t2, t1);                        // t1 = [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2]
    fp2add(t3, t1, t2);                             // t2 = [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2] + A24minus*(X-Z)^4 - A24plus*(X+Z)^4
    fp2sqr_mont(t2, t2);                            // t2 = t2^2
    fp2mul_mont(t4, t2, Q->X);                      // X3 = 2*X*t2
    fp2sub(t3, t1, t1);                             // t1 = A24minus*(X-Z)^4 - A24plus*(X+Z)^4 - [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2]
    fp2sqr_mont(t1, t1);                            // t1 = t1^2
    fp2mul_mont(t0, t1, Q->Z);                      // Z3 = 2*Z*t1
}


void xTPLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24minus, const f2elm_t A24plus, const int e)
{ // Computes [3^e](X:Z) on Montgomery curve with projective constant via e repeated triplings.
  // Input: projective Montgomery x-coordinates P = (XP:ZP), such that xP=XP/ZP and Montgomery curve constants A24plus = A+2C and A24minus = A-2C.
  // Output: projective Montgomery x-coordinates Q <- (3^e)*P.
    int i;
        
    copy_words((digit_t*)P, (digit_t*)Q, 2*2*NWORDS_FIELD);

    for (i = 0; i < e; i++) {
        xTPL(Q, Q, A24minus, A24plus);
    }
}


void get_3_isog(const point_proj_t P, f2elm_t A24minus, f2elm_t A24plus, f2elm_t* coeff)
{ // Computes the corresponding 3-isogeny of a projective Montgomery point (X3:Z3) of order 3.
  // Input:  projective point of order three P = (X3:Z3).
  // Output: the 3-isogenous Montgomery curve with projective coefficient A/C. 
    f2elm_t t0, t1, t2, t3, t4;
    
    mp2_sub_p2(P->X, P->Z, coeff[0]);               // coeff0 = X-Z
    fp2sqr_mont(coeff[0], t0);                      // t0 = (X-Z)^2
    mp2_add(P->X, P->Z, coeff[1]);                  // coeff1 = X+Z
    fp2sqr_mont(coeff[1], t1);                      // t1 = (X+Z)^2
    mp2_add(P->X, P->X, t3);                        // t3 = 2*X
    fp2sqr_mont(t3, t3);                            // t3 = 4*X^2 
    fp2sub(t3, t0, t2);                             // t2 = 4*X^2 - (X-Z)^2 
    fp2sub(t3, t1, t3);                             // t3 = 4*X^2 - (X+Z)^2
    mp2_add(t0, t3, t4);                            // t4 = 4*X^2 - (X+Z)^2 + (X-Z)^2 
    mp2_add(t4, t4, t4);                            // t4 = 2(4*X^2 - (X+Z)^2 + (X-Z)^2) 
    mp2_add(t1, t4, t4);                            // t4 = 8*X^2 - (X+Z)^2 + 2*(X-Z)^2
    fp2mul_mont(t2, t4, A24minus);                  // A24minus = [4*X^2 - (X-Z)^2]*[8*X^2 - (X+Z)^2 + 2*(X-Z)^2]
    mp2_add(t1, t2, t4);                            // t4 = 4*X^2 + (X+Z)^2 - (X-Z)^2
    mp2_add(t4, t4, t4);                            // t4 = 2(4*X^2 + (X+Z)^2 - (X-Z)^2) 
    mp2_add(t0, t4, t4);                            // t4 = 8*X^2 + 2*(X+Z)^2 - (X-Z)^2
    fp2mul_mont(t3, t4, A24plus);                   // A24plus = [4*X^2 - (X+Z)^2]*[8*X^2 + 2*(X+Z)^2 - (X-Z)^2]
}


void eval_3_isog(point_proj_t Q, const f2elm_t* coeff)
{ // Computes the 3-isogeny R=phi(X:Z), given projective point (X3:Z3) of order 3 on a Montgomery curve and 
  // a point P with 2 coefficients in coeff (computed in the function get_3_isog()).
  // Inputs: projective points P = (X3:Z3) and Q = (X:Z).
  // Output: the projective point Q <- phi(Q) = (X3:Z3). 
    f2elm_t t0, t1, t2;

    mp2_add(Q->X, Q->Z, t0);                      // t0 = X+Z
    mp2_sub_p2(Q->X, Q->Z, t1);                   // t1 = X-Z
    fp2mul_mont(coeff[0], t0, t0);                // t0 = coeff0*(X+Z)
    fp2mul_mont(coeff[1], t1, t1);                // t1 = coeff1*(X-Z)
    mp2_add(t0, t1, t2);                          // t2 = coeff0*(X+Z) + coeff1*(X-Z)
    mp2_sub_p2(t1, t0, t0);                       // t0 = coeff1*(X-Z) - coeff0*(X+Z)
    fp2sqr_mont(t2, t2);                          // t2 = [coeff0*(X+Z) + coeff1*(X-Z)]^2
    fp2sqr_mont(t0, t0);                          // t0 = [coeff1*(X-Z) - coeff0*(X+Z)]^2
    fp2mul_mont(Q->X, t2, Q->X);                  // X3final = X*[coeff0*(X+Z) + coeff1*(X-Z)]^2        
    fp2mul_mont(Q->Z, t0, Q->Z);                  // Z3final = Z*[coeff1*(X-Z) - coeff0*(X+Z)]^2
}


void inv_3_way(f2elm_t z1, f2elm_t z2, f2elm_t z3)
{ // 3-way simultaneous inversion
  // Input:  z1,z2,z3
  // Output: 1/z1,1/z2,1/z3 (override inputs).
    f2elm_t t0, t1, t2, t3;

    fp2mul_mont(z1, z2, t0);                      // t0 = z1*z2
    fp2mul_mont(z3, t0, t1);                      // t1 = z1*z2*z3
    fp2inv_mont(t1);                              // t1 = 1/(z1*z2*z3)
    fp2mul_mont(z3, t1, t2);                      // t2 = 1/(z1*z2) 
    fp2mul_mont(t2, z2, t3);                      // t3 = 1/z1
    fp2mul_mont(t2, z1, z2);                      // z2 = 1/z2
    fp2mul_mont(t0, t1, z3);                      // z3 = 1/z3
    fp2copy(t3, z1);                              // z1 = 1/z1
}


void get_A(const f2elm_t xP, const f2elm_t xQ, const f2elm_t xR, f2elm_t A)
{ // Given the x-coordinates of P, Q, and R, returns the value A corresponding to the Montgomery curve E_A: y^2=x^3+A*x^2+x such that R=Q-P on E_A.
  // Input:  the x-coordinates xP, xQ, and xR of the points P, Q and R.
  // Output: the coefficient A corresponding to the curve E_A: y^2=x^3+A*x^2+x.
    f2elm_t t0, t1, one = {0};
    
    fpcopy((digit_t*)&Montgomery_one, one[0]);
    fp2add(xP, xQ, t1);                           // t1 = xP+xQ
    fp2mul_mont(xP, xQ, t0);                      // t0 = xP*xQ
    fp2mul_mont(xR, t1, A);                       // A = xR*t1
    fp2add(t0, A, A);                             // A = A+t0
    fp2mul_mont(t0, xR, t0);                      // t0 = t0*xR
    fp2sub(A, one, A);                            // A = A-1
    fp2add(t0, t0, t0);                           // t0 = t0+t0
    fp2add(t1, xR, t1);                           // t1 = t1+xR
    fp2add(t0, t0, t0);                           // t0 = t0+t0
    fp2sqr_mont(A, A);                            // A = A^2
    fp2inv_mont(t0);                              // t0 = 1/t0
    fp2mul_mont(A, t0, A);                        // A = A*t0
    fp2sub(A, t1, A);                             // Afinal = A-t1
}


void j_inv(const f2elm_t A, const f2elm_t C, f2elm_t jinv)
{ // Computes the j-invariant of a Montgomery curve with projective constant.
  // Input: A,C in GF(p^2).
  // Output: j=256*(A^2-3*C^2)^3/(C^4*(A^2-4*C^2)), which is the j-invariant of the Montgomery curve B*y^2=x^3+(A/C)*x^2+x or (equivalently) j-invariant of B'*y^2=C*x^3+A*x^2+C*x.
    f2elm_t t0, t1;
    
    fp2sqr_mont(A, jinv);                           // jinv = A^2        
    fp2sqr_mont(C, t1);                             // t1 = C^2
    fp2add(t1, t1, t0);                             // t0 = t1+t1
    fp2sub(jinv, t0, t0);                           // t0 = jinv-t0
    fp2sub(t0, t1, t0);                             // t0 = t0-t1
    fp2sub(t0, t1, jinv);                           // jinv = t0-t1
    fp2sqr_mont(t1, t1);                            // t1 = t1^2
    fp2mul_mont(jinv, t1, jinv);                    // jinv = jinv*t1
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2sqr_mont(t0, t1);                            // t1 = t0^2
    fp2mul_mont(t0, t1, t0);                        // t0 = t0*t1
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2inv_mont(jinv);                              // jinv = 1/jinv 
    fp2mul_mont(jinv, t0, jinv);                    // jinv = t0*jinv
}


void xDBLADD(point_proj_t P, point_proj_t Q, const f2elm_t xPQ, const f2elm_t A24)
{ // Simultaneous doubling and differential addition.
  // Input: projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, affine difference xPQ=x(P-Q) and Montgomery curve constant A24=(A+2)/4.
  // Output: projective Montgomery points P <- 2*P = (X2P:Z2P) such that x(2P)=X2P/Z2P, and Q <- P+Q = (XQP:ZQP) such that = x(Q+P)=XQP/ZQP. 
    f2elm_t t0, t1, t2;

    mp2_add(P->X, P->Z, t0);                        // t0 = XP+ZP
    mp2_sub_p2(P->X, P->Z, t1);                     // t1 = XP-ZP
    fp2sqr_mont(t0, P->X);                          // XP = (XP+ZP)^2
    mp2_sub_p2(Q->X, Q->Z, t2);                     // t2 = XQ-ZQ
    mp2_add(Q->X, Q->Z, Q->X);                      // XQ = XQ+ZQ
    fp2mul_mont(t0, t2, t0);                        // t0 = (XP+ZP)*(XQ-ZQ)
    fp2sqr_mont(t1, P->Z);                          // ZP = (XP-ZP)^2
    fp2mul_mont(t1, Q->X, t1);                      // t1 = (XP-ZP)*(XQ+ZQ)
    mp2_sub_p2(P->X, P->Z, t2);                     // t2 = (XP+ZP)^2-(XP-ZP)^2
    fp2mul_mont(P->X, P->Z, P->X);                  // XP = (XP+ZP)^2*(XP-ZP)^2
    fp2mul_mont(A24, t2, Q->X);                     // XQ = A24*[(XP+ZP)^2-(XP-ZP)^2]
    mp2_sub_p2(t0, t1, Q->Z);                       // ZQ = (XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)
    mp2_add(Q->X, P->Z, P->Z);                      // ZP = A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2
    mp2_add(t0, t1, Q->X);                          // XQ = (XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)
    fp2mul_mont(P->Z, t2, P->Z);                    // ZP = [A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2]*[(XP+ZP)^2-(XP-ZP)^2]
    fp2sqr_mont(Q->Z, Q->Z);                        // ZQ = [(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
    fp2sqr_mont(Q->X, Q->X);                        // XQ = [(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2
    fp2mul_mont(Q->Z, xPQ, Q->Z);                   // ZQ = xPQ*[(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
}


static void swap_points(point_proj_t P, point_proj_t Q, const digit_t option)
{ // Swap points.
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    digit_t temp;
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++) {
        temp = option & (P->X[0][i] ^ Q->X[0][i]);
        P->X[0][i] = temp ^ P->X[0][i]; 
        Q->X[0][i] = temp ^ Q->X[0][i];  
        temp = option & (P->X[1][i] ^ Q->X[1][i]);
        P->X[1][i] = temp ^ P->X[1][i]; 
        Q->X[1][i] = temp ^ Q->X[1][i];
        temp = option & (P->Z[0][i] ^ Q->Z[0][i]);
        P->Z[0][i] = temp ^ P->Z[0][i]; 
        Q->Z[0][i] = temp ^ Q->Z[0][i];
        temp = option & (P->Z[1][i] ^ Q->Z[1][i]);
        P->Z[1][i] = temp ^ P->Z[1][i]; 
        Q->Z[1][i] = temp ^ Q->Z[1][i]; 
    }
}


static void LADDER3PT(const f2elm_t xP, const f2elm_t xQ, const f2elm_t xPQ, const digit_t* m, const unsigned int AliceOrBob, point_proj_t R, const f2elm_t A)
{
    point_proj_t R0 = {0}, R2 = {0};
    f2elm_t A24 = {0};
    digit_t mask;
    int i, nbits, bit, swap, prevbit = 0;

    if (AliceOrBob == ALICE) {
        nbits = OALICE_BITS;
    } else {
        nbits = OBOB_BITS - 1;
    }

    // Initializing constant
    fpcopy((digit_t*)&Montgomery_one, A24[0]);
    mp2_add(A24, A24, A24);
    mp2_add(A, A24, A24);
    fp2div2(A24, A24);  
    fp2div2(A24, A24);  // A24 = (A+2)/4

    // Initializing points
    fp2copy(xQ, R0->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R0->Z);
    fp2copy(xPQ, R2->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R2->Z);
    fp2copy(xP, R->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R->Z);
    fpzero((digit_t*)(R->Z)[1]);

    // Main loop
    for (i = 0; i < nbits; i++) {
        bit = (m[i >> LOG2RADIX] >> (i & (RADIX-1))) & 1;
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(R, R2, mask);
        xDBLADD(R0, R2, R->X, A24);
        fp2mul_mont(R2->X, R->Z, R2->X);
    }
    swap = 0 ^ prevbit;
    mask = 0 - (digit_t)swap;
    swap_points(R, R2, mask);
}

#ifdef COMPRESS

static void CompletePoint(const point_proj_t P, point_full_proj_t R)
{ // Complete point on A = 0 curve
    f2elm_t xz, s2, r2, yz, invz, t0, t1, one = {0};

    fpcopy((digit_t*)&Montgomery_one, one[0]);
    fp2mul_mont(P->X, P->Z, xz);
    fpsub(P->X[0], P->Z[1], t0[0]);
    fpadd(P->X[1], P->Z[0], t0[1]);
    fpadd(P->X[0], P->Z[1], t1[0]);
    fpsub(P->X[1], P->Z[0], t1[1]);
    fp2mul_mont(t0, t1, s2);
    fp2mul_mont(xz, s2, r2);
    sqrt_Fp2(r2, yz);
    fp2copy(P->Z,invz);
    fp2inv_mont_bingcd(invz);    
    fp2mul_mont(P->X, invz, R->X);
    fp2sqr_mont(invz, t0);
    fp2mul_mont(yz, t0, R->Y);
    fp2copy(one, R->Z);
}


void CompleteMPoint(const f2elm_t A, point_proj_t P, point_full_proj_t R)
{ // Given an xz-only representation on a montgomery curve, compute its affine representation
    f2elm_t zero = {0}, one = {0}, xz, yz, s2, r2, invz, temp0, temp1;

    fpcopy((digit_t*)&Montgomery_one, one[0]);    
    if (memcmp(P->Z[0], zero,NBITS_TO_NBYTES(NBITS_FIELD)) != 0 || memcmp(P->Z[1], zero, NBITS_TO_NBYTES(NBITS_FIELD)) != 0) {
        fp2mul_mont(P->X, P->Z, xz);       // xz = x*z;
        fpsub(P->X[0], P->Z[1], temp0[0]);
        fpadd(P->X[1], P->Z[0], temp0[1]);
        fpadd(P->X[0], P->Z[1], temp1[0]);
        fpsub(P->X[1], P->Z[0], temp1[1]);        
        fp2mul_mont(temp0, temp1, s2);     // s2 = (x + i*z)*(x - i*z);
        fp2mul_mont(A, xz, temp0);
        fp2add(temp0, s2, temp1);
        fp2mul_mont(xz, temp1, r2);        // r2 = xz*(A*xz + s2);
        sqrt_Fp2(r2, yz);
        fp2copy(P->Z, invz);
        fp2inv_mont_bingcd(invz);        
        fp2mul_mont(P->X, invz, R->X);
        fp2sqr_mont(invz, temp0);
        fp2mul_mont(yz, temp0, R->Y);      // R = EM![x*invz, yz*invz^2];
        fp2copy(one, R->Z);
    } else {
        fp2copy(zero, R->X);
        fp2copy(one, R->Y); 
        fp2copy(zero, R->Z);               // R = EM!0;
    }
}


void Double(point_proj_t P, point_proj_t Q, f2elm_t A24, const int k)
{ // Doubling of a Montgomery point in projective coordinates (X:Z) over affine curve coefficient A. 
  // Input: projective Montgomery x-coordinates P = (X1:Z1), where x1=X1/Z1 and Montgomery curve constants (A+2)/4.
  // Output: projective Montgomery x-coordinates Q = 2*P = (X2:Z2). 
    f2elm_t temp, a, b, c, aa, bb;    
    
    fp2copy(P->X, Q->X);
    fp2copy(P->Z, Q->Z);
    
    for (int j = 0; j < k; j++) {
        fp2add(Q->X, Q->Z, a);
        fp2sub(Q->X, Q->Z, b);
        fp2sqr_mont(a, aa);
        fp2sqr_mont(b, bb);
        fp2sub(aa, bb, c);
        fp2mul_mont(aa, bb, Q->X);
        fp2mul_mont(A24, c, temp);
        fp2add(temp, bb, temp);
        fp2mul_mont(c, temp, Q->Z);
    }
}


void xTPL_fast(const point_proj_t P, point_proj_t Q, const f2elm_t A2)
{ // Montgomery curve (E: y^2 = x^3 + A*x^2 + x) x-only tripling at a cost 5M + 6S + 9A = 27p + 61a.
  // Input : projective Montgomery x-coordinates P = (X:Z), where x=X/Z and Montgomery curve constant A/2. 
  // Output: projective Montgomery x-coordinates Q = 3*P = (X3:Z3).
       f2elm_t t1, t2, t3, t4;
       
       fp2sqr_mont(P->X, t1);        // t1 = x^2
       fp2sqr_mont(P->Z, t2);        // t2 = z^2
       fp2add(t1, t2, t3);           // t3 = t1 + t2
       fp2add(P->X, P->Z, t4);       // t4 = x + z
       fp2sqr_mont(t4, t4);          // t4 = t4^2
       fp2sub(t4, t3, t4);           // t4 = t4 - t3
       fp2mul_mont(A2, t4, t4);      // t4 = t4*A2
       fp2add(t3, t4, t4);           // t4 = t4 + t3
       fp2sub(t1, t2, t3);           // t3 = t1 - t2
       fp2sqr_mont(t3, t3);          // t3 = t3^2
       fp2mul_mont(t1, t4, t1);      // t1 = t1*t4
       fp2shl(t1, 2, t1);            // t1 = 4*t1
       fp2sub(t1, t3, t1);           // t1 = t1 - t3
       fp2sqr_mont(t1, t1);          // t1 = t1^2
       fp2mul_mont(t2, t4, t2);      // t2 = t2*t4
       fp2shl(t2, 2, t2);            // t2 = 4*t2
       fp2sub(t2, t3, t2);           // t2 = t2 - t3
       fp2sqr_mont(t2, t2);          // t2 = t2^2
       fp2mul_mont(P->X, t2, Q->X);  // x = x*t2
       fp2mul_mont(P->Z, t1, Q->Z);  // z = z*t1    
}


void xTPLe_fast(point_proj_t P, point_proj_t Q, const f2elm_t A2, int e)
{ // Computes [3^e](X:Z) on Montgomery curve with projective constant via e repeated triplings. e triplings in E costs k*(5M + 6S + 9A)
  // Input: projective Montgomery x-coordinates P = (X:Z), where x=X/Z, Montgomery curve constant A2 = A/2 and the number of triplings e.
  // Output: projective Montgomery x-coordinates Q <- [3^e]P.    
    point_proj_t T;

    copy_words((digit_t*)P, (digit_t*)T, 2*2*NWORDS_FIELD);
    for (int j = 0; j < e; j++) { 
        xTPL_fast(T, T, A2);
    }
    copy_words((digit_t*)T, (digit_t*)Q, 2*2*NWORDS_FIELD);
}


void ADD(const point_full_proj_t P, const f2elm_t QX, const f2elm_t QY, const f2elm_t QZ, const f2elm_t A, point_full_proj_t R)
{ // General addition.
  // Input: projective Montgomery points P=(XP:YP:ZP) and Q=(XQ:YQ:ZQ).
  // Output: projective Montgomery point R <- P+Q = (XQP:YQP:ZQP). 
    f2elm_t t0 = {0}, t1 = {0}, t2 = {0}, t3 = {0}, t4 = {0}, t5 = {0}, t6 = {0}, t7 = {0};

    fp2mul_mont(QX, P->Z, t0);            // t0 = x2*Z1    
    fp2mul_mont(P->X, QZ, t1);            // t1 = X1*z2    
    fp2add(t0, t1, t2);                   // t2 = t0 + t1
    fp2sub(t1, t0, t3);                   // t3 = t1 - t0
    fp2mul_mont(QX, P->X, t0);            // t0 = x2*X1    
    fp2mul_mont(P->Z, QZ, t1);            // t1 = Z1*z2
    fp2add(t0, t1, t4);                   // t4 = t0 + t1
    fp2mul_mont(t0, A, t0);               // t0 = t0*A
    fp2mul_mont(QY, P->Y, t5);            // t5 = y2*Y1
    fp2sub(t0, t5, t0);                   // t0 = t0 - t5
    fp2mul_mont(t0, t1, t0);              // t0 = t0*t1
    fp2add(t0, t0, t0);                   // t0 = t0 + t0
    fp2mul_mont(t2, t4, t5);              // t5 = t2*t4
    fp2add(t5, t0, t5);                   // t5 = t5 + t0
    fp2sqr_mont(P->X, t0);                // t0 = X1 ^ 2
    fp2sqr_mont(P->Z, t6);                // t6 = Z1 ^ 2
    fp2add(t0, t6, t0);                   // t0 = t0 + t6    
    fp2add(t1, t1, t1);                   // t1 = t1 + t1
    fp2mul_mont(QY, P->X, t7);            // t7 = y2*X1
    fp2mul_mont(QX, P->Y, t6);            // t6 = x2*Y1
    fp2sub(t7, t6, t7);                   // t7 = t7 - t6
    fp2mul_mont(t1, t7, t1);              // t1 = t1*t7
    fp2mul_mont(A, t2, t7);               // t7 = A*t2
    fp2add(t7, t4, t4);                   // t4 = t4 + t7
    fp2mul_mont(t1, t4, t4);              // t4 = t1*t4
    fp2mul_mont(QY, QZ, t1);              // t1 = y2*z2
    fp2mul_mont(t0, t1, t0);              // t0 = t0*t1
    fp2sqr_mont(QZ, t1);                  // t1 = z2 ^ 2
    fp2sqr_mont(QX, t6);                  // t6 = x2 ^ 2
    fp2add(t1, t6, t1);                   // t1 = t1 + t6
    fp2mul_mont(P->Z, P->Y, t6);          // t6 = Z1*Y1
    fp2mul_mont(t1, t6, t1);              // t1 = t1*t6
    fp2sub(t0, t1, t0);                   // t0 = t0 - t1
    fp2mul_mont(t2, t0, t0);              // t0 = t2*t0
    fp2mul_mont(t5, t3, R->X);            // X3 = t5*t3    
    fp2add(t4, t0, R->Y);                 // Y3 = t4 + t0
    fp2sqr_mont(t3, t0);                  // t0 = t3 ^ 2
    fp2mul_mont(t3, t0, R->Z);            // Z3 = t3*t0
}


void Mont_ladder(const f2elm_t x, const digit_t* m, point_proj_t P, point_proj_t Q, const f2elm_t A24, const unsigned int order_bits, const unsigned int order_fullbits)
{ // The Montgomery ladder
  // Inputs: the affine x-coordinate of a point P on E: B*y^2=x^3+A*x^2+x, 
  //         scalar m
  //         curve constant A24 = (A+2)/4
  //         order_bits = subgroup order bitlength
  //         order_fullbits = smallest multiple of 32 larger than the order bitlength
  // Output: P = m*(x:1)
    unsigned int bit = 0, owords = NBITS_TO_NWORDS(order_fullbits);
    digit_t mask, scalar[NWORDS_ORDER];
    int i;
    
    // Initializing with the points (1:0) and (x:1)
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)P->X[0]);
    fpzero(P->X[1]);
    fp2zero(P->Z);
    
    fp2copy(x, Q->X);    
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)Q->Z[0]);    
    fpzero(Q->Z[1]);

    for (i = NWORDS_ORDER-1; i >= 0; i--) {
        scalar[i] = m[i];
    }
    
    for (i = order_fullbits-order_bits; i > 0; i--) {
        mp_shiftl1(scalar, owords);
    }    
    
    for (i = order_bits; i > 0; i--) {
        bit = (unsigned int)(scalar[owords-1] >> (RADIX-1));
        mp_shiftl1(scalar, owords);
        mask = 0-(digit_t)bit;

        swap_points(P, Q, mask);        
        xDBLADD(P, Q, x, A24);                     // If bit=0 then P <- 2*P and Q <- P+Q, 
        swap_points(P, Q, mask);                   // else if bit=1 then Q <- 2*Q and P <- P+Q
    }
}


void mont_twodim_scalarmult(digit_t* a, const point_t R, const point_t S, const f2elm_t A, const f2elm_t A24, point_full_proj_t P, const unsigned int order_bits)
{ // Computes P = R + [a]S  
    point_proj_t P0 = {0}, P1 = {0};
    point_full_proj_t P2 = {0};
    f2elm_t one = {0};    

    fpcopy((digit_t*)&Montgomery_one, one[0]);
    Mont_ladder(S->x, a, P0, P1, A24, order_bits, MAXBITS_ORDER);    
    recover_os(P0->X, P0->Z, P1->X, P1->Z, S->x, S->y, A, P2->X, P2->Y, P2->Z);     
    ADD(P2, R->x, R->y, one, A, P);       
}


void xDBLADD_proj(point_proj_t P, point_proj_t Q, const f2elm_t XPQ, const f2elm_t ZPQ, const f2elm_t A24)
{ // Simultaneous doubling and differential addition.
  // Input: projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, affine difference xPQ=x(P-Q) and Montgomery curve constant A24=(A+2)/4.
  // Output: projective Montgomery points P <- 2*P = (X2P:Z2P) such that x(2P)=X2P/Z2P, and Q <- P+Q = (XQP:ZQP) such that = x(Q+P)=XQP/ZQP. 
    f2elm_t t0, t1, t2;

    fp2add(P->X, P->Z, t0);                         // t0 = XP+ZP
    fp2sub(P->X, P->Z, t1);                         // t1 = XP-ZP    
    fp2sqr_mont(t0, P->X);                          // XP = (XP+ZP)^2    
    fp2sub(Q->X, Q->Z, t2);                         // t2 = XQ-ZQ
    fp2correction(t2);    
    fp2add(Q->X, Q->Z, Q->X);                       // XQ = XQ+ZQ    
    fp2mul_mont(t0, t2, t0);                        // t0 = (XP+ZP)*(XQ-ZQ)    
    fp2sqr_mont(t1, P->Z);                          // ZP = (XP-ZP)^2    
    fp2mul_mont(t1, Q->X, t1);                      // t1 = (XP-ZP)*(XQ+ZQ)    
    fp2sub(P->X, P->Z, t2);                         // t2 = (XP+ZP)^2-(XP-ZP)^2    
    fp2mul_mont(P->X, P->Z, P->X);                  // XP = (XP+ZP)^2*(XP-ZP)^2    
    fp2mul_mont(t2, A24, Q->X);                     // XQ = A24*[(XP+ZP)^2-(XP-ZP)^2]    
    fp2sub(t0, t1, Q->Z);                           // ZQ = (XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)    
    fp2add(Q->X, P->Z, P->Z);                       // ZP = A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2    
    fp2add(t0, t1, Q->X);                           // XQ = (XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)    
    fp2mul_mont(P->Z, t2, P->Z);                    // ZP = [A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2]*[(XP+ZP)^2-(XP-ZP)^2]    
    fp2sqr_mont(Q->Z, Q->Z);                        // ZQ = [(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2    
    fp2sqr_mont(Q->X, Q->X);                        // XQ = [(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2    
    fp2mul_mont(Q->X, ZPQ, Q->X);                   // XQ = ZPQ*[(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2    
    fp2mul_mont(Q->Z, XPQ, Q->Z);                   // ZQ = XPQ*[(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2          
}


void xDBL_e(const point_proj_t P, point_proj_t Q, const f2elm_t A24, const int e)
{ // Doubling of a Montgomery point in projective coordinates (X:Z) over affine curve coefficient A. 
  // Input: projective Montgomery x-coordinates P = (X1:Z1), where x1=X1/Z1 and Montgomery curve constants (A+2)/4.
  // Output: projective Montgomery x-coordinates Q = 2*P = (X2:Z2). 
    f2elm_t temp, a, b, c, aa, bb;    
    
    fp2copy(P->X,Q->X);
    fp2copy(P->Z,Q->Z);
    
    for (int j = 0; j < e; j++) {
        fp2add(Q->X, Q->Z, a);           // a = xQ + zQ
        fp2sub(Q->X, Q->Z, b);           // b = xQ - zQ
        fp2sqr_mont(a, aa);              //aa = (xQ + zQ)^2
        fp2sqr_mont(b, bb);              //bb = (xQ - zQ)^2
        fp2sub(aa, bb, c);               // c = (xQ + zQ)^2 - (xQ - zQ)^2
        fp2mul_mont(aa, bb, Q->X);       // xQ = (xQ + zQ)^2 * (xQ - zQ)^2
        fp2mul_mont(A24, c, temp);       // temp = A24 * ((xQ + zQ)^2 - (xQ - zQ)^2)
        fp2add(temp, bb, temp);          // temp = A24 * ((xQ + zQ)^2 - (xQ - zQ)^2) + (xQ - zQ)^2
        fp2mul_mont(c, temp, Q->Z);      // temp =  (A24 * ((xQ + zQ)^2 - (xQ - zQ)^2) + (xQ - zQ)^2) * ((xQ + zQ)^2 - (xQ - zQ)^2)
    }
}


void Ladder(const point_proj_t P, const digit_t* m, const f2elm_t A, const unsigned int order_bits, point_proj_t R) 
{
    point_proj_t R0, R1;
    f2elm_t A24 = {0};
    unsigned int bit = 0;
    digit_t mask;
    int j, swap, prevbit = 0;    
        
    fpcopy((digit_t*)&Montgomery_one, A24[0]);
    fpadd(A24[0], A24[0], A24[0]);
    fp2add(A, A24, A24);
    fp2div2(A24, A24);  
    fp2div2(A24, A24);  // A24 = (A+2)/4          

    j = order_bits - 1;
    bit = (m[j >> LOG2RADIX] >> (j & (RADIX-1))) & 1;
    while (bit == 0) {
        j--;
        bit = (m[j >> LOG2RADIX] >> (j & (RADIX-1))) & 1;
    }

    // R0 <- P, R1 <- 2P
    fp2copy(P->X, R0->X);
    fp2copy(P->Z, R0->Z);
    xDBL_e(P, R1, A24, 1);    
    
    // Main loop
    for (int i = j - 1;  i >= 0; i--) {
        bit = (m[i >> LOG2RADIX] >> (i & (RADIX-1))) & 1;
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(R0, R1, mask);
        xDBLADD_proj(R0, R1, P->X, P->Z, A24);
    }
    swap = 0 ^ prevbit;
    mask = 0 - (digit_t)swap;
    swap_points(R0, R1, mask);    
    
    fp2copy(R0->X, R->X);
    fp2copy(R0->Z, R->Z);
}

#endif


/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: ephemeral supersingular isogeny Diffie-Hellman key exchange (SIDH)
*********************************************************************************************/ 



static void init_basis(digit_t *gen, f2elm_t XP, f2elm_t XQ, f2elm_t XR)
{ // Initialization of basis points
    
    fpcopy(gen,                  XP[0]);
    fpcopy(gen +   NWORDS_FIELD, XP[1]);
    fpcopy(gen + 2*NWORDS_FIELD, XQ[0]);
    fpcopy(gen + 3*NWORDS_FIELD, XQ[1]);
    fpcopy(gen + 4*NWORDS_FIELD, XR[0]);
    fpcopy(gen + 5*NWORDS_FIELD, XR[1]);
}


void random_mod_order_A(unsigned char* random_digits)
{  // Generation of Alice's secret key  
   // Outputs random value in [0, 2^eA - 1]

    randombytes(random_digits, SECRETKEY_A_BYTES);
    random_digits[SECRETKEY_A_BYTES-1] &= MASK_ALICE;    // Masking last byte 
}


void random_mod_order_B(unsigned char* random_digits)
{  // Generation of Bob's secret key  
   // Outputs random value in [0, 2^Floor(Log(2, oB)) - 1]

    randombytes(random_digits, SECRETKEY_B_BYTES);
    random_digits[SECRETKEY_B_BYTES-1] &= MASK_BOB;     // Masking last byte 
}


int EphemeralKeyGeneration_A(const unsigned char* PrivateKeyA, unsigned char* PublicKeyA)
{ // Alice's ephemeral public key generation
  // Input:  a private key PrivateKeyA in the range [0, 2^eA - 1]. 
  // Output: the public key PublicKeyA consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R, phiP = {0}, phiQ = {0}, phiR = {0}, pts[MAX_INT_POINTS_ALICE];
    f2elm_t XPA, XQA, XRA, coeff[3], A24plus = {0}, C24 = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;
    digit_t SecretKeyA[NWORDS_ORDER] = {0};

    // Initialize basis points
    init_basis((digit_t*)A_gen, XPA, XQA, XRA);
    init_basis((digit_t*)B_gen, phiP->X, phiQ->X, phiR->X);
    fpcopy((digit_t*)&Montgomery_one, (phiP->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiQ->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiR->Z)[0]);

    // Initialize constants: A24plus = A+2C, C24 = 4C, where A=6, C=1
    fpcopy((digit_t*)&Montgomery_one, A24plus[0]);
    mp2_add(A24plus, A24plus, A24plus);
    mp2_add(A24plus, A24plus, C24);
    mp2_add(A24plus, C24, A);
    mp2_add(C24, C24, A24plus);

    // Retrieve kernel point
    decode_to_digits(PrivateKeyA, SecretKeyA, SECRETKEY_A_BYTES, NWORDS_ORDER);
    LADDER3PT(XPA, XQA, XRA, SecretKeyA, ALICE, R, A);       

#if (OALICE_BITS % 2 == 1)
    point_proj_t S;

    xDBLe(R, S, A24plus, C24, (int)(OALICE_BITS-1));
    get_2_isog(S, A24plus, C24); 
    eval_2_isog(phiP, S); 
    eval_2_isog(phiQ, S); 
    eval_2_isog(phiR, S);
    eval_2_isog(R, S);
#endif

    // Traverse tree
    index = 0;        
    for (row = 1; row < MAX_Alice; row++) {
        while (index < MAX_Alice-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Alice[ii++];
            xDBLe(R, R, A24plus, C24, (int)(2*m));
            index += m;
        }
        get_4_isog(R, A24plus, C24, coeff);        

        for (i = 0; i < npts; i++) {
            eval_4_isog(pts[i], coeff);
        }
        eval_4_isog(phiP, coeff);
        eval_4_isog(phiQ, coeff);
        eval_4_isog(phiR, coeff);

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }

    get_4_isog(R, A24plus, C24, coeff); 
    eval_4_isog(phiP, coeff);
    eval_4_isog(phiQ, coeff);
    eval_4_isog(phiR, coeff);

    inv_3_way(phiP->Z, phiQ->Z, phiR->Z);
    fp2mul_mont(phiP->X, phiP->Z, phiP->X);
    fp2mul_mont(phiQ->X, phiQ->Z, phiQ->X);
    fp2mul_mont(phiR->X, phiR->Z, phiR->X);
                
    // Format public key                   
    fp2_encode(phiP->X, PublicKeyA);
    fp2_encode(phiQ->X, PublicKeyA + FP2_ENCODED_BYTES);
    fp2_encode(phiR->X, PublicKeyA + 2*FP2_ENCODED_BYTES);

    return 0;
}


int EphemeralKeyGeneration_B(const unsigned char* PrivateKeyB, unsigned char* PublicKeyB)
{ // Bob's ephemeral public key generation
  // Input:  a private key PrivateKeyB in the range [0, 2^Floor(Log(2,oB)) - 1]. 
  // Output: the public key PublicKeyB consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R, phiP = {0}, phiQ = {0}, phiR = {0}, pts[MAX_INT_POINTS_BOB];
    f2elm_t XPB, XQB, XRB, coeff[3], A24plus = {0}, A24minus = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;
    digit_t SecretKeyB[NWORDS_ORDER] = {0};

    // Initialize basis points
    init_basis((digit_t*)B_gen, XPB, XQB, XRB);
    init_basis((digit_t*)A_gen, phiP->X, phiQ->X, phiR->X);
    fpcopy((digit_t*)&Montgomery_one, (phiP->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiQ->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiR->Z)[0]);

    // Initialize constants: A24minus = A-2C, A24plus = A+2C, where A=6, C=1
    fpcopy((digit_t*)&Montgomery_one, A24plus[0]);
    mp2_add(A24plus, A24plus, A24plus);
    mp2_add(A24plus, A24plus, A24minus);
    mp2_add(A24plus, A24minus, A);
    mp2_add(A24minus, A24minus, A24plus);

    // Retrieve kernel point
    decode_to_digits(PrivateKeyB, SecretKeyB, SECRETKEY_B_BYTES, NWORDS_ORDER);
    LADDER3PT(XPB, XQB, XRB, SecretKeyB, BOB, R, A);
    
    // Traverse tree
    index = 0;  
    for (row = 1; row < MAX_Bob; row++) {
        while (index < MAX_Bob-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Bob[ii++];
            xTPLe(R, R, A24minus, A24plus, (int)m);
            index += m;
        } 
        get_3_isog(R, A24minus, A24plus, coeff);

        for (i = 0; i < npts; i++) {
            eval_3_isog(pts[i], coeff);
        }     
        eval_3_isog(phiP, coeff);
        eval_3_isog(phiQ, coeff);
        eval_3_isog(phiR, coeff);

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }
    
    get_3_isog(R, A24minus, A24plus, coeff);
    eval_3_isog(phiP, coeff);
    eval_3_isog(phiQ, coeff);
    eval_3_isog(phiR, coeff);

    inv_3_way(phiP->Z, phiQ->Z, phiR->Z);
    fp2mul_mont(phiP->X, phiP->Z, phiP->X);
    fp2mul_mont(phiQ->X, phiQ->Z, phiQ->X);
    fp2mul_mont(phiR->X, phiR->Z, phiR->X);

    // Format public key
    fp2_encode(phiP->X, PublicKeyB);
    fp2_encode(phiQ->X, PublicKeyB + FP2_ENCODED_BYTES);
    fp2_encode(phiR->X, PublicKeyB + 2*FP2_ENCODED_BYTES);

    return 0;
}


int EphemeralSecretAgreement_A(const unsigned char* PrivateKeyA, const unsigned char* PublicKeyB, unsigned char* SharedSecretA)
{ // Alice's ephemeral shared secret computation
  // It produces a shared secret key SharedSecretA using her secret key PrivateKeyA and Bob's public key PublicKeyB
  // Inputs: Alice's PrivateKeyA is an integer in the range [0, oA-1]. 
  //         Bob's PublicKeyB consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
  // Output: a shared secret SharedSecretA that consists of one element in GF(p^2) encoded by removing leading 0 bytes.  
    point_proj_t R, pts[MAX_INT_POINTS_ALICE];
    f2elm_t coeff[3], PKB[3], jinv;
    f2elm_t A24plus = {0}, C24 = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;
    digit_t SecretKeyA[NWORDS_ORDER] = {0};
      
    // Initialize images of Bob's basis
    fp2_decode(PublicKeyB, PKB[0]);
    fp2_decode(PublicKeyB + FP2_ENCODED_BYTES, PKB[1]);
    fp2_decode(PublicKeyB + 2*FP2_ENCODED_BYTES, PKB[2]);

    // Initialize constants: A24plus = A+2C, C24 = 4C, where C=1
    get_A(PKB[0], PKB[1], PKB[2], A);
    mp_add((digit_t*)&Montgomery_one, (digit_t*)&Montgomery_one, C24[0], NWORDS_FIELD);
    mp2_add(A, C24, A24plus);
    mp_add(C24[0], C24[0], C24[0], NWORDS_FIELD);

    // Retrieve kernel point
    decode_to_digits(PrivateKeyA, SecretKeyA, SECRETKEY_A_BYTES, NWORDS_ORDER);
    LADDER3PT(PKB[0], PKB[1], PKB[2], SecretKeyA, ALICE, R, A);    

#if (OALICE_BITS % 2 == 1)
    point_proj_t S;

    xDBLe(R, S, A24plus, C24, (int)(OALICE_BITS-1));
    get_2_isog(S, A24plus, C24);
    eval_2_isog(R, S);
#endif

    // Traverse tree
    index = 0;        
    for (row = 1; row < MAX_Alice; row++) {
        while (index < MAX_Alice-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Alice[ii++];
            xDBLe(R, R, A24plus, C24, (int)(2*m));
            index += m;
        }
        get_4_isog(R, A24plus, C24, coeff);        

        for (i = 0; i < npts; i++) {
            eval_4_isog(pts[i], coeff);
        }

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }

    get_4_isog(R, A24plus, C24, coeff); 
    mp2_add(A24plus, A24plus, A24plus);                                                
    fp2sub(A24plus, C24, A24plus); 
    fp2add(A24plus, A24plus, A24plus);                    
    j_inv(A24plus, C24, jinv);
    fp2_encode(jinv, SharedSecretA);    // Format shared secret

    return 0;
}


int EphemeralSecretAgreement_B(const unsigned char* PrivateKeyB, const unsigned char* PublicKeyA, unsigned char* SharedSecretB)
{ // Bob's ephemeral shared secret computation
  // It produces a shared secret key SharedSecretB using his secret key PrivateKeyB and Alice's public key PublicKeyA
  // Inputs: Bob's PrivateKeyB is an integer in the range [0, 2^Floor(Log(2,oB)) - 1]. 
  //         Alice's PublicKeyA consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
  // Output: a shared secret SharedSecretB that consists of one element in GF(p^2) encoded by removing leading 0 bytes.  
    point_proj_t R, pts[MAX_INT_POINTS_BOB];
    f2elm_t coeff[3], PKB[3], jinv;
    f2elm_t A24plus = {0}, A24minus = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;
    digit_t SecretKeyB[NWORDS_ORDER] = {0};
      
    // Initialize images of Alice's basis
    fp2_decode(PublicKeyA, PKB[0]);
    fp2_decode(PublicKeyA + FP2_ENCODED_BYTES, PKB[1]);
    fp2_decode(PublicKeyA + 2*FP2_ENCODED_BYTES, PKB[2]);

    // Initialize constants: A24plus = A+2C, A24minus = A-2C, where C=1
    get_A(PKB[0], PKB[1], PKB[2], A);
    mp_add((digit_t*)&Montgomery_one, (digit_t*)&Montgomery_one, A24minus[0], NWORDS_FIELD);
    mp2_add(A, A24minus, A24plus);
    mp2_sub_p2(A, A24minus, A24minus);

    // Retrieve kernel point
    decode_to_digits(PrivateKeyB, SecretKeyB, SECRETKEY_B_BYTES, NWORDS_ORDER);
    LADDER3PT(PKB[0], PKB[1], PKB[2], SecretKeyB, BOB, R, A);
    
    // Traverse tree
    index = 0;  
    for (row = 1; row < MAX_Bob; row++) {
        while (index < MAX_Bob-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Bob[ii++];
            xTPLe(R, R, A24minus, A24plus, (int)m);
            index += m;
        }
        get_3_isog(R, A24minus, A24plus, coeff);

        for (i = 0; i < npts; i++) {
            eval_3_isog(pts[i], coeff);
        } 

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }
     
    get_3_isog(R, A24minus, A24plus, coeff);    
    fp2add(A24plus, A24minus, A);                 
    fp2add(A, A, A);
    fp2sub(A24plus, A24minus, A24plus);                   
    j_inv(A, A24plus, jinv);
    fp2_encode(jinv, SharedSecretB);    // Format shared secret

    return 0;
}


/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: supersingular isogeny key encapsulation (SIKE) protocol
*********************************************************************************************/ 



int crypto_kem_keypair(unsigned char *pk, unsigned char *sk)
{ // SIKE's key generation
  // Outputs: secret key sk (CRYPTO_SECRETKEYBYTES = MSG_BYTES + SECRETKEY_B_BYTES + CRYPTO_PUBLICKEYBYTES bytes)
  //          public key pk (CRYPTO_PUBLICKEYBYTES bytes) 

    // Generate lower portion of secret key sk <- s||SK
    randombytes(sk, MSG_BYTES);
    random_mod_order_B(sk + MSG_BYTES);

    // Generate public key pk
    EphemeralKeyGeneration_B(sk + MSG_BYTES, pk);

    // Append public key pk to secret key sk
    memcpy(&sk[MSG_BYTES + SECRETKEY_B_BYTES], pk, CRYPTO_PUBLICKEYBYTES);

    return 0;
}


int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk)
{ // SIKE's encapsulation
  // Input:   public key pk         (CRYPTO_PUBLICKEYBYTES bytes)
  // Outputs: shared secret ss      (CRYPTO_BYTES bytes)
  //          ciphertext message ct (CRYPTO_CIPHERTEXTBYTES = CRYPTO_PUBLICKEYBYTES + MSG_BYTES bytes)
    unsigned char ephemeralsk[SECRETKEY_A_BYTES];
    unsigned char jinvariant[FP2_ENCODED_BYTES];
    unsigned char h[MSG_BYTES];
    unsigned char temp[CRYPTO_CIPHERTEXTBYTES+MSG_BYTES];

    // Generate ephemeralsk <- G(m||pk) mod oA 
    randombytes(temp, MSG_BYTES);
    memcpy(&temp[MSG_BYTES], pk, CRYPTO_PUBLICKEYBYTES);
    shake256(ephemeralsk, SECRETKEY_A_BYTES, temp, CRYPTO_PUBLICKEYBYTES+MSG_BYTES);
    ephemeralsk[SECRETKEY_A_BYTES - 1] &= MASK_ALICE;

    // Encrypt
    EphemeralKeyGeneration_A(ephemeralsk, ct);
    EphemeralSecretAgreement_A(ephemeralsk, pk, jinvariant);
    shake256(h, MSG_BYTES, jinvariant, FP2_ENCODED_BYTES);
    for (int i = 0; i < MSG_BYTES; i++) {
        ct[i + CRYPTO_PUBLICKEYBYTES] = temp[i] ^ h[i];
    }

    // Generate shared secret ss <- H(m||ct)
    memcpy(&temp[MSG_BYTES], ct, CRYPTO_CIPHERTEXTBYTES);
    shake256(ss, CRYPTO_BYTES, temp, CRYPTO_CIPHERTEXTBYTES+MSG_BYTES);

    return 0;
}


int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk)
{ // SIKE's decapsulation
  // Input:   secret key sk         (CRYPTO_SECRETKEYBYTES = MSG_BYTES + SECRETKEY_B_BYTES + CRYPTO_PUBLICKEYBYTES bytes)
  //          ciphertext message ct (CRYPTO_CIPHERTEXTBYTES = CRYPTO_PUBLICKEYBYTES + MSG_BYTES bytes) 
  // Outputs: shared secret ss      (CRYPTO_BYTES bytes)
    unsigned char ephemeralsk_[SECRETKEY_A_BYTES];
    unsigned char jinvariant_[FP2_ENCODED_BYTES];
    unsigned char h_[MSG_BYTES];
    unsigned char c0_[CRYPTO_PUBLICKEYBYTES];
    unsigned char temp[CRYPTO_CIPHERTEXTBYTES+MSG_BYTES];

    // Decrypt
    EphemeralSecretAgreement_B(sk + MSG_BYTES, ct, jinvariant_);
    shake256(h_, MSG_BYTES, jinvariant_, FP2_ENCODED_BYTES);
    for (int i = 0; i < MSG_BYTES; i++) {
        temp[i] = ct[i + CRYPTO_PUBLICKEYBYTES] ^ h_[i];
    }

    // Generate ephemeralsk_ <- G(m||pk) mod oA
    memcpy(&temp[MSG_BYTES], &sk[MSG_BYTES + SECRETKEY_B_BYTES], CRYPTO_PUBLICKEYBYTES);
    shake256(ephemeralsk_, SECRETKEY_A_BYTES, temp, CRYPTO_PUBLICKEYBYTES+MSG_BYTES);
    ephemeralsk_[SECRETKEY_A_BYTES - 1] &= MASK_ALICE;
    
    // Generate shared secret ss <- H(m||ct), or output ss <- H(s||ct) in case of ct verification failure
    EphemeralKeyGeneration_A(ephemeralsk_, c0_);
    // If selector = 0 then do ss = H(m||ct), else if selector = -1 load s to do ss = H(s||ct)
    int8_t selector = ct_compare(c0_, ct, CRYPTO_PUBLICKEYBYTES);
    ct_cmov(temp, sk, MSG_BYTES, selector);
    memcpy(&temp[MSG_BYTES], ct, CRYPTO_CIPHERTEXTBYTES);
    shake256(ss, CRYPTO_BYTES, temp, CRYPTO_CIPHERTEXTBYTES+MSG_BYTES);

    return 0;
}
