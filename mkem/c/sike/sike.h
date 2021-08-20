/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: API header file for P434
*********************************************************************************************/  

#ifndef P434_API_H
#define P434_API_H


/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: configuration file and platform-dependent macros
*********************************************************************************************/  


#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>


// Definition of operating system

#define OS_WIN       1
#define OS_NIX       2

#if defined(__WINDOWS__)        // Microsoft Windows OS
    #define OS_TARGET OS_WIN
#elif defined(__NIX__) || defined(__APPLE__)          // Unix-like operative systems
    #define OS_TARGET OS_NIX
#else
    #error -- "Unsupported OS"
#endif




// Definition of compiler

#define COMPILER_VC      1
#define COMPILER_GCC     2
#define COMPILER_CLANG   3

#if defined(_MSC_VER)           // Microsoft Visual C compiler
    #define COMPILER COMPILER_VC
#elif defined(__GNUC__)         // GNU GCC compiler
    #define COMPILER COMPILER_GCC   
#elif defined(__clang__)        // Clang compiler
    #define COMPILER COMPILER_CLANG
#else
    #error -- "Unsupported COMPILER"
#endif


// Definition of the targeted architecture and basic data types
    
#define TARGET_AMD64        1
#define TARGET_x86          2
#define TARGET_S390X        3
#define TARGET_ARM          4
#define TARGET_ARM64        5

#if defined(_AMD64_)
    #define TARGET TARGET_AMD64
    #define RADIX           64
    #define LOG2RADIX       6  
    typedef uint64_t        digit_t;        // Unsigned 64-bit digit
    typedef uint32_t        hdigit_t;       // Unsigned 32-bit digit
#elif defined(_X86_)
    #define TARGET TARGET_x86
    #define RADIX           32
    #define LOG2RADIX       5  
    typedef uint32_t        digit_t;        // Unsigned 32-bit digit
    typedef uint16_t        hdigit_t;       // Unsigned 16-bit digit  
#elif defined(_S390X_)
    #define TARGET TARGET_S390X
    #define RADIX           64
    #define LOG2RADIX       6  
    typedef uint64_t        digit_t;        // Unsigned 64-bit digit
    typedef uint32_t        hdigit_t;       // Unsigned 32-bit digit
#elif defined(_ARM_)
    #define TARGET TARGET_ARM
    #define RADIX           32
    #define LOG2RADIX       5  
    typedef uint32_t        digit_t;        // Unsigned 32-bit digit
    typedef uint16_t        hdigit_t;       // Unsigned 16-bit digit  
#elif defined(_ARM64_) || defined(__aarch64__)
    #define TARGET TARGET_ARM64
    #define RADIX           64
    #define LOG2RADIX       6  
    typedef uint64_t        digit_t;        // Unsigned 64-bit digit
    typedef uint32_t        hdigit_t;       // Unsigned 32-bit digit
#else
    #error -- "Unsupported ARCHITECTURE"
#endif

#define RADIX64             64


// Selection of generic, portable implementation

#define GENERIC_IMPLEMENTATION


// Extended datatype support

typedef uint64_t uint128_t[2];
    

// Macro definitions

#define NBITS_TO_NBYTES(nbits)      (((nbits)+7)/8)                                          // Conversion macro from number of bits to number of bytes
#define NBITS_TO_NWORDS(nbits)      (((nbits)+(sizeof(digit_t)*8)-1)/(sizeof(digit_t)*8))    // Conversion macro from number of bits to number of computer words
#define NBYTES_TO_NWORDS(nbytes)    (((nbytes)+sizeof(digit_t)-1)/sizeof(digit_t))           // Conversion macro from number of bytes to number of computer words

// Macro to avoid compiler warnings when detecting unreferenced parameters
#define UNREFERENCED_PARAMETER(PAR) ((void)(PAR))


// Macros for endianness
// 32-bit byte swap
#if (COMPILER == COMPILER_GCC || COMPILER == COMPILER_CLANG)
    #define BSWAP32(i) __builtin_bswap32((i))
#else
    #define BSWAP32(i) ((((i) >> 24) & 0xff) | (((i) >> 8) & 0xff00) | (((i) & 0xff00) << 8) | ((i) << 24))
#endif

// 64-bit byte swap
#if (COMPILER == COMPILER_GCC || COMPILER == COMPILER_CLANG)
    #define BSWAP64(i) __builtin_bswap64((i))
#else
    #define BSWAP64(i) ((BSWAP32((i) >> 32) & 0xffffffff) | (BSWAP32(i) << 32))
#endif

#if RADIX == 32
    #define BSWAP_DIGIT(i) BSWAP32((i))
#elif RADIX == 64
    #define BSWAP_DIGIT(i) BSWAP64((i))
#endif

// Host to little endian, little endian to host
#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
    #define _BIG_ENDIAN_
    #define HTOLE_64(i) BSWAP64((i))
    #define LETOH_64(i) BSWAP64((i))
#else
    #define _LITTLE_ENDIAN_
    #define HTOLE_64(i) (i)
    #define LETOH_64(i) (i)
#endif


/********************** Constant-time unsigned comparisons ***********************/

// The following functions return 1 (TRUE) if condition is true, 0 (FALSE) otherwise

static __inline unsigned int is_digit_nonzero_ct(digit_t x)
{ // Is x != 0?
    return (unsigned int)((x | (0-x)) >> (RADIX-1));
}

static __inline unsigned int is_digit_zero_ct(digit_t x)
{ // Is x = 0?
    return (unsigned int)(1 ^ is_digit_nonzero_ct(x));
}

static __inline unsigned int is_digit_lessthan_ct(digit_t x, digit_t y)
{ // Is x < y?
    return (unsigned int)((x ^ ((x ^ y) | ((x - y) ^ y))) >> (RADIX-1)); 
}


/********************** Macros for platform-dependent operations **********************/

#if defined(GENERIC_IMPLEMENTATION) || (TARGET == TARGET_ARM)

// Digit multiplication
#define MUL(multiplier, multiplicand, hi, lo)                                                     \
    digit_x_digit((multiplier), (multiplicand), &(lo));
    
// Digit addition with carry
#define ADDC(carryIn, addend1, addend2, carryOut, sumOut)                                         \
    { digit_t tempReg = (addend1) + (digit_t)(carryIn);                                           \
    (sumOut) = (addend2) + tempReg;                                                               \
    (carryOut) = (is_digit_lessthan_ct(tempReg, (digit_t)(carryIn)) | is_digit_lessthan_ct((sumOut), tempReg)); }

// Digit subtraction with borrow
#define SUBC(borrowIn, minuend, subtrahend, borrowOut, differenceOut)                             \
    { digit_t tempReg = (minuend) - (subtrahend);                                                 \
    unsigned int borrowReg = (is_digit_lessthan_ct((minuend), (subtrahend)) | ((borrowIn) & is_digit_zero_ct(tempReg)));  \
    (differenceOut) = tempReg - (digit_t)(borrowIn);                                              \
    (borrowOut) = borrowReg; }
    
// Shift right with flexible datatype
#define SHIFTR(highIn, lowIn, shift, shiftOut, DigitSize)                                         \
    (shiftOut) = ((lowIn) >> (shift)) ^ ((highIn) << (DigitSize - (shift)));
    
// Shift left with flexible datatype
#define SHIFTL(highIn, lowIn, shift, shiftOut, DigitSize)                                         \
    (shiftOut) = ((highIn) << (shift)) ^ ((lowIn) >> (DigitSize - (shift)));

#elif (TARGET == TARGET_AMD64 && OS_TARGET == OS_WIN)

// Digit multiplication
#define MUL(multiplier, multiplicand, hi, lo)                                                     \
    (lo) = _umul128((multiplier), (multiplicand), (hi));                

// Digit addition with carry
#define ADDC(carryIn, addend1, addend2, carryOut, sumOut)                                         \
    (carryOut) = _addcarry_u64((carryIn), (addend1), (addend2), &(sumOut));

// Digit subtraction with borrow
#define SUBC(borrowIn, minuend, subtrahend, borrowOut, differenceOut)                             \
    (borrowOut) = _subborrow_u64((borrowIn), (minuend), (subtrahend), &(differenceOut));

// Digit shift right
#define SHIFTR(highIn, lowIn, shift, shiftOut, DigitSize)                                         \
    (shiftOut) = __shiftright128((lowIn), (highIn), (shift));

// Digit shift left
#define SHIFTL(highIn, lowIn, shift, shiftOut, DigitSize)                                         \
    (shiftOut) = __shiftleft128((lowIn), (highIn), (shift));

// 64x64-bit multiplication
#define MUL128(multiplier, multiplicand, product)                                                 \
    (product)[0] = _umul128((multiplier), (multiplicand), &(product)[1]);

// 128-bit addition with output carry
#define ADC128(addend1, addend2, carry, addition)                                                 \
    (carry) = _addcarry_u64(0, (addend1)[0], (addend2)[0], &(addition)[0]);                       \
    (carry) = _addcarry_u64((carry), (addend1)[1], (addend2)[1], &(addition)[1]); 

#define MULADD128(multiplier, multiplicand, addend, carry, result);                               \
    { uint128_t product;                                                                          \
      MUL128(multiplier, multiplicand, product);                                                  \
      ADC128(addend, product, carry, result); }

#elif ((TARGET == TARGET_AMD64 || TARGET == TARGET_ARM64) && OS_TARGET == OS_NIX)

// Digit multiplication
#define MUL(multiplier, multiplicand, hi, lo)                                                     \
    { uint128_t tempReg = (uint128_t)(multiplier) * (uint128_t)(multiplicand);                    \
    *(hi) = (digit_t)(tempReg >> RADIX);                                                          \
    (lo) = (digit_t)tempReg; }

// Digit addition with carry
#define ADDC(carryIn, addend1, addend2, carryOut, sumOut)                                         \
    { uint128_t tempReg = (uint128_t)(addend1) + (uint128_t)(addend2) + (uint128_t)(carryIn);     \
    (carryOut) = (digit_t)(tempReg >> RADIX);                                                     \
    (sumOut) = (digit_t)tempReg; }  
    
// Digit subtraction with borrow
#define SUBC(borrowIn, minuend, subtrahend, borrowOut, differenceOut)                             \
    { uint128_t tempReg = (uint128_t)(minuend) - (uint128_t)(subtrahend) - (uint128_t)(borrowIn); \
    (borrowOut) = (digit_t)(tempReg >> (sizeof(uint128_t)*8 - 1));                                \
    (differenceOut) = (digit_t)tempReg; }

// Digit shift right
#define SHIFTR(highIn, lowIn, shift, shiftOut, DigitSize)                                         \
    (shiftOut) = ((lowIn) >> (shift)) ^ ((highIn) << (RADIX - (shift)));

// Digit shift left
#define SHIFTL(highIn, lowIn, shift, shiftOut, DigitSize)                                         \
    (shiftOut) = ((highIn) << (shift)) ^ ((lowIn) >> (RADIX - (shift)));

#endif



/*********************** Key encapsulation mechanism API ***********************/

#define CRYPTO_SECRETKEYBYTES     374    // MSG_BYTES + SECRETKEY_B_BYTES + CRYPTO_PUBLICKEYBYTES bytes
#define CRYPTO_PUBLICKEYBYTES     330
#define CRYPTO_BYTES               16
#define CRYPTO_CIPHERTEXTBYTES    346    // CRYPTO_PUBLICKEYBYTES + MSG_BYTES bytes  

// Algorithm name
#define CRYPTO_ALGNAME "SIKEp434"  

// SIKE's key generation
// It produces a private key sk and computes the public key pk.
// Outputs: secret key sk (CRYPTO_SECRETKEYBYTES = 374 bytes)
//          public key pk (CRYPTO_PUBLICKEYBYTES = 330 bytes) 
int crypto_kem_keypair(unsigned char *pk, unsigned char *sk);

// SIKE's encapsulation
// Input:   public key pk         (CRYPTO_PUBLICKEYBYTES = 330 bytes)
// Outputs: shared secret ss      (CRYPTO_BYTES = 16 bytes)
//          ciphertext message ct (CRYPTO_CIPHERTEXTBYTES = 346 bytes)
int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk);

// SIKE's decapsulation
// Input:   secret key sk         (CRYPTO_SECRETKEYBYTES = 374 bytes)
//          ciphertext message ct (CRYPTO_CIPHERTEXTBYTES = 346 bytes) 
// Outputs: shared secret ss      (CRYPTO_BYTES = 16 bytes)
int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk);


// Encoding of keys for KEM-based isogeny system "SIKEp434" (wire format):
// ----------------------------------------------------------------------
// Elements over GF(p434) are encoded in 55 octets in little endian format (i.e., the least significant octet is located in the lowest memory address). 
// Elements (a+b*i) over GF(p434^2), where a and b are defined over GF(p434), are encoded as {a, b}, with a in the lowest memory portion.
//
// Private keys sk consist of the concatenation of a 16-byte random value, a value in the range [0, 2^217-1] and the public key pk. In the SIKE API, 
// private keys are encoded in 374 octets in little endian format. 
// Public keys pk consist of 3 elements in GF(p434^2). In the SIKE API, pk is encoded in 330 octets. 
// Ciphertexts ct consist of the concatenation of a public key value and a 16-byte value. In the SIKE API, ct is encoded in 330 + 16 = 346 octets.  
// Shared keys ss consist of a value of 16 octets.


/*********************** Ephemeral key exchange API ***********************/

#define SIDH_SECRETKEYBYTES_A    27
#define SIDH_SECRETKEYBYTES_B    28
#define SIDH_PUBLICKEYBYTES     330
#define SIDH_BYTES              110

// SECURITY NOTE: SIDH supports ephemeral Diffie-Hellman key exchange. It is NOT secure to use it with static keys.
// See "On the Security of Supersingular Isogeny Cryptosystems", S.D. Galbraith, C. Petit, B. Shani and Y.B. Ti, in ASIACRYPT 2016, 2016.
// Extended version available at: http://eprint.iacr.org/2016/859  

// Generation of Alice's secret key 
// Outputs random value in [0, 2^216 - 1] to be used as Alice's private key
void random_mod_order_A(unsigned char* random_digits);

// Generation of Bob's secret key 
// Outputs random value in [0, 2^Floor(Log(2,3^137)) - 1] to be used as Bob's private key
void random_mod_order_B(unsigned char* random_digits);

// Alice's ephemeral public key generation
// Input:  a private key PrivateKeyA in the range [0, 2^216 - 1], stored in 27 bytes. 
// Output: the public key PublicKeyA consisting of 3 GF(p434^2) elements encoded in 330 bytes.
int EphemeralKeyGeneration_A(const unsigned char* PrivateKeyA, unsigned char* PublicKeyA);

// Bob's ephemeral key-pair generation
// It produces a private key PrivateKeyB and computes the public key PublicKeyB.
// The private key is an integer in the range [0, 2^Floor(Log(2,3^137)) - 1], stored in 28 bytes. 
// The public key consists of 3 GF(p434^2) elements encoded in 330 bytes.
int EphemeralKeyGeneration_B(const unsigned char* PrivateKeyB, unsigned char* PublicKeyB);

// Alice's ephemeral shared secret computation
// It produces a shared secret key SharedSecretA using her secret key PrivateKeyA and Bob's public key PublicKeyB
// Inputs: Alice's PrivateKeyA is an integer in the range [0, 2^216 - 1], stored in 27 bytes. 
//         Bob's PublicKeyB consists of 3 GF(p434^2) elements encoded in 330 bytes.
// Output: a shared secret SharedSecretA that consists of one element in GF(p434^2) encoded in 110 bytes.
int EphemeralSecretAgreement_A(const unsigned char* PrivateKeyA, const unsigned char* PublicKeyB, unsigned char* SharedSecretA);

// Bob's ephemeral shared secret computation
// It produces a shared secret key SharedSecretB using his secret key PrivateKeyB and Alice's public key PublicKeyA
// Inputs: Bob's PrivateKeyB is an integer in the range [0, 2^Floor(Log(2,3^137)) - 1], stored in 28 bytes. 
//         Alice's PublicKeyA consists of 3 GF(p434^2) elements encoded in 330 bytes.
// Output: a shared secret SharedSecretB that consists of one element in GF(p434^2) encoded in 110 bytes.
int EphemeralSecretAgreement_B(const unsigned char* PrivateKeyB, const unsigned char* PublicKeyA, unsigned char* SharedSecretB);


// Encoding of keys for KEX-based isogeny system "SIDHp434" (wire format):
// ----------------------------------------------------------------------
// Elements over GF(p434) are encoded in 55 octets in little endian format (i.e., the least significant octet is located in the lowest memory address). 
// Elements (a+b*i) over GF(p434^2), where a and b are defined over GF(p434), are encoded as {a, b}, with a in the lowest memory portion.
//
// Private keys PrivateKeyA and PrivateKeyB can have values in the range [0, 2^216-1] and [0, 2^Floor(Log(2,3^137)) - 1], resp. In the SIDH API, 
// Alice's and Bob's private keys are encoded in 27 and 28 octets, resp., in little endian format. 
// Public keys PublicKeyA and PublicKeyB consist of 3 elements in GF(p434^2). In the SIDH API, they are encoded in 330 octets. 
// Shared keys SharedSecretA and SharedSecretB consist of one element in GF(p434^2). In the SIDH API, they are encoded in 110 octets.


#endif
