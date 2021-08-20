#include <stdint.h>
#include "params.h"
#include "cbd.h"

/*************************************************
* Name:        load32_littleendian
*
* Description: load 4 bytes into a 32-bit integer
*              in little-endian order
*
* Arguments:   - const uint8_t *x: pointer to input byte array
*
* Returns 32-bit unsigned integer loaded from x
**************************************************/
static uint32_t load32_littleendian(const uint8_t x[4])
{
  uint32_t r;
  r  = (uint32_t)x[0];
  r |= (uint32_t)x[1] << 8;
  r |= (uint32_t)x[2] << 16;
  r |= (uint32_t)x[3] << 24;
  return r;
}

/*************************************************
* Name:        load24_littleendian
*
* Description: load 3 bytes into a 32-bit integer
*              in little-endian order.
*              This function is only needed for Kyber-512
*
* Arguments:   - const uint8_t *x: pointer to input byte array
*
* Returns 32-bit unsigned integer loaded from x (most significant byte is zero)
**************************************************/
static uint32_t load24_littleendian(const uint8_t x[3])
{
  uint32_t r;
  r  = (uint32_t)x[0];
  r |= (uint32_t)x[1] << 8;
  r |= (uint32_t)x[2] << 16;
  return r;
}

#if KYBER_ETA1 == 6 || KYBER_ETA2 == 6
// Load 5 bytes into a 64-bit integer little endian.
static uint64_t load48_littleendian(const uint8_t x[6])
{
  return ((uint64_t)x[0]) | ((uint64_t)x[1] << 8)  | ((uint64_t)x[2] << 16) |
    ((uint64_t)x[3] << 24) | ((uint64_t)x[4] << 32) | ((uint64_t)x[5] << 40);
}

// Given an array of uniformly random bytes, compute polynomial with
// coefficients distributed according to a centered binomial distribution
// with parameter eta=6.
static void cbd6(poly *r, const uint8_t buf[12*KYBER_N/8])
{
  // The distribution at hand is exactly the same as that
  // of (a₁ + a₂ + a₃ + a₄ + a₅ + a₆) - (b₁ + b₂ + b₃ + b₄ + b₅ + b₆)
  // with a_i,b_i~U(1).  Thus we need 12 bits per coefficient, thus 384 bytes
  // of input entropy.
  unsigned int i, j;

  uint64_t t, d;
  int16_t a, b;

  // We compute two coefficients at a time.
  for (i = 0; i < KYBER_N/4; i++) {
    // We interpret t as a₁ + 2a₂ + 4a₃ + … + 64b₁ + 128b₂ + …
    t = load48_littleendian(buf+3*i);
    d = t & 0x41041041041; // a₁ + 64b₁ + …
    d += (t >> 1) & 0x41041041041; // a₁+a₂ + 64(b₁+b₂) + …
    d += (t >> 2) & 0x41041041041; // a₁+a₂+a₃ + 64(b₁+b₂+b₃) + …
    d += (t >> 3) & 0x41041041041; // a₁+…+a₄ + 64(b₁+…+b₄) + …
    d += (t >> 4) & 0x41041041041; // a₁+…+a₅ + 64(b₁+…+b₅) + …
    d += (t >> 5) & 0x41041041041; // a₁+…+a₆ + 64(b₁+…+b₆) + …

    for (j = 0; j < 4; j++) {
      a = d & 7; // a₁ + … + a₆
      d >>= 6;
      b = d & 7; // b₁ + … + b₆
      d >>= 6;

      r->coeffs[4*i+j] = a - b;
    }
  }
}
#endif


/*************************************************
* Name:        cbd2
*
* Description: Given an array of uniformly random bytes, compute
*              polynomial with coefficients distributed according to
*              a centered binomial distribution with parameter eta=2
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *buf: pointer to input byte array
**************************************************/
static void cbd2(poly *r, const uint8_t buf[2*KYBER_N/4])
{
  unsigned int i,j;
  uint32_t t,d;
  int16_t a,b;

  for(i=0;i<KYBER_N/8;i++) {
    t  = load32_littleendian(buf+4*i);
    d  = t & 0x55555555;
    d += (t>>1) & 0x55555555;

    for(j=0;j<8;j++) {
      a = (d >> (4*j+0)) & 0x3;
      b = (d >> (4*j+2)) & 0x3;
      r->coeffs[8*i+j] = a - b;
    }
  }
}

/*************************************************
* Name:        cbd3
*
* Description: Given an array of uniformly random bytes, compute
*              polynomial with coefficients distributed according to
*              a centered binomial distribution with parameter eta=3.
*              This function is only needed for Kyber-512
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *buf: pointer to input byte array
**************************************************/
static void cbd3(poly *r, const uint8_t buf[3*KYBER_N/4])
{
  unsigned int i,j;
  uint32_t t,d;
  int16_t a,b;

  for(i=0;i<KYBER_N/4;i++) {
    t  = load24_littleendian(buf+3*i);
    d  = t & 0x00249249;
    d += (t>>1) & 0x00249249;
    d += (t>>2) & 0x00249249;

    for(j=0;j<4;j++) {
      a = (d >> (6*j+0)) & 0x7;
      b = (d >> (6*j+3)) & 0x7;
      r->coeffs[4*i+j] = a - b;
    }
  }
}

void poly_cbd_eta1(poly *r, const uint8_t buf[KYBER_ETA1*KYBER_N/4])
{
#if KYBER_ETA1 == 6
  cbd6(r, buf);
#elif KYBER_ETA1 == 3
  cbd3(r, buf);
#elif KYBER_ETA1 == 2
  cbd2(r, buf);
#else
#error "This implementation requires eta1 in {2, 3, 6}"
#endif
}

void poly_cbd_eta2(poly *r, const uint8_t buf[KYBER_ETA2*KYBER_N/4])
{
#if KYBER_ETA2 == 6
  cbd6(r, buf);
#elif KYBER_ETA2 == 3
  cbd3(r, buf);
#elif KYBER_ETA2 == 2
  cbd2(r, buf);
#else
#error "This implementation requires eta1 in {2, 3, 6}"
#endif
}
