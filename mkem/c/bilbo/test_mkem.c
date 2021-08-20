#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "mkem.h"
#include "randombytes.h"

#define NTESTS 1000

static int test_keys()
{
  uint8_t group[MKEM_GROUPBYTES];

  uint8_t pk1[MKEM_PUBLICKEYBYTES];
  uint8_t pk2[MKEM_PUBLICKEYBYTES];

  uint8_t sk1[MKEM_SECRETKEYBYTES];
  uint8_t sk2[MKEM_SECRETKEYBYTES];

  uint8_t cti[MKEM_CTIBYTES];
  uint8_t ctd1[MKEM_CTDBYTES];
  uint8_t ctd2[MKEM_CTDBYTES];

  const uint8_t* pks[2] = {pk1, pk2};
  uint8_t* ctds[2] = {ctd1, ctd2};

  uint8_t ss_a[MKEM_SSBYTES];
  uint8_t ss_b1[MKEM_SSBYTES];
  uint8_t ss_b2[MKEM_SSBYTES];

  mkem_group(group);

  mkem_keypair(pk1, sk1, group);
  mkem_keypair(pk2, sk2, group);

  mkem_enc(ctds, cti, ss_a, group, pks, 2);

  mkem_dec(ss_b1, cti, ctd1, group, pk1, sk1);
  mkem_dec(ss_b2, cti, ctd2, group, pk2, sk2);

  if(memcmp(ss_a, ss_b1, MKEM_SSBYTES) ||
          memcmp(ss_b2, ss_b1, MKEM_SSBYTES)) {
    printf("ERROR keys\n");
    return 1;
  }

  return 0;
}


int main(void)
{
  unsigned int i;
  int r;

  for(i=0;i<NTESTS;i++) {
    r  = test_keys();
    if(r)
      return 1;
  }

  printf("MKEM_CTIBYTES       %d\n", MKEM_CTIBYTES);
  printf("MKEM_CTDBYTES       %d\n", MKEM_CTDBYTES);
  printf("MKEM_PUBLICKEYBYTES %d\n", MKEM_PUBLICKEYBYTES);
  printf("MKEM_SECRETKEYBYTES %d\n", MKEM_SECRETKEYBYTES);

  return 0;
}
