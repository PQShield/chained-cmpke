#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <locale.h>

#include "mkem.h"
#include "cpucycles.h"
#include "speed_print.h"

#define NTESTS 1000
#define NRECIPS 100

uint64_t t[NTESTS];

int main()
{
  unsigned int i;
  uint8_t ss[MKEM_SSBYTES];
  uint8_t mpks[MKEM_PUBLICKEYBYTES*NRECIPS];
  uint8_t msks[MKEM_SECRETKEYBYTES*NRECIPS];
  uint8_t mcti[MKEM_CTIBYTES];
  uint8_t mctd[MKEM_CTDBYTES*NRECIPS];
  uint8_t mgrp[MKEM_GROUPBYTES];
  uint8_t* mpctd[NRECIPS];
  const uint8_t* mppks[NRECIPS];

  setlocale(LC_NUMERIC, "");
  init_cpucycles();

  mkem_group(mgrp);
  for(i=0;i<NRECIPS;i++) {
      mppks[i] = mpks + MKEM_PUBLICKEYBYTES*i;
      mpctd[i] = mctd + MKEM_CTDBYTES*i;
      mkem_keypair(
        mpks + MKEM_PUBLICKEYBYTES*i,
        msks + MKEM_SECRETKEYBYTES*i, mgrp
      );
  }

  for(i=0;i<NTESTS;i++) {
    t[i] = cpucycles();
    mkem_keypair(mpks, msks, mgrp);
  }
  print_results("mkem_keypair: ", t, NTESTS);

  for(i=0;i<NTESTS;i++) {
    t[i] = cpucycles();
    mkem_enc(mpctd, mcti, ss, mgrp, mppks, 1);
  }
  print_results("mkem_enc (x1): ", t, NTESTS);

  for(i=0;i<NTESTS;i++) {
    t[i] = cpucycles();
    mkem_enc(mpctd, mcti, ss, mgrp, mppks, 10);
  }
  print_results("mkem_enc (x10): ", t, NTESTS);

  for(i=0;i<NTESTS;i++) {
    t[i] = cpucycles();
    mkem_enc(mpctd, mcti, ss, mgrp, mppks, NRECIPS);
  }
  print_results("mkem_enc (x100): ", t, NTESTS);

  for(i=0;i<NTESTS;i++) {
    t[i] = cpucycles();
    mkem_dec(ss, mcti, mctd, mgrp, mpks, msks);
  }
  print_results("mkem_dec: ", t, NTESTS);

  return 0;
}
