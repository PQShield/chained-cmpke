/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: benchmarking/testing isogeny-based key encapsulation mechanism
*********************************************************************************************/ 

#include <stdio.h>
#include <string.h>

#include "sike.h"

#define TEST_LOOPS        10


static int cryptotest_kem()
{ // Testing KEM
    unsigned int i;
    unsigned char sk[CRYPTO_SECRETKEYBYTES] = {0};
    unsigned char pk[CRYPTO_PUBLICKEYBYTES] = {0};
    unsigned char ct[CRYPTO_CIPHERTEXTBYTES] = {0};
    unsigned char ss[CRYPTO_BYTES] = {0};
    unsigned char ss_[CRYPTO_BYTES] = {0};
    bool passed = true;

    printf("\n\nTESTING ISOGENY-BASED KEY ENCAPSULATION MECHANISM SIKEp434\n");
    printf("--------------------------------------------------------------------------------------------------------\n\n");

    for (i = 0; i < TEST_LOOPS; i++) 
    {
        crypto_kem_keypair(pk, sk);
        crypto_kem_enc(ct, ss, pk);
        crypto_kem_dec(ss_, ct, sk);
        
        if (memcmp(ss, ss_, CRYPTO_BYTES) != 0) {
            passed = false;
            break;
        }
    }

    if (passed == true) printf("  KEM tests .................................................... PASSED");
    else { printf("  KEM tests ... FAILED"); printf("\n"); return 1; }
    printf("\n"); 

    return 0;
}


int main()
{
    int Status = 0;
    
    Status = cryptotest_kem();             // Test key encapsulation mechanism
    if (Status != 0) {
        printf("\n\n   Error detected: KEM_ERROR_SHARED_KEY \n\n");
        return 1;
    }

    return Status;
}
