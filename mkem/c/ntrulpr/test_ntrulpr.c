#include "namespace.h"
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "ntrulpr.h"
#include "randombytes.h"

#define NTESTS 1000

static int test_ntt2(void)
{
    int16_t ag[p];
    int16_t ag2[p];
    int8_t a[p];
    int i = 0;


    for (i = 0; i < p; i++) {
        ag[i] = q12;
        ag2[i] = q12;
        a[i] = 0;
    }

    for (i = 0; i < 200; i++) {
        a[i] = 1;
    }

    mult_small_ntt(ag, ag, a);
    mult_small(ag2, a);

    for (i = 0; i < p; i++) {
        if (ag[i] != ag2[i]) {
            printf("Mismatch at index: %d %d %d\n", i, ag[i], ag2[i]);
            return 1;
        }
    }
    return 0;
}

static int test_ntt(void)
{
    int16_t ag[p];
    int16_t ag2[p];
    int8_t a[p];
    int i = 0;
    uint8_t group[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};


    generator(ag, group);
    generator(ag2, group);
    short_random(a);
    mult_small_ntt(ag, ag, a);
    mult_small(ag2, a);

    for (i = 0; i < p; i++) {
        if (ag[i] != ag2[i]) {
            printf("Mismatch at index: %d\n", i);
            return 1;
        }
    }
    return 0;
}

static int test_rounding_encoding(void)
{
    int16_t v[p];
    int16_t v2[p];
    uint8_t buf[1076];
    int16_t c = -3939;

    do {
        for (int i = 0; i < p; i++) {
            v[i] = c;
            if (c != 3939)
                c += 3;
        }
        crypto_encode_rounded(buf, v);
        crypto_decode_rounded(v2, buf);
        for (int i = 0; i < p; i++) {
            if (v[i] != v2[i]) {
                printf("%d %d %d", i, v[i], v2[i]);
                return 1;
            }
        }
    } while (c != 3939);
    return 0;
}

int main(void)
{
    if(test_rounding_encoding())
        return 1;
    if(test_ntt())
        return 2;
    if(test_ntt2())
        return 3;

    return 0;
}
