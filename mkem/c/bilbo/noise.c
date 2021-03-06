/********************************************************************************************
* FrodoKEM: Learning with Errors Key Encapsulation
*
* Abstract: noise sampling functions
*********************************************************************************************/

#include <stdint.h>

#include "api.h"
#include "common.h"
#include "params.h"

static const uint32_t CDF_TABLE[CDF_TABLE_LEN] = CDF_TABLE_DATA;

void PQCLEAN_FRODOKEM640SHAKE_OPT_sample_n(uint16_t* out, uint32_t *in, size_t n) {
    // Fills vector s with n samples from the noise distribution which requires 32 bits to sample.
    // The distribution is specified by its CDF.
    // Input: pseudo-random values (4*n bytes) passed in s.
    // Output: samples
    size_t i;
    unsigned int j;

    for (i = 0; i < n; ++i) {
        uint16_t sample = 0;
        uint32_t prnd = in[i] >> 1;    // Drop the least significant bit
        uint32_t sign = in[i] & 0x1;    // Pick the least significant bit

        // No need to compare with the last value.
        for (j = 0; j < (unsigned int)(CDF_TABLE_LEN - 1); j++) {
            // Constant time comparison: 1 if CDF_TABLE[j] < s, 0 otherwise. Uses the fact that CDF_TABLE[j] and s fit in 31 bits.
            sample += (uint32_t)(CDF_TABLE[j] - prnd) >> 31;
        }
        // Assuming that sign is either 0 or 1, flips sample iff sign = 1
        out[i] = ((-sign) ^ sample) + sign;
    }
}
