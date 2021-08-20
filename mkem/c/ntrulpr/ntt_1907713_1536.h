#pragma once

#include <stdint.h>

// NTT for the polynomial ring GF(1907713)[x]/<x+1>

// Performs in-place forward NTT.
void ntt_1907713_1536(int32_t pp[1536]);

// Performs in-place reverse NTT and multiplies by the
// Montgomery factor R.
void ntt_1907713_1536_inverse(int32_t pp[1536]);

// Set p to the "pointwise" product of p, q and R⁻¹.
//
// Thus, this is multiplication in the NTT and Montgomery
// representtation.
void ntt_1907713_1536_mul(int32_t pp[1536], int32_t qq[1536]);
