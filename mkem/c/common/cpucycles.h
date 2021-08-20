#ifndef CPUCYCLES_H
#define CPUCYCLES_H

#include <stdint.h>

#if defined(__APPLE__) && defined(__aarch64__)

uint64_t cpucycles(void);
void init_cpucycles(void);

#elif defined(USE_RDPMC)  /* Needs echo 2 > /sys/devices/cpu/rdpmc */

static inline void init_cpucycles() { }
static inline uint64_t cpucycles(void) {
  const uint32_t ecx = (1U << 30) + 1;
  uint64_t result;

  __asm__ volatile ("rdpmc; shlq $32,%%rdx; orq %%rdx,%%rax"
    : "=a" (result) : "c" (ecx) : "rdx");

  return result;
}

#else

static inline void init_cpucycles() { }
static inline uint64_t cpucycles(void) {
  uint64_t result;

  __asm__ volatile ("rdtsc; shlq $32,%%rdx; orq %%rdx,%%rax"
    : "=a" (result) : : "%rdx");

  return result;
}

#endif

uint64_t cpucycles_overhead(void);

#endif
