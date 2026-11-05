// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// common.hpp — GMP declarations, timers and a seeded generator shared by the library, tests and benchmarks.
#pragma once
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <x86intrin.h>
#include <gmp.h>

#if GMP_LIMB_BITS != 64 || GMP_NAIL_BITS != 0
#error "LaneBPSW requires GMP with 64-bit limbs and zero nail bits"
#endif

// GMP 6.3.0 entry points that are exported by libgmp.so.10 but not documented: mpz_millerrabin is declared in
// gmp.h, the others in GMP's internal gmp-impl.h. mpz_stronglucas (and mpz_lucas_mod) are used by the tests only.
extern "C" {
int __gmpz_millerrabin(mpz_srcptr n, int reps);            // BPSW (strong base 2 + strong Lucas) + (reps-24) MR
int __gmpz_stronglucas(mpz_srcptr n, mpz_ptr x, mpz_ptr y);
mp_limb_t __gmpn_mod_1s_4p(mp_srcptr ap, mp_size_t n, mp_limb_t b, const mp_limb_t cps[7]);
void __gmpn_mod_1s_4p_cps(mp_limb_t cps[7], mp_limb_t b);
}
#define mpz_millerrabin __gmpz_millerrabin

// 1 if the CPU can run this build: AVX-512F and AVX-512 IFMA, or AVX-512F only for EMULATE_IFMA=1.
static inline int lanebpsw_cpu_ok() {
    __builtin_cpu_init();
#ifdef LANEBPSW_EMULATE_IFMA
    return __builtin_cpu_supports("avx512f") != 0;
#else
    return __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512ifma");
#endif
}
static inline int lanebpsw_gmp_630() { return strcmp(gmp_version, "6.3.0") == 0; }
// tests and benchmarks: exit code 77 (skipped) on a CPU that cannot run this build
static inline void require_cpu() {
    if (!lanebpsw_cpu_ok()) { printf("SKIP: this CPU lacks AVX-512F or AVX-512 IFMA (EMULATE_IFMA=1 needs AVX-512F only)\n"); exit(77); }
}

static inline double now_s() { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
static inline double cpu_s() { struct timespec t; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
static inline uint64_t rdtsc_now() { unsigned aux; return __rdtscp(&aux); }

// splitmix64 — deterministic seeds for reproducible inputs
static inline uint64_t splitmix64(uint64_t &s) {
    uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
// random k-bit number with top bit set (and optionally top two bits), deterministic from seed
static inline void rand_kbit(mpz_t r, int k, uint64_t &seed, bool top_two = false) {
    int nl = (k + 63) / 64;
    mp_limb_t *d = mpz_limbs_write(r, nl);
    for (int i = 0; i < nl; i++) d[i] = splitmix64(seed);
    int topbits = k - 64 * (nl - 1);
    if (topbits < 64) d[nl - 1] &= ((mp_limb_t)1 << topbits) - 1;
    d[nl - 1] |= (mp_limb_t)1 << (topbits - 1);
    if (top_two && k >= 2) {
        if (topbits >= 2) d[nl - 1] |= (mp_limb_t)1 << (topbits - 2);
        else d[nl - 2] |= (mp_limb_t)1 << 63;
    }
    mpz_limbs_finish(r, nl);
}
