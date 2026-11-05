// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// test_paths.cpp — targeted tests of the paths that random inputs rarely or never reach (paper, Section 6.2):
//   (1) fallback of the lane-packed test to GMP: a Jacobi symbol 0 before the first -1, no -1 below 256 (squares),
//       n < 2^64, even n, and n of more than 4102 bits (no kernel);
//   (2) next prime across block ends (the output does not depend on the block length, so tiny blocks force many
//       block changes and short last batches) and across powers of two (batches of mixed bit length);
//   (3) next prime for n + 1 < 2^201, which is passed to GMP, and just above that limit;
//   (4) the chiavdf variant above 2049 bits (lanes against the scalar version, also on strong base-2 pseudoprimes),
//       and the scalar version above 4102 bits against GMP's Baillie-PSW test.
// Every comparison is exact; the program prints a summary and returns 1 on any mismatch.
#include "lanebpsw/chiabpsw8.hpp"
#include "lanebpsw/nextprime_fast.hpp"
#include <vector>

static int g_checks = 0, g_bad = 0;
static void expect(bool ok, const char *what, const mpz_t n) {
    g_checks++;
    if (!ok) { g_bad++; gmp_printf("MISMATCH %s: n = %Zx\n", what, n); }
}

static mpz_t g_small;   // product of the odd primes up to 14951
static bool no_small_factor(const mpz_t n) { mpz_t g; mpz_init(g); mpz_gcd(g, n, g_small); bool r = mpz_cmp_ui(g, 1) == 0; mpz_clear(g); return r; }

static int lane_test(const mpz_t n) { BpswFn f = bpsw_for_bits((int)mpz_sizeinbase(n, 2)); return f ? f(n) : mpz_millerrabin(n, 25); }

int main() {
    require_cpu();
    uint64_t seed = 20261031;
    mpz_t n, m, p, q; mpz_init(n); mpz_init(m); mpz_init(p); mpz_init(q);
    mpz_init_set_ui(g_small, 1);
    for (unsigned a = 3; a <= 14951; a += 2) { bool pr = true; for (unsigned b = 3; b * b <= a; b += 2) if (a % b == 0) { pr = false; break; } if (pr) mpz_mul_ui(g_small, g_small, a); }

    // (1) fallback of the lane-packed test
    uint64_t fb0 = g_bpsw.fallback;
    int jac0 = 0;
    for (int k : {130, 256, 1024}) {
        // a symbol 0 before the first -1: |D| = 7, 11, 13, 15 (via 3 | n), 17 (fixed list) and 19, 23, 29 (search loop);
        // |D| = 25 cannot give 0 here, since (n/25) = (n/5)^2 and (n/5) = +1
        for (unsigned d : {7u, 11u, 13u, 15u, 17u, 19u, 23u, 29u}) {
            for (int found = 0, tries = 0; found < 4 && tries < 100000; tries++) {
                rand_kbit(m, k - 5, seed); mpz_mul_ui(n, m, d == 15 ? 3 : d); mpz_setbit(n, 0);
                bool ok = true;   // all earlier |D| of GMP's sequence give +1
                for (unsigned e : {5u, 7u, 11u, 13u, 15u, 17u, 19u, 23u, 25u}) { if (e >= d) break; if (mpz_kronecker_ui(n, e) != 1) { ok = false; break; } }
                if (!ok || mpz_kronecker_ui(n, d) != 0) continue;
                found++; jac0++;
                expect(lane_test(n) == mpz_millerrabin(n, 25), "fallback, Jacobi 0", n);
            }
        }
        // squares: no D with (n/|D|) = -1
        for (int i = 0; i < 5; i++) { rand_kbit(m, k / 2, seed); mpz_setbit(m, 0); mpz_mul(n, m, m); expect(lane_test(n) == mpz_millerrabin(n, 25), "fallback, square", n); }
        // even numbers
        for (int i = 0; i < 5; i++) { rand_kbit(n, k, seed); mpz_clrbit(n, 0); expect(lane_test(n) == mpz_millerrabin(n, 25), "fallback, even", n); }
    }
    // n < 2^64 and n with more than 4102 bits
    for (int i = 0; i < 200; i++) {
        rand_kbit(n, 40 + i % 25, seed); mpz_setbit(n, 0); expect(lane_test(n) == mpz_millerrabin(n, 25), "n < 2^64", n);
        mpz_nextprime(p, n); expect(lane_test(p) == mpz_millerrabin(p, 25), "prime < 2^64", p);
    }
    expect(bpsw_for_bits(4102) != nullptr && bpsw_for_bits(4103) == nullptr, "kernel limit 4102 bits", n);
    rand_kbit(n, 4103, seed); mpz_nextprime(p, n);
    expect(lane_test(p) == mpz_millerrabin(p, 25) && lane_test(p) == 1, "4103-bit prime via GMP", p);
    uint64_t fallbacks = g_bpsw.fallback - fb0;
    printf("(1) fallback: %d inputs with a Jacobi symbol 0, squares, even, < 2^64 and > 4102 bits; %lu fallbacks taken\n", jac0, (unsigned long)fallbacks);

    // (2) next prime: tiny blocks and powers of two
    NPStats st; g_final_lanes = true;
    long np_checks = 0;
    for (int k : {256, 512, 1024}) {
        int reps = k == 1024 ? 6 : 20;
        for (uint32_t L : {1u, 2u, 3u, 5u, 8u, 9u, 17u}) {
            for (int i = 0; i < reps; i++) {
                rand_kbit(n, k, seed);
                uint64_t b0 = st.blocks;
                nextprime_fast(p, n, 0, L, st);
                mpz_nextprime(q, n);
                expect(mpz_cmp(p, q) == 0, "next prime with tiny blocks", n);
                np_checks++;
                (void)b0;
            }
        }
    }
    printf("(2a) next prime, block lengths 1..17: %ld queries, %lu blocks in all\n", np_checks, (unsigned long)st.blocks);
    uint64_t mixed0 = st.scalar_fallback;
    for (int k : {256, 300, 512, 1024}) {
        for (int j = 1; j <= 40; j++) {
            mpz_set_ui(n, 0); mpz_setbit(n, k); mpz_sub_ui(n, n, 2 * (unsigned long)j + (j & 1));   // just below 2^k
            nextprime_fast(p, n, 0, 0, st); mpz_nextprime(q, n);
            expect(mpz_cmp(p, q) == 0, "next prime across 2^k", n);
            np_checks++;
        }
    }
    printf("(2b) next prime below powers of two: 160 queries, %lu batches of mixed bit length\n", (unsigned long)(st.scalar_fallback - mixed0));

    // (3) next prime around the 201-bit limit
    int small = 0;
    for (int k = 2; k <= 210; k++) {
        for (int i = 0; i < 3; i++) {
            rand_kbit(n, k, seed);
            nextprime_fast(p, n, 0, 0, st); mpz_nextprime(q, n);
            expect(mpz_cmp(p, q) == 0, "next prime, small n", n); small++;
        }
        mpz_set_ui(n, 0); mpz_setbit(n, k); mpz_sub_ui(n, n, 1);   // n + 1 = 2^k
        nextprime_fast(p, n, 0, 0, st); mpz_nextprime(q, n);
        expect(mpz_cmp(p, q) == 0, "next prime, n = 2^k - 1", n); small++;
    }
    printf("(3) next prime for 2..210 bits: %d queries\n", small);

    // (4) chiavdf variant above 2049 bits, and the scalar version above 4102 bits
    g_pp.build();
    int chia = 0;
    for (int k : {2304, 3072, 4096, 4102}) {
        ChiaFn f = chia_for_bits(k);
        for (int i = 0; i < 2; i++) {   // primes
            do { rand_kbit(n, k, seed); mpz_nextprime(p, n); } while ((int)mpz_sizeinbase(p, 2) != k);
            int a = f(p), b = chia_bpsw_ref(p); expect(a == b && a != 0, "chiavdf lanes = scalar, prime", p); chia++;
        }
        for (int i = 0; i < 6; i++) {   // odd numbers without factors up to 14951: they reach the strong test
            do { rand_kbit(n, k, seed); mpz_setbit(n, 0); } while (!no_small_factor(n));
            int a = f(n), b = chia_bpsw_ref(n); expect(a == b, "chiavdf lanes = scalar, odd", n); chia++;
        }
    }
    // strong base-2 pseudoprimes and Mersenne primes between 2049 and 4102 bits: they reach the Lucas-V test
    mpz_set_ui(n, 0); mpz_setbit(n, 4096); mpz_add_ui(n, n, 1);   // F12
    { ChiaFn f = chia_for_bits(4097); g_chia_skip_td = true; int a = f(n); g_chia_skip_td = false;
      int b2 = chia_sprp2_ref(n); expect(b2 == 1 && a == 0, "F12: strong base-2 pseudoprime, rejected by Lucas-V", n); chia++; }
    int mers = 0;
    for (unsigned e = 2051; e < 4102; e += 2) {
        bool pr = true; for (unsigned b = 3; b * b <= e; b += 2) if (e % b == 0) { pr = false; break; }
        if (!pr) continue;
        if (++mers % 8 != 1 && e != 2203 && e != 2281 && e != 3217) continue;   // every eighth exponent, and the Mersenne primes
        mpz_set_ui(n, 0); mpz_setbit(n, e); mpz_sub_ui(n, n, 1);
        ChiaFn f = chia_for_bits((int)e);
        int a = f(n), b = chia_bpsw_ref(n), g = mpz_millerrabin(n, 25) != 0;
        expect(a == b && (a != 0) == g, "Mersenne number: lanes = scalar = GMP", n); chia++;
    }
    for (unsigned e : {4253u, 4423u, 4127u, 4157u}) {   // scalar path: Mersenne primes M4253, M4423 and composites M4127, M4157
        mpz_set_ui(n, 0); mpz_setbit(n, e); mpz_sub_ui(n, n, 1);
        expect(chia_for_bits((int)e) == nullptr, "no chiavdf kernel above 4102 bits", n);
        int a = chia_bpsw_ref(n) != 0, g = mpz_millerrabin(n, 25) != 0;
        expect(a == g, "scalar chiavdf version = GMP above 4102 bits", n); chia++;
    }
    for (int i = 0; i < 4; i++) {
        do { rand_kbit(n, 4200, seed); mpz_setbit(n, 0); } while (!no_small_factor(n));
        int a = chia_bpsw_ref(n) != 0, g = mpz_millerrabin(n, 25) != 0;
        expect(a == g, "scalar chiavdf version = GMP, 4200-bit odd", n); chia++;
    }
    printf("(4) chiavdf variant above 2049 bits and scalar version above 4102 bits: %d inputs\n", chia);

    printf("ALL: checks %d, MISMATCHES %d\n", g_checks, g_bad);
    return g_bad != 0;
}
