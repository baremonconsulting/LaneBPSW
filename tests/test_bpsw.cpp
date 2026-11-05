// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// test_bpsw.cpp — exactness and speed of the lane-packed BPSW (bpsw8.hpp) against GMP 6.3 mpz_millerrabin(n, 25).
//   (1) intermediate Lucas values vs GMP's own mpz_lucas_mod (different formulas: U-chain vs our V-ladder)
//   (2) sub-predicates vs GMP stronglucas / reference strong tests
//   (3) full predicate on random primes, random composites without small factors, composite Fermat numbers
//       F7..F11 and composite Mersenne numbers (both strong base-2 pseudoprimes), prime Mersenne numbers,
//       Chernick Carmichael numbers (Fermat pseudoprimes to every coprime base)
//   (4) timing on primes (the case that occurs once per next-prime query)
#include "lanebpsw/bpsw8.hpp"
#include <vector>
#include <algorithm>
extern "C" int __gmpz_lucas_mod(mpz_ptr V, mpz_ptr Qk, long Q, mp_bitcnt_t b0, mpz_srcptr n, mpz_ptr T1, mpz_ptr T2);

static int sprp_ref(const mpz_t n, const mpz_t b) {
    mpz_t nm1, d, y; mpz_init(nm1); mpz_init(d); mpz_init(y);
    mpz_sub_ui(nm1, n, 1); unsigned long s = mpz_scan1(nm1, 0); mpz_tdiv_q_2exp(d, nm1, s);
    mpz_powm(y, b, d, n);
    int ok = 0;
    if (mpz_cmp_ui(y, 1) == 0 || mpz_cmp(y, nm1) == 0) ok = 1;
    for (unsigned long i = 1; i < s && !ok; i++) { mpz_powm_ui(y, y, 2, n); if (mpz_cmp(y, nm1) == 0) ok = 1; }
    mpz_clear(nm1); mpz_clear(d); mpz_clear(y);
    return ok;
}

static mpz_t g_prim;   // product of odd primes <= 255
static bool no_small_factor(const mpz_t n) { mpz_t g; mpz_init(g); mpz_gcd(g, n, g_prim); bool r = mpz_cmp_ui(g, 1) == 0; mpz_clear(g); return r; }

static int g_bad = 0, g_cnt = 0, g_fb = 0;
static void check_one(const mpz_t n, const char *tag, bool verbose = false) {
    int k = (int)mpz_sizeinbase(n, 2);
    BpswFn f = bpsw_for_bits(k);
    if (!f) return;
    uint64_t fb0 = g_bpsw.fallback;
    g_bdbg.on = true;
    int got = f(n);
    g_bdbg.on = false;
    bool fell = g_bpsw.fallback != fb0;
    g_fb += fell;
    int ref = mpz_millerrabin(n, 25) != 0;
    g_cnt++;
    if (got != ref) { g_bad++; gmp_printf("MISMATCH %s k=%d got=%d ref=%d n=%Zx\n", tag, k, got, ref, n); }
    if (!fell) {
        // sub-predicates: base-2 strong test and the extra MR round vs reference; Lucas vs GMP stronglucas
        mpz_t two, V, Qk; mpz_init_set_ui(two, 2); mpz_init(V); mpz_init(Qk);
        mpz_realloc2(V, 128 * (mpz_size(n) + 4)); mpz_realloc2(Qk, 128 * (mpz_size(n) + 4));   // GMP uses V's limbs as scratch (D=5 path)
        int r2 = sprp_ref(n, two), rm = sprp_ref(n, g_bdbg.a), rl = __gmpz_stronglucas(n, V, Qk) != 0;
        if (r2 != g_bdbg.base2 || rm != g_bdbg.mr || rl != g_bdbg.lucas) {
            g_bad++; gmp_printf("SUBMISMATCH %s k=%d b2 %d/%d mr %d/%d luc %d/%d D=%u\n", tag, k, g_bdbg.base2, r2, g_bdbg.mr, rm,
                               g_bdbg.lucas, rl, g_bdbg.absD);
        }
        // Lucas values vs GMP mpz_lucas_mod (U-chain): when it returns 0 and b0 > 1 it leaves V = V_d, Qk = Q^d
        mpz_t T1, T2; mpz_init(T1); mpz_init(T2);
        mpz_t np1; mpz_init(np1); mpz_add_ui(np1, n, 1); unsigned long b0 = mpz_scan1(np1, 0);
        int early = __gmpz_lucas_mod(V, Qk, g_bdbg.Q, b0, n, T1, T2);
        if (!early && b0 > 1) {
            mpz_mod(V, V, n); mpz_mod(Qk, Qk, n);
            if (mpz_cmp(V, g_bdbg.vd) != 0 || mpz_cmp(Qk, g_bdbg.qd) != 0) { g_bad++; printf("LUCAS VALUE MISMATCH %s k=%d\n", tag, k); }
        }
        if (verbose) printf("  %s k=%d D=%u Q=%ld base2=%d lucas=%d mr=%d -> %d (GMP %d)\n", tag, k, g_bdbg.absD, g_bdbg.Q,
                            g_bdbg.base2, g_bdbg.lucas, g_bdbg.mr, got, ref);
        mpz_clear(T1); mpz_clear(T2); mpz_clear(np1); mpz_clear(two); mpz_clear(V); mpz_clear(Qk);
    } else if (verbose) printf("  %s k=%d fallback -> %d (GMP %d)\n", tag, k, got, ref);
}

int main(int argc, char **argv) {
    require_cpu();
    const bool timing = !(argc > 1 && strcmp(argv[1], "--checks-only") == 0);
    mpz_init_set_ui(g_prim, 1);
    for (unsigned p = 3; p < 256; p += 2) { bool pr = true; for (unsigned q = 3; q * q <= p; q += 2) if (p % q == 0) pr = false; if (pr) mpz_mul_ui(g_prim, g_prim, p); }
    uint64_t seed = 20261028;
    mpz_t n; mpz_init(n);
    // (3a) random primes and random composites without small factors
    int sizes[] = {72, 130, 200, 256, 300, 384, 512, 768, 1024, 1536, 2048, 3072, 4096};
    for (int k : sizes) {
        int np = k <= 1024 ? 300 : (k <= 2048 ? 60 : 12), nc = np;
        for (int i = 0; i < np; i++) {
            do { rand_kbit(n, k, seed); mpz_nextprime(n, n); } while ((int)mpz_sizeinbase(n, 2) != k);
            check_one(n, "prime");
        }
        for (int i = 0; i < nc; i++) {
            do { rand_kbit(n, k, seed); mpz_setbit(n, 0); } while (!no_small_factor(n));
            // force long 2-adic tails sometimes
            if (i % 5 == 1) { mpz_tdiv_q_2exp(n, n, 20); mpz_mul_2exp(n, n, 20); mpz_add_ui(n, n, 1); }      // n = 1 mod 2^20
            if (i % 5 == 2) { mpz_tdiv_q_2exp(n, n, 20); mpz_mul_2exp(n, n, 20); mpz_sub_ui(n, n, 1); }      // n = -1 mod 2^20
            if (!no_small_factor(n)) continue;
            check_one(n, "composite");
        }
        printf("k=%4d done: checks %d, mismatches %d, fallbacks %d\n", k, g_cnt, g_bad, g_fb); fflush(stdout);
    }
    // primes with long tails: n = 1 + c*2^e and n = c*2^e - 1
    for (int k : {256, 1024}) for (int e : {30, 60, 100}) for (int t = 0; t < 20; t++) {
        mpz_t c; mpz_init(c); rand_kbit(c, k - e, seed);
        for (int sgn = 0; sgn < 2; sgn++) {
            for (;;) { mpz_mul_2exp(n, c, e); if (sgn) mpz_sub_ui(n, n, 1); else mpz_add_ui(n, n, 1);
                       if (mpz_probab_prime_p(n, 30)) break; mpz_add_ui(c, c, 1); }
            check_one(n, sgn ? "prime_m1tail" : "prime_p1tail");
        }
        mpz_clear(c);
    }
    printf("tail primes done: checks %d, mismatches %d, fallbacks %d\n", g_cnt, g_bad, g_fb);
    // (3b) Fermat numbers F7..F11 (composite, strong pseudoprimes to base 2)
    for (int m = 7; m <= 11; m++) { mpz_set_ui(n, 0); mpz_setbit(n, 1UL << m); mpz_add_ui(n, n, 1); check_one(n, "Fermat F_m", true); }
    // (3c) Mersenne numbers 2^p - 1, p prime in [131, 1300]: composite ones are strong base-2 pseudoprimes, d_L = 1
    int mers_c = 0, mers_p = 0;
    for (unsigned p = 131; p < 1300; p += 2) {
        bool pr = true; for (unsigned q = 3; q * q <= p; q += 2) if (p % q == 0) pr = false;
        if (!pr) continue;
        mpz_set_ui(n, 0); mpz_setbit(n, p); mpz_sub_ui(n, n, 1);
        int before = g_bad; check_one(n, "Mersenne", false);
        if (mpz_millerrabin(n, 25)) mers_p++; else mers_c++;
        if (g_bad != before) printf("  (Mersenne p=%u)\n", p);
    }
    printf("Mersenne numbers: %d composite, %d prime (M521, M607, M1279 expected prime)\n", mers_c, mers_p);
    // (3d) Chernick Carmichael numbers (6k+1)(12k+1)(18k+1)
    int found = 0;
    for (int bits : {64, 100}) {
        mpz_t kk, a, b, c; mpz_init(kk); mpz_init(a); mpz_init(b); mpz_init(c);
        rand_kbit(kk, bits, seed);
        for (int it = 0; it < 3000000 && found < (bits == 64 ? 3 : 6); it++) {
            mpz_add_ui(kk, kk, 1);
            mpz_mul_ui(a, kk, 6); mpz_add_ui(a, a, 1);
            if (!mpz_probab_prime_p(a, 1)) continue;
            mpz_mul_ui(b, kk, 12); mpz_add_ui(b, b, 1);
            if (!mpz_probab_prime_p(b, 1)) continue;
            mpz_mul_ui(c, kk, 18); mpz_add_ui(c, c, 1);
            if (!mpz_probab_prime_p(c, 1)) continue;
            mpz_mul(n, a, b); mpz_mul(n, n, c);
            found++;
            check_one(n, "Carmichael", true);
        }
        mpz_clear(kk); mpz_clear(a); mpz_clear(b); mpz_clear(c);
    }
    printf("ALL: checks %d, MISMATCHES %d, fallbacks %d\n", g_cnt, g_bad, g_fb);
    if (!timing) return g_bad != 0;
    // (4) timing on primes: lanes vs GMP, interleaved, min and median
    printf("\nTiming on primes (the final check of a next-prime query):\n");
    for (int k : {256, 512, 768, 1024, 1536, 2048, 3072, 4096}) {
        int cnt = k <= 1024 ? 100 : (k <= 2048 ? 40 : 12);
        std::vector<double> tl, tg;
        BpswFn f = bpsw_for_bits(k);
        for (int i = 0; i < cnt; i++) {
            do { rand_kbit(n, k, seed); mpz_nextprime(n, n); } while ((int)mpz_sizeinbase(n, 2) != k);
            for (int rep = 0; rep < 3; rep++) {
                double t0 = cpu_s(); volatile int r1 = f(n); double t1 = cpu_s(); volatile int r2 = mpz_millerrabin(n, 25); double t2 = cpu_s();
                (void)r1; (void)r2; tl.push_back(t1 - t0); tg.push_back(t2 - t1);
            }
        }
        std::sort(tl.begin(), tl.end()); std::sort(tg.begin(), tg.end());
        double ml = 0, mg = 0; for (double x : tl) ml += x; for (double x : tg) mg += x; ml /= tl.size(); mg /= tg.size();
        printf("k=%4d: lanes median %.1f us (mean %.1f) | GMP mpz_millerrabin(n,25) median %.1f us (mean %.1f) | ratio(median) %.2fx\n",
               k, 1e6 * tl[tl.size() / 2], 1e6 * ml, 1e6 * tg[tg.size() / 2], 1e6 * mg, tg[tg.size() / 2] / tl[tl.size() / 2]);
        fflush(stdout);
    }
    return g_bad != 0;
}
