// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// test_kernels.cpp — the 8-lane kernels against GMP and scalar arithmetic, for every instantiated digit count:
//   (1) the two IFMA operations (or their emulation, make EMULATE_IFMA=1) against 128-bit products, with factors
//       that have bits above bit 51 (ignored by the instructions);
//   (2) Montgomery multiplication: r = xy R^{-1} mod N and r < 2N for x, y < fN, eight different moduli per call,
//       including the edges 0 and fN - 1, with f = 4 for 52 ND >= k + 4 (filter) and f = 6 for 52 ND >= k + 6 (test);
//   (3) the dedicated squaring: bit for bit equal to the multiplication of x by itself;
//   (4) the 8-lane base-2 Fermat test against mpz_powm on random odd numbers and primes of the smallest and largest
//       bit length of each kernel.
// Prints a summary and returns 1 on any mismatch.
#include "lanebpsw/fermat8.hpp"
#include "lanebpsw/lanes.hpp"
#include <vector>

static int g_checks = 0, g_bad = 0;
static void expect(bool ok, const char *what, int nd, int k) {
    g_checks++;
    if (!ok && g_bad++ < 20) printf("MISMATCH %s: ND=%d k=%d\n", what, nd, k);
}

static void test_ifma_ops(uint64_t &seed) {
    alignas(64) uint64_t a[8], b[8], c[8], lo[8], hi[8];
    const uint64_t M = MASK52;
    for (int t = 0; t < 20000; t++) {
        for (int l = 0; l < 8; l++) {
            a[l] = splitmix64(seed); b[l] = splitmix64(seed); c[l] = splitmix64(seed);
            if (t % 4 == 1) { a[l] |= M; b[l] |= M; }                 // largest 52-bit factors
            if (t % 4 == 2) { a[l] &= M; b[l] = M; c[l] = ~0ULL; }    // wrap-around of the 64-bit sum
        }
        __m512i va = _mm512_load_si512(a), vb = _mm512_load_si512(b), vc = _mm512_load_si512(c);
        _mm512_store_si512(lo, mlo(vc, va, vb));
        _mm512_store_si512(hi, mhi(vc, va, vb));
        for (int l = 0; l < 8; l++) {
            unsigned __int128 p = (unsigned __int128)(a[l] & M) * (b[l] & M);
            bool ok = lo[l] == c[l] + (uint64_t)(p & M) && hi[l] == c[l] + (uint64_t)(p >> 52);
            expect(ok, "IFMA operation", 0, 52);
        }
    }
}

template <int ND>
static void test_mont(int k, int trials, uint64_t &seed, unsigned f = 4) {
    mpz_t mods[8], x[8], y[8], r, t, R, Rinv;
    for (int l = 0; l < 8; l++) { mpz_init(mods[l]); mpz_init(x[l]); mpz_init(y[l]); }
    mpz_init(r); mpz_init(t); mpz_init(R); mpz_init(Rinv);
    mpz_setbit(R, 52 * ND);
    for (int tr = 0; tr < trials; tr++) {
        for (int l = 0; l < 8; l++) { rand_kbit(mods[l], k, seed); mpz_setbit(mods[l], 0); }
        LaneMods<ND> L; lanes_setup<ND>(L, mods, false, false);
        alignas(64) uint64_t xd[8][ND], yd[8][ND], rd[8][ND];
        for (int l = 0; l < 8; l++) {
            mpz_mul_ui(t, mods[l], f);                                    // fN
            rand_kbit(x[l], 52 * ND - 1, seed); mpz_mod(x[l], x[l], t);
            rand_kbit(y[l], 52 * ND - 1, seed); mpz_mod(y[l], y[l], t);
            if (tr % 5 == 0) { mpz_sub_ui(x[l], t, 1); mpz_sub_ui(y[l], t, 1); }   // fN - 1
            if (tr % 7 == 1 && l < 4) mpz_set_ui(x[l], 0);
            mpz_to_d52(x[l], xd[l], ND); mpz_to_d52(y[l], yd[l], ND);
        }
        __m512i vx[ND], vy[ND], vr[ND], vs[ND], vq[ND];
        load_lanes<ND>(vx, xd); load_lanes<ND>(vy, yd);
        ammx<ND>(vr, vx, vy, L.c);
        ammx<ND>(vq, vx, vx, L.c);
        store_lanes<ND>(rd, vr);
        if (f == 4) {   // the dedicated squaring is used by the filter only (inputs < 4N)
            if constexpr (ND <= 24) ams<ND>(vs, vx, L.c); else ams_rolled<ND>(vs, vx, L.c);
            bool same = true;
            for (int j = 0; j < ND; j++) same &= _mm512_cmpeq_epi64_mask(vq[j], vs[j]) == 0xFF;
            expect(same, "squaring = multiplication", ND, k);
        }
        for (int l = 0; l < 8; l++) {
            bool digits_ok = true;
            for (int j = 0; j < ND; j++) digits_ok &= rd[l][j] <= MASK52;
            uint64_t limbs[160]; int nl = (52 * ND + 63) / 64;
            d52_to_limbs(rd[l], ND, limbs, nl);
            mpz_import(r, nl, -1, 8, 0, 0, limbs);
            mpz_invert(Rinv, R, mods[l]);
            mpz_mul(t, x[l], y[l]); mpz_mul(t, t, Rinv); mpz_sub(t, t, r);
            mpz_mul_2exp(Rinv, mods[l], 1);                               // 2N
            expect(digits_ok && mpz_divisible_p(t, mods[l]) && mpz_cmp(r, Rinv) < 0, "Montgomery product", ND, k);
        }
    }
    for (int l = 0; l < 8; l++) { mpz_clear(mods[l]); mpz_clear(x[l]); mpz_clear(y[l]); }
    mpz_clear(r); mpz_clear(t); mpz_clear(R); mpz_clear(Rinv);
}

template <int ND>
static void test_fermat(int k, int batches, uint64_t &seed) {
    int nl = (k + 63) / 64;
    mpz_t v[8], e, r, two; for (int l = 0; l < 8; l++) mpz_init(v[l]);
    mpz_init(e); mpz_init(r); mpz_init_set_ui(two, 2);
    for (int b = 0; b < batches; b++) {
        for (int l = 0; l < 8; l++) {
            do {
                rand_kbit(v[l], k, seed); mpz_setbit(v[l], 0);
                if ((l + b) % 3 == 0) mpz_nextprime(v[l], v[l]);
            } while ((int)mpz_sizeinbase(v[l], 2) != k);
        }
        const mp_limb_t *cand[8];
        for (int l = 0; l < 8; l++) cand[l] = mpz_limbs_read(v[l]);
        uint8_t m = fermat2_8<ND>(cand, nl, k);
        for (int l = 0; l < 8; l++) {
            mpz_sub_ui(e, v[l], 1); mpz_powm(r, two, e, v[l]);
            expect(((m >> l) & 1) == (mpz_cmp_ui(r, 1) == 0), "8-lane Fermat test", ND, k);
        }
    }
    for (int l = 0; l < 8; l++) mpz_clear(v[l]);
    mpz_clear(e); mpz_clear(r); mpz_clear(two);
}

// k from the smallest to the largest bit length that the kernel serves (52 ND >= k + 4 in the filter)
template <int ND>
static void test_nd(int kmin, int trials, int batches, uint64_t &seed) {
    int kmax = 52 * ND - 4;
    for (int k : {kmin, (kmin + kmax) / 2, kmax}) test_mont<ND>(k, trials, seed);
    if (ND > 5) test_mont<ND>(kmax - 2, trials, seed, 6);   // inputs < 6N of the lane-packed test (52 ND >= k + 6)
    for (int k : {kmin, kmax}) test_fermat<ND>(k, batches, seed);
    printf("ND=%2d: %4d to %4d bits\n", ND, kmin, kmax);
}

int main() {
    require_cpu();
    uint64_t seed = 20261031;
#ifdef LANEBPSW_EMULATE_IFMA
    printf("IFMA emulated with AVX-512F (EMULATE_IFMA=1)\n");
#endif
    test_ifma_ops(seed);
    test_nd<5>(65, 200, 3, seed);
    test_nd<6>(257, 200, 3, seed);
    test_nd<10>(309, 150, 3, seed);
    test_nd<15>(517, 100, 2, seed);
    test_nd<20>(777, 80, 2, seed);
    test_nd<30>(1037, 40, 1, seed);
    test_nd<40>(1557, 25, 1, seed);
    test_nd<60>(2077, 12, 1, seed);
    test_nd<79>(3117, 8, 1, seed);
    printf("ALL: checks %d, MISMATCHES %d\n", g_checks, g_bad);
    return g_bad != 0;
}
