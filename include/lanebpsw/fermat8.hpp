// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// fermat8.hpp — 8-lane base-2 Fermat test, one candidate (and modulus) per lane: the filter of the searches.
#pragma once
#include "common.hpp"
#include "ifma.hpp"

// Per-lane setup: digits of N, k0, and Montgomery representation of 2 (2R mod N).
// cand[l] points to nl limbs; all candidates must have exactly k bits (k = bit length).
template <int ND>
struct Fermat8Setup {
    alignas(64) uint64_t nd[8][ND];
    alignas(64) uint64_t twor[8][ND];
    alignas(64) uint64_t k0[8];
};

// scratch mpz objects reused across calls (single-threaded lab code)
struct Scratch {
    mpz_t t, n, r;
    Scratch() { mpz_init(t); mpz_init(n); mpz_init(r); }
    ~Scratch() { mpz_clear(t); mpz_clear(n); mpz_clear(r); }
};
// Function-local construction prevents AVX-targeted translation units from running
// compiler-generated AVX-512 global initializers while the shared object is loaded.
static Scratch &g_scr_instance() { static Scratch v; return v; }
#define g_scr (g_scr_instance())

// 2R mod N computed with GMP (R = 2^(52*ND)).
template <int ND>
static inline void setup_lane(Fermat8Setup<ND> &s, int l, const mp_limb_t *limbs, int nl) {
    limbs_to_d52(limbs, nl, s.nd[l], ND);
    s.k0[l] = neg_inv52(limbs[0]);
    mpz_t nview;
    mpz_roinit_n(nview, limbs, nl);
    mpz_set_ui(g_scr.t, 0);
    mpz_setbit(g_scr.t, 52 * ND + 1);       // 2R
    mpz_tdiv_r(g_scr.r, g_scr.t, nview);    // 2R mod N
    uint64_t tmp[2 * 128] = {0};
    size_t cnt = 0;
    mpz_export(tmp, &cnt, -1, 8, 0, 0, g_scr.r);
    limbs_to_d52(tmp, nl, s.twor[l], ND);
}

// Lockstep left-to-right binary exponentiation 2^(N_l - 1) mod N_l for 8 lanes.
// All N_l have bit length k. Returns bit l set iff 2^(N_l-1) == 1 (mod N_l).
template <int ND>
static uint8_t fermat2_run(const Fermat8Setup<ND> &s, const mp_limb_t *const cand[8], int k) {
    Mont8<ND> c;
    load_lanes<ND>(c.N, s.nd);
    c.N[ND] = _mm512_setzero_si512();
    c.k0 = _mm512_load_si512((const void *)s.k0);
    __m512i x[ND];
    load_lanes<ND>(x, s.twor);
    // exponent E = N - 1: bits b>=1 equal those of N, bit 0 is 0 (N odd)
    for (int b = k - 2; b >= 1; b--) {
        amsx<ND>(x, x, c);
        int li = b >> 6, sh = b & 63;
        unsigned m = 0;
        for (int l = 0; l < 8; l++) m |= (unsigned)((cand[l][li] >> sh) & 1) << l;
        dbl_mask<ND>(x, (__mmask8)m);          // always executed: no branch on secret exponent bits
    }
    amsx<ND>(x, x, c);   // bit 0 of E is 0: square only
    // convert out of Montgomery form: amm(x, 1) is fully reduced in [1, N-1]
    __m512i one[ND];
    one[0] = _mm512_set1_epi64(1);
    for (int j = 1; j < ND; j++) one[j] = _mm512_setzero_si512();
    ammx<ND>(x, x, one, c);
    __mmask8 ok = _mm512_cmpeq_epi64_mask(x[0], _mm512_set1_epi64(1));
    for (int j = 1; j < ND; j++) ok &= _mm512_cmpeq_epi64_mask(x[j], _mm512_setzero_si512());
    return (uint8_t)ok;
}

template <int ND>
static uint8_t fermat2_8(const mp_limb_t *const cand[8], int nl, int k) {
    Fermat8Setup<ND> s;
    for (int l = 0; l < 8; l++) setup_lane<ND>(s, l, cand[l], nl);
    return fermat2_run<ND>(s, cand, k);
}
