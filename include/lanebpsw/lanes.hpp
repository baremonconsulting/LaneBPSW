// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// lanes.hpp — set-up of the 8-lane Montgomery context from GMP integers: digits of the modulus per lane,
// -N^{-1} mod 2^52 and optionally R mod N (Montgomery form of 1).
#pragma once
#include "common.hpp"
#include "ifma.hpp"

template <int ND>
struct LaneMods {
    alignas(64) Mont8<ND> c;
    alignas(64) __m512i rm[ND];     // R mod N   (Montgomery 1), fully reduced
    alignas(64) __m512i nm1[ND];    // N - 1, fully reduced digits
    bool have_rm = false;
};

static inline void mpz_to_d52(const mpz_t v, uint64_t *dig, int ND) {
    int nl = (int)mpz_size(v);
    limbs_to_d52(mpz_limbs_read(v), nl, dig, ND);
}

// moduli given as mpz (odd, all < 2^(52*ND-4)); if same==true, mod[0] is broadcast to all lanes
template <int ND>
static inline void lanes_setup(LaneMods<ND> &L, const mpz_t *mod, bool same, bool want_rm) {
    alignas(64) uint64_t nd[8][ND], rmd[8][ND], k0[8];
    mpz_t t; mpz_init(t);
    for (int l = 0; l < 8; l++) {
        const mpz_t &n = mod[same ? 0 : l];
        if (same && l > 0) { memcpy(nd[l], nd[0], sizeof nd[0]); memcpy(rmd[l], rmd[0], sizeof rmd[0]); k0[l] = k0[0]; continue; }
        mpz_to_d52(n, nd[l], ND);
        k0[l] = neg_inv52(mpz_getlimbn(n, 0));
        if (want_rm) {
            mpz_set_ui(t, 0); mpz_setbit(t, 52 * ND); mpz_tdiv_r(t, t, n);
            mpz_to_d52(t, rmd[l], ND);
        }
    }
    load_lanes<ND>(L.c.N, nd);
    L.c.N[ND] = _mm512_setzero_si512();
    L.c.k0 = _mm512_load_si512((const void *)k0);
    for (int j = 0; j < ND; j++) L.nm1[j] = L.c.N[j];
    L.nm1[0] = _mm512_sub_epi64(L.nm1[0], _mm512_set1_epi64(1));   // N odd: no borrow
    if (want_rm) { load_lanes<ND>(L.rm, rmd); L.have_rm = true; }
    mpz_clear(t);
}
