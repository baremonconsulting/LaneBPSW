// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// ifma.hpp — 8-lane AVX-512 IFMA Montgomery arithmetic in radix 2^52, vertical layout: vector j holds digit j
// of eight numbers, one per lane; each lane may have its own modulus N_l.
//
// The filter uses x,y < 4N and R > 16N. Lucas callers use x,y < 6N and kernel selection gives R > 64N.
// Montgomery reduction is below N + 36N^2/R < 25N/16 < 2N. Caller-specific state and signed-normalisation
// invariants are documented in docs/ALGORITHM_INVARIANTS.md.
// Digits passed to vpmadd52* must be < 2^52 (the instruction reads only the low 52 bits),
// accumulators may be unnormalized (< 2^64).
#pragma once
#include <immintrin.h>
#include <stdint.h>
#include <string.h>

#define IFMA_INLINE __attribute__((always_inline)) inline

#ifndef LANEBPSW_EMULATE_IFMA
static IFMA_INLINE __m512i mlo(__m512i acc, __m512i a, __m512i b) { return _mm512_madd52lo_epu64(acc, a, b); }
static IFMA_INLINE __m512i mhi(__m512i acc, __m512i a, __m512i b) { return _mm512_madd52hi_epu64(acc, a, b); }
#else
// Exact emulation of VPMADD52LUQ / VPMADD52HUQ with AVX-512F (make EMULATE_IFMA=1), so that the tests can run on
// CPUs with AVX-512F but without IFMA. Much slower: timings of such a build mean nothing. The 52-bit factors are
// split into 26-bit halves, a = a1 2^26 + a0, so that ab = a1 b1 2^52 + (a0 b1 + a1 b0) 2^26 + a0 b0.
static IFMA_INLINE void mul52_emu(__m512i a, __m512i b, __m512i &lo, __m512i &hi) {
    const __m512i m26 = _mm512_set1_epi64((1LL << 26) - 1), m52 = _mm512_set1_epi64((1LL << 52) - 1);
    a = _mm512_and_si512(a, m52); b = _mm512_and_si512(b, m52);          // the instructions read the low 52 bits
    const __m512i a0 = _mm512_and_si512(a, m26), a1 = _mm512_srli_epi64(a, 26);
    const __m512i b0 = _mm512_and_si512(b, m26), b1 = _mm512_srli_epi64(b, 26);
    const __m512i mid = _mm512_add_epi64(_mm512_mul_epu32(a0, b1), _mm512_mul_epu32(a1, b0));                  // < 2^53
    const __m512i low = _mm512_add_epi64(_mm512_mul_epu32(a0, b0), _mm512_slli_epi64(_mm512_and_si512(mid, m26), 26)); // < 2^53
    lo = _mm512_and_si512(low, m52);
    hi = _mm512_add_epi64(_mm512_add_epi64(_mm512_mul_epu32(a1, b1), _mm512_srli_epi64(mid, 26)), _mm512_srli_epi64(low, 52));
}
static IFMA_INLINE __m512i mlo(__m512i acc, __m512i a, __m512i b) { __m512i l, h; mul52_emu(a, b, l, h); return _mm512_add_epi64(acc, l); }
static IFMA_INLINE __m512i mhi(__m512i acc, __m512i a, __m512i b) { __m512i l, h; mul52_emu(a, b, l, h); return _mm512_add_epi64(acc, h); }
#endif

static const uint64_t MASK52 = (1ULL << 52) - 1;

// digits needed for k-bit moduli so that 52*ND >= k+4
static inline int nd_for_bits(int k) {
    return k > 0 ? (int)(((int64_t)k + 55) / 52) : 0;
}

// ---- conversion helpers (scalar) ----
// limbs (64-bit, little endian) -> ND digits of 52 bits
static inline void limbs_to_d52(const uint64_t *limbs, int nl, uint64_t *dig, int ND) {
    for (int j = 0; j < ND; j++) {
        int bit = 52 * j, li = bit >> 6, sh = bit & 63;
        uint64_t v = 0;
        if (li < nl) {
            v = limbs[li] >> sh;
            if (sh > 12 && li + 1 < nl) v |= limbs[li + 1] << (64 - sh);
        }
        dig[j] = v & MASK52;
    }
}
static inline void d52_to_limbs(const uint64_t *dig, int ND, uint64_t *limbs, int nl) {
    memset(limbs, 0, sizeof(uint64_t) * nl);
    for (int j = 0; j < ND; j++) {
        int bit = 52 * j, li = bit >> 6, sh = bit & 63;
        if (li < nl) limbs[li] |= dig[j] << sh;
        if (sh > 12 && li + 1 < nl) limbs[li + 1] |= dig[j] >> (64 - sh);
    }
}
// -N^{-1} mod 2^52 from the low limb (N odd)
static inline uint64_t neg_inv52(uint64_t n0) {
    uint64_t x = n0;                    // correct to 3 bits (n*n = 1 mod 8)
    for (int i = 0; i < 5; i++) x *= 2 - n0 * x;   // 6,12,24,48,96 bits
    return (0 - x) & MASK52;
}

// ---- 8-lane Montgomery context ----
template <int ND>
struct alignas(64) Mont8 {
    __m512i N[ND + 1];   // N[ND] = 0 (padding used by the fused loop)
    __m512i k0;
};

// Load per-lane digit arrays dig[l][ND] into vertical vectors.
template <int ND>
static inline void load_lanes(__m512i *v, const uint64_t dig[8][ND]) {
    alignas(64) uint64_t tmp[8];
    for (int j = 0; j < ND; j++) {
        for (int l = 0; l < 8; l++) tmp[l] = dig[l][j];
        v[j] = _mm512_load_si512((const void *)tmp);
    }
}
template <int ND>
static inline void store_lanes(uint64_t dig[8][ND], const __m512i *v) {
    alignas(64) uint64_t tmp[8];
    for (int j = 0; j < ND; j++) {
        _mm512_store_si512((void *)tmp, v[j]);
        for (int l = 0; l < 8; l++) dig[l][j] = tmp[l];
    }
}

// ---- Montgomery multiplication, fused CIOS (one pass per outer digit) ----
// r = a*b*R^{-1} mod N (lane-wise), r < 2N if a,b < 4N.  r may alias a or b.
template <int ND>
static IFMA_INLINE void amm(__m512i *r, const __m512i *a, const __m512i *b, const Mont8<ND> &c) {
    const __m512i zero = _mm512_setzero_si512();
    const __m512i mask = _mm512_set1_epi64(MASK52);
    __m512i acc[ND + 1];
#pragma GCC unroll 128
    for (int j = 0; j <= ND; j++) acc[j] = zero;
#pragma GCC unroll 128
    for (int i = 0; i < ND; i++) {
        const __m512i ai = a[i];
        __m512i t0 = mlo(acc[0], ai, b[0]);
        const __m512i m = mlo(zero, t0, c.k0);
        t0 = mlo(t0, m, c.N[0]);
        const __m512i carry = _mm512_srli_epi64(t0, 52);
        // new acc[j] = acc[j+1] + hi(ai*b_j) + lo(ai*b_{j+1}) + hi(m*n_j) + lo(m*n_{j+1})
#pragma GCC unroll 128
        for (int j = 0; j < ND - 1; j++) {
            __m512i t = acc[j + 1];
            t = mhi(t, ai, b[j]);
            t = mlo(t, ai, b[j + 1]);
            t = mhi(t, m, c.N[j]);
            t = mlo(t, m, c.N[j + 1]);
            acc[j] = t;
        }
        acc[ND - 1] = mhi(mhi(acc[ND], ai, b[ND - 1]), m, c.N[ND - 1]);
        acc[0] = _mm512_add_epi64(acc[0], carry);
        acc[ND] = zero;
    }
    __m512i carry = zero;
#pragma GCC unroll 128
    for (int j = 0; j < ND; j++) {
        __m512i t = _mm512_add_epi64(acc[j], carry);
        r[j] = _mm512_and_si512(t, mask);
        carry = _mm512_srli_epi64(t, 52);
    }
}

// Same, with a rolled outer loop (for large ND, keeps code size in L1i).
template <int ND>
static inline void amm_rolled(__m512i *r, const __m512i *a, const __m512i *b, const Mont8<ND> &c) {
    const __m512i zero = _mm512_setzero_si512();
    const __m512i mask = _mm512_set1_epi64(MASK52);
    alignas(64) __m512i acc[ND + 1];
    for (int j = 0; j <= ND; j++) acc[j] = zero;
    for (int i = 0; i < ND; i++) {
        const __m512i ai = a[i];
        __m512i t0 = mlo(acc[0], ai, b[0]);
        const __m512i m = mlo(zero, t0, c.k0);
        t0 = mlo(t0, m, c.N[0]);
        const __m512i carry = _mm512_srli_epi64(t0, 52);
#pragma GCC unroll 8
        for (int j = 0; j < ND - 1; j++) {
            __m512i t = acc[j + 1];
            t = mhi(t, ai, b[j]);
            t = mlo(t, ai, b[j + 1]);
            t = mhi(t, m, c.N[j]);
            t = mlo(t, m, c.N[j + 1]);
            acc[j] = t;
        }
        acc[ND - 1] = mhi(mhi(acc[ND], ai, b[ND - 1]), m, c.N[ND - 1]);
        acc[0] = _mm512_add_epi64(acc[0], carry);
        acc[ND] = zero;
    }
    __m512i carry = zero;
    for (int j = 0; j < ND; j++) {
        __m512i t = _mm512_add_epi64(acc[j], carry);
        r[j] = _mm512_and_si512(t, mask);
        carry = _mm512_srli_epi64(t, 52);
    }
}

template <int ND>
static IFMA_INLINE void ammx(__m512i *r, const __m512i *a, const __m512i *b, const Mont8<ND> &c) {
    if constexpr (ND <= 24) amm<ND>(r, a, b, c);
    else amm_rolled<ND>(r, a, b, c);
}

// ---- Dedicated Montgomery squaring (symmetric products), r = a^2 * R^{-1} mod N, r < 2N if a < 4N ----
// Row-wise (CIOS-like) with shifting accumulators. Row i adds a_i^2 at position 2i and 2*a_i*a_j (j>i) at
// position i+j, using precomputed digits of 2A: E_j = (2 a_j) mod 2^52 (no carry-in) for the first off-diagonal
// digit of a row and D_j = E_j | (a_{j-1} >> 51) for the following ones (digits of 2*A_{>i} are exactly
// E_{i+1}, D_{i+2}, ..., D_{ND-1}; the top digit a_{ND-1} >> 51 is 0 because a < 4N < R/4).
// All product contributions to position i come from rows <= i, so m_i can be formed at the start of row i.
// The result equals amm(a, a) bit for bit (same T = A^2, same M). IFMA count 3*ND^2 + 2*ND vs 4*ND^2 + ND
// (1,240 vs 1,620 at ND = 20, counted in the object code).
template <int ND>
static IFMA_INLINE void ams(__m512i *r, const __m512i *a, const Mont8<ND> &c) {
    const __m512i zero = _mm512_setzero_si512();
    const __m512i mask = _mm512_set1_epi64(MASK52);
    __m512i E[ND], D[ND];
#pragma GCC unroll 128
    for (int j = 0; j < ND; j++) {
        E[j] = _mm512_and_si512(_mm512_slli_epi64(a[j], 1), mask);
        D[j] = j ? _mm512_or_si512(E[j], _mm512_srli_epi64(a[j - 1], 51)) : E[0];
    }
    __m512i acc[ND + 1];
#pragma GCC unroll 128
    for (int j = 0; j <= ND; j++) acc[j] = zero;
#pragma GCC unroll 128
    for (int i = 0; i < ND; i++) {
        const __m512i ai = a[i];
        acc[i] = mlo(acc[i], ai, ai);
        acc[i + 1] = mhi(acc[i + 1], ai, ai);
        if (i + 1 < ND) {
            acc[i + 1] = mlo(acc[i + 1], ai, E[i + 1]);
            acc[i + 2] = mhi(acc[i + 2], ai, E[i + 1]);
        }
#pragma GCC unroll 128
        for (int j = i + 2; j < ND; j++) {
            acc[j] = mlo(acc[j], ai, D[j]);
            acc[j + 1] = mhi(acc[j + 1], ai, D[j]);
        }
        const __m512i m = mlo(zero, acc[0], c.k0);
        __m512i t0 = mlo(acc[0], m, c.N[0]);
        const __m512i carry = _mm512_srli_epi64(t0, 52);
        // new acc[j] = acc[j+1] + hi(m n_j) + lo(m n_{j+1})
#pragma GCC unroll 128
        for (int j = 0; j < ND - 1; j++) acc[j] = mlo(mhi(acc[j + 1], m, c.N[j]), m, c.N[j + 1]);
        acc[ND - 1] = mhi(acc[ND], m, c.N[ND - 1]);
        acc[0] = _mm512_add_epi64(acc[0], carry);
        acc[ND] = zero;
    }
    __m512i carry = zero;
#pragma GCC unroll 128
    for (int j = 0; j < ND; j++) {
        __m512i t = _mm512_add_epi64(acc[j], carry);
        r[j] = _mm512_and_si512(t, mask);
        carry = _mm512_srli_epi64(t, 52);
    }
}

// Same with a rolled outer loop (large ND): one fused pass per row, like amm_rolled. Row i (positions shifted by i):
// new acc[j] = acc[j+1] + hi(m n_j) + lo(m n_{j+1}) + product terms landing at position i+j+1, namely
//   j = i-1: lo(a_i^2);  j = i: hi(a_i^2) + lo(a_i E_{i+1});  j = i+1: hi(a_i E_{i+1}) + lo(a_i D_{i+2});
//   j >= i+2: hi(a_i D_j) + lo(a_i D_{j+1})   (with E_ND = D_ND = 0).
// Row 0 adds lo(a_0^2) to acc[0] before m_0 is formed (position 0).
template <int ND>
static inline void ams_rolled(__m512i *r, const __m512i *a, const Mont8<ND> &c) {
    const __m512i zero = _mm512_setzero_si512();
    const __m512i mask = _mm512_set1_epi64(MASK52);
    alignas(64) __m512i E[ND + 2], D[ND + 2], acc[ND + 1];
    for (int j = 0; j < ND; j++) {
        E[j] = _mm512_and_si512(_mm512_slli_epi64(a[j], 1), mask);
        D[j] = j ? _mm512_or_si512(E[j], _mm512_srli_epi64(a[j - 1], 51)) : E[0];
    }
    E[ND] = E[ND + 1] = D[ND] = D[ND + 1] = zero;
    for (int j = 0; j <= ND; j++) acc[j] = zero;
    acc[0] = mlo(acc[0], a[0], a[0]);
    for (int i = 0; i < ND; i++) {
        const __m512i ai = a[i];
        const __m512i m = mlo(zero, acc[0], c.k0);
        const __m512i t0 = mlo(acc[0], m, c.N[0]);
        const __m512i carry = _mm512_srli_epi64(t0, 52);
        int j = 0;
        // j < i-1: reduction terms only
#pragma GCC unroll 8
        for (; j < i - 1; j++) acc[j] = mlo(mhi(acc[j + 1], m, c.N[j]), m, c.N[j + 1]);
        if (i >= 1) {   // j = i-1: + lo(a_i^2)
            acc[j] = mlo(mlo(mhi(acc[j + 1], m, c.N[j]), m, c.N[j + 1]), ai, ai); j++;
        }
        // j = i: + hi(a_i^2) + lo(a_i E_{i+1})
        acc[j] = mlo(mhi(mlo(mhi(acc[j + 1], m, c.N[j]), m, c.N[j + 1]), ai, ai), ai, E[i + 1]); j++;
        if (j < ND) {   // j = i+1: + hi(a_i E_{i+1}) + lo(a_i D_{i+2})
            acc[j] = mlo(mhi(mlo(mhi(acc[j + 1], m, c.N[j]), m, c.N[j + 1]), ai, E[i + 1]), ai, D[i + 2]); j++;
        }
#pragma GCC unroll 8
        for (; j < ND; j++)   // j >= i+2: + hi(a_i D_j) + lo(a_i D_{j+1})
            acc[j] = mlo(mhi(mlo(mhi(acc[j + 1], m, c.N[j]), m, c.N[j + 1]), ai, D[j]), ai, D[j + 1]);
        acc[0] = _mm512_add_epi64(acc[0], carry);
        acc[ND] = zero;
    }
    __m512i carry = zero;
    for (int j = 0; j < ND; j++) {
        __m512i t = _mm512_add_epi64(acc[j], carry);
        r[j] = _mm512_and_si512(t, mask);
        carry = _mm512_srli_epi64(t, 52);
    }
}

// Squaring dispatcher: dedicated squaring if USE_AMS is defined, otherwise the general multiplication.
template <int ND>
static IFMA_INLINE void amsx(__m512i *r, const __m512i *a, const Mont8<ND> &c) {
#ifdef USE_AMS
    if constexpr (ND <= 24) ams<ND>(r, a, c);
    else ams_rolled<ND>(r, a, c);
#else
    ammx<ND>(r, a, a, c);
#endif
}

// Lane-masked doubling: x <- 2x in lanes where bit is set, then renormalize digits (ripple carry).
// Executed unconditionally by the callers (the same instructions for every exponent bit).
template <int ND>
static IFMA_INLINE void dbl_mask(__m512i *x, __mmask8 k) {
    const __m512i mask = _mm512_set1_epi64(MASK52);
    __m512i carry = _mm512_setzero_si512();
#pragma GCC unroll 128
    for (int j = 0; j < ND; j++) {
        __m512i t = _mm512_mask_add_epi64(x[j], k, x[j], x[j]);
        t = _mm512_add_epi64(t, carry);
        x[j] = _mm512_and_si512(t, mask);
        carry = _mm512_srli_epi64(t, 52);
    }
}
