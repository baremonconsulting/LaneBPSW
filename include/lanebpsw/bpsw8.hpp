// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// bpsw8.hpp — LaneBPSW: GMP 6.3.0's mpz_millerrabin(n, 25) for one odd n >= 2^64, with all chains of the test
// evaluated together (paper, Sections 2-4):
// the three independent exponentiation chains of GMP 6.3's test run simultaneously in the 8 lanes of the
// vertical IFMA kernel (all lanes use the same modulus n):
//   lanes 0-3 : strong Lucas test, Montgomery-ladder form (V_k, V_{k+1}, Q^k, Q^{k+1}) -> 4 independent
//               products per bit of d_L = (n+1)/2^{s_L}
//   lane  4   : strong base-2 test, left-to-right (square + masked doubling) over d_2 = (n-1)/2^{s_2}
//   lanes 5-6 : extra Miller–Rabin round with GMP's fixed-seed base a, right-to-left (z <- z^2, acc <- acc*z)
//   lane  7   : idle
// One lockstep step = one 8-lane Montgomery multiplication, ~max(log2 d_L, log2 d_2) steps in total.
// The short squaring tails (s-1 steps) and the final comparisons are done with GMP on exported values.
//
// Exactness: the result equals GMP's predicate SPRP_2(n) && SLPRP_{D(n)}(n) && SPRP_{a(n)}(n), where
//   D(n): first |D| in 5,7,11,13,15,17,19,23,25,29,31,35,... with Jacobi (n/|D|) = -1 (GMP 6.3 stronglucas.c),
//         D = |D| if |D| = 1 mod 4 else -|D|, P = 1, Q = (1-D)/4  (D = 5 is GMP's mpn_strongfibo path: same predicate)
//   a(n): 3 + mpz_urandomm(default-seeded gmp_randstate, (n-5)/2)   (GMP 6.3 millerrabin.c)
// Any irregular case (Jacobi 0, D search beyond 255, n < 2^64) is delegated to mpz_millerrabin itself.
#pragma once
#include "lanes.hpp"

struct BpswStats { uint64_t calls = 0, fallback = 0; };
static BpswStats g_bpsw;

// GMP 6.3 D selection (absolute value), 0 = irregular -> caller falls back to GMP
static unsigned gmp_lucas_absD(const mpz_t n) {
    static const unsigned first[6] = {5, 7, 11, 13, 15, 17};
    for (unsigned d : first) {
        int j = mpz_kronecker_ui(n, d);
        if (j == -1) return d;
        if (j == 0) return 0;
    }
    unsigned D = 19, Ddiff = 4;
    while (D < 256) {
        int j = mpz_kronecker_ui(n, D);
        if (j == -1) return D;
        if (j == 0) return 0;
        D += Ddiff; Ddiff = 6 - Ddiff;
    }
    return 0;   // includes perfect squares (no D exists)
}

static inline int mpz_bit_or0(const mp_limb_t *p, int nl, long b) {
    if (b < 0) return 0;
    long li = b >> 6;
    if (li >= nl) return 0;
    return (int)((p[li] >> (b & 63)) & 1);
}

static inline void d52_to_mpz(mpz_t z, const uint64_t *dig, int ND) {
    int nl = (52 * ND + 63) / 64;
    uint64_t limbs[160];
    d52_to_limbs(dig, ND, limbs, nl);
    mpz_import(z, nl, -1, 8, 0, 0, limbs);
}

// strong-test tail with GMP: y = b^d mod n given (fully reduced), s = v2(n-1)
static int sprp_tail(mpz_t y, const mpz_t n, const mpz_t nm1, unsigned long s) {
    if (mpz_cmp_ui(y, 1) == 0 || mpz_cmp(y, nm1) == 0) return 1;
    for (unsigned long i = 1; i < s; i++) {
        mpz_mul(y, y, y); mpz_mod(y, y, n);
        if (mpz_cmp(y, nm1) == 0) return 1;
    }
    return 0;
}

// Montgomery representation v*R mod n (v may be negative)
static inline void to_mont_mpz(mpz_t out, const mpz_t v, const mpz_t n, int ND) {
    mpz_mul_2exp(out, v, 52 * ND);
    mpz_mod(out, out, n);
}

struct BpswScratch {
    mpz_t nm1, np1, d2, dL, a, nm, t, u, w, vd, vd1, qd, x2, xm;
    gmp_randstate_t rs0;
    BpswScratch() {
        mpz_init(nm1); mpz_init(np1); mpz_init(d2); mpz_init(dL); mpz_init(a); mpz_init(nm); mpz_init(t);
        mpz_init(u); mpz_init(w); mpz_init(vd); mpz_init(vd1); mpz_init(qd); mpz_init(x2); mpz_init(xm);
    }
};
static BpswScratch &g_bs_instance() { static BpswScratch v; return v; }
#define g_bs (g_bs_instance())

// Optional instrumentation: last computed Lucas values (normal form) for testing
struct BpswDebug { bool on = false; int base2 = 0, lucas = 0, mr = 0; mpz_t vd, vd1, qd, x2, xm, a; long Q = 0; unsigned absD = 0;
                   BpswDebug() { mpz_init(vd); mpz_init(vd1); mpz_init(qd); mpz_init(x2); mpz_init(xm); mpz_init(a); } };
static BpswDebug &g_bdbg_instance() { static BpswDebug v; return v; }
#define g_bdbg (g_bdbg_instance())

template <int ND>
static int bpsw_lanes(const mpz_t n) {
    g_bpsw.calls++;
    const int k = (int)mpz_sizeinbase(n, 2);
    if (k <= 64 || 52 * ND < k + 6 || mpz_even_p(n)) { g_bpsw.fallback++; return mpz_millerrabin(n, 25); }
    unsigned absD = gmp_lucas_absD(n);
    if (absD == 0) { g_bpsw.fallback++; return mpz_millerrabin(n, 25); }
    const long Dsg = (absD & 2) ? -(long)absD : (long)absD;
    const long Q = (1 - Dsg) / 4;
    BpswScratch &S = g_bs;
    // exponents
    mpz_sub_ui(S.nm1, n, 1);
    unsigned long s2 = mpz_scan1(S.nm1, 0);
    mpz_tdiv_q_2exp(S.d2, S.nm1, s2);
    mpz_add_ui(S.np1, n, 1);
    unsigned long sL = mpz_scan1(S.np1, 0);
    mpz_tdiv_q_2exp(S.dL, S.np1, sL);
    // GMP's extra-round base: 3 + urandomm(default state, (n-5)/2)
    mpz_tdiv_q_2exp(S.nm, n, 1); mpz_sub_ui(S.nm, S.nm, 2);
    {
        gmp_randstate_t rs; gmp_randinit_default(rs);
        mpz_urandomm(S.a, rs, S.nm); mpz_add_ui(S.a, S.a, 3);
        gmp_randclear(rs);
    }
    // lanes: same modulus
    LaneMods<ND> L;
    { mpz_t mods[1]; mpz_init_set(mods[0], n); lanes_setup<ND>(L, mods, true, true); mpz_clear(mods[0]); }
    // initial Montgomery values per lane
    alignas(64) uint64_t init[8][ND];
    mpz_set_ui(S.t, 0); mpz_setbit(S.t, 52 * ND); mpz_mod(S.t, S.t, n);          // 1~  (R mod n)
    mpz_to_d52(S.t, init[1], ND);                                               // V_1 = P = 1
    mpz_to_d52(S.t, init[2], ND);                                               // Q^0 = 1
    mpz_to_d52(S.t, init[4], ND);                                               // x = 1 (base 2)
    mpz_to_d52(S.t, init[6], ND);                                               // acc = 1
    mpz_to_d52(S.t, init[7], ND);                                               // idle
    mpz_mul_2exp(S.u, S.t, 1); mpz_mod(S.u, S.u, n); mpz_to_d52(S.u, init[0], ND);   // V_0 = 2
    mpz_set_si(S.u, Q); to_mont_mpz(S.w, S.u, n, ND); mpz_to_d52(S.w, init[3], ND);  // Q^1 = Q
    to_mont_mpz(S.w, S.a, n, ND); mpz_to_d52(S.w, init[5], ND);                      // z = a
    __m512i st[ND], A[ND], Bv[ND], r[ND];
    load_lanes<ND>(st, init);
    const mp_limb_t *pL = mpz_limbs_read(S.dL); const int nL = (int)mpz_size(S.dL);
    const mp_limb_t *p2 = mpz_limbs_read(S.d2); const int n2 = (int)mpz_size(S.d2);
    const long T = (long)std::max(mpz_sizeinbase(S.dL, 2), mpz_sizeinbase(S.d2, 2));
    const __m512i idxA0 = _mm512_setr_epi64(0, 0, 2, 2, 4, 5, 6, 7), idxA1 = _mm512_setr_epi64(1, 0, 3, 2, 4, 5, 6, 7);
    const __m512i idxB0 = _mm512_setr_epi64(0, 1, 2, 3, 4, 5, 5, 7), idxB1 = _mm512_setr_epi64(1, 1, 3, 3, 4, 5, 5, 7);
    const __m512i idxR0 = _mm512_setr_epi64(0, 1, 2, 3, 4, 5, 6, 7), idxR1 = _mm512_setr_epi64(1, 0, 3, 2, 4, 5, 6, 7);
    const __m512i idxS0 = _mm512_setr_epi64(2, 2, 0, 0, 0, 0, 0, 0), idxS1 = _mm512_setr_epi64(2, 3, 0, 0, 0, 0, 0, 0);
    const __m512i shS0 = _mm512_setr_epi64(1, 0, 0, 0, 0, 0, 0, 0), shS1 = _mm512_setr_epi64(0, 1, 0, 0, 0, 0, 0, 0);
    const __m512i shN0 = _mm512_setr_epi64(2, 1, 0, 0, 0, 0, 0, 0), shN1 = _mm512_setr_epi64(1, 2, 0, 0, 0, 0, 0, 0);
    const __m512i mask52 = _mm512_set1_epi64(MASK52);
    for (long i = 0; i < T; i++) {
        const int bL = mpz_bit_or0(pL, nL, T - 1 - i);   // Lucas ladder bit (MSB first)
        const int b2 = mpz_bit_or0(p2, n2, T - 1 - i);   // base-2 left-to-right bit
        const int bm = mpz_bit_or0(p2, n2, i);           // MR right-to-left bit
        const __mmask8 mb = (__mmask8)(0 - bL);          // all-ones iff bL (branch-free selection)
        const __m512i iA = _mm512_mask_blend_epi64(mb, idxA0, idxA1), iB = _mm512_mask_blend_epi64(mb, idxB0, idxB1);
        const __m512i iR = _mm512_mask_blend_epi64(mb, idxR0, idxR1), iS = _mm512_mask_blend_epi64(mb, idxS0, idxS1);
        const __m512i shS = _mm512_mask_blend_epi64(mb, shS0, shS1), shN = _mm512_mask_blend_epi64(mb, shN0, shN1);
        for (int j = 0; j < ND; j++) { A[j] = _mm512_permutexvar_epi64(iA, st[j]); Bv[j] = _mm512_permutexvar_epi64(iB, st[j]); }
        ammx<ND>(r, A, Bv, L.c);
        const __mmask8 keep = (__mmask8)(0x80 | ((bm ^ 1) << 6));   // lane 7 idle; lane 6 keeps acc unless bit set
        __m512i carry = _mm512_setzero_si512();
        for (int j = 0; j < ND; j++) {
            __m512i v = _mm512_permutexvar_epi64(iR, r[j]);
            __m512i sb = _mm512_maskz_permutexvar_epi64(0x03, iS, st[j]);
            sb = _mm512_sllv_epi64(sb, shS);
            __m512i an = _mm512_maskz_sllv_epi64(0x03, L.c.N[j], shN);
            v = _mm512_add_epi64(_mm512_sub_epi64(v, sb), an);
            v = _mm512_mask_mov_epi64(v, keep, st[j]);
            v = _mm512_add_epi64(v, carry);                 // signed ripple normalization
            st[j] = _mm512_and_si512(v, mask52);
            carry = _mm512_srai_epi64(v, 52);
        }
        dbl_mask<ND>(st, (__mmask8)(b2 << 4));
    }
    // out of Montgomery form: values in [0, n]
    __m512i one[ND];
    one[0] = _mm512_set1_epi64(1);
    for (int j = 1; j < ND; j++) one[j] = _mm512_setzero_si512();
    ammx<ND>(r, st, one, L.c);
    alignas(64) uint64_t out[8][ND];
    store_lanes<ND>(out, r);
    d52_to_mpz(S.vd, out[0], ND); mpz_mod(S.vd, S.vd, n);
    d52_to_mpz(S.vd1, out[1], ND); mpz_mod(S.vd1, S.vd1, n);
    d52_to_mpz(S.qd, out[2], ND); mpz_mod(S.qd, S.qd, n);
    d52_to_mpz(S.x2, out[4], ND); mpz_mod(S.x2, S.x2, n);
    d52_to_mpz(S.xm, out[6], ND); mpz_mod(S.xm, S.xm, n);
    if (g_bdbg.on) {
        mpz_set(g_bdbg.vd, S.vd); mpz_set(g_bdbg.vd1, S.vd1); mpz_set(g_bdbg.qd, S.qd);
        mpz_set(g_bdbg.x2, S.x2); mpz_set(g_bdbg.xm, S.xm); mpz_set(g_bdbg.a, S.a); g_bdbg.Q = Q; g_bdbg.absD = absD;
    }
    // base 2 and MR tails
    int ok2 = sprp_tail(S.x2, n, S.nm1, s2);
    int okm = sprp_tail(S.xm, n, S.nm1, s2);
    // Lucas: U_d == 0  <=>  2 V_{d+1} - P V_d == 0 (gcd(D, n) = 1);  V_d == 0;  V_{d 2^r} == 0, r < s_L
    int okL = 0;
    mpz_mul_2exp(S.t, S.vd1, 1); mpz_sub(S.t, S.t, S.vd); mpz_mod(S.t, S.t, n);
    if (mpz_sgn(S.t) == 0 || mpz_sgn(S.vd) == 0) okL = 1;
    else {
        for (unsigned long rr = 1; rr < sL; rr++) {
            mpz_mul(S.t, S.vd, S.vd); mpz_submul_ui(S.t, S.qd, 2); mpz_mod(S.vd, S.t, n);
            if (mpz_sgn(S.vd) == 0) { okL = 1; break; }
            mpz_mul(S.t, S.qd, S.qd); mpz_mod(S.qd, S.t, n);
        }
    }
    if (g_bdbg.on) { g_bdbg.base2 = ok2; g_bdbg.lucas = okL; g_bdbg.mr = okm; }
    return (ok2 && okL && okm) ? 1 : 0;
}

typedef int (*BpswFn)(const mpz_t n);
static BpswFn bpsw_for_bits(int k) {
    if (k <= 0) return nullptr;
    int64_t need = ((int64_t)k + 57) / 52;
    if (need <= 6) return bpsw_lanes<6>;
    if (need <= 10) return bpsw_lanes<10>;
    if (need <= 15) return bpsw_lanes<15>;
    if (need <= 20) return bpsw_lanes<20>;
    if (need <= 30) return bpsw_lanes<30>;
    if (need <= 40) return bpsw_lanes<40>;
    if (need <= 60) return bpsw_lanes<60>;
    if (need <= 79) return bpsw_lanes<79>;
    return nullptr;
}
