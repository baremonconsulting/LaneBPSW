// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// chiabpsw8.hpp — the acceptance predicate of Chia's chiavdf `is_prime_bpsw` (src/primetest.h; the same for n > 14951 in release 1.2.0 and at commit f41299e):
//   (1) gcd trial division with `pprods` (greedy 64-bit products of consecutive odd primes 3..14951, 379 entries);
//       for operands with fewer than 16 limbs only the first 5 + (379-5)*limbs^3/4096 products are used;
//   (2) strong base-2 Miller–Rabin;
//   (3) Lucas-V probable-prime test V_{n+1} == 2Q (mod n) with Baillie–Fiori–Wagstaff "Method A*" parameters:
//       first D in 5,-7,9,-11,13,... (|D| < 1000) with Jacobi (D/n) = -1; D = 5 -> P = Q = 5, else P = 1, Q = (1-D)/4.
// Evaluated with intra-test lane packing on the 8-lane IFMA kernel: lanes 0-3 Lucas ladder (V_k, V_{k+1}, Q^k, Q^{k+1})
// over e = n+1, lane 4 base-2 left-to-right, lane 5 the product P*Q^k (needed in the same step for V_{2k+1}), lanes 6-7 idle.
// A scalar GMP reference of the same predicate is provided for cross-checks.
#pragma once
#include "bpsw8.hpp"
#include <vector>

struct PProds {
    std::vector<uint64_t> v;
    unsigned maxp = 0;
    void build() {
        if (!v.empty()) return;
        std::vector<char> comp(15001, 0);   // chiavdf: odd primes 3..14951
        unsigned __int128 prod = 1;
        uint64_t cur = 1;
        for (unsigned p = 3; p <= 14951; p += 2) {
            if (comp[p]) continue;
            for (unsigned q = p * p; q <= 15000; q += 2 * p) comp[q] = 1;
            prod = (unsigned __int128)cur * p;
            if (prod >> 64) { v.push_back(cur); cur = p; } else cur = (uint64_t)prod;
            maxp = p;
        }
        v.push_back(cur);
    }
};
static PProds &g_pp_instance() { static PProds v; return v; }
#define g_pp (g_pp_instance())
static bool g_chia_skip_td = false;   // tests only: evaluate the Lucas-V part on numbers with small factors

static int chia_npprods(int limbs) {
    int n = (int)g_pp.v.size();
    const int min_pprods = 5;
    if (limbs < 16) {
        int size_cube = limbs * limbs * limbs;
        n = min_pprods + (n - min_pprods) * size_cube / (16 * 16 * 16);
    }
    return n;
}
// chiavdf trial division semantics: false iff gcd(n, pprods[i]) != 1 for some used i.
// Reference (as chiavdf): mpn_gcd_1 per product (~180 ns each). Fast equivalent: r = n mod pprods[i] (mpn_mod_1, ~37 ns),
// then p | n  <=>  p | r, tested for every prime p of the product by r*p^{-1} mod 2^64 <= floor((2^64-1)/p).
static bool chia_td_gcd(const mpz_t n) {
    int limbs = (int)mpz_size(n), np = chia_npprods(limbs);
    const mp_limb_t *d = mpz_limbs_read(n);
    for (int i = 0; i < np; i++) if (mpn_gcd_1(d, limbs, g_pp.v[i]) != 1) return false;
    return true;
}
struct PPDiv { std::vector<uint32_t> start; std::vector<uint64_t> inv, lim; };
static PPDiv &g_ppd_instance() { static PPDiv v; return v; }
#define g_ppd (g_ppd_instance())
static void ppd_build() {
    if (!g_ppd.start.empty()) return;
    g_pp.build();
    std::vector<char> comp(14952, 0);
    std::vector<uint32_t> ps;
    for (unsigned p = 3; p <= 14951; p += 2) { if (comp[p]) continue; for (unsigned q = p * p; q <= 14951; q += 2 * p) comp[q] = 1; ps.push_back(p); }
    size_t k = 0;
    for (size_t i = 0; i < g_pp.v.size(); i++) {
        g_ppd.start.push_back((uint32_t)k);
        unsigned __int128 prod = 1;
        while (k < ps.size()) {
            unsigned __int128 t = prod * ps[k];
            if (t > g_pp.v[i] || g_pp.v[i] % ps[k] != 0) break;
            prod = t;
            uint64_t p = ps[k], x = p;
            for (int it = 0; it < 6; it++) x *= 2 - p * x;      // p^{-1} mod 2^64
            g_ppd.inv.push_back(x); g_ppd.lim.push_back(~0ULL / p);
            k++;
        }
    }
    g_ppd.start.push_back((uint32_t)k);
}
static bool chia_td(const mpz_t n) {
    ppd_build();
    int limbs = (int)mpz_size(n), np = chia_npprods(limbs);
    const mp_limb_t *d = mpz_limbs_read(n);
    for (int i = 0; i < np; i++) {
        uint64_t r = mpn_mod_1(d, limbs, g_pp.v[i]);
        for (uint32_t j = g_ppd.start[i]; j < g_ppd.start[i + 1]; j++) if (r * g_ppd.inv[j] <= g_ppd.lim[j]) return false;
    }
    return true;
}
static bool chia_find_pq(const mpz_t n, long &P, long &Q) {
    for (long ad = 5; ad < 1000; ad += 2) {
        long d = (ad % 4 == 1) ? ad : -ad;
        if (mpz_si_kronecker(d, n) == -1) {
            if (d == 5) { P = 5; Q = 5; } else { P = 1; Q = (1 - d) / 4; }
            return true;
        }
    }
    return false;
}

// ---- scalar GMP reference (independent straightforward code) ----
static bool chia_sprp2_ref(const mpz_t n) {
    mpz_t nm1, d, y, two; mpz_init(nm1); mpz_init(d); mpz_init(y); mpz_init_set_ui(two, 2);
    mpz_sub_ui(nm1, n, 1); unsigned long s = mpz_scan1(nm1, 0); mpz_tdiv_q_2exp(d, nm1, s);
    mpz_powm(y, two, d, n);
    bool ok = mpz_cmp_ui(y, 1) == 0 || mpz_cmp(y, nm1) == 0;
    for (unsigned long r = 1; r < s && !ok; r++) { mpz_mul(y, y, y); mpz_mod(y, y, n); if (mpz_cmp(y, nm1) == 0) ok = true; }
    mpz_clear(nm1); mpz_clear(d); mpz_clear(y); mpz_clear(two);
    return ok;
}
// V_e mod n by the classic left-to-right U/V doubling method (independent of the lane ladder):
//   U_2k = U_k V_k, V_2k = V_k^2 - 2Q^k;  U_{2k+1} = (P U_2k + V_2k)/2, V_{2k+1} = (D U_2k + P V_2k)/2  (mod n)
static void lucas_v_ref(mpz_t V, const mpz_t e, long P, long Q, const mpz_t n) {
    long D = P * P - 4 * Q;
    mpz_t u, v, qk, t, inv2; mpz_init_set_ui(u, 1); mpz_init_set_si(v, P); mpz_init_set_si(qk, Q); mpz_init(t); mpz_init(inv2);
    mpz_mod(v, v, n); mpz_mod(qk, qk, n);
    mpz_set_ui(inv2, 2); mpz_invert(inv2, inv2, n);
    for (long b = (long)mpz_sizeinbase(e, 2) - 2; b >= 0; b--) {
        mpz_mul(u, u, v); mpz_mod(u, u, n);                                  // U_2k
        mpz_mul(t, v, v); mpz_submul_ui(t, qk, 2); mpz_mod(v, t, n);         // V_2k
        mpz_mul(qk, qk, qk); mpz_mod(qk, qk, n);                             // Q^2k
        if (mpz_tstbit(e, b)) {
            mpz_t u2; mpz_init(u2);
            mpz_mul_si(u2, u, P); mpz_add(u2, u2, v); mpz_mul(u2, u2, inv2); mpz_mod(u2, u2, n);      // U_{2k+1}
            mpz_mul_si(t, u, D); mpz_t pv; mpz_init(pv); mpz_mul_si(pv, v, P); mpz_add(t, t, pv); mpz_clear(pv);
            mpz_mul(t, t, inv2); mpz_mod(v, t, n);                                                    // V_{2k+1}
            mpz_set(u, u2); mpz_clear(u2);
            mpz_mul_si(qk, qk, Q); mpz_mod(qk, qk, n);                                                // Q^{2k+1}
        }
    }
    mpz_set(V, v);
    mpz_clear(u); mpz_clear(v); mpz_clear(qk); mpz_clear(t); mpz_clear(inv2);
}
// chiavdf 1.2.1 and later do not apply the small-prime product trial division
// to n <= 14951. Answer this finite domain exactly before the common predicate.
static int chia_small(const mpz_t n) {
    if (mpz_cmp_ui(n, 2) < 0) return 0;
    if (mpz_cmp_ui(n, 8) < 0) return (0xac >> mpz_get_ui(n)) & 1;
    if (mpz_even_p(n)) return 0;
    if (!chia_sprp2_ref(n)) return 0;
    long P, Q;
    if (!chia_find_pq(n, P, Q)) return 0;
    mpz_t e, V; mpz_init(e); mpz_init(V);
    mpz_add_ui(e, n, 1);
    lucas_v_ref(V, e, P, Q, n);
    if (Q < 0 && mpz_sgn(V) != 0) mpz_sub(V, V, n);
    int ok = mpz_cmp_si(V, 2 * Q) == 0;
    mpz_clear(e); mpz_clear(V);
    return ok;
}
static int chia_bpsw_ref(const mpz_t n) {
    if (mpz_cmp_ui(n, 14951) <= 0) return chia_small(n) ? 2 : 0;
    if (mpz_even_p(n)) return mpz_cmp_ui(n, 2) == 0 ? 2 : 0;
    if (!chia_td(n)) return 0;
    if (!chia_sprp2_ref(n)) return 0;
    long P, Q; if (!chia_find_pq(n, P, Q)) return 0;
    mpz_t e, V, tq; mpz_init(e); mpz_init(V); mpz_init(tq);
    mpz_add_ui(e, n, 1); lucas_v_ref(V, e, P, Q, n);
    mpz_set_si(tq, 2 * Q); mpz_mod(tq, tq, n);
    int ok = mpz_cmp(V, tq) == 0;
    mpz_clear(e); mpz_clear(V); mpz_clear(tq);
    return ok ? 2 : 0;
}

// ---- lane-packed version (strong base-2 + Lucas-V); trial division and parameter search with GMP ----
struct ChiaDebug { int sprp2 = -1, lucv = -1; long P = 0, Q = 0; };
static ChiaDebug g_cdbg;

template <int ND>
static int chia_bpsw_lanes(const mpz_t n) {
    if (!g_chia_skip_td && mpz_cmp_ui(n, 14951) <= 0) return chia_small(n) ? 2 : 0;
    if (mpz_even_p(n)) return mpz_cmp_ui(n, 2) == 0 ? 2 : 0;
    if (!g_chia_skip_td && !chia_td(n)) return 0;
    long P, Q;
    if (!chia_find_pq(n, P, Q)) return 0;     // chiavdf: no D below 1000 -> not prime (checked after base 2; same result)
    const int k = (int)mpz_sizeinbase(n, 2);
    if (52 * ND < k + 6 || k < 3) return chia_bpsw_ref(n);
    BpswScratch &S = g_bs;
    mpz_sub_ui(S.nm1, n, 1);
    unsigned long s2 = mpz_scan1(S.nm1, 0);
    mpz_tdiv_q_2exp(S.d2, S.nm1, s2);
    mpz_add_ui(S.np1, n, 1);                   // Lucas exponent e = n + 1 (full)
    LaneMods<ND> L;
    { mpz_t mods[1]; mpz_init_set(mods[0], n); lanes_setup<ND>(L, mods, true, true); mpz_clear(mods[0]); }
    alignas(64) uint64_t init[8][ND];
    mpz_set_ui(S.t, 0); mpz_setbit(S.t, 52 * ND); mpz_mod(S.t, S.t, n);            // 1~
    mpz_to_d52(S.t, init[2], ND); mpz_to_d52(S.t, init[4], ND);
    mpz_to_d52(S.t, init[6], ND); mpz_to_d52(S.t, init[7], ND);
    mpz_mul_2exp(S.u, S.t, 1); mpz_mod(S.u, S.u, n); mpz_to_d52(S.u, init[0], ND); // V_0 = 2
    mpz_set_si(S.u, P); to_mont_mpz(S.w, S.u, n, ND); mpz_to_d52(S.w, init[1], ND); // V_1 = P
    mpz_to_d52(S.w, init[5], ND);                                                     // lane 5 constant: P~
    mpz_set_si(S.u, Q); to_mont_mpz(S.w, S.u, n, ND); mpz_to_d52(S.w, init[3], ND); // Q^1
    __m512i st[ND], A[ND], Bv[ND], r[ND];
    load_lanes<ND>(st, init);
    const mp_limb_t *pe = mpz_limbs_read(S.np1); const int ne = (int)mpz_size(S.np1);
    const mp_limb_t *p2 = mpz_limbs_read(S.d2); const int n2 = (int)mpz_size(S.d2);
    const long T = (long)mpz_sizeinbase(S.np1, 2);
    const __m512i idxA0 = _mm512_setr_epi64(0, 0, 2, 2, 4, 2, 6, 7), idxA1 = _mm512_setr_epi64(1, 0, 3, 2, 4, 2, 6, 7);
    const __m512i idxB0 = _mm512_setr_epi64(0, 1, 2, 3, 4, 5, 6, 7), idxB1 = _mm512_setr_epi64(1, 1, 3, 3, 4, 5, 6, 7);
    const __m512i idxR0 = _mm512_setr_epi64(0, 1, 2, 3, 4, 5, 6, 7), idxR1 = _mm512_setr_epi64(1, 0, 3, 2, 4, 5, 6, 7);
    const __m512i idx23 = _mm512_setr_epi64(2, 3, 0, 0, 0, 0, 0, 0), idx55 = _mm512_set1_epi64(5);
    const __m512i mask52 = _mm512_set1_epi64(MASK52);
    for (long i = 0; i < T; i++) {
        const int bL = mpz_bit_or0(pe, ne, T - 1 - i);
        const int b2 = mpz_bit_or0(p2, n2, T - 1 - i);
        const __mmask8 mb = (__mmask8)(0 - bL);
        const __m512i iA = _mm512_mask_blend_epi64(mb, idxA0, idxA1), iB = _mm512_mask_blend_epi64(mb, idxB0, idxB1);
        const __m512i iR = _mm512_mask_blend_epi64(mb, idxR0, idxR1);
        const __mmask8 mst = (__mmask8)(1u << bL), mr = (__mmask8)(2u >> bL);   // lane taking 2*S_{2|3} / lane taking r5
        for (int j = 0; j < ND; j++) { A[j] = _mm512_permutexvar_epi64(iA, st[j]); Bv[j] = _mm512_permutexvar_epi64(iB, st[j]); }
        ammx<ND>(r, A, Bv, L.c);
        __m512i carry = _mm512_setzero_si512();
        for (int j = 0; j < ND; j++) {
            __m512i v = _mm512_permutexvar_epi64(iR, r[j]);
            __m512i s1 = _mm512_slli_epi64(_mm512_maskz_permutexvar_epi64(mst, idx23, st[j]), 1);
            __m512i s5 = _mm512_maskz_permutexvar_epi64(mr, idx55, r[j]);
            __m512i an = _mm512_or_si512(_mm512_maskz_slli_epi64(mst, L.c.N[j], 2), _mm512_maskz_slli_epi64(mr, L.c.N[j], 1));
            v = _mm512_add_epi64(_mm512_sub_epi64(_mm512_sub_epi64(v, s1), s5), an);
            v = _mm512_mask_mov_epi64(v, 0xE0, st[j]);        // lanes 5-7 constant
            v = _mm512_add_epi64(v, carry);
            st[j] = _mm512_and_si512(v, mask52);
            carry = _mm512_srai_epi64(v, 52);
        }
        dbl_mask<ND>(st, (__mmask8)(b2 << 4));
    }
    __m512i one[ND];
    one[0] = _mm512_set1_epi64(1);
    for (int j = 1; j < ND; j++) one[j] = _mm512_setzero_si512();
    ammx<ND>(r, st, one, L.c);
    alignas(64) uint64_t out[8][ND];
    store_lanes<ND>(out, r);
    d52_to_mpz(S.vd, out[0], ND); mpz_mod(S.vd, S.vd, n);       // V_{n+1}
    d52_to_mpz(S.x2, out[4], ND); mpz_mod(S.x2, S.x2, n);       // 2^{d2}
    int ok2 = sprp_tail(S.x2, n, S.nm1, s2);
    mpz_set_si(S.t, 2 * Q); mpz_mod(S.t, S.t, n);
    int okv = mpz_cmp(S.vd, S.t) == 0;
    g_cdbg.sprp2 = ok2; g_cdbg.lucv = okv; g_cdbg.P = P; g_cdbg.Q = Q;
    return (ok2 && okv) ? 2 : 0;
}
typedef int (*ChiaFn)(const mpz_t n);
static ChiaFn chia_for_bits(int k) {
    if (k <= 0) return nullptr;
    int64_t need = ((int64_t)k + 57) / 52;
    if (need <= 6) return chia_bpsw_lanes<6>;
    if (need <= 10) return chia_bpsw_lanes<10>;
    if (need <= 15) return chia_bpsw_lanes<15>;
    if (need <= 20) return chia_bpsw_lanes<20>;
    if (need <= 30) return chia_bpsw_lanes<30>;
    if (need <= 40) return chia_bpsw_lanes<40>;
    if (need <= 60) return chia_bpsw_lanes<60>;
    if (need <= 79) return chia_bpsw_lanes<79>;
    return nullptr;
}
