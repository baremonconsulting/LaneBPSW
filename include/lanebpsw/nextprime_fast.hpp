// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// nextprime_fast.hpp — next prime with the output of GMP 6.3.0's mpz_nextprime for every input (paper, Section 3.3):
// returns the first survivor of GMP's sieve above n that mpz_millerrabin(p, 25) accepts (BPSW + 1 fixed-seed MR).
// Sieve survivors are tested speculatively 8 at a time (base-2 Fermat in IFMA lanes) and confirmed
// in ascending order: a lane is accepted only if all lower survivors were rejected.
// Correctness: Fermat-2 failure  =>  strong-2 failure  =>  mpz_millerrabin()==0, so skipping such
// candidates never changes the answer; every Fermat-passing candidate is re-checked by GMP's own test.
#pragma once
#include "fermat8.hpp"
#include "sieve.hpp"
#include "bpsw8.hpp"

// final confirmation: false = GMP mpz_millerrabin(p,25); true = lane-packed BPSW (bpsw8.hpp, GMP-identical predicate)
static bool g_final_lanes = false;
static inline int final_check(const mpz_t p) {
    if (g_final_lanes) { BpswFn f = bpsw_for_bits((int)mpz_sizeinbase(p, 2)); if (f) return f(p); }
    return mpz_millerrabin(p, 25);
}

struct NPStats {
    uint64_t queries = 0, blocks = 0, batches = 0, lanes = 0, waste = 0, finals = 0, scalar_fallback = 0;
    double t_res = 0, t_mark = 0, t_test = 0, t_final = 0;
};

static PrimeTable &g_pt_instance() { static PrimeTable v; return v; }
#define g_pt (g_pt_instance())

typedef uint8_t (*Fermat8Fn)(const mp_limb_t *const cand[8], int nl, int k);
static Fermat8Fn fermat_for_bits(int kc) {
    int nd = nd_for_bits(kc);
    if (nd <= 5) return fermat2_8<5>;
    if (nd <= 6) return fermat2_8<6>;     // 257-308 bits, e.g. chiavdf's 264-bit proof prime
    if (nd <= 10) return fermat2_8<10>;
    if (nd <= 15) return fermat2_8<15>;
    if (nd <= 20) return fermat2_8<20>;
    if (nd <= 30) return fermat2_8<30>;
    if (nd <= 40) return fermat2_8<40>;
    if (nd <= 60) return fermat2_8<60>;
    if (nd <= 79) return fermat2_8<79>;
    return nullptr;
}

// GMP 6.3 mpz/nextprime.c: sieve_limit = isqrt(floor(nbits^5 / 124^2)) for nbits < 12818, else 150000001,
// where nbits = bit length of n+1; findnext() then sieves with the odd primes 3..limit EXCEPT the largest
// one (its loop runs prime_limit = #primes in [5,limit] times starting from 3).
static uint32_t gmp_sieve_limit(size_t nbits) {
    if (nbits >= 12818) return 150000001;
    mpz_t t; mpz_init(t);
    mpz_ui_pow_ui(t, nbits, 5); mpz_tdiv_q_ui(t, t, 124 * 124); mpz_sqrt(t, t);
    uint32_t r = (uint32_t)mpz_get_ui(t); mpz_clear(t);
    return r;
}

// B: sieve bound (0 = exactly GMP's prime set => candidate sets identical to mpz_nextprime),
// L: odd slots per block (GMP uses 5*nbits)
static void nextprime_fast(mpz_t out, const mpz_t n, uint32_t B, uint32_t L, NPStats &st, bool timing = false) {
    st.queries++;
    if (mpz_sizeinbase(n, 2) < 64 || mpz_cmp_ui(n, 2) < 0) { mpz_nextprime(out, n); return; }
    bool gmp_exact = (B == 0);
    if (gmp_exact) {
        mpz_t np1; mpz_init(np1); mpz_add_ui(np1, n, 1);
        size_t nbits = mpz_sizeinbase(np1, 2);
        mpz_clear(np1);
        // GMP sieves with the first nbits/2 odd primes when nbits/2 <= 100 (nextprime.c, findnext); that
        // candidate set is not reproduced here, so such inputs go to GMP itself.
        if (nbits <= 201) { mpz_nextprime(out, n); return; }
        B = gmp_sieve_limit(nbits);
        if (L == 0) L = 5 * (uint32_t)nbits;
    }
    if (g_pt.B < B) g_pt.build(B);
    mpz_t A, c[8];
    mpz_init(A);
    for (int l = 0; l < 8; l++) mpz_init(c[l]);
    mpz_add_ui(A, n, mpz_odd_p(n) ? 2 : 1);
    std::vector<uint32_t> res(g_pt.count_upto(B));
    std::vector<uint8_t> comp(L);
    std::vector<uint32_t> surv;
    surv.reserve(L);
    bool found = false;
    while (!found) {
        st.blocks++;
        double t0 = timing ? cpu_s() : 0;
        int nl = mpz_size(A);
        size_t np = residues(mpz_limbs_read(A), nl, g_pt, B, res.data());
        if (gmp_exact && np > 0) np--;           // GMP never uses the largest prime <= limit
        double t1 = timing ? cpu_s() : 0;
        mark_plain(comp.data(), L, res.data(), np, g_pt, false, 0);
        surv.clear();
        for (uint32_t i = 0; i < L; i++) if (!comp[i]) surv.push_back(i);
        double t2 = timing ? cpu_s() : 0;
        if (timing) { st.t_res += t1 - t0; st.t_mark += t2 - t1; }
        for (size_t s = 0; s < surv.size() && !found; s += 8) {
            int cnt = (int)std::min<size_t>(8, surv.size() - s);
            for (int l = 0; l < 8; l++) mpz_add_ui(c[l], A, 2 * (uint64_t)surv[s + std::min(l, cnt - 1)]);
            int kc = (int)mpz_sizeinbase(c[0], 2);
            bool same = (int)mpz_sizeinbase(c[7], 2) == kc;
            double t3 = timing ? cpu_s() : 0;
            uint8_t mask;
            Fermat8Fn f = same ? fermat_for_bits(kc) : nullptr;
            if (f) {
                const mp_limb_t *cp[8];
                for (int l = 0; l < 8; l++) cp[l] = mpz_limbs_read(c[l]);
                mask = f(cp, (int)mpz_size(c[0]), kc);
                st.batches++;
            } else {          // mixed bit lengths (crossing 2^k) or unsupported size: scalar path
                mask = 0xFF;  // let the exact GMP test decide every lane
                st.scalar_fallback++;
            }
            double t4 = timing ? cpu_s() : 0;
            if (timing) st.t_test += t4 - t3;
            for (int l = 0; l < cnt; l++) {
                st.lanes++;
                if ((mask >> l) & 1) {
                    st.finals++;
                    double t5 = timing ? cpu_s() : 0;
                    int ok = final_check(c[l]);
                    if (timing) st.t_final += cpu_s() - t5;
                    if (ok) { mpz_set(out, c[l]); found = true; st.waste += cnt - 1 - l; break; }
                }
            }
        }
        mpz_add_ui(A, A, 2 * (uint64_t)L);
    }
    mpz_clear(A);
    for (int l = 0; l < 8; l++) mpz_clear(c[l]);
}
