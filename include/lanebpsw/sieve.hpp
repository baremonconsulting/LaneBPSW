// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// sieve.hpp — table of odd primes, residues of a big number modulo all of them (GMP's mpn_mod_1s_4p on
// products of primes below 2^62), and marking of the multiples in a block of odd numbers.
#pragma once
#include "common.hpp"
#include <vector>
#include <algorithm>

struct PrimeGroup {
    uint64_t prod;        // product of primes p[first..end) (< 2^62)
    uint64_t prod_shifted;
    uint64_t cps[7];
    uint32_t first, end;
};

struct PrimeTable {
    uint32_t B = 0;                 // all odd primes 3..B
    std::vector<uint32_t> p;
    std::vector<uint64_t> fm;       // Lemire fastmod constant: floor(2^64 / p) + 1
    std::vector<PrimeGroup> g;
    void build(uint32_t bound) {
        B = bound;
        std::vector<uint8_t> c(bound + 1, 0);
        p.clear();
        for (uint64_t i = 3; i <= bound; i += 2) {
            if (!c[i]) {
                p.push_back((uint32_t)i);
                for (uint64_t j = i * i; j <= bound; j += 2 * i) c[j] = 1;
            }
        }
        fm.resize(p.size());
        for (size_t i = 0; i < p.size(); i++) fm[i] = UINT64_C(0xFFFFFFFFFFFFFFFF) / p[i] + 1;
        g.clear();
        for (uint32_t i = 0; i < p.size();) {
            PrimeGroup G;
            G.prod = 1; G.first = i;
            while (i < p.size() && G.prod <= (UINT64_C(1) << 62) / p[i]) G.prod *= p[i++];
            G.end = i;
            __gmpn_mod_1s_4p_cps(G.cps, G.prod);
            G.prod_shifted = G.prod << G.cps[1];
            g.push_back(G);
        }
    }
    size_t count_upto(uint32_t b) const { return std::upper_bound(p.begin(), p.end(), b) - p.begin(); }
};

// residues res[i] = A mod p[i] for all primes p[i] <= bound (A given by limbs).
static inline size_t residues(const mp_limb_t *a, int nl, const PrimeTable &T, uint32_t bound, uint32_t *res) {
    size_t np = T.count_upto(bound);
    for (const PrimeGroup &G : T.g) {
        if (G.first >= np) break;
        uint64_t R = __gmpn_mod_1s_4p(a, nl, G.prod_shifted, G.cps);
        uint32_t e = std::min<uint32_t>(G.end, (uint32_t)np);
        for (uint32_t i = G.first; i < e; i++) res[i] = (uint32_t)(R % T.p[i]);
    }
    return np;
}

// Plain search: slots i=0..L-1 represent odd A+2i. comp[i]=1 if divisible by some p<=bound
// (and A+2i != p). Returns number of primes used.
static inline void mark_plain(uint8_t *comp, uint32_t L, const uint32_t *res, size_t np, const PrimeTable &T,
                              bool A_small, uint64_t A_lo) {
    memset(comp, 0, L);
    for (size_t k = 0; k < np; k++) {
        uint32_t q = T.p[k], r = res[k];
        // A + 2i == 0 (mod q)  ->  i == (q - r) * (q+1)/2 (mod q)
        uint32_t d = r ? q - r : 0;
        uint64_t i0 = (d & 1) ? ((uint64_t)d + q) >> 1 : d >> 1;
        if (A_small && A_lo + 2 * i0 == q) i0 += q;   // do not cross out q itself
        for (uint64_t i = i0; i < L; i += q) comp[i] = 1;
    }
}
