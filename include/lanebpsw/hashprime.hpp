// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// hashprime.hpp — HashPrime(seed, length, bitmask) of Chia's chiavdf (src/proof_common.h, release 1.2.0 and commit
// f41299e): hash-and-check prime generation. SHA-256 over an incrementing big-endian "sprout", concatenated to `length`
// bits, mask bits and bit 0 set, accepted if is_prime_bpsw accepts. Two implementations with the same output:
//   hp_scalar  every candidate through the scalar implementation of chiavdf's predicate (chia_bpsw_ref)
//   hp_lanes   candidates in counter order -> chiavdf's trial division -> 8-lane base-2 Fermat filter -> the
//              Fermat-passing candidates confirmed in order with the lane-packed predicate (chiabpsw8.hpp)
// Exactness: a candidate that fails the Fermat test fails the strong base-2 test (C2), so chiavdf rejects it too; the
// candidates are examined in counter order (paper, Sections 3.3 and 4).
#pragma once
#ifndef OPENSSL_SUPPRESS_DEPRECATED
#define OPENSSL_SUPPRESS_DEPRECATED       // SHA256_Init/Update/Final: deprecated in OpenSSL 3, but much faster per call
#endif
#include "chiabpsw8.hpp"
#include "nextprime_fast.hpp"     // fermat_for_bits
#include "lanebpsw.h"
#include <openssl/sha.h>
#include <vector>
#include <string>

static void hp_next_candidate(mpz_t p, std::vector<uint8_t> &sprout, int length_bits, const int *mask, int nmask,
                              std::vector<uint8_t> &blob) {
    const size_t nbytes = (size_t)length_bits / 8;
    blob.clear();
    while (blob.size() < nbytes) {
        for (size_t i = sprout.size(); i-- > 0;) { sprout[i]++; if (sprout[i]) break; }
        uint8_t h[32];
        SHA256_CTX c; SHA256_Init(&c); SHA256_Update(&c, sprout.data(), sprout.size()); SHA256_Final(h, &c);   // 70 ns per call against 430-470 ns for the one-shot SHA256() of OpenSSL 3.0.13 (bench/bench_sha.c)
        size_t take = std::min<size_t>(32, nbytes - blob.size());
        blob.insert(blob.end(), h, h + take);
    }
    mpz_import(p, nbytes, 1, 1, 1, 0, blob.data());
    for (int i = 0; i < nmask; i++) mpz_setbit(p, mask[i]);
    mpz_setbit(p, 0);
}

struct HPStats { long cand = 0, td_pass = 0, batches = 0, finals = 0; };

static int hp_scalar(mpz_t out, const uint8_t *seed, int seedlen, int bits, const int *mask, int nmask, HPStats &st) {
    std::vector<uint8_t> initial(seed, seed + seedlen), sprout(initial), blob;
    mpz_t p; mpz_init(p);
    for (;;) {
        hp_next_candidate(p, sprout, bits, mask, nmask, blob);
        st.cand++;
        if (chia_bpsw_ref(p)) { mpz_set(out, p); mpz_clear(p); return LANEBPSW_OK; }
        if (sprout == initial) { mpz_clear(p); return LANEBPSW_SEARCH_EXHAUSTED; }
    }
}

static bool hp_confirm_survivors(mpz_t out, mpz_t q[8], int nq, Fermat8Fn ferm, int bits,
                                 HPStats &st, bool full_batch) {
    uint8_t m = 0xFF;
    bool same = full_batch;
    for (int l = 1; l < nq; l++)
        same &= mpz_size(q[l]) == mpz_size(q[0]) && mpz_sizeinbase(q[l], 2) == mpz_sizeinbase(q[0], 2);
    if (full_batch && ferm && same && (int)mpz_sizeinbase(q[0], 2) == bits) {
        const mp_limb_t *cp[8];
        for (int l = 0; l < 8; l++) cp[l] = mpz_limbs_read(q[l]);
        m = ferm(cp, (int)mpz_size(q[0]), bits);
        st.batches++;
    }
    for (int l = 0; l < nq; l++) {
        if (full_batch && !((m >> l) & 1)) continue;
        st.finals++;
        size_t qbits = mpz_sizeinbase(q[l], 2);
        ChiaFn candidate_fn = qbits <= 4102 ? chia_for_bits((int)qbits) : nullptr;
        int ok = candidate_fn ? candidate_fn(q[l]) : chia_bpsw_ref(q[l]);
        if (ok) { mpz_set(out, q[l]); return true; }
    }
    return false;
}

static int hp_lanes(mpz_t out, const uint8_t *seed, int seedlen, int bits, const int *mask, int nmask, HPStats &st) {
    std::vector<uint8_t> initial(seed, seed + seedlen), sprout(initial), blob;
    Fermat8Fn ferm = fermat_for_bits(bits);
    mpz_t q[8]; for (int l = 0; l < 8; l++) mpz_init(q[l]);
    int nq = 0;
    for (;;) {
        hp_next_candidate(q[nq], sprout, bits, mask, nmask, blob);
        st.cand++;
        bool orbit_end = sprout == initial;
        // chiavdf 1.2.1+ deliberately skips product trial division for this
        // finite domain because n may itself be one of the small primes.
        if (mpz_cmp_ui(q[nq], 14951) <= 0 || chia_td(q[nq])) {
            st.td_pass++;
            nq++;
        }
        if (nq == 8) {
            if (hp_confirm_survivors(out, q, nq, ferm, bits, st, true)) {
                for (int l = 0; l < 8; l++) mpz_clear(q[l]);
                return LANEBPSW_OK;
            }
            nq = 0;
        }
        if (orbit_end) {
            // Trial-division survivors that do not fill a complete SIMD batch still
            // have to be confirmed in original candidate order before exhaustion.
            bool found = nq > 0 && hp_confirm_survivors(out, q, nq, ferm, bits, st, false);
            for (int l = 0; l < 8; l++) mpz_clear(q[l]);
            return found ? LANEBPSW_OK : LANEBPSW_SEARCH_EXHAUSTED;
        }
    }
}
