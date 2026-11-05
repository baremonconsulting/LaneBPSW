// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// lanebpsw.cpp — liblanebpsw.so: the C interface of lanebpsw.h, and the functions with string arguments used by the
// Python scripts (hp_run, hp_predicate, hp_td_check, hp_info).
#include "lanebpsw/hashprime.hpp"
#include "lanebpsw.h"

extern "C" {

int lanebpsw_millerrabin25(mpz_srcptr n) {
    if (!lanebpsw_gmp_630()) return mpz_millerrabin(n, 25);
    if (mpz_sgn(n) <= 0) return mpz_millerrabin(n, 25);
    size_t bits = mpz_sizeinbase(n, 2);
    if (bits > 4102) return mpz_millerrabin(n, 25);
    BpswFn f = bpsw_for_bits((int)bits);
    return f ? f(n) : mpz_millerrabin(n, 25);
}

void lanebpsw_nextprime(mpz_ptr rop, mpz_srcptr n) {
    if (!lanebpsw_gmp_630()) { mpz_nextprime(rop, n); return; }
    if (mpz_sizeinbase(n, 2) > 4100) { mpz_nextprime(rop, n); return; }   // beyond the kernels
    static NPStats st;                                                     // statistics only
    g_final_lanes = true;
    nextprime_fast(rop, n, 0, 0, st);
}

int lanebpsw_chiavdf_is_prime(mpz_srcptr n) {
    g_pp.build();
    if (mpz_cmp_ui(n, 14951) <= 0) return chia_small(n);
    size_t bits = mpz_sizeinbase(n, 2);
    ChiaFn f = bits <= 4102 ? chia_for_bits((int)bits) : nullptr;
    return (f ? f(n) : chia_bpsw_ref(n)) != 0;
}

int lanebpsw_hashprime(mpz_ptr rop, const uint8_t *seed, size_t seedlen, int bits, const int *mask, int nmask) {
    if (bits <= 0 || bits > LANEBPSW_MAX_HASHPRIME_BITS || bits % 8 != 0 || nmask < 0 ||
        nmask > LANEBPSW_MAX_MASK_ENTRIES || (nmask > 0 && !mask) || !seed || seedlen == 0 ||
        seedlen > LANEBPSW_MAX_SEED_BYTES)
        return LANEBPSW_INVALID_ARGUMENT;
    for (int i = 0; i < nmask; i++) if (mask[i] < 0 || mask[i] >= bits) return LANEBPSW_INVALID_ARGUMENT;
    g_pp.build();
    HPStats st;
    mpz_t result; mpz_init(result);
    int status = hp_lanes(result, seed, (int)seedlen, bits, mask, nmask, st);
    if (status == LANEBPSW_OK) mpz_set(rop, result);
    mpz_clear(result);
    return status;
}

int lanebpsw_create_discriminant(mpz_ptr rop, const uint8_t *seed, size_t seedlen, int bits) {
    // Validate before evaluating bits - 1 (INT_MIN would otherwise overflow).
    if (bits <= 0 || bits % 8 != 0) return LANEBPSW_INVALID_ARGUMENT;
    const int mask[4] = {0, 1, 2, bits - 1};
    int r = lanebpsw_hashprime(rop, seed, seedlen, bits, mask, 4);
    if (r == 0) mpz_neg(rop, rop);
    return r;
}

} // extern "C"

// ---- functions used by the Python scripts (bench/hp_bench3.py, tests/hp_check.py, scripts/run_hp.sh) ----
extern "C" long hp_run(int mode, const uint8_t *seed, int seedlen, int bits, const int *mask, int nmask,
                       char *outhex, int outcap, double *seconds) {
    if (!seed || seedlen <= 0 || seedlen > LANEBPSW_MAX_SEED_BYTES || bits <= 0 ||
        bits > LANEBPSW_MAX_HASHPRIME_BITS || bits % 8 != 0 || nmask < 0 ||
        nmask > LANEBPSW_MAX_MASK_ENTRIES || (nmask > 0 && !mask) || !outhex || outcap <= 0 || !seconds)
        return LANEBPSW_INVALID_ARGUMENT;
    for (int i = 0; i < nmask; i++) if (mask[i] < 0 || mask[i] >= bits) return LANEBPSW_INVALID_ARGUMENT;
    g_pp.build();
    HPStats st;
    mpz_t out; mpz_init(out);
    double t0 = cpu_s();
    int status = mode == 0 ? hp_scalar(out, seed, seedlen, bits, mask, nmask, st)
                           : hp_lanes(out, seed, seedlen, bits, mask, nmask, st);
    *seconds = cpu_s() - t0;
    if (status != LANEBPSW_OK) { outhex[0] = '\0'; mpz_clear(out); return status; }
    std::string s(mpz_sizeinbase(out, 16) + 2, '\0');
    mpz_get_str(&s[0], 16, out);
    snprintf(outhex, outcap, "%s", s.c_str());
    mpz_clear(out);
    return st.cand;
}
extern "C" int hp_predicate(int mode, const char *hex) {   // for cross-checks from Python
    g_pp.build();
    mpz_t n; mpz_init_set_str(n, hex, 16);
    int r;
    if (mode == 0) r = chia_bpsw_ref(n);
    else { size_t bits = mpz_sizeinbase(n, 2); ChiaFn f = bits <= 4102 ? chia_for_bits((int)bits) : nullptr; r = f ? f(n) : chia_bpsw_ref(n); }
    mpz_clear(n);
    return r;
}
extern "C" long hp_td_check(int bits, int count, unsigned long long seed) {   // fast TD == gcd TD on random odd numbers
    g_pp.build(); ppd_build();
    mpz_t n; mpz_init(n); uint64_t sd = seed; long bad = 0, pass = 0;
    for (int i = 0; i < count; i++) { rand_kbit(n, bits, sd); mpz_setbit(n, 0); bool a = chia_td_gcd(n), b = chia_td(n); bad += a != b; pass += a; }
    mpz_clear(n);
    return bad * 1000000 + pass;
}
extern "C" int hp_info(char *buf, int cap) {
    g_pp.build();
    return snprintf(buf, cap, "pprods=%zu maxp=%u first=%016llx second=%016llx third=%016llx last=%016llx", g_pp.v.size(), g_pp.maxp,
                    (unsigned long long)g_pp.v[0], (unsigned long long)g_pp.v[1], (unsigned long long)g_pp.v[2],
                    (unsigned long long)g_pp.v.back());
}
