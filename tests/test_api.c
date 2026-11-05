/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2026 Jan Suchánek, Baremon s.r.o.
 * test_api.c — the C interface of liblanebpsw.so against the reference functions:
 *   lanebpsw_millerrabin25     = mpz_millerrabin(n, 25)   1 to 3000, random odd and even numbers of 2 to 4200 bits, primes
 *   lanebpsw_nextprime         = mpz_nextprime            0, 1, small, negative and random n of 1 to 2100 bits, aliasing
 *   lanebpsw_chiavdf_is_prime  = primality for n <= 70000 (what chiavdf 1.2.1+ returns there), and = the scalar
 *                                version (hp_predicate) on random odd numbers up to 2048 bits
 *   lanebpsw_create_discriminant / lanebpsw_hashprime: consistency with the scalar HashPrime (hp_run)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "lanebpsw.h"

/* functions of liblanebpsw.so with string arguments (see src/lanebpsw.cpp) */
long hp_run(int mode, const uint8_t *seed, int seedlen, int bits, const int *mask, int nmask, char *outhex, int outcap, double *seconds);
int hp_predicate(int mode, const char *hex);

static int checks = 0, bad = 0;
static void expect(int ok, const char *what, mpz_srcptr n) {
    checks++;
    if (!ok) { bad++; gmp_printf("MISMATCH %s: n = %Zd\n", what, n); }
}

int main(void) {
    if (!lanebpsw_cpu_supported()) { printf("SKIP: this CPU lacks AVX-512F or AVX-512 IFMA\n"); return 77; }
    if (strcmp(gmp_version, "6.3.0") != 0) printf("WARNING: GMP %s; the outputs are meant to equal those of GMP 6.3.0\n", gmp_version);
    gmp_randstate_t rs; gmp_randinit_mt(rs); gmp_randseed_ui(rs, 20261031UL);
    mpz_t n, p, q; mpz_init(n); mpz_init(p); mpz_init(q);

    /* millerrabin25 */
    for (unsigned long v = 1; v <= 3000; v++) {
        mpz_set_ui(n, v);
        expect(lanebpsw_millerrabin25(n) == mpz_millerrabin(n, 25), "millerrabin25, small", n);
    }
    for (int k = 2; k <= 4200; k += (k < 300 ? 7 : 97)) {
        for (int i = 0; i < 4; i++) {
            mpz_urandomb(n, rs, k); mpz_setbit(n, k - 1); if (i < 3) mpz_setbit(n, 0);
            if (mpz_cmp_ui(n, 5) < 0) mpz_set_ui(n, 5);
            expect(lanebpsw_millerrabin25(n) == mpz_millerrabin(n, 25), "millerrabin25, random", n);
        }
        if (k >= 64 && k <= 2100) {
            mpz_urandomb(n, rs, k); mpz_setbit(n, k - 1); mpz_nextprime(p, n);
            expect(lanebpsw_millerrabin25(p) == mpz_millerrabin(p, 25), "millerrabin25, prime", p);
        }
    }

    /* nextprime */
    long vals[] = {-100, -1, 0, 1, 2, 3, 4, 10, 310242, 310243, 1000000};
    for (unsigned i = 0; i < sizeof vals / sizeof vals[0]; i++) {
        mpz_set_si(n, vals[i]); lanebpsw_nextprime(p, n); mpz_nextprime(q, n);
        expect(mpz_cmp(p, q) == 0, "nextprime, small", n);
    }
    for (int k = 1; k <= 2100; k += (k < 260 ? 3 : 131)) {
        mpz_urandomb(n, rs, k); mpz_setbit(n, k - 1);
        lanebpsw_nextprime(p, n); mpz_nextprime(q, n);
        expect(mpz_cmp(p, q) == 0, "nextprime, random", n);
    }
    mpz_urandomb(n, rs, 512); mpz_setbit(n, 511); mpz_nextprime(q, n);
    lanebpsw_nextprime(n, n);   /* aliasing */
    expect(mpz_cmp(n, q) == 0, "nextprime, rop = n", n);

    /* chiavdf predicate: primality up to 70000, then against the scalar version */
    for (unsigned long v = 0; v <= 70000; v++) {
        mpz_set_ui(n, v);
        int prime = v >= 2 && mpz_probab_prime_p(n, 30) > 0;
        expect(lanebpsw_chiavdf_is_prime(n) == prime, "chiavdf_is_prime, small", n);
    }
    char hex[1200];
    for (int k = 20; k <= 2048; k += 29) {
        for (int i = 0; i < 6; i++) {
            mpz_urandomb(n, rs, k); mpz_setbit(n, k - 1); mpz_setbit(n, 0);
            if (i == 0) mpz_nextprime(n, n);
            if (mpz_cmp_ui(n, 14951) <= 0) continue;
            mpz_get_str(hex, 16, n);
            expect(lanebpsw_chiavdf_is_prime(n) == (hp_predicate(0, hex) != 0), "chiavdf_is_prime vs scalar", n);
        }
    }

    /* HashPrime and CreateDiscriminant against the scalar HashPrime */
    for (int bits = 264; bits <= 1024; bits += 760) {
        for (int s = 0; s < 8; s++) {
            uint8_t seed[32]; for (int j = 0; j < 32; j++) seed[j] = (uint8_t)(17 * s + 31 * j + bits);
            int mask_b[1] = {bits - 1}, mask_d[4] = {0, 1, 2, bits - 1};
            double t;
            lanebpsw_hashprime(p, seed, 32, bits, mask_b, 1);
            hp_run(0, seed, 32, bits, mask_b, 1, hex, sizeof hex, &t); mpz_set_str(q, hex, 16);
            expect(mpz_cmp(p, q) == 0, "hashprime vs scalar", p);
            lanebpsw_create_discriminant(p, seed, 32, bits);
            hp_run(0, seed, 32, bits, mask_d, 4, hex, sizeof hex, &t); mpz_set_str(q, hex, 16); mpz_neg(q, q);
            expect(mpz_cmp(p, q) == 0 && (int)mpz_sizeinbase(p, 2) == bits && mpz_fdiv_ui(p, 8) == 1, "create_discriminant", p);
        }
    }
    int m0[1] = {0};
    expect(lanebpsw_hashprime(p, (const uint8_t *)"x", 1, 8, NULL, 0) == LANEBPSW_OK && mpz_cmp_ui(p, 89) == 0,
           "hashprime small-domain chiavdf 1.2.1+ semantics", p);
    mpz_set_ui(p, 12345);
    expect(lanebpsw_hashprime(p, (const uint8_t *)"x", 1, 260, m0, 1) == LANEBPSW_INVALID_ARGUMENT &&
           mpz_cmp_ui(p, 12345) == 0, "hashprime rejects bits % 8 != 0 and preserves output", p);
    expect(lanebpsw_hashprime(p, (const uint8_t *)"x", 0, 264, m0, 1) == LANEBPSW_INVALID_ARGUMENT &&
           mpz_cmp_ui(p, 12345) == 0, "hashprime rejects an empty seed and preserves output", p);
    int mout[1] = {264};
    expect(lanebpsw_hashprime(p, (const uint8_t *)"x", 1, 264, mout, 1) == LANEBPSW_INVALID_ARGUMENT,
           "hashprime rejects mask outside bits", p);
    expect(lanebpsw_create_discriminant(p, (const uint8_t *)"x", 1, 0) == LANEBPSW_INVALID_ARGUMENT,
           "create_discriminant rejects bits = 0", p);
    expect(lanebpsw_create_discriminant(p, (const uint8_t *)"x", 1, INT_MIN) == LANEBPSW_INVALID_ARGUMENT,
           "create_discriminant rejects INT_MIN", p);
    expect(lanebpsw_create_discriminant(p, (const uint8_t *)"x", 1, INT_MAX - 7) == LANEBPSW_INVALID_ARGUMENT,
           "create_discriminant rejects excessive aligned bit count", p);
    expect(lanebpsw_create_discriminant(p, (const uint8_t *)"x", 1, LANEBPSW_MAX_HASHPRIME_BITS + 8) ==
           LANEBPSW_INVALID_ARGUMENT, "create_discriminant rejects its documented resource limit", p);

    /* A one-byte sprout and 4096-bit candidate consume 16 hashes per candidate: exactly 16 boundary states. */
    uint8_t zero_seed[1] = {0};
    mpz_set_ui(p, 12345);
    expect(lanebpsw_create_discriminant(p, zero_seed, 1, 4096) == LANEBPSW_SEARCH_EXHAUSTED &&
           mpz_cmp_ui(p, 12345) == 0, "finite HashPrime orbit exhausts and preserves output", p);
    int all_bits[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    expect(lanebpsw_hashprime(p, zero_seed, 1, 8, all_bits, 8) == LANEBPSW_SEARCH_EXHAUSTED &&
           mpz_cmp_ui(p, 12345) == 0, "impossible all-bits mask exhausts", p);

    printf("ALL: checks %d, MISMATCHES %d\n", checks, bad);
    mpz_clear(n); mpz_clear(p); mpz_clear(q); gmp_randclear(rs);
    return bad != 0;
}
