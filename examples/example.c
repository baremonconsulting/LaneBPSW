/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2026 Jan Suchánek, Baremon s.r.o.
 * example.c — the C interface: next prime and the discriminant of chiavdf. Built by make as build/example.
 */
#include <stdio.h>
#include <gmp.h>
#include <lanebpsw.h>

int main(void) {
    mpz_t n, p;
    mpz_init_set_str(n, "1000000000000000000000000000000000000000000000000000000000000000000000000000000", 10);
    mpz_init(p);
    int supported = lanebpsw_cpu_supported();
    if (supported)
        lanebpsw_nextprime(p, n);      /* the same result as mpz_nextprime(p, n) */
    else
        mpz_nextprime(p, n);
    gmp_printf("%Zd\n", p);

    if (supported) {
        unsigned char seed[32] = {1, 2, 3};
        lanebpsw_create_discriminant(p, seed, sizeof seed, 1024); /* = CreateDiscriminant(seed, 1024) of chiavdf */
        gmp_printf("%Zx\n", p);
    } else {
        puts("CreateDiscriminant skipped: AVX-512 IFMA is unavailable.");
    }
    mpz_clear(n); mpz_clear(p);
    return 0;
}
