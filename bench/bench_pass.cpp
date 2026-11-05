// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// bench_pass.cpp — cost of one lane-packed BPSW test relative to (a) one 8-lane base-2 Fermat batch of the same
// size and (b) scalar GMP: mpz_millerrabin(n, 25) and one base-2 Fermat test mpz_powm(2, n-1, n).
// Inputs are primes (the case of the final check). Calls are interleaved; medians over all calls.
// usage: bench_pass [reps]
#include "lanebpsw/nextprime_fast.hpp"
#include <algorithm>
#include <vector>
int main(int argc, char **argv) {
    require_cpu();
    int reps = argc > 1 ? atoi(argv[1]) : 5;
    uint64_t seed = 20261029;
    mpz_t n, two, e, y, c[8]; mpz_init(n); mpz_init_set_ui(two, 2); mpz_init(e); mpz_init(y);
    for (int l = 0; l < 8; l++) mpz_init(c[l]);
    printf("bits  primes  lane-BPSW_us  fermat8_us  GMP-MR25_us  GMP-Fermat_us  BPSW/fermat8  GMP-MR25/GMP-Fermat  GMP-MR25/lane-BPSW\n");
    for (int k : {256, 512, 1024, 2048, 4096}) {
        int cnt = k <= 1024 ? 40 : (k <= 2048 ? 16 : 6);
        BpswFn f = bpsw_for_bits(k);
        Fermat8Fn g = fermat_for_bits(k);
        std::vector<double> ta, tb, tc, td;
        for (int i = 0; i < cnt; i++) {
            for (int l = 0; l < 8; l++) {
                do { rand_kbit(c[l], k, seed); mpz_nextprime(c[l], c[l]); } while ((int)mpz_sizeinbase(c[l], 2) != k);
            }
            mpz_set(n, c[0]);
            const mp_limb_t *cp[8];
            for (int l = 0; l < 8; l++) cp[l] = mpz_limbs_read(c[l]);
            for (int r = 0; r < reps; r++) {
                double t0 = cpu_s(); volatile int a = f(n);
                double t1 = cpu_s(); volatile uint8_t b = g(cp, (int)mpz_size(c[0]), k);
                double t2 = cpu_s(); volatile int m = mpz_millerrabin(n, 25);
                double t3 = cpu_s(); mpz_sub_ui(e, n, 1); mpz_powm(y, two, e, n);
                double t4 = cpu_s();
                (void)a; (void)b; (void)m;
                ta.push_back(t1 - t0); tb.push_back(t2 - t1); tc.push_back(t3 - t2); td.push_back(t4 - t3);
            }
        }
        auto med = [](std::vector<double> v) { std::sort(v.begin(), v.end()); return v[v.size() / 2]; };
        double A = med(ta), B = med(tb), C = med(tc), D = med(td);
        printf("%4d  %6d  %12.1f  %10.1f  %11.1f  %13.1f  %12.2f  %19.2f  %18.2f\n", k, cnt, 1e6 * A, 1e6 * B, 1e6 * C, 1e6 * D,
               A / B, C / D, C / A);
        fflush(stdout);
    }
    return 0;
}
