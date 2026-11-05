// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// bench_np3.cpp — 3-way paired benchmark on identical random N (order rotated per query):
//   GMP 6.3 mpz_nextprime | lanes + GMP final check (previous version) | lanes + lane-packed BPSW final check
// All three outputs must be identical.  usage: bench_np3 bits queries seed [csv]
#include "lanebpsw/nextprime_fast.hpp"
#include <algorithm>
int main(int argc, char **argv) {
    require_cpu();
    if (argc < 4) { fprintf(stderr, "usage: bench_np3 bits queries seed [csv]\n"); return 2; }
    int bits = atoi(argv[1]); int Q = atoi(argv[2]); uint64_t seed = strtoull(argv[3], 0, 10);
    const char *csv = argc > 4 ? argv[4] : nullptr;
    g_pt.build(gmp_sieve_limit(bits + 1) + 1000);
    FILE *f = csv ? fopen(csv, "w") : nullptr;
    if (f) fprintf(f, "i,bits,gap,t_gmp,t_lanes_gmpfinal,t_lanes_lanefinal\n");
    mpz_t n, p0, p1, p2; mpz_init(n); mpz_init(p0); mpz_init(p1); mpz_init(p2);
    std::vector<double> t0v, t1v, t2v; NPStats st1, st2; int bad = 0;
    for (int i = 0; i < Q; i++) {
        rand_kbit(n, bits, seed);
        double t[3];
        for (int r = 0; r < 3; r++) {
            int which = (i + r) % 3;
            double a = cpu_s();
            if (which == 0) mpz_nextprime(p0, n);
            else if (which == 1) { g_final_lanes = false; nextprime_fast(p1, n, 0, 0, st1); }
            else { g_final_lanes = true; nextprime_fast(p2, n, 0, 0, st2); }
            t[which] = cpu_s() - a;
        }
        if (mpz_cmp(p0, p1) || mpz_cmp(p0, p2)) bad++;
        t0v.push_back(t[0]); t1v.push_back(t[1]); t2v.push_back(t[2]);
        if (f) { mpz_sub(p1, p0, n); fprintf(f, "%d,%d,%lu,%.9f,%.9f,%.9f\n", i, bits, mpz_get_ui(p1), t[0], t[1], t[2]); }
    }
    auto mean = [](const std::vector<double> &v) { double s = 0; for (double x : v) s += x; return s / v.size(); };
    auto med = [](std::vector<double> v) { std::sort(v.begin(), v.end()); return v[v.size() / 2]; };
    double m0 = mean(t0v), m1 = mean(t1v), m2 = mean(t2v);
    std::vector<double> r02(Q), r12(Q); for (int i = 0; i < Q; i++) { r02[i] = t0v[i] / t2v[i]; r12[i] = t1v[i] / t2v[i]; }
    printf("bits=%d Q=%d mismatches=%d | mean ms: GMP %.4f  lanes+GMPfinal %.4f  lanes+BPSWlanes %.4f | ratio of means vs GMP: old %.2fx new %.2fx | new vs old %.3fx | median per-query ratio vs GMP %.2fx, new vs old %.3fx | bpsw fallbacks %lu/%lu\n",
           bits, Q, bad, 1e3 * m0, 1e3 * m1, 1e3 * m2, m0 / m1, m0 / m2, m1 / m2, med(r02), med(r12), (unsigned long)g_bpsw.fallback, (unsigned long)g_bpsw.calls);
    if (f) fclose(f);
    return bad != 0;
}
