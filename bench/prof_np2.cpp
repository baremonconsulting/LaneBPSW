// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// prof_np2.cpp — time breakdown of the next-prime search (residues, marking, filter, final test) per query.
// usage: prof_np2 bits queries seed B Lmult lanes   (B = 0: GMP's sieve set and block length; lanes = 1: lane-packed final test)
#include "lanebpsw/nextprime_fast.hpp"
int main(int argc, char **argv) {
    require_cpu();
    if (argc < 6) { fprintf(stderr, "usage: prof_np2 bits queries seed B Lmult [lanes]\n"); return 2; }
    int bits = atoi(argv[1]); int Q = atoi(argv[2]); uint64_t seed = strtoull(argv[3], 0, 10);
    uint32_t B = (uint32_t)atol(argv[4]); double Lm = atof(argv[5]);
    uint32_t L = B == 0 ? 0 : (uint32_t)(Lm * bits);
    g_pt.build(std::max<uint32_t>(B ? B : gmp_sieve_limit(bits + 1) + 1000, 1000));
    mpz_t n, p; mpz_init(n); mpz_init(p);
    NPStats st; g_final_lanes = argc > 6 && atoi(argv[6]) == 1;
    double t0 = cpu_s();
    for (int i = 0; i < Q; i++) { rand_kbit(n, bits, seed); nextprime_fast(p, n, B, L, st, true); }
    double tt = cpu_s() - t0;
    printf("bits=%d B=%u L=%u Q=%d  total %.3f ms/q | res %.3f mark %.3f test %.3f final %.3f other %.3f (ms/q) | lanes/q %.1f finals/q %.2f\n",
        bits, B, L, Q, 1e3*tt/Q, 1e3*st.t_res/Q, 1e3*st.t_mark/Q, 1e3*st.t_test/Q, 1e3*st.t_final/Q,
        1e3*(tt - st.t_res - st.t_mark - st.t_test - st.t_final)/Q, (double)st.lanes/Q, (double)st.finals/Q);
}
