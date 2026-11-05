// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// cyc_ams.cpp — TSC ticks per 8-lane Montgomery multiplication (ammx) and dedicated squaring (ams), min of 15 trials.
#include "lanebpsw/lanes.hpp"
#include <algorithm>
#include <vector>
template <int ND> void run(int k) {
    uint64_t seed = 5; mpz_t mods[8]; for (int l = 0; l < 8; l++) { mpz_init(mods[l]); rand_kbit(mods[l], k, seed); mpz_setbit(mods[l], 0); }
    LaneMods<ND> L; lanes_setup<ND>(L, mods, false, true);
    __m512i x[ND]; for (int j = 0; j < ND; j++) x[j] = L.rm[j];
    const int reps = 2000; std::vector<double> ca, cs;
    for (int t = 0; t < 15; t++) {
        uint64_t t0 = rdtsc_now(); for (int i = 0; i < reps; i++) ammx<ND>(x, x, x, L.c); uint64_t t1 = rdtsc_now();
        for (int i = 0; i < reps; i++) { if constexpr (ND <= 24) ams<ND>(x, x, L.c); else ams_rolled<ND>(x, x, L.c); } uint64_t t2 = rdtsc_now();
        ca.push_back((double)(t1 - t0) / reps); cs.push_back((double)(t2 - t1) / reps);
    }
    alignas(64) uint64_t sd[8][ND]; store_lanes<ND>(sd, x); volatile uint64_t s = sd[0][0]; (void)s;
    double a = *std::min_element(ca.begin(), ca.end()), b = *std::min_element(cs.begin(), cs.end());
    printf("ND=%2d: TSC ticks per 8-lane op: amm %.0f (IFMA %d -> %.2f IFMA/tick)  ams %.0f (IFMA %d -> %.2f/tick)  ratio %.3f\n",
           ND, a, 4*ND*ND+ND, (4.0*ND*ND+ND)/a, b, 3*ND*ND+2*ND, (3.0*ND*ND+2*ND)/b, a/b);
}
int main() { require_cpu(); run<5>(256); run<10>(512); run<20>(1024); run<40>(2048); run<79>(4096); }
