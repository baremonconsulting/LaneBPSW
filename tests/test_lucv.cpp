// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// test_lucv.cpp — the Lucas-V part of the chiavdf variant (Method A*) on the five composites below 10^15 that pass
// the Lucas-V test [Baillie, Fiori and Wagstaff 2021, Table 2]: the lanes must accept them like the scalar reference,
// and the strong base-2 test (C2) must reject them. Also two primes and one composite as controls.
// Trial division is switched off so that the small composites reach the Lucas-V test. Returns 1 on any mismatch.
#include "lanebpsw/chiabpsw8.hpp"
int main() {
    require_cpu();
    g_pp.build(); g_chia_skip_td = true;
    struct Case { const char *n; int lucv, sprp2; };   // expected results
    const Case cases[] = {{"913", 1, 0}, {"150267335403", 1, 0}, {"430558874533", 1, 0}, {"14760229232131", 1, 0},
                          {"936916995253453", 1, 0}, {"1000003", 1, 1}, {"999983", 1, 1}, {"1000001", 0, 0}};
    mpz_t n, e, V, tq; mpz_init(n); mpz_init(e); mpz_init(V); mpz_init(tq);
    int bad = 0, cnt = 0;
    for (const Case &c : cases) {
        mpz_set_str(n, c.n, 10);
        long P, Q; bool f = chia_find_pq(n, P, Q);
        mpz_add_ui(e, n, 1); lucas_v_ref(V, e, P, Q, n); mpz_set_si(tq, 2 * Q); mpz_mod(tq, tq, n);
        int ref = f && mpz_cmp(V, tq) == 0;
        g_cdbg = ChiaDebug();
        chia_bpsw_lanes<6>(n);
        bool ok = f && ref == c.lucv && g_cdbg.lucv == ref && g_cdbg.sprp2 == c.sprp2;
        cnt++; bad += !ok;
        printf("n=%-16s prime=%d P=%ld Q=%ld  V-check ref=%d lanes=%d  sprp2(lanes)=%d%s\n", c.n, mpz_probab_prime_p(n, 30) > 0,
               P, Q, ref, g_cdbg.lucv, g_cdbg.sprp2, ok ? "" : "  MISMATCH");
    }
    printf("ALL: checks %d, MISMATCHES %d\n", cnt, bad);
    mpz_clear(n); mpz_clear(e); mpz_clear(V); mpz_clear(tq);
    return bad != 0;
}
