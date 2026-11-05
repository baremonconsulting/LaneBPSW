// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// Finite-orbit and ordered-partial-batch checks for the internal HashPrime implementations.
#include "lanebpsw/hashprime.hpp"
#include <climits>
#include <cstdio>

static int checks = 0, bad = 0;
static void expect(bool ok, const char *what) {
    checks++;
    if (!ok) { bad++; std::fprintf(stderr, "MISMATCH %s\n", what); }
}

int main() {
    require_cpu();
    expect(nd_for_bits(INT_MIN) == 0 && nd_for_bits(0) == 0, "non-positive selector boundary");
    expect(nd_for_bits(INT_MAX) > 0, "INT_MAX selector arithmetic is defined");
    expect(bpsw_for_bits(4102) != nullptr && bpsw_for_bits(4103) == nullptr,
           "BPSW kernel transition");
    expect(chia_for_bits(4102) != nullptr && chia_for_bits(4103) == nullptr,
           "chiavdf kernel transition");

    mpz_t scalar, lanes; mpz_init(scalar); mpz_init(lanes);
    bool saw_counter_wrap = false;
    for (unsigned s = 0; s < 256; s++) {
        uint8_t seed = (uint8_t)s;
        HPStats a, b;
        int ra = hp_scalar(scalar, &seed, 1, 8, nullptr, 0, a);
        int rb = hp_lanes(lanes, &seed, 1, 8, nullptr, 0, b);
        expect(ra == rb && (ra != LANEBPSW_OK || mpz_cmp(scalar, lanes) == 0),
               "small finite orbit preserves scalar first accepted result");
        saw_counter_wrap |= a.cand > (long)(255 - s);
    }
    expect(saw_counter_wrap, "successful one-byte search crosses counter wrap");
    mpz_t partial[8]; for (int i = 0; i < 8; i++) mpz_init(partial[i]);
    mpz_set_ui(partial[0], 89); mpz_set_ui(partial[1], 91);
    HPStats partial_stats; mpz_set_ui(lanes, 12345);
    expect(hp_confirm_survivors(lanes, partial, 2, nullptr, 8, partial_stats, false) &&
           mpz_cmp_ui(lanes, 89) == 0 && partial_stats.finals == 1,
           "partial survivor batch confirms the first accepted candidate in order");
    for (int i = 0; i < 8; i++) mpz_clear(partial[i]);

    uint8_t seed = 0;
    int all_bits[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    HPStats a, b;
    mpz_set_ui(scalar, 12345); mpz_set_ui(lanes, 12345);
    expect(hp_scalar(scalar, &seed, 1, 8, all_bits, 8, a) == LANEBPSW_SEARCH_EXHAUSTED,
           "scalar finite orbit exhaustion");
    expect(hp_lanes(lanes, &seed, 1, 8, all_bits, 8, b) == LANEBPSW_SEARCH_EXHAUSTED,
           "lane finite orbit exhaustion");
    expect(mpz_cmp_ui(scalar, 12345) == 0 && mpz_cmp_ui(lanes, 12345) == 0,
           "internal exhausted searches preserve output");
    mpz_clear(scalar); mpz_clear(lanes);
    std::printf("ALL: checks %d, MISMATCHES %d\n", checks, bad);
    return bad != 0;
}
