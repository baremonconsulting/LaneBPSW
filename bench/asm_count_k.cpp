// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// asm_count_k.cpp — noinline wrappers of the ND = 20 and ND = 40 kernels, for instruction counts with objdump.
#include "lanebpsw/ifma.hpp"
__attribute__((noinline)) void f_amm20(__m512i *r, const __m512i *a, const Mont8<20> &c) { amm<20>(r, a, a, c); }
__attribute__((noinline)) void f_ams20(__m512i *r, const __m512i *a, const Mont8<20> &c) { ams<20>(r, a, c); }
__attribute__((noinline)) void f_amm40(__m512i *r, const __m512i *a, const Mont8<40> &c) { amm_rolled<40>(r, a, a, c); }
__attribute__((noinline)) void f_ams40(__m512i *r, const __m512i *a, const Mont8<40> &c) { ams_rolled<40>(r, a, c); }
