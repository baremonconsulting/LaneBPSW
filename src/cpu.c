// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// Baseline-ISA CPU dispatch entry point. This file must not be compiled with
// the AVX-512 ARCH flags used for the computational translation unit.
#include "lanebpsw.h"

int lanebpsw_cpu_supported(void) {
    __builtin_cpu_init();
#ifdef LANEBPSW_EMULATE_IFMA
    return __builtin_cpu_supports("avx512f") != 0;
#else
    return __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512ifma");
#endif
}
