#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Jan Suchánek, Baremon s.r.o.
# Next prime with four builds (g++ / clang, with and without the dedicated squaring), three interleaved repetitions
# with new inputs (Section 6.3). Summarise with: python3 -I scripts/ab_table.py <output directory>/ab_builds.log
# usage: scripts/ab_builds.sh [output directory]    environment: CPU (default 1), GMP_INCLUDE, GMP_LIB
set -eu
cd "$(dirname "$0")/.."
OUT=${1:-results/$(date +%F)}; CPU=${CPU:-1}
mkdir -p "$OUT"
for v in "gcc_base g++ 0" "gcc_ams g++ 1" "clang_base clang++ 0" "clang_ams clang++ 1"; do
  set -- $v
  make -s BUILD=build/$1 CXX=$2 SQUARING=$3 build/$1/bench_np3
done
log="$OUT/ab_builds.log"
{ g++ --version | head -1; clang++ --version | head -1; } > "$log"
for rep in 1 2 3; do
  for spec in "256 1500" "1024 300" "2048 60" "4096 16"; do
    set -- $spec
    for v in gcc_base gcc_ams clang_base clang_ams; do
      echo -n "$v rep$rep " >> "$log"
      taskset -c $CPU build/$v/bench_np3 $1 $2 $((100 + rep)) >> "$log"
    done
  done
done
echo AB_DONE >> "$log"
python3 -I scripts/ab_table.py "$log"
