#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Jan Suchánek, Baremon s.r.o.
# Main measurements of the paper with the default build (clang, dedicated squaring), one run at a time on one CPU:
#   test_bpsw.log    5201 checks against GMP and the timing of one test (Table 3)
#   np3.log, np3_*.csv   next prime against mpz_nextprime (Table 4)
#   prof_np2.log     time breakdown of the search (Section 6.3); pairs of lines: GMP / lane-packed final test
#   bench_pass.log   both tests in units of a base-2 Fermat test (Section 6.3)
#   hp_check.log     chiavdf predicate, lanes against the scalar version (Section 6.2; needs gmpy2)
# usage: scripts/run_all.sh [output directory]    environment: CPU (default 1), PY (python with gmpy2), BUILD,
#        and the make variables (GMP_INCLUDE, GMP_LIB, ...)
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=${1:-results/$(date +%F)}; CPU=${CPU:-1}; PY=${PY:-python3}; B=${BUILD:-build}
export LANEBPSW_LIB=${LANEBPSW_LIB:-$B/liblanebpsw.so}
make -s BUILD="$B" all || exit 1
mkdir -p "$OUT"
{ lscpu | grep -E 'Model name|CPU family|Model:'; grep -o -m1 avx512ifma /proc/cpuinfo; uname -r; ${CXX:-clang++} --version | head -1; } > "$OUT/env.txt" 2>&1
echo "== test_bpsw"; taskset -c $CPU $B/test_bpsw > "$OUT/test_bpsw.log" 2>&1; echo "EXIT $?" >> "$OUT/test_bpsw.log"
echo "== next prime"
: > "$OUT/np3.log"
for spec in "256 3000" "512 1500" "768 800" "1024 600" "1536 200" "2048 150" "3072 60" "4096 40"; do
  set -- $spec
  taskset -c $CPU $B/bench_np3 $1 $2 2026 "$OUT/np3_$1.csv" >> "$OUT/np3.log" 2>&1
done
echo "== profile"
: > "$OUT/prof_np2.log"
for spec in "256 2000" "512 1000" "1024 400" "2048 100"; do
  set -- $spec
  for fl in 0 1; do taskset -c $CPU $B/prof_np2 $1 $2 7 0 0 $fl >> "$OUT/prof_np2.log" 2>&1; done
done
echo "== Fermat units"; taskset -c $CPU $B/bench_pass 5 > "$OUT/bench_pass.log" 2>&1
echo "== chiavdf predicate"; taskset -c $CPU $PY -I tests/hp_check.py > "$OUT/hp_check.log" 2>&1
echo "RUN_ALL_DONE ($OUT)"
