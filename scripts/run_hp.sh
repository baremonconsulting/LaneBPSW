#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Jan Suchánek, Baremon s.r.o.
# HashPrime against chiavdf (Table 5): chiavdf with its own GMP and with the system GMP, each in its own process,
# alternating over four parts of the seeds; then the trial-division check against chiavdf's gcd form (Section 6.2).
# Needs the chiavdf wheel (pip install chiavdf==1.2.0) in the Python of PY; run scripts/make_chiavdf_sysgmp.sh first.
# usage: scripts/run_hp.sh [output directory]    environment: CPU (default 1), PY, BUILD
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=${1:-results/$(date +%F)}; CPU=${CPU:-1}; PY=${PY:-python3}; B=${BUILD:-build}
export LANEBPSW_LIB=${LANEBPSW_LIB:-$B/liblanebpsw.so} CHIAVDF_SYSGMP=${CHIAVDF_SYSGMP:-$B/chiavdf_sysgmp}
make -s BUILD="$B" lib || exit 1
mkdir -p "$OUT"
out="$OUT/hashprime2.jsonl"; : > "$out"
for spec in "264 1000 11" "1024 300 12"; do
  set -- $spec
  for c in 0 1 2 3; do
    for v in shipped sysgmp; do taskset -c $CPU $PY -I bench/hp_bench3.py $v $1 $2 $3 $c 4 >> "$out"; done
  done
done
$PY -I bench/hp_agg.py "$out" > "$OUT/hashprime2_summary.json"
taskset -c $CPU $PY -I bench/hp_td_check.py > "$OUT/td_check.log"
cat "$OUT/hashprime2_summary.json" "$OUT/td_check.log"
