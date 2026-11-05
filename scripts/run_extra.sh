#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Jan Suchánek, Baremon s.r.o.
# Kernel measurements and side checks (Sections 5, 6.2, 6.3, 7):
#   ams_cycles.log       TSC ticks of the 8-lane multiplication and squaring, g++ and clang
#   asm_counts_nd20.log  IFMA instructions and vector moves of the unrolled l = 20 kernels
#   test_lucv.log        the five Lucas-V composites of Baillie, Fiori and Wagstaff
#   gmp_builds.log       system GMP against the GMP bundled with chiavdf (needs PY with chiavdf)
#   sha_api.log          OpenSSL one-shot SHA256() against SHA256_Init/Update/Final
# usage: scripts/run_extra.sh [output directory]    environment: CPU (default 1), PY, GMP_INCLUDE, GMP_LIB
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=${1:-results/$(date +%F)}; CPU=${CPU:-1}; PY=${PY:-python3}
mkdir -p "$OUT"
make -s BUILD=build/gcc_ams CXX=g++ build/gcc_ams/cyc_ams
make -s BUILD=build build/cyc_ams build/test_lucv build/bench_sha build/bench_gmp_builds
{ echo "## g++"; taskset -c $CPU build/gcc_ams/cyc_ams; echo "## clang"; taskset -c $CPU build/cyc_ams; } > "$OUT/ams_cycles.log" 2>&1
taskset -c $CPU build/test_lucv > "$OUT/test_lucv.log" 2>&1
INC=""
if [ -n "${GMP_INCLUDE:-}" ]; then INC="-I$GMP_INCLUDE"; fi
for c in g++ clang++; do
  $c -O3 -march=native -std=c++17 -Iinclude $INC -c -o build/asm_$c.o bench/asm_count_k.cpp
  for f in f_amm20 f_ams20; do
    objdump -d --no-show-raw-insn build/asm_$c.o | awk -v f="$f" '/^[0-9a-f]+ </ {p = index($0, f) > 0; next} /^$/ {p=0} p' > build/asm_$f.txt
    echo "$c $f: total $(grep -c . build/asm_$f.txt) ifma $(grep -c vpmadd52 build/asm_$f.txt) vmovdq $(grep -c vmovdq build/asm_$f.txt) vmovdq-mem $(grep vmovdq build/asm_$f.txt | grep -c '(')"
  done
done > "$OUT/asm_counts_nd20.log" 2>&1
SYS=$(readlink -f "$(ldconfig -p | awk '/libgmp\.so\.10 / {print $NF; exit}')")
WHEEL=$($PY -c "import chiavdf, os, glob; d = os.path.dirname(os.path.dirname(chiavdf.__file__)); print(glob.glob(os.path.join(d, 'chiavdf.libs', 'libgmp-*'))[0])" 2>/dev/null || true)
if [ -n "$WHEEL" ]; then
  { echo "# A = $SYS, B = $WHEEL; ratios are time(A)/time(B)"; taskset -c $CPU build/bench_gmp_builds "$SYS" "$WHEEL"; } > "$OUT/gmp_builds.log" 2>&1
fi
taskset -c $CPU build/bench_sha > "$OUT/sha_api.log" 2>&1
echo "EXTRA_DONE ($OUT)"
