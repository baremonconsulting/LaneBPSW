#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Jan Suchánek, Baremon s.r.o.
# A copy of the installed chiavdf package whose bundled GMP is replaced by the system GMP (Table 5, "system GMP").
# The chiavdf wheel loads chiavdf.libs/libgmp-<hash>.so.10.5.0; the copy gets the system libgmp.so.10 under that name.
# environment: PY (python with chiavdf), CHIAVDF_SYSGMP (target, default $BUILD/chiavdf_sysgmp), SYSGMP_LIB (optional)
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-python3}
DEST=${CHIAVDF_SYSGMP:-${BUILD:-build}/chiavdf_sysgmp}
SITE=$($PY -c "import chiavdf, os; print(os.path.dirname(os.path.dirname(chiavdf.__file__)))")
LIB=${SYSGMP_LIB:-$(ldconfig -p | awk '/libgmp\.so\.10 / {print $NF; exit}')}
LIB=$(readlink -f "$LIB")
rm -rf "$DEST"; mkdir -p "$DEST"
cp -a "$SITE/chiavdf" "$SITE/chiavdf.libs" "$DEST/"
n=0; for f in "$DEST"/chiavdf.libs/libgmp-*.so*; do cp "$LIB" "$f"; n=$((n + 1)); done
[ $n -eq 1 ] || { echo "expected one bundled libgmp, found $n" >&2; exit 1; }
echo "$DEST: chiavdf from $SITE with $LIB"
