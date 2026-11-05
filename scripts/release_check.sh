#!/bin/bash
# Release gate for a supported native AVX-512 IFMA host.
set -euo pipefail
cd "$(dirname "$0")/.."

if [ -n "${BUILD_ROOT:-}" ]; then BROOT=$BUILD_ROOT; remove_broot=0; else BROOT=$(mktemp -d); remove_broot=1; fi
release_tmp=$(mktemp -d)
trap 'rm -rf "$release_tmp"; if [ "$remove_broot" = 1 ]; then rm -rf "$BROOT"; fi' EXIT
if ! grep -qw avx512f /proc/cpuinfo || ! grep -qw avx512ifma /proc/cpuinfo; then
  echo "RELEASE_CHECK_FAIL: native AVX-512F and AVX-512 IFMA are required" >&2
  exit 77
fi

gmp_version=$(python3 - <<'PY'
import ctypes, ctypes.util
lib = ctypes.CDLL(ctypes.util.find_library("gmp"))
print(ctypes.c_char_p.in_dll(lib, "__gmp_version").value.decode())
PY
)
test "$gmp_version" = 6.3.0 || { echo "RELEASE_CHECK_FAIL: GMP 6.3.0 required, found $gmp_version" >&2; exit 1; }

for compiler in clang gcc; do
  if [ "$compiler" = clang ]; then CC=clang; CXX=clang++; else CC=gcc; CXX=g++; fi
  command -v "$CC" >/dev/null && command -v "$CXX" >/dev/null || {
    echo "RELEASE_CHECK_FAIL: $CC/$CXX are required" >&2; exit 1;
  }
  B="$BROOT/$compiler"
  make -j"$(nproc)" BUILD="$B" CC="$CC" CXX="$CXX" all
  "$B"/test_api | tee "$release_tmp/test_api-$compiler.log"
  grep -q 'MISMATCHES 0' "$release_tmp/test_api-$compiler.log"
  make BUILD="$B" CC="$CC" CXX="$CXX" test
  python3 - "$B/liblanebpsw.so" <<'PY'
import ctypes, pathlib, sys
lib = ctypes.CDLL(str(pathlib.Path(sys.argv[1]).resolve()))
if lib.lanebpsw_cpu_supported() != 1:
    raise SystemExit("native CPU support probe failed")
PY
done

echo "RELEASE_CHECK_PASS: native GCC and Clang; GMP $gmp_version"
