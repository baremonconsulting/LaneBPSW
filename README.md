# LaneBPSW

Lane-packed Baillie–PSW testing for output-compatible prime generation with AVX-512 IFMA.

GMP's `mpz_nextprime` and `HashPrime` of chiavdf, the library of Chia's verifiable delay function, return the first
candidate that passes a fixed test of Baillie–PSW type. LaneBPSW evaluates all parts of one such test of one integer,
a strong probable-prime test to base 2, a strong Lucas test and, in GMP, one more Miller–Rabin round, together in the
eight lanes of one AVX-512 IFMA Montgomery multiplication. The parameters, residues and final checks are those of the
original. The implementation is designed to reproduce the named GMP 6.3.0 and chiavdf results in the documented
domains; the supplied tests found no mismatch, but equivalence is not formally proved. The method, the equivalence
argument and the measurements belong in the accompanying preprint (see [Citation](#citation)); this repository holds
the code and raw data.

The implementation contracts and remaining proof obligations are summarized in
[`docs/ALGORITHM_INVARIANTS.md`](docs/ALGORITHM_INVARIANTS.md); upstream identities
used for comparison are listed in [`docs/REFERENCE_SOURCES.md`](docs/REFERENCE_SOURCES.md).

## Functions

| function (`include/lanebpsw.h`) | returns the value of | lane-packed for | otherwise |
|---|---|---|---|
| `lanebpsw_millerrabin25(n)` | `mpz_millerrabin(n, 25)`, GMP 6.3.0 | odd regular n, 2^64 ≤ n < 2^4102 | `mpz_millerrabin` |
| `lanebpsw_nextprime(rop, n)` | `mpz_nextprime(rop, n)`, GMP 6.3.0 | 2^201 ≤ n + 1, n < 2^4100 | `mpz_nextprime` |
| `lanebpsw_chiavdf_is_prime(n)` | `is_prime_bpsw(n)`, chiavdf 1.2.0 to 1.2.2, for n > 14 951; primality for smaller n, as in 1.2.1 and later | 14 951 < n < 2^4102 | scalar implementation |
| `lanebpsw_hashprime(rop, seed, len, bits, mask, nmask)` | `HashPrime(seed, bits, mask)`, chiavdf 1.2.1–1.2.2 semantics for small candidates and the common 1.2.x predicate above 14 951 | candidates below 2^4102 | scalar implementation |
| `lanebpsw_create_discriminant(rop, seed, len, bits)` | the corresponding `CreateDiscriminant(seed, bits)` | as above | |

An odd n is regular if no Jacobi symbol (n/m) is 0 before the first −1 in GMP's sequence m = 5, 7, 11, 13, 15, 17,
19, 23, 25, … and the first −1 occurs for m ≤ 255 (paper, Section 2). The searches test eight candidates at a time
with a base-2 Fermat test and confirm the survivors in the original order.
HashPrime returns `LANEBPSW_SEARCH_EXHAUSTED` when the finite sprout orbit has no accepted candidate; invalid or
exhausted calls leave the output unchanged. Public resource limits and status values are in `include/lanebpsw.h`.

## Requirements

* x86-64 CPU with AVX-512F and AVX-512 IFMA. Call `lanebpsw_cpu_supported()` before any computational API; the
  shared library and this baseline-ISA check can be loaded safely on an unsupported x86-64 CPU.
* Linux, GNU make, clang 18 (default) or g++ 13, C++17.
* GMP 6.3.0 with 64-bit limbs and zero nail bits, as a shared library, with `gmp.h`. Other limb configurations are
  rejected at compile time. If another runtime GMP version is linked, the GMP-facing functions use the reference
  GMP implementation rather than the lane-packed path; `test_api` reports the version mismatch.
* OpenSSL 3 (`libcrypto` and headers), for SHA-256 in HashPrime.
* For the Python scripts: Python 3 with gmpy2, and chiavdf 1.2.0 from PyPI for the HashPrime benchmark.

On Ubuntu 24.04: `sudo apt install clang make libgmp-dev libssl-dev`.

## Build and test

```sh
make              # build/liblanebpsw.so, tests, benchmarks, build/example
make test         # correctness tests, a few minutes
make test-full    # the same, test_bpsw also with the timing of Table 3
```

Variables: `CXX=g++ CC=gcc`; `ARCH` (default `-march=native`; `ARCH="-march=x86-64-v4 -mavx512ifma"` builds for any
CPU with AVX-512 IFMA); `SQUARING=0` (no dedicated squaring); `GMP_INCLUDE=<dir>`,
`GMP_LIB=<path>/libgmp.so.10`, and `GMP_RPATH=<path>` if GMP is not in the default paths; `BUILD=<dir>`. Supplying an
absolute `GMP_LIB` affects link time but does not by itself force the dynamic loader to use that file; use
`GMP_RPATH`, an appropriate loader configuration or `LD_LIBRARY_PATH`, and verify the loaded library with `ldd`.
On a machine without AVX-512 IFMA, `-march=native` does not enable the IFMA instructions and the build fails; set
`ARCH` to build for another machine.

Before a release, run `scripts/release_check.sh` on a native AVX-512 IFMA host. It refuses unsupported hardware,
requires GMP 6.3.0, builds from scratch with both GCC and Clang, runs the complete functional suite and checks the
public CPU probe. The manual `native release qualification` workflow targets a correspondingly labelled self-hosted
runner; a generic hosted runner that skips computation is not release evidence.

On a CPU with AVX-512F but without IFMA, `make test BUILD=build-emu EMULATE_IFMA=1` runs the tests with an exact
AVX-512F emulation of the two IFMA instructions. It checks correctness only; its timings mean nothing. The tests
exit with code 77 on a CPU that cannot run the build.

| test | checks |
|---|---|
| `test_kernels` | the IFMA operations, Montgomery multiplication and dedicated squaring against GMP for every kernel size; the 8-lane Fermat test against `mpz_powm` |
| `test_bpsw` | 5201 inputs against `mpz_millerrabin(n, 25)`, including Fermat and Mersenne numbers and Carmichael numbers (paper, Section 6.2) |
| `test_lucv` | the Lucas-V composites of Baillie, Fiori and Wagstaff: accepted by the Lucas-V part, rejected by the base-2 test |
| `test_paths` | paths that random inputs rarely reach: the fallback to GMP, short sieve blocks, batches across 2^k, next prime below 2^201, the chiavdf variant above 2049 bits and the scalar version above 4102 bits |
| `test_api` | the C interface, resource boundaries and finite-orbit status against GMP and scalar implementations |
| `test_hashprime` | scalar/lane equality across small finite orbits, counter wrap, partial batches and selector boundaries |
| `tests/test_hp_measurement_gate.py` | exact archived summary plus semantic fault injections for measurement gates |
| `tests/hp_check.py` | the chiavdf predicate, lane-packed against scalar, 2105 inputs (needs gmpy2) |

## Use from C

```c
#include <gmp.h>
#include <lanebpsw.h>

mpz_t n, p;
/* ... */
if (lanebpsw_cpu_supported())
    lanebpsw_nextprime(p, n);      /* the same result as mpz_nextprime(p, n) */
else
    mpz_nextprime(p, n);
```

```sh
cc -O2 -Iinclude my.c -Lbuild -llanebpsw -lgmp -Wl,-rpath,"$PWD/build"
```

`examples/example.c` is a complete program. The functions keep internal state and are not thread-safe.
`liblanebpsw.so` also exports `hp_run`, `hp_predicate`, `hp_td_check` and `hp_info`, which the Python scripts call
through ctypes; they are not part of the interface.

## Layout

| path | content | paper |
|---|---|---|
| `include/lanebpsw.h`, `src/lanebpsw.cpp` | C interface, `liblanebpsw.so` | |
| `include/lanebpsw/ifma.hpp` | 8-lane Montgomery multiplication and dedicated squaring, radix 2^52 | 5 |
| `include/lanebpsw/bpsw8.hpp` | lane-packed `mpz_millerrabin(n, 25)` | 3.1, 4 |
| `include/lanebpsw/chiabpsw8.hpp` | chiavdf variant, exact trial division, scalar implementation of (C1)–(C3) | 3.2 |
| `include/lanebpsw/fermat8.hpp` | 8-lane base-2 Fermat test, one modulus per lane | 3.3 |
| `include/lanebpsw/sieve.hpp`, `nextprime_fast.hpp` | GMP's sieve, next prime | 3.3 |
| `include/lanebpsw/hashprime.hpp` | HashPrime, scalar and lane-packed | 3.3 |
| `include/lanebpsw/lanes.hpp`, `common.hpp` | lane set-up, GMP declarations, timers, seeded generator | |
| `tests/`, `examples/` | tests, example | 6.2 |
| `bench/`, `scripts/` | benchmarks and the scripts that ran them | 5–7 |
| `results/` | raw logs of the paper ([results/README.md](results/README.md)) | 6, 7 |

## Measurement evidence

Two distinct environments are included and must not be combined into one universal performance claim:

* `results/2026-10-29/`: original Intel Xeon Emerald Rapids VM measurements used by the draft tables. The result
  README identifies which current API and test files were added later.
* `results/2026-11-04-azure-amd/`: native validation on one Azure VM exposing AMD EPYC 9V74 and AVX-512 IFMA. It
  contains the environment, 268 291 pre-fix functional checks with no mismatch and all raw benchmark CSV files.
  The reviewed post-fix suite increased the count to 268 294 checks with no mismatch; its summary is supplied with
  the prepublication release artefact.

On the Azure/AMD VM, the measured ratio of mean warm thread CPU time for repeated `nextprime` after initialisation
was 3.96–4.82× versus GMP 6.3.0 across the eight tested sizes from 256 to 4096 bits. This is not a cold-start,
wall-clock, cross-machine or universal speed guarantee. The standalone final test and its incremental contribution
have different ratios; see the raw logs and report.

## Reproducing the measurements

| paper | command | log in `results/2026-10-29/` |
|---|---|---|
| Table 3, Section 6.2 | `taskset -c 1 build/test_bpsw` | `test_bpsw.log` |
| Table 4 | `scripts/run_all.sh` | `np3.log`, `np3_*.csv` |
| Section 6.3: profile, Fermat units | `scripts/run_all.sh` | `prof_np2.log`, `bench_pass.log` |
| Section 6.2: chiavdf variant | `scripts/run_all.sh` | `hp_check.log` |
| Table 5, Section 6.2: trial division | `scripts/make_chiavdf_sysgmp.sh`, `scripts/run_hp.sh` | `hashprime2*.json*`, `td_check.log` |
| Section 6.3: compilers | `scripts/ab_builds.sh` | `ab_builds.log` |
| Sections 5, 6.3: kernels, SHA-256, GMP builds | `scripts/run_extra.sh` | `ams_cycles.log`, `asm_counts_nd20.log`, `sha_api.log`, `gmp_builds.log`, `test_lucv.log` |
| Section 7: other instance | | `results/2026-10-28/` |

The scripts pin each run to vCPU 1 (`CPU=...` to change), use thread CPU time and write to `results/<date>/`.
The main next-prime benchmark prepares its prime table before timed queries, so it measures warm repeated operation.
Run times drift by 10 to 30 % between runs and instances; compare ratios of interleaved runs.

## Limitations

* The equality of outputs rests on details of GMP 6.3.0 with 64-bit limbs (sieve set and bound, block length, the
  sequence of D, the extra base, the threshold 35 · 2^46) and of chiavdf 1.2.0 to 1.2.2 (trial-division products,
  bound on D, construction of the candidates). It is tested, not proved. A release that changes any of these details
  needs the same change here.
* HashPrime masks must contain only bit positions in `[0, bits)`. For candidates through 14951 the implementation
  follows chiavdf 1.2.1 and later, which skip the small-prime product trial division in that finite domain.
  A finite seed orbit can contain no accepted candidate; this is reported as `LANEBPSW_SEARCH_EXHAUSTED`, without
  silently choosing a different seed.
* For public inputs only: the code is not constant-time.
* Not thread-safe.
* Uses GMP functions that `libgmp.so` exports but GMP does not document: `mpz_millerrabin`, `mpn_mod_1s_4p` and
  `mpn_mod_1s_4p_cps`, and in the tests `mpz_stronglucas` and `mpz_lucas_mod`.
* Uses `SHA256_Init`, `SHA256_Update` and `SHA256_Final`, which OpenSSL 3 deprecates; they are faster per call than
  `SHA256()`.
* The included performance evidence comes from documented virtualised Intel and AMD environments. It is not a
  guarantee for other processors, VM types, cold starts or workloads.

## Citation

Jan Suchánek, *LaneBPSW: Lane-packed Baillie–PSW testing for bit-exact prime generation with AVX-512 IFMA*, preprint,
2026. See `CITATION.cff`.

## Licence

Apache License 2.0; see `LICENSE` and `NOTICE`. GMP and OpenSSL are linked, not included.

## AI assistance

The tests, the scripts and the measurements were produced with OpenAI Codex under the author's direction.
