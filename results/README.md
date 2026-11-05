# Results

The directories are separate evidence sets. Do not merge their absolute timings or describe them as repetitions on
one machine. Unless explicitly stated otherwise, benchmark times are warm thread CPU times, not process cold-start
or wall-clock latency.

Raw output of the measurements in the paper. Machine (`2026-10-29/env.txt`): virtual machine with two vCPUs of an
Intel Xeon at a nominal 2.1 GHz (CPUID family 6, model 207, Emerald Rapids), Linux 6.18, Ubuntu 24.04, GMP 6.3.0
(Ubuntu package, configured with `--disable-fat`), clang 18.1.3 and g++ 13.3, OpenSSL 3.0.13, chiavdf 1.2.0 wheel
from PyPI. Every run was pinned to vCPU 1; times are thread CPU times.

The programs were built from the same code before it was arranged into this repository: the kernels, the tests, the
searches and HashPrime are those in `include/lanebpsw/`. Added afterwards, and not used for these logs: the C interface
(`src/lanebpsw.cpp`), `tests/test_kernels.cpp`, `tests/test_paths.cpp`, `tests/test_api.c`, the IFMA emulation
(`EMULATE_IFMA=1`) and the CPU check at the start of each program; `test_lucv` now also prints a summary line.

Labels in the next-prime logs: `GMP` is `mpz_nextprime`; `lanes+GMPfinal` (ratio "old") is the search with the
8-lane filter and `mpz_millerrabin(p, 25)` as confirming test, the column "filter + GMP" of Table 4; `lanes+BPSWlanes`
(ratio "new") is the search with the lane-packed test as confirming test.

## 2026-10-29: Tables 3 to 5, Sections 5 and 6

| file | content | paper |
|---|---|---|
| `env.txt` | CPU, kernel, compiler, chiavdf version | 6.1 |
| `test_bpsw.log` | 5201 checks against `mpz_millerrabin(n, 25)`, then the time of one test of a prime | Table 3, 6.2 |
| `np3.log`, `np3_<bits>.csv` | next prime, three implementations on the same inputs; per-query times in the CSV files | Table 4, 6.2 |
| `prof_np2.log` | time per query: residues, marking, filter, final test; per size one line with `mpz_millerrabin` and one with the lane-packed test as final test | 6.3 |
| `bench_pass.log` | both tests in units of a base-2 Fermat test of the same implementation | 6.3 |
| `ab_builds.log` | next prime with four builds (g++ and clang, with and without the dedicated squaring), three repetitions; then `old-binary` (built on 28 October) against `new-binary` at 2048 bits | 6.3, 7 |
| `hashprime2.jsonl`, `hashprime2_summary.json` | HashPrime per seed: chiavdf with its own and with the system GMP, our scalar and lane-packed versions; summary with means, medians and ratios | Table 5, 6.2 |
| `td_check.log` | trial division with multiplicative inverses against chiavdf's gcd form, 4 × 50 000 random odd numbers | 6.2 |
| `hp_check.log` | chiavdf predicate, lane-packed against the scalar version, 2105 inputs | 6.2 |
| `test_lucv.log` | the five Lucas-V composites of Baillie, Fiori and Wagstaff | 6.2 |
| `gmp_builds.log` | system GMP against the GMP bundled with the chiavdf wheel | 6.3 |
| `ams_cycles.log` | time-stamp counter ticks of the 8-lane multiplication and squaring, g++ and clang | 5, 7 |
| `asm_counts_nd20.log` | IFMA instructions and vector moves of the unrolled ℓ = 20 kernels in the object code | 5 |
| `sha_api.log` | OpenSSL `SHA256()` against `SHA256_Init/Update/Final` | 5 |

## 2026-10-28: Section 7

Another instance of the same machine type, one day earlier, with the code of that day.

| file | content | paper |
|---|---|---|
| `ab_ams_np.log` | next prime with the four builds, three repetitions | 7 |
| `ams_cycles_gcc_vs_clang.log` | ticks of the kernels; 1117 ticks for the multiplication at ℓ = 20 with clang | 7 |

This version of `cyc_ams` printed 3ℓ² + 3ℓ as the IFMA count of the squaring; the object code of the current kernels
has 3ℓ² + 2ℓ (`2026-10-29/asm_counts_nd20.log`).

## 2026-11-04: independent native Azure/AMD validation

`2026-11-04-azure-amd/` contains a native build and validation from one Azure VM with four vCPUs exposing an
AMD EPYC 9V74, AVX-512F and AVX-512 IFMA, Ubuntu 24.04.5, clang 18.1.3 and system GMP 6.3.0. The environment and
source hashes are in `environment.txt`; `00_REPORT.md` states the scope and limitations.

The archived pre-fix validation ran 268 291 functional comparisons without mismatch. Its main benchmark ran 6 350
paired queries across 256–4096 bits without mismatch and measured a 3.96–4.82× ratio of mean warm thread CPU time
for repeated `nextprime` versus GMP 6.3.0. These figures describe this VM and workload only. A subsequent reviewed
fix release ran 268 294 comparisons without mismatch; its release report belongs with that source artefact.

The `run_all.status` file is retained as raw provenance, but the individual logs and CSV files are the evidence:
the earlier orchestration script did not reliably propagate every child failure. Current scripts use fail-fast
status propagation.

## Running the measurements again

```sh
scripts/run_all.sh results/$(date +%F)       # Tables 3 and 4, Section 6.3, hp_check (needs gmpy2)
scripts/run_extra.sh results/$(date +%F)     # Section 5: kernel ticks, instruction counts, SHA-256; Lucas-V; GMP builds
PY=/path/to/python scripts/make_chiavdf_sysgmp.sh
PY=/path/to/python scripts/run_hp.sh results/$(date +%F)    # Table 5 and the trial-division check
scripts/ab_builds.sh results/$(date +%F)     # compilers and squaring (Section 6.3)
```

`PY` must be a Python with gmpy2 and with chiavdf 1.2.0 (`pip install gmpy2 chiavdf==1.2.0`). The environment
variable `CPU` selects the vCPU (default 1). Run times drift by 10 to 30 % between runs and between instances, so
compare ratios of interleaved runs, not absolute times.
