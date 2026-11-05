# LaneBPSW remediation evidence

Tested source commit: `aaac0b0587f70ca18b5444bf62d9faf482eadf59`  
Host: Azure `PGTESTING`; full identity is in `environment-aaac0b0.txt`.

## Results

* Native release qualification passed from a clean archive under Clang 18.1.3 and GCC 13.3.0 with GMP 6.3.0.
  Both compilers executed the complete computational suite. See `release-check-aaac0b0.log`.
* The GCC UBSan boundary build completed with no diagnostic; `test_hashprime` reported 265 checks and `test_api`
  reported 73966 checks, both with zero mismatch. See `ubsan-aaac0b0.log`.
* The exact AVX-512F emulation build completed the full suite with zero mismatch and the CPU probe returned the
  documented exact value. See `emulation-aaac0b0.log`.
* Independent exact REDC (49680 lane checks), matrix Lucas/power (1478 checks), pinned-source verification and
  chiavdf 1.2.2 differential oracle (15105 checks) all passed with zero mismatch. See
  `independent-aaac0b0.log`.
* The original 2600 HashPrime records reproduce the original summary exactly. Fault injections for a wrong value,
  missing row, duplicate row, changed digest, malformed input and trial-division disagreement are rejected.

The source-tree identity is recorded in `source-tree-aaac0b0.sha256`. This evidence commit adds only these logs and
does not change the tested implementation.

## Scope and open gates

No full performance campaign was rerun. The remediation adds a finite-orbit boundary comparison to the HashPrime
candidate loop, so the archived pre-fix performance data must not be presented as a post-fix remeasurement. All
published performance claims remain tied to their named historical environments.

A genuinely AVX-512F-unsupported VM was not available during this run; baseline shared-library loading remains
covered by the generic CI job, but that is not native computational evidence. The release workflow now requires a
self-hosted AVX-512 IFMA runner and executes both compiler builds.

The manuscript/preprint link is still an external publication input. `CITATION.cff` must not be changed from the
repository placeholder until the author supplies the versioned manuscript or stable private URL. No external
review, publication, CRoCS contact or merge was performed.
