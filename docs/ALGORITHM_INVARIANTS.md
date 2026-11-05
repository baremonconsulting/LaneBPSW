# Algorithm and implementation invariants

This note records the contracts that connect the reference predicates to the lane-packed implementation. It is a
review map, not a claim of a complete formal proof.

## Lane mapping and ordered searches

Each AVX-512 vector holds one radix-2^52 digit for eight independent logical lanes. A lane never contributes state
to another lane. `nextprime_fast` and HashPrime construct candidates in reference order. The base-2 Fermat stage is
only a speculative rejection filter: failure implies failure of the strong base-2 condition, while every survivor
is evaluated by the complete predicate. Survivors are confirmed in original candidate order, including a final
batch of fewer than eight candidates.

HashPrime increments the big-endian seed (the sprout) once for each SHA-256 block. Orbit exhaustion is tested only
at candidate boundaries, after all blocks for that candidate and after its surviving partial batch has been
confirmed. This preserves every successful reference output and returns `LANEBPSW_SEARCH_EXHAUSTED` when a finite
orbit has no accepted candidate.

## Acceptance predicates and parameter selection

The GMP-facing path uses the same strong base-2 test, Lucas parameter sequence, strong Lucas test, additional
Miller–Rabin base and fallback domains documented for GMP 6.3.0. Irregular Jacobi selection, small inputs and inputs
beyond the kernel domain fall back to GMP.

The chiavdf path uses its base-2 strong probable-prime condition and Method A* Lucas-V condition. `chia_find_pq`
selects D, P and Q in the reference order and rejects when the reference search bound is reached. Trial division is
exact; for candidates through 14951 the 1.2.1+ semantics skip product trial division.

## Montgomery and state bounds

The Montgomery radix is `R = 2^(52*ND)`. Digits passed to VPMADD52 are below 2^52. In the Fermat filter, states
satisfy `x,y < 4N` and the selected kernel has `R > 16N`. Lucas callers can reach `x,y < 6N`; their kernel selection
gives `R > 64N`. Thus a reduction is bounded by `N + 36N^2/R < N + 36N/64 = 25N/16 < 2N`.

The radix conversion and carry pass restore digits below 2^52. Signed Lucas normalisation adds or subtracts the
modulus before reuse, maintaining the caller bounds. Intermediate unsigned accumulators and the emulated 52-bit
products are covered by the independent REDC tests and UBSan boundary builds. These arguments establish the stated
headroom; they do not constitute a formal equivalence proof for every instruction sequence.

## Kernel and fallback domains

Selectors use checked wide arithmetic. Lane kernels cover the documented limits through 4102 bits; larger supported
public inputs use scalar/reference paths. HashPrime separately enforces the resource limits in `lanebpsw.h`.
Computational entry points require an exact `1` from `lanebpsw_cpu_supported()`. Public inputs only are supported:
execution is not constant-time, and internal mutable tables make the library not thread-safe.

## Remaining proof obligations

The tests compare kernels, REDC, Lucas recurrences, predicates and generated outputs with independent or upstream
references over the documented vectors. They are strong empirical evidence, but exhaustive equivalence, all
possible compiler transformations and all microarchitectural behaviours remain outside the current proof scope.
