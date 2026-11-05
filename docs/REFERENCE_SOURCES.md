# Pinned reference sources

Release qualification uses the following upstream identities. A review bundle must archive the named files and
record their SHA-256 values before comparison.

* GMP 6.3.0 release source: `https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz`
* GMP 6.3.0 mirror commit `32d04bb5400b9d6b6cd98f7756e5f65b09c1d0c9`:
  `mpz/stronglucas.c` blob `0bf1ce019291d6d498ca8f2af38a93847b7e53a4` and
  `mpz/millerrabin.c` blob `98c4d6af39ce6f31c30d6286632847acbbb214e8`.
* chiavdf 1.2.1 commit `e3e1f042aa2bb3bf18f149782978660d7c9ce49f`: `src/primetest.h` blob
  `b82c7e85f4b8fd35570dddf2bdc1221e8366fed8`.
* chiavdf 1.2.2 commit `dc84ae53799cd285af3e6fda85181161f9a8d0c0`: `src/primetest.h` blob
  `b82c7e85f4b8fd35570dddf2bdc1221e8366fed8`, `src/pprods.h` blob
  `d821bc85a6b8bf63d574424aa2eb82f70932ef99`, and `src/checked_cast.h` blob
  `57ddaac344a87b7f8261012be05ad1d241fca51f`.
* Historical chiavdf comparison point: commit prefix `f41299e`, as identified beside the corresponding headers.

The exact downloaded commit IDs and file hashes belong in generated review evidence. Moving branch names or an
unverified package cache are not acceptable reference identities.
