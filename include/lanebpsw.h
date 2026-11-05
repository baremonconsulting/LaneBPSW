/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2026 Jan Suchánek, Baremon s.r.o.
 *
 * lanebpsw.h — C interface of liblanebpsw.so.
 *
 * The functions return exactly what the named reference functions return; they are faster where the AVX-512 IFMA
 * kernels apply and call the reference (or a scalar version) elsewhere. The CPU must support AVX-512F and
 * AVX-512 IFMA (check with lanebpsw_cpu_supported() before the first call). The functions keep internal state and
 * are not thread-safe. They are meant for public inputs: the code is not constant-time.
 */
#ifndef LANEBPSW_H
#define LANEBPSW_H

#include <stddef.h>
#include <stdint.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

/* 1 if this CPU supports AVX-512F and AVX-512 IFMA (AVX-512F only in a build with EMULATE_IFMA=1), otherwise 0. */
int lanebpsw_cpu_supported(void);

enum {
    LANEBPSW_OK = 0,
    LANEBPSW_INVALID_ARGUMENT = -1,
    LANEBPSW_SEARCH_EXHAUSTED = -2
};

/* Resource limits for HashPrime inputs. The scalar path remains available above the 4102-bit lane-kernel limit. */
#define LANEBPSW_MAX_HASHPRIME_BITS 1048576
#define LANEBPSW_MAX_SEED_BYTES 4096
#define LANEBPSW_MAX_MASK_ENTRIES 1048576

/* The value of GMP 6.3.0's mpz_millerrabin(n, 25) for n > 0. Odd n with 2^64 <= n < 2^4102 for which GMP's search
 * for the Lucas parameter D is regular use the lane-packed test; all other n are passed to mpz_millerrabin itself.
 * Like mpz_millerrabin itself, it returns 0 for n = 2. */
int lanebpsw_millerrabin25(mpz_srcptr n);

/* The result of GMP 6.3.0's mpz_nextprime(rop, n) for every n; rop may alias n. Inputs with n + 1 < 2^201 or with
 * more than 4100 bits are passed to mpz_nextprime itself. */
void lanebpsw_nextprime(mpz_ptr rop, mpz_srcptr n);

/* 1 if is_prime_bpsw of chiavdf accepts n, otherwise 0. For n > 14951 this is the function of chiavdf 1.2.0 to 1.2.2;
 * for smaller n it is primality, as in chiavdf 1.2.1 and later (1.2.0 rejected the primes up to 14951). */
int lanebpsw_chiavdf_is_prime(mpz_srcptr n);

/* HashPrime(seed, bits, mask) of chiavdf: rop = the same integer. seedlen must be in
 * [1, LANEBPSW_MAX_SEED_BYTES], bits a positive multiple of 8 no larger than LANEBPSW_MAX_HASHPRIME_BITS, nmask in
 * [0, LANEBPSW_MAX_MASK_ENTRIES], and every mask entry in [0, bits). Small candidates use the chiavdf 1.2.1+ rule
 * that skips product trial division through 14951. Returns LANEBPSW_OK, LANEBPSW_INVALID_ARGUMENT, or
 * LANEBPSW_SEARCH_EXHAUSTED if the finite sprout orbit contains no accepted candidate. On error rop is unchanged. */
int lanebpsw_hashprime(mpz_ptr rop, const uint8_t *seed, size_t seedlen, int bits, const int *mask, int nmask);

/* CreateDiscriminant(seed, bits) of chiavdf: rop = -HashPrime(seed, bits, {0, 1, 2, bits - 1}). It uses the same
 * status and output-preservation contract as lanebpsw_hashprime. */
int lanebpsw_create_discriminant(mpz_ptr rop, const uint8_t *seed, size_t seedlen, int bits);

#ifdef __cplusplus
}
#endif

#endif /* LANEBPSW_H */
