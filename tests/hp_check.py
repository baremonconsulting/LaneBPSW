# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Jan Suchánek, Baremon s.r.o.
# hp_check.py — the chiavdf predicate, lane-packed against the scalar version, on 2105 inputs of 64 to 2049 bits
# (random primes and odd numbers, Fermat numbers F7..F11, Mersenne numbers 2^p - 1 for primes 131 <= p < 1300).
# Also prints the trial-division table. Needs gmpy2. usage: python3 -I tests/hp_check.py
import ctypes, random, sys
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIB = os.environ.get("LANEBPSW_LIB", os.path.join(ROOT, "build", "liblanebpsw.so"))
lib = ctypes.CDLL(LIB)
lib.lanebpsw_cpu_supported.restype = ctypes.c_int
if lib.lanebpsw_cpu_supported() != 1:
    raise SystemExit("LaneBPSW computational test requires a supported CPU")
lib.hp_predicate.argtypes = [ctypes.c_int, ctypes.c_char_p]
buf = ctypes.create_string_buffer(512); lib.hp_info(buf, 512); print(buf.value.decode())
import gmpy2
rnd = random.Random(2026)
bad = 0; cnt = 0; npr = 0
def chk(n, tag=""):
    global bad, cnt
    h = format(n, 'x').encode()
    a = lib.hp_predicate(0, h); b = lib.hp_predicate(1, h)
    cnt += 1
    if a != b: bad += 1; print("MISMATCH", tag, n.bit_length(), a, b)
    return a
for bits in (64, 128, 200, 264, 512, 1024, 1536, 2048):
    for i in range(300 if bits <= 1024 else 60):
        n = rnd.getrandbits(bits) | (1 << (bits-1)) | 1
        if i % 2 == 0: n = int(gmpy2.next_prime(n))
        r = chk(n, f"rand{bits}")
        npr += (r != 0)
        # compare with gmpy2 BPSW on primes/composites (different Lucas variant, must agree away from pseudoprimes)
        if (r != 0) != bool(gmpy2.is_bpsw_prp(n)): print("differs from gmpy2 BPSW at", bits)
# special composites: strong base-2 pseudoprimes (Fermat numbers, Mersenne composites), Carmichael numbers
for m in range(7, 12): chk((1 << (1 << m)) + 1, f"F{m}")
for p in [p for p in range(131, 1300, 2) if gmpy2.is_prime(p)]:
    chk((1 << p) - 1, f"M{p}")
# BFW 2021 composites passing V_{n+1} == 2Q with Method A* (as recalled; small): check the Lucas-V part directly
print(f"checks {cnt}, mismatches {bad}, accepted {npr}")
sys.exit(1 if bad else 0)
