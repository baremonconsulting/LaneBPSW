# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Jan Suchánek, Baremon s.r.o.
# hp_bench3.py — HashPrime / CreateDiscriminant benchmark (paper, Table 5), one chiavdf variant per process:
#   shipped : the installed chiavdf wheel (its bundled GMP)
#   sysgmp  : a copy of the wheel whose bundled libgmp is replaced by the system GMP (scripts/make_chiavdf_sysgmp.sh,
#             directory given by CHIAVDF_SYSGMP, default build/chiavdf_sysgmp)
# against our scalar and lane versions (liblanebpsw.so, system GMP), rotating order on the same seeds.
# Prints one JSON line per seed. usage: python3 -I bench/hp_bench3.py variant bits N seedval chunk nchunks
import ctypes, sys, time, random, json, hashlib
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIB = os.environ.get("LANEBPSW_LIB", os.path.join(ROOT, "build", "liblanebpsw.so"))
SYSGMP = os.environ.get("CHIAVDF_SYSGMP", os.path.join(ROOT, "build", "chiavdf_sysgmp"))
if len(sys.argv) != 7:
    raise SystemExit("usage: hp_bench3.py variant bits N seedval chunk nchunks")
variant, bits, N, seedval, chunk, nchunks = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
if variant not in ("shipped", "sysgmp") or bits <= 0 or N <= 0 or nchunks <= 0 or not 0 <= chunk < nchunks:
    raise SystemExit("invalid benchmark arguments")
if variant == "sysgmp":
    sys.path.insert(0, SYSGMP)
import chiavdf
gmp_maps = sorted({l.split()[-1] for l in open("/proc/self/maps") if "libgmp-" in l})
if len(gmp_maps) != 1 or (gmp_maps[0].startswith(os.path.realpath(SYSGMP)) != (variant == "sysgmp")):
    raise SystemExit(f"unexpected GMP mapping: {gmp_maps}")
lib = ctypes.CDLL(LIB)
lib.lanebpsw_cpu_supported.restype = ctypes.c_int
if lib.lanebpsw_cpu_supported() != 1:
    raise SystemExit("LaneBPSW computational benchmark requires a supported CPU")
lib.hp_run.restype = ctypes.c_long
lib.hp_run.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
                       ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
rnd = random.Random(seedval)
mask = (ctypes.c_int * 4)(0, 1, 2, bits - 1)
buf = ctypes.create_string_buffer(2048); sec = ctypes.c_double()
for i in range(N):
    seed = bytes(rnd.getrandbits(8) for _ in range(32))
    if i % nchunks != chunk: continue
    res, t, cand = {}, {}, {}
    for r in range(3):
        w = (i + r) % 3
        if w == 0:
            t0 = time.thread_time(); d = chiavdf.create_discriminant(seed, bits); t1 = time.thread_time()
            res[0] = -int(d, 16); t["chiavdf"] = t1 - t0
        else:
            c = lib.hp_run(w - 1, seed, len(seed), bits, mask, 4, buf, 2048, ctypes.byref(sec))
            if c < 0:
                raise SystemExit(f"HashPrime failed: variant={variant} bits={bits} seed_index={i} mode={w} status={c}")
            res[w] = int(buf.value.decode(), 16); t["scalar" if w == 1 else "lanes"] = sec.value; cand[w] = c
    row = {"variant": variant, "bits": bits, "i": i, "ok": len(set(res.values())) == 1,
           "digest": hashlib.sha256(str(res[0]).encode()).hexdigest()[:16], "t": t,
           "cand_scalar": cand[1], "cand_lanes": cand[2]}
    print(json.dumps(row), flush=True)
    if not row["ok"]:
        raise SystemExit(f"HashPrime output mismatch: variant={variant} bits={bits} seed_index={i} digests="
                         + ",".join(hashlib.sha256(str(res[w]).encode()).hexdigest()[:16] for w in range(3)))
