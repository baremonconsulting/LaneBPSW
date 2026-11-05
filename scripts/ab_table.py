# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Jan Suchánek, Baremon s.r.o.
# Table of the build comparison (results/ab_builds.log): ratio of mean times GMP / Algorithm (lane BPSW final),
# per build and size: mean over the repetitions, with min and max; plus the same-day old/new binary check.
import re, sys, statistics
fn = sys.argv[1] if len(sys.argv) > 1 else "results/ab_builds.log"
rx = re.compile(r"^(\S+) rep(\d+) bits=(\d+) Q=(\d+) mismatches=(\d+) \| mean ms: GMP ([\d.]+)\s+lanes\+GMPfinal ([\d.]+)\s+lanes\+BPSWlanes ([\d.]+)")
rows = {}
bad = 0; runs = 0; distinct = set()
for line in open(fn):
    m = rx.match(line)
    if not m: continue
    v, rep, bits, q, mm, g, b, c = m.groups()
    bits, q, mm, rep = int(bits), int(q), int(mm), int(rep)
    g, b, c = float(g), float(b), float(c)
    rows.setdefault((v, bits), []).append((rep, g / c, g, c, b / c))
    bad += mm; runs += q
    distinct.add((v.split('-')[0] if 'binary' in v else 'ab', rep, bits, q))
print(f"query runs {runs}, mismatches {bad}")
builds = ["gcc_base", "gcc_ams", "clang_base", "clang_ams"]
for bits in (256, 1024, 2048, 4096):
    out = [f"{bits:5d}"]
    for v in builds:
        r = [x[1] for x in rows.get((v, bits), [])]
        if r: out.append(f"{v}: {statistics.mean(r):.2f} ({min(r):.2f}-{max(r):.2f})")
    print(" | ".join(out))
    for v in builds:
        r = rows.get((v, bits), [])
        if r: print(f"      {v}: GMP ms {[round(x[2],4) for x in r]}  alg ms {[round(x[3],4) for x in r]}  filter/alg {[round(x[4],3) for x in r]}")
for v in ("old-binary", "new-binary"):
    r = rows.get((v, 2048), [])
    if r: print(v, "2048:", [(x[0], round(x[1], 2), round(x[2], 2), round(x[3], 3)) for x in r])
