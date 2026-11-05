# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Jan Suchánek, Baremon s.r.o.
"""Validate and aggregate the complete HashPrime measurement domain."""

import json
import math
import statistics as st
import string
import sys

EXPECTED = {264: 1000, 1024: 300}
VARIANTS = ("shipped", "sysgmp")
TIMERS = ("chiavdf", "scalar", "lanes")


def fail(message):
    raise SystemExit("HashPrime measurement invalid: " + message)


def load_rows(path):
    rows = []
    with open(path, encoding="utf-8") as source:
        for line_number, line in enumerate(source, 1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError as exc:
                fail(f"malformed JSON at line {line_number}: {exc}")
            if not isinstance(row, dict):
                fail(f"line {line_number} is not a JSON object")
            row["_line"] = line_number
            rows.append(row)
    return rows


def validate(rows):
    expected_total = sum(count * len(VARIANTS) for count in EXPECTED.values())
    if len(rows) != expected_total:
        fail(f"expected {expected_total} rows, found {len(rows)}")
    indexed = {}
    for row in rows:
        line = row["_line"]
        try:
            bits, variant, index = row["bits"], row["variant"], row["i"]
            ok, digest, timings = row["ok"], row["digest"], row["t"]
            cand_scalar, cand_lanes = row["cand_scalar"], row["cand_lanes"]
        except KeyError as exc:
            fail(f"missing field {exc.args[0]!r} at line {line}")
        if type(bits) is not int or bits not in EXPECTED:
            fail(f"unexpected bits at line {line}: {bits!r}")
        if variant not in VARIANTS or type(index) is not int or not 0 <= index < EXPECTED[bits]:
            fail(f"unexpected domain key at line {line}: bits={bits!r} variant={variant!r} i={index!r}")
        key = (bits, variant, index)
        if key in indexed:
            fail(f"duplicate key {key} at lines {indexed[key]['_line']} and {line}")
        if ok is not True:
            fail(f"output mismatch at line {line}: bits={bits} variant={variant} i={index}")
        if (not isinstance(digest, str) or len(digest) != 16 or
                any(character not in string.hexdigits for character in digest)):
            fail(f"invalid digest at line {line}")
        if not isinstance(timings, dict) or set(timings) != set(TIMERS):
            fail(f"invalid timing fields at line {line}")
        if any(type(timings[name]) not in (int, float) or not math.isfinite(timings[name]) or timings[name] <= 0
               for name in TIMERS):
            fail(f"invalid timing value at line {line}")
        if type(cand_scalar) is not int or type(cand_lanes) is not int or cand_scalar <= 0 or cand_lanes <= 0:
            fail(f"invalid candidate count at line {line}")
        indexed[key] = row
    for bits, count in EXPECTED.items():
        for index in range(count):
            left = indexed[(bits, "shipped", index)]
            right = indexed[(bits, "sysgmp", index)]
            if left["digest"] != right["digest"]:
                fail(f"cross-variant digest mismatch: bits={bits} i={index} "
                     f"shipped={left['digest']} sysgmp={right['digest']}")
    return indexed


def main():
    if len(sys.argv) != 2:
        fail("usage: hp_agg.py HASHPRIME_JSONL")
    indexed = validate(load_rows(sys.argv[1]))
    for bits, count in EXPECTED.items():
        by = {variant: [indexed[(bits, variant, index)] for index in range(count)] for variant in VARIANTS}
        cs = [row["cand_scalar"] for row in by["shipped"]]
        out = {"bits": bits, "seeds": count, "all_outputs_equal_within_runs": True,
               "shipped_equals_sysgmp": True,
               "candidates_to_prime_mean": round(st.mean(cs), 1),
               "candidates_to_prime_median": st.median(cs)}
        for variant in VARIANTS:
            tc = [row["t"]["chiavdf"] for row in by[variant]]
            ts = [row["t"]["scalar"] for row in by[variant]]
            tl = [row["t"]["lanes"] for row in by[variant]]
            out[variant] = {
                "chiavdf_mean_ms": round(1e3 * st.mean(tc), 4),
                "chiavdf_median_ms": round(1e3 * st.median(tc), 4),
                "scalar_mean_ms": round(1e3 * st.mean(ts), 4),
                "scalar_median_ms": round(1e3 * st.median(ts), 4),
                "lanes_mean_ms": round(1e3 * st.mean(tl), 4),
                "lanes_median_ms": round(1e3 * st.median(tl), 4),
                "chiavdf_over_lanes": round(st.mean(tc) / st.mean(tl), 3),
                "chiavdf_over_lanes_median_per_seed": round(st.median([a / b for a, b in zip(tc, tl)]), 3),
                "chiavdf_over_scalar": round(st.mean(tc) / st.mean(ts), 3),
                "scalar_over_lanes": round(st.mean(ts) / st.mean(tl), 3),
                "max_ms": {"chiavdf": round(1e3 * max(tc), 3), "lanes": round(1e3 * max(tl), 3)}}
        print(json.dumps(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
