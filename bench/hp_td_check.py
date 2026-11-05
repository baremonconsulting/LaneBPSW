# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Jan Suchánek, Baremon s.r.o.
"""Run the exact trial-division cross-check and fail on any disagreement."""

import ctypes
import os


def validate_result(bits, count, result):
    disagreements, passing = divmod(result, 1000000)
    if disagreements:
        raise RuntimeError(f"trial-division disagreement at bits={bits}: {disagreements}")
    return f"bits={bits} inputs={count} disagreements={disagreements} passing_TD={passing}"


def main():
    lib = ctypes.CDLL(os.environ["LANEBPSW_LIB"])
    lib.hp_td_check.restype = ctypes.c_long
    lib.hp_td_check.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_ulonglong]
    lib.lanebpsw_cpu_supported.restype = ctypes.c_int
    if lib.lanebpsw_cpu_supported() != 1:
        raise SystemExit("trial-division check requires a supported CPU")
    for bits in (264, 512, 1024, 2048):
        count = 50000
        try:
            print(validate_result(bits, count, lib.hp_td_check(bits, count, 7 + bits)))
        except RuntimeError as exc:
            raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
