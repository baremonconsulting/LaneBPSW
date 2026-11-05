# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Jan Suchánek, Baremon s.r.o.
"""Fault-injection tests for the HashPrime measurement acceptance gate."""

import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/2026-10-29/hashprime2.jsonl"
SUMMARY = ROOT / "results/2026-10-29/hashprime2_summary.json"


def command(path):
    return subprocess.run([sys.executable, str(ROOT / "bench/hp_agg.py"), str(path)],
                          text=True, capture_output=True)


def write_rows(path, rows):
    path.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")


def main():
    rows = [json.loads(line) for line in SOURCE.read_text(encoding="utf-8").splitlines() if line]
    with tempfile.TemporaryDirectory() as directory:
        temp = Path(directory) / "rows.jsonl"
        valid = command(SOURCE)
        if valid.returncode or valid.stdout != SUMMARY.read_text(encoding="utf-8"):
            raise SystemExit("original 2600-row archive did not reproduce its exact summary")

        mutations = []
        wrong = copy.deepcopy(rows); wrong[0]["ok"] = False; mutations.append(("wrong value", wrong))
        mutations.append(("missing pair", rows[1:]))
        mutations.append(("duplicate", rows + [copy.deepcopy(rows[0])]))
        digest = copy.deepcopy(rows); digest[1300]["digest"] = "0" * 16; mutations.append(("changed digest", digest))
        for name, changed in mutations:
            write_rows(temp, changed)
            result = command(temp)
            if result.returncode == 0:
                raise SystemExit(f"fault injection unexpectedly passed: {name}")
        temp.write_text("{not-json}\n", encoding="utf-8")
        if command(temp).returncode == 0:
            raise SystemExit("malformed input unexpectedly passed")

    spec = importlib.util.spec_from_file_location("hp_td_check", ROOT / "bench/hp_td_check.py")
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    if "disagreements=0" not in module.validate_result(264, 10, 7):
        raise SystemExit("valid TD result rejected")
    try:
        module.validate_result(264, 10, 1000007)
    except RuntimeError:
        pass
    else:
        raise SystemExit("TD disagreement unexpectedly passed")
    print("HashPrime measurement gate: original archive PASS; six fault injections rejected")


if __name__ == "__main__":
    main()
