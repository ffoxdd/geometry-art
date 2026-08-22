#!/usr/bin/env python3
"""Time generate_globe over a fixed set of configurations.

Iteration counts are reproducible because every run is seeded, so a
change that claims to reduce work can be checked against a recorded
run rather than against an impression.
"""

import argparse
import re
import subprocess
import sys
import time

CONFIGURATIONS = [
    ("constant", 200),
    ("quadratic", 200),
    ("quadratic", 500),
    ("noise", 200),
]

SEED = 7
RESULT_PATTERN = re.compile(
    r"Final\s+(\w+) after (\d+) outer / (\d+) inner iterations, "
    r"capacity RMS ([0-9.e+-]+)"
)


def main():
    arguments = parse_arguments()
    print(
        f"{'field':<10} {'sites':>6} {'outcome':>9} {'outer':>6} "
        f"{'inner':>7} {'capacity rms':>13} {'seconds':>8}"
    )

    for field, sites in CONFIGURATIONS:
        run(arguments, field, sites)


def run(arguments, field, sites):
    command = [
        arguments.binary,
        "--render", "false",
        "--density-field", field,
        "--points", str(sites),
        "--seed", str(SEED),
        "--inner-solver", arguments.inner_solver,
        "--warm-start", arguments.warm_start,
        "--output-dir", arguments.output_dir,
    ]

    started = time.monotonic()
    completed = subprocess.run(command, capture_output=True, text=True)
    elapsed = time.monotonic() - started

    if completed.returncode != 0:
        print(f"{field:<10} {sites:>6}   failed", file=sys.stderr)
        print(completed.stderr, file=sys.stderr)
        return

    match = RESULT_PATTERN.search(completed.stdout)

    if match is None:
        print(f"{field:<10} {sites:>6}   no result line", file=sys.stderr)
        return

    outcome, outer, inner, capacity = match.groups()
    print(
        f"{field:<10} {sites:>6} {outcome:>9} {outer:>6} "
        f"{inner:>7} {capacity:>13} {elapsed:>8.1f}"
    )


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="./build-release/generate_globe")
    parser.add_argument("--inner-solver", default="lbfgs")
    parser.add_argument("--warm-start", default="lloyd")
    parser.add_argument("--output-dir", default="/tmp/globe-benchmark")
    return parser.parse_args()


if __name__ == "__main__":
    main()
