"""Check the published evidence table and replay binary subspace searches."""

import argparse
import csv
from pathlib import Path
import re
import subprocess

from verify_independent import field


ROOT = Path(__file__).resolve().parent
RECORD = "full-replay.txt"
SAVED_ROW = re.compile(r"^(\d+) (\d+) ([0-9a-f]+) GOOD is_E=1 basis= (.+)$")
RESULT = re.compile(r"^GOOD is_E=1 basis= (.+)$")


def is_f4_case(n, m):
    return n == 2 and 7 <= m <= 12


def read_cases():
    saved = {}
    for line in (ROOT / RECORD).read_text().splitlines():
        match = SAVED_ROW.fullmatch(line)
        if match:
            n, m, modulus, basis = match.groups()
            key = int(n), int(m)
            if key in saved:
                raise ValueError(f"Duplicate saved search for {key}")
            saved[key] = modulus, basis.split()

    cases = {}
    with (ROOT / "search-evidence.csv").open(newline="") as file:
        for row in csv.DictReader(file):
            n, m = int(row["n"]), int(row["m"])
            key = n, m
            modulus = row["modulus_hex"]
            basis = row["surviving_basis_hex"].split()
            if key in cases:
                raise ValueError(f"Duplicate evidence row for {key}")
            if (row["exhaustive_search"], row["record"]) != ("full replay", RECORD):
                raise ValueError(f"Wrong execution source for {key}")
            if saved.get(key) != (modulus, basis) or len(basis) != n:
                raise ValueError(f"Modulus or survivor disagrees with saved search for {key}")
            if row["separate_F4_check"] != ("yes" if is_f4_case(n, m) else ""):
                raise ValueError(f"Wrong separate membership flag for {key}")
            if is_f4_case(n, m):
                z = int(basis[1], 16)
                multiply = field(2 * m, int(modulus, 16))
                if basis[0] != "1" or z in (0, 1) or (multiply(z, z) ^ z ^ 1) != 0:
                    raise ValueError(f"Printed survivor fails the separate F4 test for {key}")
            cases[key] = modulus, basis
    if not cases or cases.keys() != saved.keys():
        raise ValueError("The evidence table and saved search cover different cases")
    print(f"Matched {len(cases)} evidence rows to the saved log and checked the six F4 entries.",
          flush=True)
    return cases


def run(name, *args):
    return subprocess.check_output([str(ROOT / name), *(str(x) for x in args)], text=True)


def replay_case(n, m, modulus, basis):
    if m == 2:
        output = run("search_quadratic", n, modulus)
    else:
        output = run("search_general", n, m, modulus)
    survivors = [line for line in output.splitlines() if line.startswith("GOOD ")]
    match = RESULT.fullmatch(survivors[0]) if len(survivors) == 1 else None
    if match is None or match[1].split() != basis:
        raise ValueError(f"Replayed survivor disagrees with the evidence table for {(n, m)}")
    if is_f4_case(n, m):
        trials = (1 << (2 * m - 1)) - 1
        if not re.search(rf"^dimension=1 .* extension_trials={trials}$", output, re.M):
            raise ValueError(f"Wrong number of normalized candidates for {(n, m)}")
        print(f"m={m}: {trials} normalized candidates tested; one survivor", flush=True)
    print(n, m, modulus, survivors[0], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--quick", action="store_true", help="Replay selected smaller searches")
    modes.add_argument("--records-only", action="store_true", help="Check the table and saved log")
    modes.add_argument("--case", nargs=2, type=int, metavar=("N", "M"),
                       help="Replay one published case")
    args = parser.parse_args()
    cases = read_cases()
    if args.records_only:
        return
    if args.case:
        key = tuple(args.case)
        if key not in cases:
            parser.error(f"no published case for {key}")
        replay_case(*key, *cases[key])
        return
    for (n, m), (modulus, basis) in sorted(cases.items()):
        if not args.quick or n * m <= 12 or is_f4_case(n, m):
            replay_case(n, m, modulus, basis)


if __name__ == "__main__":
    main()
