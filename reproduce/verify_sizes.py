"""Check the character-bound table by integer arithmetic."""
import argparse
import csv
from math import prod
from pathlib import Path
import subprocess


EXPECTED = {
    (3, 2): (19, [19]),
    (4, 2): (32, [16, 20, 24, 26, 28, 29, 31, 32]),
    (4, 3): (43, [43]),
}
EXPECTED.update({(4, k): (31 // k, []) for k in range(4, 32)})


def prime_divisor_count(n):
    count = 0
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            count += 1
            while n % divisor == 0:
                n //= divisor
        divisor += 1
    return count + (n > 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prove-primes", action="store_true",
                        help="Prove primality of every listed factor with PARI/GP")
    args = parser.parse_args()
    tested = {case: [] for case in EXPECTED}
    passing = {case: [] for case in EXPECTED}
    primes = set()
    with Path(__file__).with_name("parameter-factors.csv").open(newline="") as file:
        for row in csv.DictReader(file):
            m, k, r = (int(row[key]) for key in ("order", "relative_dimension", "exponent"))
            q = 2**r
            factors = list(map(int, row["prime_factors"].split()))
            if min(factors) <= 1 or prod(factors) != (q**(m*k) - 1) // (q - 1):
                raise ValueError(f"Invalid factorization for {(m, k, r)}")
            distinct = set(factors)
            primes.update(distinct)
            tested[m, k].append(r)
            omega_degree = prime_divisor_count(m*k)
            bound = (m*k - 2) * (2**len(distinct) - 1)
            if q > 2*omega_degree and (q + 1)**2 > q*bound**2:
                passing[m, k].append(r)
    for case, (last, expected) in EXPECTED.items():
        if tested[case] != list(range(1, last + 1)) or passing[case] != expected:
            raise ValueError(f"Parameter table disagrees for {case}")
        if expected:
            print(f"Order {case[0]}, relative dimension {case[1]}: passing exponents {expected}")
    first_order_four_dimension = min(k*r for (m, k), exponents in passing.items()
                                     if m == 4 for r in exponents)
    if first_order_four_dimension != 32:
        raise ValueError("Smallest order-four binary output dimension disagrees")
    print("No order-four case with binary output dimension below 32 meets the bound.")
    if args.prove_primes:
        vector = "[" + ",".join(map(str, sorted(primes))) + "]"
        command = f"v={vector}; print(vecmin(apply(isprime,v))); quit()\n"
        result = subprocess.run(["gp", "-fq"], input=command, text=True,
                                capture_output=True, check=True)
        if result.stdout.strip() != "1":
            raise ValueError("PARI/GP primality check failed")
        print(f"Proved primality of {len(primes)} distinct factors with PARI/GP.")


if __name__ == "__main__":
    main()
