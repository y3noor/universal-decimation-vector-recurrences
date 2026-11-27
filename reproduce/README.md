# Replaying the recurrence checks

Run from the repository root. See the [root README](../README.md)
for the Python and compiler requirements. The records and size
checks use Python.

```sh
sh reproduce/replay.sh quick
sh reproduce/replay.sh full
sh reproduce/replay.sh records
sh reproduce/replay.sh case 3 3
python3 reproduce/verify_sizes.py
python3 reproduce/verify_sizes.py --prove-primes
```

Quick reruns the feedback enumeration, the independent small searches,
and 18 normalized subspace searches. These cover the cases through ambient
binary degree 12 and the six cases with output dimension two and orders 7
through 12. Quick also checks the cubic example and the size inequalities.
Full runs the same checks and reruns all 31 normalized subspace searches.

Records compares the 31 rows of [search-evidence.csv](search-evidence.csv)
with the moduli and surviving bases in [full-replay.txt](full-replay.txt).
It also checks six degree-two subfield memberships by polynomial arithmetic.
Case checks the records and reruns one search. Give the output dimension
first and the order second.

The size check uses integers to check the products in
[parameter-factors.csv](parameter-factors.csv) and both inequalities.
The optional `--prove-primes` check uses `gp` from PARI/GP to prove
primality of the factors. The factors come from Table 2- of the
[Cunningham Project](https://homes.cerias.purdue.edu/~ssw/cun/).

Quick and full build `search_quadratic` and `search_general` and write
`results.json` and `independent_verification.json` in `reproduce/`.
Case builds the search program for the selected case. The saved log also
records the feedback enumeration and independent small searches.
