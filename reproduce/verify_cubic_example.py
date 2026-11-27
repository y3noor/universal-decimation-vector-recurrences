"""Verify the order-four cubic intermediate-field example."""
from math import gcd

from verify_independent import field, rank

multiply = field(12, 0x1053)
powers = [1]
for _ in range(4095):
    powers.append(multiply(powers[-1], 2))
assert powers[-1] == 1
assert set(powers[:-1]) == set(range(1, 4096))

beta, z, gamma = (powers[e] for e in (91, 273, 2449))
assert multiply(multiply(beta, beta), beta) == z
assert 4095 // gcd(273, 4095) == 15
assert all(multiply(multiply(x, x), x) != z for x in powers[0:4095:273])
assert gamma == beta ^ 1 and gcd(2449, 4095) == 1
assert rank([1, z, multiply(z, z)]) == 3
gamma2 = multiply(gamma, gamma)
assert multiply(gamma2, gamma) ^ gamma2 ^ gamma ^ 1 ^ z == 0
print("Verified cubic example: beta=0x%03x, z=0x%03x, gamma=0x%03x." % (beta,z,gamma))
