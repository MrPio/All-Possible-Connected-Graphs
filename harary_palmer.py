import sys
import time
from collections import Counter
from math import comb

known = {}

sys.set_int_max_str_digits(100000)


def harary_palmer(k) -> int:
    if k in known.keys():
        return known[k]
    known[k] = pow(2, comb(k, 2)) - (
        sum(comb(k, j) * j * pow(2, comb(k - j, 2)) * harary_palmer(j) for j in range(1, k))) // k
    return known[k]


for n in []:
    start_time = time.time()
    print(harary_palmer(n) % 10)
    # print(f"Elapsed time: {(time.time() - start_time)} sec")

print([(pow(2, comb(n, 2)) - harary_palmer(n)) % 10 for n in range(3, 401)])
