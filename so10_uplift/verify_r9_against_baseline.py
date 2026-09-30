"""
Separate verification step: compare the independently derived r^8_n, r^9_n
(data/r9.json from derive_r9.py) with the baseline Sage pipeline
(baseline_SO8_construct.factor_susy_d4 and baseline_SO8_to_SO9.greedy_so8_to_so9).
Also checks, with Sage's own branching, that r^9_n restricts to r^8_n.
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from sage.all import WeylCharacterRing
from baseline_SO8_construct import fock_series, factor_susy_d4
from baseline_SO8_to_SO9 import greedy_so8_to_so9

D4 = WeylCharacterRing("D4", style="coroots")
B4 = WeylCharacterRing("B4", style="coroots")


def char(ring, entries):
    c = ring(0)
    for e in entries:
        c += e["mult"] * ring(*e["dynkin"])
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "data", "r9.json"))
    ap.add_argument("--max-level", type=int, default=10)
    args = ap.parse_args()
    data = json.load(open(args.data))["levels"]
    N = min(args.max_level, max(int(k) for k in data))

    F = fock_series(N)
    for n in range(1, N + 1):
        r8 = char(D4, data[str(n)]["r8_D4"])
        r9 = char(B4, data[str(n)]["r9_B4"])
        base8 = factor_susy_d4(F[n])
        base9 = greedy_so8_to_so9(base8)
        assert r8 == base8, f"r^8_{n} differs from baseline"
        assert r9 == base9, f"r^9_{n} differs from baseline"
        assert r9.branch(D4) == r8, f"Sage branching of r^9_{n} != r^8_{n}"
        print(f"level {n}: r^8, r^9 agree with baseline; Sage branching B4->D4 ok")
    print("ALL OK")


if __name__ == "__main__":
    main()
