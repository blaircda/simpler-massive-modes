"""
Separate verification of data/R10.json (never imported by the derivation):

 1. interlacing branching rule (uplift.interlacing) agrees with Sage's own
    D5 -> B4 branching on all D5 irreps with l1 <= 3 (both chiralities);
 2. stored R9_B4 equals the cumulative sum of r^9 in data/r9.json;
 3. Sage branching of every stored uplift (canonical + extremal; both
    chirality assignments for canonical) reproduces R^9_n, n <= --sage-max;
 4. canonical R^10_n equals the hardcoded baseline for n <= 10, exactly
    (Dynkin labels, including the chirality convention a5 >= a4).
"""

import argparse
import json
import os
import sys
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from sage.all import WeylCharacterRing
from uplift import interlacing
from weights import dynkin_B4, dynkin_D5, is_dom_D
from baseline_hardcoded_level10_SO10 import get_so10_character

D5 = WeylCharacterRing("D5", style="coroots")
B4 = WeylCharacterRing("B4", style="coroots")


def d5_char(entries, flip=False):
    c = D5(0)
    for e in entries:
        lam = list(e["ambient2"])
        if flip:
            lam[4] = -lam[4]
        c += e["mult"] * D5(*dynkin_D5(lam))
    return c


def b4_char(entries):
    c = B4(0)
    for e in entries:
        c += e["mult"] * B4(*e["dynkin"])
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r9", default=os.path.join(HERE, "data", "r9.json"))
    ap.add_argument("--R10", default=os.path.join(HERE, "data", "R10.json"))
    ap.add_argument("--sage-max", type=int, default=12)
    args = ap.parse_args()
    r9 = json.load(open(args.r9))["levels"]
    R10 = json.load(open(args.R10))["levels"]
    N = max(int(k) for k in R10)

    # 1. interlacing vs Sage
    count = 0
    for lam in product(range(0, 7), repeat=5):
        lam = (lam[0], lam[1], lam[2], lam[3], lam[4] - 3)
        if not is_dom_D(lam) or len({a % 2 for a in lam}) != 1 or lam[0] > 6:
            continue
        mine = B4(0)
        for mu in interlacing(lam):
            mine += B4(*dynkin_B4(mu))
        assert D5(*dynkin_D5(lam)).branch(B4, rule="default") == mine, lam
        count += 1
    print(f"1. interlacing rule == Sage branching on {count} D5 irreps")

    # 2. cumulative sums
    acc = {}
    for n in range(1, N + 1):
        for e in r9[str(n)]["r9_B4"]:
            acc[tuple(e["ambient2"])] = acc.get(tuple(e["ambient2"]), 0) + e["mult"]
        stored = {tuple(e["ambient2"]): e["mult"] for e in R10[str(n)]["R9_B4"]}
        assert stored == acc, f"R9 at level {n} is not the cumulative sum"
    print(f"2. R9_B4 = cumulative sum of r^9, levels 1..{N}")

    # 3. Sage branching of stored uplifts
    for n in range(1, min(N, args.sage_max) + 1):
        e = R10[str(n)]
        target = b4_char(e["R9_B4"])
        sols = [e["canonical_R10_D5"]] + list(e["extremal_R10_D5"].values())
        for s in sols:
            assert d5_char(s).branch(B4, rule="default") == target, f"level {n}"
        assert d5_char(e["canonical_R10_D5"], flip=True).branch(B4, rule="default") == target
        print(f"3. level {n}: Sage branching of {len(sols)} stored uplifts (+ conjugate) == R^9_{n}")

    # 4. baseline comparison
    for n in range(1, min(N, 10) + 1):
        assert d5_char(R10[str(n)]["canonical_R10_D5"]) == get_so10_character(n), \
            f"canonical R^10_{n} != hardcoded baseline"
    print(f"4. canonical R^10_n == hardcoded baseline for n = 1..{min(N, 10)}")
    print("ALL OK")


if __name__ == "__main__":
    main()
