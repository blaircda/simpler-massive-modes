"""
Verification step: compare the BBMS cumulative lifts (data/bbms_graded.json)
with the class-level uplifts of so10_uplift/ (data/R10.json canonical,
data/R10_all_uplifts.json for n <= 10).

Reports, per level:
 * whether BBMS, merged to conjugation classes, equals the canonical uplift;
 * whether it is among the enumerated non-negative uplifts;
 * the chirality-odd part of BBMS (lambda with l5 != 0, minus conjugate),
   i.e. the information the class-level pipeline cannot see.
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "so10_uplift"))
from weights import dynkin_D5


def classes(entries):
    out = {}
    for e in entries:
        lam = tuple(e["ambient2"])
        lam = lam[:4] + (abs(lam[4]),)
        out[lam] = out.get(lam, 0) + e["mult"]
    return {k: v for k, v in out.items() if v}


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "data", "bbms_graded.json")
    bb = json.load(open(path))
    can = json.load(open(os.path.join(ROOT, "so10_uplift", "data", "R10.json")))["levels"]
    allu = json.load(open(os.path.join(ROOT, "so10_uplift", "data", "R10_all_uplifts.json")))["levels"]
    for n, rec in sorted(bb["levels"].items(), key=lambda t: int(t[0])):
        if "R_cumulative_D5" not in rec:
            print(f"level {n}: no BBMS solution stored")
            continue
        B = classes(rec["R_cumulative_D5"])
        C = classes(can[n]["canonical_R10_D5"])
        member = None
        if n in allu:
            cl = [tuple(c) for c in allu[n]["classes_ambient2"]]
            member = any({cl[i]: k for i, k in enumerate(u) if k} == B
                         for u in allu[n]["uplifts"])
        diff = {dynkin_D5(k): B.get(k, 0) - C.get(k, 0) for k in set(B) | set(C)
                if B.get(k, 0) != C.get(k, 0)}
        odd = {}
        for e in rec["R_cumulative_D5"]:
            lam = tuple(e["ambient2"])
            if lam[4] > 0:
                odd[dynkin_D5(lam)] = odd.get(dynkin_D5(lam), 0) + e["mult"]
            elif lam[4] < 0:
                key = dynkin_D5(lam[:4] + (-lam[4],))
                odd[key] = odd.get(key, 0) - e["mult"]
        odd = {k: v for k, v in odd.items() if v}
        print(f"level {n}: unique_in_domain {rec['unique_in_domain']}, nonneg "
              f"{rec['cumulative_nonneg']}, classes == canonical {not diff}, "
              f"among enumerated uplifts {member}")
        if diff:
            print(f"   BBMS - canonical (classes, Dynkin a5>=a4): {diff}")
        print(f"   chirality-odd part (V_lam - V_lambar, a5>a4 labels): {odd}")


if __name__ == "__main__":
    main()
