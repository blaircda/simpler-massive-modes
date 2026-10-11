"""
Per-level diagnostics of a stored BBMS graded lift (default: the unitarity-
window dataset).  Reads data/*.json only.

Per level n:
  gaps      : histogram of Delta - b(lambda) over the graded terms of X_n
              (+ for positive, - for negative coefficients); max gap
  R_n       : number of D5 irreps (with multiplicity) N_n, compared with the
              uplift invariant sum_{a>=b} m_n(a,a,b,b) of R^9_n
              (uplift_structure S1, recomputed here from R_n|B4)
  margin    : smallest multiplicity in R_n (chirality-resolved) and how many
              irreps sit at it; and the same for conjugation classes
              (V_lambda + V_lambdabar merged), which is what existence needs
              (uplift_structure S4)
  unique    : nullity 0 in the window

Output: results/level_stats.txt
"""

import json
import os
import sys
from collections import Counter
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "so10_uplift"))
from uplift import restrict


def bound(dyn):
    k, l, m, p, q = dyn
    return 1 + k + 2 * l + 3 * m + Fraction(5, 2) * p + Fraction(3, 2) * q


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "data", "bbms_unitarity_g8.json")
    levels = json.load(open(path))["levels"]
    lines = []
    for n in sorted(map(int, levels)):
        rec = levels[str(n)]
        if "R_cumulative_D5" not in rec:
            lines.append(f"level {n:2d}: no solution stored")
            continue
        gaps = Counter()
        for e in rec["graded"]:
            g = Fraction(e["Delta2"], 2) - bound(e["dynkin"])
            assert g.denominator == 1
            gaps[(int(g), "+" if e["mult"] > 0 else "-")] += 1
        R = {tuple(e["ambient2"]): e["mult"] for e in rec["R_cumulative_D5"]}
        tot = sum(R.values())
        R9 = restrict(R)
        inv = sum(m for mu, m in R9.items() if mu[0] == mu[1] and mu[2] == mu[3])
        mn = min(R.values())
        nmin = sum(1 for v in R.values() if v == mn)
        C = {}
        for lam, v in R.items():
            key = lam[:4] + (abs(lam[4]),)
            C[key] = C.get(key, 0) + v
        C = {k: v for k, v in C.items() if v}
        cmn = min(C.values())
        negres = sum(1 for v in R.values() if v < 0)
        hist = " ".join(f"{g}{s}:{c}" for (g, s), c in sorted(gaps.items()))
        lines.append(f"level {n:2d}: unique {rec['unique_in_domain']}, max gap "
                     f"{max(g for g, _ in gaps)}; N_n {tot} (S1 invariant {inv}, "
                     f"{'ok' if tot == inv else 'MISMATCH'}); #classes {len(R)}; "
                     f"min mult {mn} (x{nmin}), #negative {negres}; classes: min {cmn} "
                     f"{'NONNEG' if cmn >= 0 else 'NEGATIVE'}; gaps {hist}")
    out = "\n".join(lines) + "\n"
    print(out, end="")
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", "level_stats.txt"), "w") as f:
        f.write(out)


if __name__ == "__main__":
    main()
