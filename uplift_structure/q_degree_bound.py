"""
S18 (model): the degree bound behind "the starts equal Qhat on the upper half of every ramp".

Lemma (proof in STATUS.md, S18).  Give a monomial q^D y^e the excess eps = 2e - D (additive under
products).  Every monomial of O = (1 - qy) Q has eps <= 0, and eps = 0 only for the single term
-q^2 y with weight 0: the T_2-difference on the vacuum line.  The proof counts the largest excess
factor by factor in
    O = q^2 y^-2 ground(y/q) det_7(1 - (q/y) g) (1 - qy)(P - 1 + y) / (1 + R):
prefactor -6, ground states <= 3, Koszul <= 0, 1/(1+R) <= 0, bracket <= 3 (only via the level-0 term
-q y^2; otherwise <= 2, from the four level-1 charge-1/2 fermions at +1/2 each).

This script checks the lemma weight by weight: it builds V exactly (q_spin.build, mod two primes)
to depth Dmax, multiplies by (1 - qy), and verifies that every (D, e) slice with eps >= 0 vanishes
except the vacuum term.  It also reports the largest excess attained by the tensor and spinor parts
(the proof gives -1 and -3/2).  No spectrum data are read.

Usage: q_degree_bound.py [--dmax 18]   Output: results/q_degree_bound.txt
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import q_spin
from common import RESULTS

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dmax", type=int, default=18)
    args = ap.parse_args()
    Dmax = args.dmax
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    log(f"O = (1 - qy) Q built exactly to depth D = {Dmax}, weight by weight; excess eps = 2e - D")
    for p in (q_spin.P1, q_spin.P2):
        S = q_spin.build(Dmax, p, log)
        S.binom(2, 2, (0, 0, 0), -1)                          # (1 - q y): depth 1, charge 1
        viol, vac, top = [], None, {0: None, 1: None}
        for t, Zt in enumerate(S.Z):
            for i, A in enumerate(Zt):
                if not A.any():
                    continue
                tq, sq = t + 2, S.smin[t] + 2 * i - 4          # Q = q y^-2 V: doubled depth and charge
                e2 = 2 * sq - tq                              # 2 eps
                if (tq, sq) == (4, 2):                        # D = 2, e = 1: the vacuum term
                    b = S.B[t]
                    rest = A.copy()
                    rest[b, b, b] = 0
                    vac = (int(A[b, b, b]) == p - 1) and not rest.any()
                    continue
                if e2 >= 0:
                    viol.append((tq / 2, sq / 2))
                par = tq % 2
                top[par] = e2 if top[par] is None else max(top[par], e2)
        log(f"  mod {p}: slices with eps >= 0 other than the vacuum term: {len(viol)}"
            + (f" (e.g. D, e = {viol[:4]})" if viol else "")
            + f"; vacuum term (D, e) = (2, 1) is exactly -1 at weight 0: {vac}"
            + f"; largest eps attained: tensor {top[0] / 2 if top[0] is not None else '-'}, "
              f"spinor {top[1] / 2 if top[1] is not None else '-'}")
        del S
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "q_degree_bound.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
