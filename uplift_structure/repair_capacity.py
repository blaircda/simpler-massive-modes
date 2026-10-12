"""
S13 (part 3): is the onset of a repair level an aggregate capacity effect?

Along the max-height Regge induction (n = 16..N), for the family of lifts of
rho_n whose non-full classes sit at levels L <= K = K_L(n) (repair_levels.py),
the least and the largest possible number of non-full classes at each level
L <= K (ILP), and their ratio.  A ratio near 1 would mean the level's total
capacity is nearly used up.

Usage: repair_capacity.py [--data path]   Output: results/repair_capacity.txt
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import repair_levels as RL
from common import RESULTS, SO10
from regge_stress import T, load_R9, residual

sys.path.insert(0, SO10)
from uplift import candidates, height, interlacing, solve_ilp

HERE = os.path.dirname(os.path.abspath(__file__))


def extreme(m, K, Lt, maximize):
    from sage.all import MixedIntegerLinearProgram
    cands = [l for l in candidates(m) if RL.lev(l)[1] == 0 or RL.lev(l)[0] <= K]
    cover = {mu: [] for mu in m}
    for j, lam in enumerate(cands):
        for mu in interlacing(lam):
            cover[mu].append(j)
    p = MixedIntegerLinearProgram(maximization=maximize, solver="GLPK")
    x = p.new_variable(integer=True, nonnegative=True)
    for mu, js in cover.items():
        p.add_constraint(p.sum(x[j] for j in js) == m[mu])
    p.set_objective(p.sum(x[j] for j, lam in enumerate(cands) if RL.lev(lam)[1] > 0 and RL.lev(lam)[0] == Lt))
    return int(round(p.solve()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "data", "r9_to40.json"))
    args = ap.parse_args()
    N, r9, R9 = load_R9(args.data)
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    A = {}
    for n in range(1, N + 1):
        rho = residual(r9[n], A)
        if n >= 16:
            K = 0
            import half_b
            while not half_b.ilp_family(rho, lambda l, KK=K: RL.lev(l)[1] == 0 or RL.lev(l)[0] <= KK):
                K += 1
            row = []
            for Lt in range(1, K + 1):
                lo, hi = extreme(rho, K, Lt, False), extreme(rho, K, Lt, True)
                row.append(f"L={Lt}: min {lo} max {hi} ratio {lo / hi if hi else 0:.3f}")
            log(f"n={n} (K_L={K}): " + " | ".join(row) + f"  [{time.time() - t0:.0f}s]")
        E = solve_ilp(rho, lambda l: height(l[:4]), maximize=True) if rho else {}
        A = T(A)
        for l, c in E.items():
            A[l] = A.get(l, 0) + c
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "repair_capacity.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
