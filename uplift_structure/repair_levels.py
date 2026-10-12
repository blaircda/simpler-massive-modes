"""
S13: mapping the repair levels of half (b).

A class lambda = (l1..l5) is "full" when l5 is minimal (floor = 0 for tensors,
1/2 for spinors).  A non-full class skips the bottom J = l5 - floor levels of
mu4 and reaches up to level L = l4 - floor.  For every residual rho_n of the
Regge induction (max-height and min-height rules) the script reports:

  K_L(n): least K such that rho_n lifts with non-full classes only at L <= K;
  K_J(n): least J such that rho_n lifts with every class having l5 - floor <= J;
  the least number of non-full classes at the top needed level L = K_L(n),
  and the (L, J) types of the non-full classes in that minimal lift.

The onset n_K = min{n : K_L(n) >= K} is the "repair level K" onset
(S12: n_1 = 6, n_2 = 34 for the max-height rule).

Usage: repair_levels.py [--data path] [--nmin n] [--rules max-height,min-height]   Output: results/repair_levels[_<N>].txt
"""

import argparse
import os
import sys
import time
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import half_b
from common import RESULTS, SO10
from regge_stress import T, load_R9, residual

sys.path.insert(0, SO10)
from uplift import candidates, height, interlacing, restrict, solve_ilp

HERE = os.path.dirname(os.path.abspath(__file__))


def lev(l):
    par = l[0] % 2
    return (l[3] - par) // 2, (l[4] - par) // 2          # (L, J)


def min_top(m, K):
    """Lift with non-full classes only at L <= K, minimising those at L = K; returns (count, Counter of (L,J))."""
    from sage.all import MixedIntegerLinearProgram
    cands = [l for l in candidates(m) if lev(l)[1] == 0 or lev(l)[0] <= K]
    cover = {mu: [] for mu in m}
    for j, lam in enumerate(cands):
        for mu in interlacing(lam):
            cover[mu].append(j)
    p = MixedIntegerLinearProgram(maximization=False, solver="GLPK")
    x = p.new_variable(integer=True, nonnegative=True)
    for mu, js in cover.items():
        p.add_constraint(p.sum(x[j] for j in js) == m[mu])
    p.set_objective(p.sum(x[j] for j, lam in enumerate(cands) if lev(lam)[1] > 0 and lev(lam)[0] == K))
    val = p.solve()
    xv = p.get_values(x)
    sol = {cands[j]: int(round(v)) for j, v in xv.items() if round(v)}
    assert restrict(sol) == m
    types = Counter()
    for l, c in sol.items():
        if lev(l)[1] > 0:
            types[lev(l)] += c
    return int(round(val)), types


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "data", "r9_to40.json"))
    ap.add_argument("--nmin", type=int, default=2)
    ap.add_argument("--rules", default="max-height,min-height")
    args = ap.parse_args()
    N, r9, R9 = load_R9(args.data)
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    log(f"data {os.path.basename(args.data)}: n <= {N}")
    for rule in args.rules.split(","):
        log(f"[{rule} induction]  n: K_L K_J | min #non-full at L=K_L | non-full (L,J) types in that lift")
        A, KL, onsets = {}, 0, {}
        for n in range(1, N + 1):
            rho = residual(r9[n], A)
            if rho and n >= args.nmin:
                KL = 0
                while not half_b.ilp_family(rho, lambda l, K=KL: lev(l)[1] == 0 or lev(l)[0] <= K):
                    KL += 1
                for K in range(1, KL + 1):
                    onsets.setdefault(K, n)
                KJ = 0
                while not half_b.ilp_family(rho, lambda l, J=KJ: lev(l)[1] <= J):
                    KJ += 1
                cnt, types = min_top(rho, KL) if KL else (0, Counter())
                log(f"  n={n:2d}: K_L={KL} K_J={KJ} | {cnt} | {dict(sorted(types.items()))}  [{time.time() - t0:.0f}s]")
            E = solve_ilp(rho, lambda l: height(l[:4]), maximize=(rule == "max-height")) if rho else {}
            A = T(A)
            for l, c in E.items():
                A[l] = A.get(l, 0) + c
        log(f"  onsets of repair levels (first n >= {args.nmin} with K_L >= K): {onsets}")
    os.makedirs(RESULTS, exist_ok=True)
    name = "repair_levels.txt" if N == 40 else f"repair_levels_{N}.txt"
    with open(os.path.join(RESULTS, name), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
