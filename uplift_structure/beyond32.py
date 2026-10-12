"""
S12/S11 check beyond n = 32 with data/r9_to40.json (from so10_uplift/derive_r9.py
--max-level 40; identical to data/r9_to32.json for n <= 32).

  (1) the Regge induction to n = 40 under the max-height and min-height rules,
      with half (a) (rho_n >= 0) and half (b) (rho_n liftable) checked at every
      step and the uplift verified by restriction;
  (2) along the max-height induction, n >= 30: does rho_n have an l5-extremal
      lift (global ILP over classes with l5 in {floor, l4}), and does the
      no-lookahead level recursion of half_b.py succeed?
  (3) n = 34..36: the minimal number of classes with intermediate l5 that a
      lift of rho_n must use (ILP), and which classes they are;
  (4) the S9 column monotonicities r9_n(nu + w) >= r9_(n-d)(nu) for all nu,
      n <= 40: w1 at +1 (Regge monotonicity, the input to half (a)), w2 at +3,
      w3 at +6, 2w4 at +10; and w1 at +2 (violations expected only at the
      gap [n-2,0,0,0]).

Usage: beyond32.py [--data path]      Output: results/beyond32.txt
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import half_b
from common import RESULTS, SO10
from regge_stress import T, load_R9, residual

sys.path.insert(0, SO10)
from uplift import candidates, height, interlacing, restrict, solve_ilp
from weights import dim_irrep, dynkin_D5

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = lambda l: l[4] in (l[3] % 2, l[3])


def min_intermediate(m):
    from sage.all import MixedIntegerLinearProgram
    cands = candidates(m)
    cover = {mu: [] for mu in m}
    for j, lam in enumerate(cands):
        for mu in interlacing(lam):
            cover[mu].append(j)
    p = MixedIntegerLinearProgram(maximization=False, solver="GLPK")
    x = p.new_variable(integer=True, nonnegative=True)
    for mu, js in cover.items():
        p.add_constraint(p.sum(x[j] for j in js) == m[mu])
    p.set_objective(p.sum(x[j] for j, lam in enumerate(cands) if not EXT(lam)))
    val = p.solve()
    xv = p.get_values(x)
    sol = {cands[j]: int(round(v)) for j, v in xv.items() if round(v)}
    assert restrict(sol) == m
    return int(round(val)), {l: c for l, c in sol.items() if not EXT(l)}


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
    log(f"data {os.path.basename(args.data)}: n <= {N}")
    for rule in ("max-height", "min-height"):
        A, last = {}, None
        for n in range(1, N + 1):
            rho = residual(r9[n], A)
            if any(v < 0 for v in rho.values()):
                last = f"half (a) fails at n={n}"
                break
            E = solve_ilp(rho, lambda l: height(l[:4]), maximize=(rule == "max-height")) if rho else {}
            if E is None:
                last = f"half (b) fails at n={n}"
                break
            if rule == "max-height" and n >= 30:
                ext = half_b.ilp_family(rho, EXT)
                rec = half_b.level_lift(rho, "ilp-min") is not None
                log(f"  (2) n={n}: rho_n l5-extremal lift {'exists' if ext else 'does NOT exist'}; "
                    f"no-lookahead level recursion {'succeeds' if rec else 'fails'}")
                if 34 <= n <= 36:
                    k, ne = min_intermediate(rho)
                    log(f"  (3) n={n}: at least {k} intermediate-l5 classes needed, e.g. " + ", ".join(
                        f"{c}x{dim_irrep(l, 'D', 5)}{list(dynkin_D5(l))}" for l, c in sorted(ne.items())[:4]))
            A = T(A)
            for l, c in E.items():
                A[l] = A.get(l, 0) + c
            assert restrict(A) == R9[n]
        log(f"(1) Regge induction [{rule}]: {last or f'runs to n={N} (both halves at every step)'}; "
            f"|A_N| = {sum(A.values())} [{time.time() - t0:.0f}s]")
    for name, w, d in (("w1", (2, 0, 0, 0), 1), ("w2", (2, 2, 0, 0), 3), ("w3", (2, 2, 2, 0), 6),
                       ("2w4", (2, 2, 2, 2), 10), ("w1", (2, 0, 0, 0), 2)):
        viol, tot, gap = 0, 0, 0
        for n in range(d + 1, N + 1):
            for nu, v in r9[n - d].items():
                tot += 1
                if r9[n].get(tuple(x + y for x, y in zip(nu, w)), 0) < v:
                    viol += 1
                    gap += nu[1:] == (0, 0, 0) and nu[0] == 2 * (n - d - 1)
        log(f"(4) {name} at +{d}: {viol} violations of {tot} (n <= {N}); {gap} of them at the top [n-d-1,0,0,0]")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "beyond32.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
