"""
S5 step 2(a): is one uplift simultaneously extremal for the top content?

For an uplift A of R = R9_{n-1}, top(A)(nu) = sum of A(lambda) over lambda
with nu in box(lambda) and l1 = nu1 (= R(nu) - Res T^{-1}(A)(nu)).
Conjecture (a) at mu = nu + e1 for ALL uplifts  <=>  r9_n(mu) >= max_A top(A)(nu).

This script computes, per level and per weight nu in supp R:
  * M(nu)  = LP maximum of top(A)(nu) over the uplift polytope (an upper
             bound on the integer maximum);
  * top(A)(nu) for the min-height and max-height (ILP) uplifts of R;
and reports whether the min-height uplift attains M(nu) at every nu
(simultaneous extremality).  If it does, conjecture (a) for all uplifts is
equivalent to the single inequality r9_n >= S(top A_min).

It also re-checks (a) in that form: r9_n(nu + e1) >= M(nu) for all nu.

Usage: regge_extremal.py [nmin nmax]      Output: results/regge_extremal.txt
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS
from regge_stress import DEFAULT_EXT, load_R9

from uplift import candidates, height, interlacing, solve_ilp


def top_all(A):
    out = {}
    for lam, c in A.items():
        for mu in interlacing(lam):
            if mu[0] == lam[0]:
                out[mu] = out.get(mu, 0) + c
    return out


def main():
    from sage.all import MixedIntegerLinearProgram
    nmin, nmax = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (2, 24)
    N, r9, R9 = load_R9(DEFAULT_EXT)
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    for n in range(nmin, min(nmax, N) + 1):
        t0 = time.time()
        R = R9[n - 1]
        cands = candidates(R)
        cover = {mu: [] for mu in R}
        topof = {}
        for j, lam in enumerate(cands):
            for mu in interlacing(lam):
                cover[mu].append(j)
                if mu[0] == lam[0]:
                    topof.setdefault(mu, []).append(j)
        p = MixedIntegerLinearProgram(maximization=True, solver="GLPK")
        x = p.new_variable(nonnegative=True)
        for mu, js in cover.items():
            p.add_constraint(p.sum(x[j] for j in js) == R[mu])
        p.set_objective(p.sum(x[j] for j in range(len(cands))))
        b = p.get_backend()
        cols = [next(iter(x[j].dict())) for j in range(len(cands))]
        ncols = b.ncols()
        M = {}
        for nu in R:
            coeffs = [0] * ncols
            for j in topof.get(nu, []):
                coeffs[cols[j]] = 1
            b.set_objective(coeffs)
            b.solve()
            M[nu] = b.get_objective_value()
        Amin = solve_ilp(R, lambda l: height(l[:4]), maximize=False)
        Amax = solve_ilp(R, lambda l: height(l[:4]), maximize=True)
        tmin, tmax = top_all(Amin), top_all(Amax)
        miss_min = [nu for nu in R if tmin.get(nu, 0) < M[nu] - 1e-7]
        miss_max = [nu for nu in R if tmax.get(nu, 0) < M[nu] - 1e-7]
        frac = [nu for nu in R if abs(M[nu] - round(M[nu])) > 1e-7]
        viol = [nu for nu in R if r9[n].get((nu[0] + 2,) + nu[1:], 0) < M[nu] - 1e-7]
        log(f"n={n:2d}: {len(R)} weights; min-height attains max top at all weights: {not miss_min}"
            f" (misses {len(miss_min)}); max-height misses {len(miss_max)}; fractional LP maxima {len(frac)};"
            f" (a) violated by LP max at {len(viol)} weights [{time.time() - t0:.0f}s]")
        if miss_min:
            ex = miss_min[:3]
            log("      e.g. " + "; ".join(f"nu={nu}: top(Amin)={tmin.get(nu, 0)}, LP max={M[nu]:.2f}" for nu in ex))
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "regge_extremal.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
