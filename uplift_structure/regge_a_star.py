"""
S5 step 2(a): an uplift-independent sufficient condition for half (a).

Lemma (proved in STATUS S5): for every B4 weight nu, with
    nu~ = (nu1+1, nu1+1, nu3, nu4)   (doubled: +2 on the first two coordinates' common value),
the functional y_nu = delta_nu - delta_{nu+e1} + delta_{nu~} satisfies
y_nu . Res V_lambda >= [top slice of lambda contains nu] for every D5 class lambda.
Hence for every non-negative uplift A of R = R9_{n-1}:
    top(A)(nu) <= B(nu) := R(nu) - R(nu+e1) + R(nu~),
and half (a) of the Regge continuation conjecture follows from
    (a*)   R9_n(nu+e1) >= R9_{n-1}(nu) + R9_{n-1}(nu~)     for all nu.

This script (1) checks the lemma's pointwise inequality on all classes of a
large box, (2) checks (a*) for n <= 32, and (3) for n <= NLP compares B(nu)
with the LP maximum of top(A)(nu) (is the certificate optimal?).

Output: results/regge_a_star.txt
"""

import os
import sys
import time
from itertools import product

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS
from regge_stress import DEFAULT_EXT, load_R9

from uplift import candidates, interlacing

NLP = 20


def tilde(nu):
    return (nu[0] + 2, nu[0] + 2, nu[2], nu[3])


def e1(nu):
    return (nu[0] + 2,) + nu[1:]


def check_lemma(L=7):
    """Pointwise: y_nu . Res V_lambda >= [top(lambda) contains nu], all classes l1 <= L, all nu."""
    bad = 0
    lams = []
    for par in (0, 1):
        for l in product(range(L, -1, -1), repeat=4):
            if list(l) != sorted(l, reverse=True):
                continue
            l = [2 * a + par for a in l]
            for l5 in range(l[3], -1, -2):
                lams.append(tuple(l) + (l5,))
    for lam in lams:
        box = set(interlacing(lam))
        nus = set(box) | {(m[0] - 2,) + m[1:] for m in box if m[0] - 2 >= m[1]}
        for nu in nus:
            lhs = (nu in box) - (e1(nu) in box) + (tilde(nu) in box)
            rhs = 1 if (nu in box and nu[0] == lam[0]) else 0
            bad += lhs < rhs
    return len(lams), bad


def main():
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    nl, bad = check_lemma()
    log(f"(1) lemma pointwise check on {nl} classes (l1 <= 7): violations {bad}")

    N, r9, R9 = load_R9(DEFAULT_EXT)
    log(f"(2) (a*): R9_n(nu+e1) >= R9_(n-1)(nu) + R9_(n-1)(nu~), n <= {N}")
    allok = True
    for n in range(2, N + 1):
        R, Rn = R9[n - 1], R9[n]
        viol, tight = [], 0
        for nu in R:
            lhs = Rn.get(e1(nu), 0)
            rhs = R.get(nu, 0) + R.get(tilde(nu), 0)
            if lhs < rhs:
                viol.append((nu, lhs, rhs))
            elif lhs == rhs:
                tight += 1
        allok &= not viol
        if viol or n in (3, 4, 8, 12, 16, 20, 24, 28, 32):
            log(f"  n={n:2d}: weights {len(R)}; violations {len(viol)}; equalities {tight}"
                + (f"; e.g. {viol[:3]}" if viol else ""))
    log(f"  (a*) holds for all n <= {N}: {allok}")

    from sage.all import MixedIntegerLinearProgram
    log(f"(3) optimality of the certificate: B(nu) vs LP max of top(A)(nu), n <= {NLP}")
    for n in range(2, NLP + 1):
        t0 = time.time()
        R = R9[n - 1]
        cands = candidates(R)
        cover, topof = {mu: [] for mu in R}, {}
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
        gap = []
        for nu in R:
            coeffs = [0] * b.ncols()
            for j in topof.get(nu, []):
                coeffs[cols[j]] = 1
            b.set_objective(coeffs)
            b.solve()
            M = b.get_objective_value()
            B = R.get(nu, 0) - R.get(e1(nu), 0) + R.get(tilde(nu), 0)
            assert M <= B + 1e-7, "LP max exceeds the certificate bound: lemma violated?"
            if B > M + 1e-7:
                gap.append((nu, B, M))
        if gap or n in (4, 8, 12, 16, 20):
            log(f"  n={n:2d}: certificate optimal at {len(R) - len(gap)}/{len(R)} weights"
                + (f"; e.g. {[(g[0], g[1], round(g[2], 2)) for g in gap[:3]]}" if gap else "")
                + f" [{time.time() - t0:.0f}s]")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "regge_a_star.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
