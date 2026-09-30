"""
Item 6: the irrep-count invariant and its fibre refinement.

For any uplift c of R^9_n (conjugation classes, l5 >= 0):

    N_n = sum_lambda c_lambda = sum_{a>=b} m_n(a,a,b,b),

and, finer, #{lambda in c : l2 = a, l4 = b} = m_n(a,a,b,b)  (fibre sizes).
Reason: V_lambda|_B4 contains exactly one mu with mu1 = mu2, mu3 = mu4,
namely (l2,l2,l4,l4).

Checks performed:
 (a) the functional f with sum_{mu in I(lambda)} f(mu) = 1 for all D5
     classes with l1 <= KMAX is unique and equals the indicator of
     {mu1 = mu2, mu3 = mu4};
 (b) for every stored level, the all-ones vector lies in the rational row
     space of the restriction matrix on admissible candidates;
 (c) N_n formula and fibre sizes agree with every stored uplift
     (canonical + extremal, all levels) and with every enumerated uplift
     (n <= 10).
"""

import json
import os
import sys
from itertools import product

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, R9, canonical, extremal, load_R10, max_level

from sage.all import QQ, matrix, vector
from uplift import candidates, interlacing

KMAX = 10  # doubled l1 bound for the universal-functional solve


def restriction_matrix(rows, cols):
    """Rows indexed by D5 classes, columns by B4 weights; entry 1 if mu in I(lambda)."""
    idx = {mu: i for i, mu in enumerate(cols)}
    return matrix(QQ, len(rows), len(cols),
                  {(i, idx[mu]): 1 for i, lam in enumerate(rows) for mu in interlacing(lam)},
                  sparse=True)


def fibre_sizes(sol):
    out = {}
    for lam, c in sol.items():
        out[(lam[1], lam[3])] = out.get((lam[1], lam[3]), 0) + c
    return out


def predicted_fibres(m):
    return {(mu[1], mu[3]): k for mu, k in m.items() if mu[0] == mu[1] and mu[2] == mu[3]}


def main():
    # (a) universal functional
    lams = [l for par in (0, 1) for l in product(range(par, KMAX + 1, 2), repeat=5)
            if l[0] >= l[1] >= l[2] >= l[3] >= l[4] >= 0]
    mus = sorted({mu for l in lams for mu in interlacing(l)})
    A = restriction_matrix(lams, mus)
    f = A.solve_right(vector(QQ, [1] * len(lams)))
    assert A.rank() == len(mus), "functional not unique"
    expected = [1 if (mu[0] == mu[1] and mu[2] == mu[3]) else 0 for mu in mus]
    assert list(f) == expected
    print(f"(a) unique f on {len(mus)} B4 weights (all D5 classes with l1 <= {KMAX // 2}): "
          f"f = indicator of mu1 = mu2, mu3 = mu4")

    L = load_R10()
    N = max_level(L)

    # (b) row space, per level
    for n in range(1, N + 1):
        m = R9(L, n)
        cands = candidates(m)
        B = restriction_matrix(cands, sorted(m)).transpose()
        assert B.stack(vector(QQ, [1] * len(cands))).rank() == B.rank(), f"level {n}"
    print(f"(b) all-ones in row space of restriction on candidates, n = 1..{N}")

    # (c) stored uplifts
    Ns = []
    for n in range(1, N + 1):
        m = R9(L, n)
        pred = predicted_fibres(m)
        for sol in [canonical(L, n)] + list(extremal(L, n).values()):
            assert fibre_sizes(sol) == pred, f"level {n}"
        assert L[str(n)]["num_classes_sum_c"] == sum(pred.values())
        Ns.append(sum(pred.values()))
    with open(os.path.join(DATA, "R10_all_uplifts.json")) as fh:
        allu = json.load(fh)["levels"]
    count = 0
    for n, rec in allu.items():
        m = R9(L, int(n))
        pred = predicted_fibres(m)
        classes = [tuple(c) for c in rec["classes_ambient2"]]
        for vec in rec["uplifts"]:
            assert fibre_sizes({l: k for l, k in zip(classes, vec) if k}) == pred
            count += 1
    print(f"(c) fibre sizes = m_n(a,a,b,b) for all stored uplifts (n <= {N}) "
          f"and all {count} enumerated uplifts (n <= {max(int(k) for k in allu)})")
    print("N_n =", ", ".join(map(str, Ns)))
    print("ALL OK")


if __name__ == "__main__":
    main()
