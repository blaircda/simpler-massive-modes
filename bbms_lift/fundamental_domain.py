"""
B3 conjecture test: is the unitarity window W = {0 <= Delta - b(lambda) <= 8}
a fundamental domain for the ideal (K) = ker(phi)?

Equivalently, since phi is injective on R(Spin10)[x^{1/2}]/(K) (B2): for
every Z there is a unique X supported in W with phi(X) = phi(Z).

Test: Z = V_mu x^Delta for sample mu and gaps Delta - b(mu) far inside and
outside W; X restricted to lambda with lambda_1 <= L (both classes, both
chiralities).  Uniqueness = nullity of the W-system; existence = consistency
of each right-hand side, read off from one reduced echelon form of [A | B]
over GF(p): b_i is in col(A) iff column i vanishes on all rows of the RREF
whose A-part is zero.

Usage: fundamental_domain.py L1 [L2 ...]      Output: results/fundamental_domain.txt
"""

import os
import sys
import time
from itertools import product

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from graded import HERE, bound2, phi

from sage.all import GF, matrix

P = 8388593
MUS = [(0, 0, 0, 0, 0), (2, 0, 0, 0, 0), (1, 1, 1, 1, 1), (1, 1, 1, 1, -1), (2, 2, 0, 0, 0),
       (3, 1, 1, 1, 1), (4, 0, 0, 0, 0), (3, 3, 1, 1, -1)]
GAPS = (-3, -1, 0, 5, 8, 9, 10, 11, 12, 14, 20)   # Delta - b(mu), in units of Delta


def doms(L):
    out = []
    for par in (0, 1):
        for l in product(range(L, -1, -1), repeat=4):
            if list(l) != sorted(l, reverse=True):
                continue
            l = [2 * a + par for a in l]
            for l5 in range(-l[3], l[3] + 1, 2):
                out.append(tuple(l) + (l5,))
    return out


def run(L, lines):
    t = time.time()
    cols = [(lam, D2) for lam in doms(L) for D2 in range(bound2(lam), bound2(lam) + 17, 2)]
    rows, ent = {}, {}
    for j, c in enumerate(cols):
        for k, m in phi(*c).items():
            rows.setdefault(k, len(rows))
            ent[(rows[k], j)] = m
    tests, bent = [], {}
    nA = len(cols)
    for mu in MUS:
        for g in GAPS:
            i = len(tests)
            tests.append((mu, g))
            for k, m in phi(mu, bound2(mu) + 2 * g).items():
                rows.setdefault(k, len(rows))
                bent[(rows[k], nA + i)] = m
    ent.update(bent)
    M = matrix(GF(P), len(rows), nA + len(tests), ent, sparse=False)
    E = M.echelon_form()
    piv = E.pivots()
    rankA = sum(1 for c in piv if c < nA)
    bad = set()
    for i, c in enumerate(piv):
        if c >= nA:  # row with zero A-part: the B columns it touches are inconsistent
            for jj in range(c, nA + len(tests)):
                if E[i, jj]:
                    bad.add(jj - nA)
    s = (f"L={L}: {nA} window unknowns, nullity {nA - rankA}; "
         f"{len(tests) - len(bad)}/{len(tests)} single terms V_mu x^Delta have a window representative")
    print(s, flush=True)
    lines.append(s)
    if bad:
        s = "   no representative (mu, gap): " + str(sorted(tests[i] for i in bad))
        print(s, flush=True)
        lines.append(s)
    lines.append(f"   [{time.time() - t:.0f}s]")


def main():
    lines = []
    for L in [int(a) for a in sys.argv[1:]] or (4, 5):
        run(L, lines)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", "fundamental_domain.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
