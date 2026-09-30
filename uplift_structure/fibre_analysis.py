"""
Item 7: coupling between fibres.

Fibre F(a,b) = conjugation classes with l2 = a, l4 = b.  Its restriction
lands in the box mu2 <= a <= mu1, mu4 <= b <= mu3.  Per level:
 (1) coupling: B4 irreps reachable from > 1 fibre, and the share of R^9_n
     they carry; connected components of the fibre graph;
 (2) kernel locality: dim of the fibre-local kernel vs the full kernel
     (restriction is injective on each fibre => fibre-local kernel = 0);
 (3) for n <= ENUM_MAX: is the set of all uplifts a Cartesian product of
     per-fibre choices?
"""

import os
import sys
from math import prod

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import R9, load_R10, max_level

from sage.all import QQ, matrix
from uplift import candidates, enumerate_uplifts, interlacing

ENUM_MAX = 10


def fib(lam):
    return (lam[1], lam[3])


def main():
    L = load_R10()
    for n in range(1, max_level(L) + 1):
        m = R9(L, n)
        cands = candidates(m)
        F = sorted({fib(l) for l in cands})

        reach = {mu: set() for mu in m}
        for l in cands:
            for mu in interlacing(l):
                reach[mu].add(fib(l))
        shared = [mu for mu in m if len(reach[mu]) > 1]
        frac = sum(m[mu] for mu in shared) / sum(m.values())

        parent = {f: f for f in F}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for fs in reach.values():
            fs = list(fs)
            for f in fs[1:]:
                parent[find(f)] = find(fs[0])
        ncomp = len({find(f) for f in F})

        mus = sorted(m)
        idx = {mu: i for i, mu in enumerate(mus)}
        A = matrix(QQ, len(mus), len(cands),
                   {(idx[mu], j): 1 for j, l in enumerate(cands) for mu in interlacing(l)},
                   sparse=True)
        kfull = len(cands) - A.rank()
        kloc = 0
        for f in F:
            js = [j for j, l in enumerate(cands) if fib(l) == f]
            kloc += len(js) - A.matrix_from_columns(js).rank()
        assert kloc == 0, f"restriction not injective on a fibre at level {n}"

        line = (f"n={n:2d} fibres={len(F):3d} components={ncomp} shared-mu={len(shared)}/{len(m)} "
                f"mass-on-shared={frac:.2f} kernel full={kfull} fibre-local={kloc}")
        if n <= ENUM_MAX:
            sols = list(enumerate_uplifts(m))
            proj = {f: set() for f in F}
            for s in sols:
                for f in F:
                    proj[f].add(tuple(sorted((l, c) for l, c in s.items() if fib(l) == f)))
            P = prod(len(v) for v in proj.values())
            line += (f" | #uplifts={len(sols)} product-of-fibre-choices={P} "
                     f"{'PRODUCT' if P == len(sols) else 'COUPLED'}")
        print(line, flush=True)


if __name__ == "__main__":
    main()
