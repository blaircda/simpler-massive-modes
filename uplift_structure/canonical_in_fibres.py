"""
Item 8 (first half): what the canonical uplift does in fibre language.

 (1) how canonical and min-height uplifts split each shared B4 irrep mu
     among the fibres of its box: marginal position of a in [mu2, mu1] and
     of b in [mu4, mu3], over non-degenerate directions (HIST_LEVELS);
 (2) label moments sum c*l_i of the canonical uplift vs their range over all
     uplifts (exact ILP), n = 6..N;  l2, l4 moments must be invariant;
 (3) uniqueness of the max sum c*l1 uplift, n = 1..N;
 (4) uniqueness of the lexicographic maximum of (sum c*l1, sum c*l3,
     sum c*l5), and agreement with the stored canonical uplift.
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import R9, canonical, extremal, load_R10, max_level

from sage.all import MixedIntegerLinearProgram
from uplift import candidates, interlacing, is_unique_optimum, restrict, solve_ilp

HIST_LEVELS = (8, 10, 12, 16, 20, 24)
LABELS = {"l1": 0, "l2": 1, "l3": 2, "l4": 3, "l5": 4}


def moment(sol, i):
    return sum(l[i] * c for l, c in sol.items())


def box_position_marginals(m, sol):
    """
    Where along each non-degenerate box direction the multiplicity of a
    shared mu goes.  a-marginal over mu with mu1 > mu2: fraction of mass
    with a = mu1 (top), a = mu2 (bottom), or strictly between.  Likewise
    b-marginal over mu with mu3 > mu4.  Degenerate directions are excluded,
    since there the corners coincide and any labelling would be arbitrary.
    """
    split = {}
    for lam, c in sol.items():
        for mu in interlacing(lam):
            d = split.setdefault(mu, {})
            d[(lam[1], lam[3])] = d.get((lam[1], lam[3]), 0) + c
    out = {}
    for name, hi_i, lo_i, pos in (("a", 0, 1, 0), ("b", 2, 3, 1)):
        h = {"top": 0, "bottom": 0, "between": 0}
        tot = 0
        for mu, k in m.items():
            if mu[hi_i] == mu[lo_i]:
                continue
            tot += k
            for f, v in split[mu].items():
                x = f[pos]
                h["top" if x == mu[hi_i] else "bottom" if x == mu[lo_i] else "between"] += v
        out[name] = {k: v / tot for k, v in h.items()} if tot else None
    return out


def lex_label_optimum(m):
    """Sequential exact ILPs: max sum c*l1, then l3, then l5.  Returns (sol, unique)."""
    cands = candidates(m)
    cover = {mu: [] for mu in m}
    for j, l in enumerate(cands):
        for mu in interlacing(l):
            cover[mu].append(j)

    def model(fixed):
        p = MixedIntegerLinearProgram(maximization=True)
        x = p.new_variable(integer=True, nonnegative=True)
        for mu, js in cover.items():
            p.add_constraint(p.sum(x[j] for j in js) == m[mu])
        for i, v in fixed:
            p.add_constraint(p.sum(l[i] * x[j] for j, l in enumerate(cands)) == v)
        return p, x

    fixed = []
    for i in (0, 2, 4):
        p, x = model(fixed)
        p.set_objective(p.sum(l[i] * x[j] for j, l in enumerate(cands)))
        fixed.append((i, int(round(p.solve()))))
    p, x = model(fixed)
    unique = True
    for j in range(len(cands)):
        p.set_objective(x[j])
        hi = p.solve()
        p.set_objective(-x[j])
        lo = -p.solve()
        if hi - lo > 0.5:
            unique = False
            break
    p.set_objective(None)
    p.solve()
    sol = {cands[j]: int(round(v)) for j, v in p.get_values(x).items() if round(v)}
    assert restrict(sol) == m
    return sol, unique


def main():
    L = load_R10()
    N = max_level(L)

    print("(1) position of the fibre along each non-degenerate box direction (fraction of mass)")
    for n in HIST_LEVELS:
        if n > N:
            continue
        m = R9(L, n)
        for name, sol in (("canonical", canonical(L, n)), ("min_height", extremal(L, n)["min_height"])):
            h = box_position_marginals(m, sol)
            print(f"  n={n:2d} {name:10s} | "
                  + " | ".join(f"{d}: " + ", ".join(f"{k} {v:.2f}" for k, v in h[d].items())
                               for d in ("a", "b") if h[d]), flush=True)

    print("(2) label moments of canonical vs [min, max] over all uplifts")
    for n in range(6, N + 1):
        m = R9(L, n)
        can = canonical(L, n)
        row = f"  n={n:2d}"
        for name, i in LABELS.items():
            f = lambda l, i=i: l[i]
            lo = moment(solve_ilp(m, f, False), i)
            hi = moment(solve_ilp(m, f, True), i)
            if name in ("l2", "l4"):
                assert lo == hi == moment(can, i), f"{name} moment not invariant at {n}"
            row += f" | {name}: {moment(can, i)} [{lo},{hi}]"
        assert moment(can, 0) == moment(solve_ilp(m, lambda l: l[0], True), 0), \
            f"canonical does not maximise sum c*l1 at level {n}"
        print(row, flush=True)

    print("(3) is the max sum c*l1 uplift unique?")
    for n in range(1, N + 1):
        u = is_unique_optimum(R9(L, n), lambda l: l[0], canonical(L, n), maximize=True)
        print(f"  n={n:2d} {u}", flush=True)

    print("(4) lex (sum c*l1, sum c*l3, sum c*l5): unique? equals canonical?")
    for n in range(1, N + 1):
        t = time.time()
        sol, unique = lex_label_optimum(R9(L, n))
        eq = sol == canonical(L, n) if unique else None
        print(f"  n={n:2d} unique={unique}" + (f" equals canonical={eq}" if unique else "")
              + f"  [{time.time() - t:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
