"""
S12: existence of uplifts of the cumulative R9_n through the l5-extremal level
recursion (S11), without the Regge induction.

All-full lift.  The SO(10) classes with minimal |l5| (0 for tensors, 1/2 for
spinors) restrict to a unitriangular basis, so R9_n has a unique lift by them:
    a(l1,l2,l3,l4) = Delta^4 R(l) = sum_{eps in {0,1}^4} (-1)^|eps| R(l + eps)
(R = multiplicity function of R9_n, zero off the dominant chamber; spinors in
coordinates shifted by 1/2).  Degenerate classes (l5 = l4) repair its
negatives level by level (S11): at level t = l4 the tops are
T_t(m) = a(m, t) + D_(t+1)(m) for m3 > t, plus boundary terms for m3 = t, with
0 <= D_t <= T_t.

The script reports, for every n <= NMAX:
  (A) where a < 0: by level l4 (conjecture: only at the bottom level);
  (B) the forced recursion D_t(m) = max(0, -a(m, t-1)) for columns with m3 >= t:
      slack T_t(m) - D_t(m) (minimum and where, d1 = n - m1), and boundary
      negatives above the bottom level;
  (C) the bottom repairs as a flow on the chains (k,c,0) (tensor c = 0, 1;
      spinor c = 0): deficit, deliverable supply, minimum prefix sum, and the
      chain's tail sums f(k) = sum_{j>=k} a(j,c,0,0) by d = n - k;
  (D) an actual uplift by the per-level ILP recursion of half_b.py, checked
      by restriction, with its irrep count.

  (E) with --bosonic D instead: the same all-full lift for the bosonic
      analogue B = C[J_inf(null cone of C^9)] of S6 (cumulative bosonic
      spectrum), degrees 4..D: where its coefficients are negative.

Usage: cumulative_lift.py [--data r9.json] [--nmax N] | --bosonic D
Output: results/cumulative_lift.txt (or results/cumulative_lift_<N>.txt,
        results/cumulative_lift_bosonic.txt)
"""

import argparse
import os
import sys
import time
from itertools import product

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, SO10
from regge_stress import DEFAULT_EXT, load_R9

sys.path.insert(0, SO10)
from uplift import restrict


def all_full(R, par):
    """a = Delta^4 R on the sector of parity par (coords x = (w - par)/2)."""
    S = {tuple((w - par) // 2 for w in k): v for k, v in R.items() if k[0] % 2 == par}
    if not S:
        return None, None

    def Rv(x):
        return S.get(x, 0) if x[0] >= x[1] >= x[2] >= x[3] >= 0 else 0

    tops = set()
    for x in S:
        for e in product((0, -1), repeat=4):
            y = tuple(p + q for p, q in zip(x, e))
            if y[0] >= y[1] >= y[2] >= y[3] >= 0:
                tops.add(y)
    a = {}
    for y in tops:
        v = sum((-1) ** sum(e) * Rv(tuple(p + q for p, q in zip(y, e))) for e in product((0, 1), repeat=4))
        if v:
            a[y] = v
    return a, S


def forced(a):
    """Forced recursion; returns per-level tops T, forced D, slacks, boundary negatives above the bottom."""
    tmax = max(y[3] for y in a)
    T = {tmax: {y[:3]: v for y, v in a.items() if y[3] == tmax}}
    D, slacks, bad_bnd = {}, [], []
    for t in range(tmax, 0, -1):
        base = {y[:3]: v for y, v in a.items() if y[3] == t - 1}
        Dt = {}
        for m, v in base.items():
            if m[2] >= t and v < 0:
                have = T[t].get(m, 0)
                slacks.append((have + v, t, m))
                Dt[m] = min(-v, max(have, 0))
        D[t] = Dt
        Tn = dict(base)
        for m, c in Dt.items():
            m1, m2, m3 = m
            for k, s in ((m, 1), ((m1, m2, t - 1), -1), ((m1, m3 - 1, t - 1), 1), ((m2 - 1, m3 - 1, t - 1), -1)):
                if k[0] >= k[1] >= k[2] >= t - 1:
                    Tn[k] = Tn.get(k, 0) + s * c
        T[t - 1] = {k: v for k, v in Tn.items() if v}
        if t - 1 > 0:
            bad_bnd += [(t - 1, k, v) for k, v in T[t - 1].items() if v < 0]
    return T, D, slacks, bad_bnd


def chain(T, D, c):
    """Bottom chain (k, c, 0): node values, repair arcs j -> k, deficit, deliverable, min prefix."""
    from sage.all import DiGraph
    T0, T1, D1 = T.get(0, {}), T.get(1, {}), D.get(1, {})
    nodes = sorted({k[0] for k in T0 if k[1] == c and k[2] == 0})
    v = {k: T0.get((k, c, 0), 0) for k in nodes}
    dem = sum(-x for x in v.values() if x < 0)
    if not dem:
        return dem, 0, None
    cap = {}
    for m, tv in T1.items():
        k, m2, m3 = m
        if m3 != c + 1 or not c <= m2 - 1 < k:
            continue
        cp = min(tv - D1.get(m, 0), T0.get((k, m2, 0), 0))
        if cp > 0:
            cap[(m2 - 1, k)] = cap.get((m2 - 1, k), 0) + cp
    G = DiGraph()
    for k, x in v.items():
        if x > 0:
            G.add_edge("s", k, x)
        if x < 0:
            G.add_edge(k, "t", 10 ** 18)
    for (j, k), cp in cap.items():
        G.add_edge(j, k, cp)
    deliver = G.flow("s", "t", value_only=True, use_edge_labels=True)
    run, minp = 0, None
    for k in nodes:
        run += v[k]
        minp = run if minp is None else min(minp, run)
    return dem, deliver, minp


def bosonic(D):
    """(E): all-full lift of B_d, d <= D (Fock series as in regge_monotone.py)."""
    import numpy as np
    from regge_monotone import adams, exp_log
    from derive_r9 import V8
    from weights import DomTable
    T4 = DomTable("D", 4, 4 * (D + 2) ** 2)
    TB = DomTable("B", 4, 4 * (D + 2) ** 2)
    b_to_d = np.array([T4.index[tuple(int(a) for a in w)] for w in TB.W] + [T4.ZERO])
    psiV = {k: adams(T4, V8, k) for k in range(1, D + 2)}
    GB = [None]
    for j in range(1, D + 2):
        g = T4.zeros()
        for k in range(1, j + 1):
            if j % k == 0:
                g += (j // k) * psiV[k]
        GB.append(g)
    FB = exp_log(T4, D + 1, GB)
    cum, out = T4.zeros(), []
    for d in range(D + 1):
        cum = cum + FB[d]
        if d >= 4:
            b = TB.decompose(cum[b_to_d])
            Bd = {tuple(int(a) for a in TB.W[i]): int(b[i]) for i in np.nonzero(b[:TB.size])[0]}
            a, _ = all_full(Bd, 0)
            neg = sorted({y[3] for y, w in a.items() if w < 0})
            out.append(f"  degree {d}: {len(Bd)} SO(9) irreps; all-full lift negative at levels {neg}, min {min(a.values())}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=DEFAULT_EXT)
    ap.add_argument("--nmax", type=int, default=None)
    ap.add_argument("--bosonic", type=int, default=None)
    args = ap.parse_args()
    if args.bosonic:
        lines = ["(E) bosonic analogue B: all-full lift by degree"] + bosonic(args.bosonic)
        print("\n".join(lines))
        os.makedirs(RESULTS, exist_ok=True)
        with open(os.path.join(RESULTS, "cumulative_lift_bosonic.txt"), "w") as f:
            f.write("\n".join(lines) + "\n")
        return
    N, r9, R9 = load_R9(args.data)
    NMAX = min(args.nmax or N, N)
    import half_b
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    log(f"data {os.path.basename(args.data)}, n <= {NMAX}")
    tails_seen = {}
    for n in range(2, NMAX + 1):
        parts = []
        for par in (0, 1):
            a, S = all_full(R9[n], par)
            if a is None:
                continue
            sec = "spinor" if par else "tensor"
            levels = sorted({y[3] for y, w in a.items() if w < 0})
            T, D, slacks, bad_bnd = forced(a)
            s = min(slacks) if slacks else None
            info = f"{sec}: a<0 at levels {levels} (min {min(a.values())})"
            info += f"; forced {len(slacks)}" + (f", min slack {s[0]} at t={s[1]} m={s[2]} d1={n - s[2][0]}" if s else "")
            info += f"; boundary negatives above bottom {len(bad_bnd)}"
            for c in ((0, 1) if par == 0 else (0,)):
                dem, deliver, minp = chain(T, D, c)
                if dem:
                    info += f"; chain c={c}: deficit {dem}, deliverable {deliver}, min prefix {minp}"
                tail, run = {}, 0
                for k in range(max((y[0] for y in a), default=0), -1, -1):
                    run += a.get((k, c, 0, 0), 0)
                    tail[n - k - par] = run
                tails_seen[(sec, c, n)] = [tail.get(d, 0) for d in range(1, 15)]
            parts.append(info)
        E = half_b.level_lift(R9[n], "ilp-min")
        ok = E is not None and restrict(E) == R9[n]
        log(f"n={n:2d}: uplift {'FOUND' if ok else 'NOT FOUND'} ({sum(E.values()) if ok else '-'} irreps) | " + " | ".join(parts))
    log("tail sums f by d = n - k (d = 1..14), largest n:")
    for (sec, c) in (("tensor", 0), ("tensor", 1), ("spinor", 0)):
        log(f"  {sec} chain c={c}: {tails_seen.get((sec, c, NMAX))}")
    log(f"[{time.time() - t0:.0f}s]")
    os.makedirs(RESULTS, exist_ok=True)
    name = "cumulative_lift.txt" if args.data == DEFAULT_EXT else f"cumulative_lift_{NMAX}.txt"
    with open(os.path.join(RESULTS, name), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
