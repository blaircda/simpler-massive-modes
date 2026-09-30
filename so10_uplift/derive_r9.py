"""
Derive the SO(9) content r^9_n of the one-sided massive IIA levels
from first principles, independently of the baseline scripts.

1. Light-cone Fock character F(q) = prod_m Sym_{q^m}(8_v) Lambda_{q^m}(8_s)
   via the exp-log recurrence
        n F_n = sum_{j=1}^n G_j F_{n-j},
        G_j = sum_{k|j} (j/k) [ psi^k(8_v) - (-1)^k psi^k(8_s) ],
   where psi^k are Adams operations (weights scaled by k).  All
   characters are stored as W(D4)-invariant dominant-weight multiplicities.

2. r^8_n = F_n / (8_v + 8_s), by exact long division on dominant weights
   ordered by the height h = (10,3,2,1), for which e1 is the unique
   highest weight of 8_v + 8_s.  Exactness is verified by re-multiplying.

3. SO(9): Spin(8) and Spin(9) share a maximal torus, so r^9_n has the same
   weights as r^8_n.  W(B4)-invariance (x4 -> -x4) is checked, then the
   B4 decomposition is the Weyl alternation sum.  Non-negativity checked.

Output: data/r9.json  (r^8_n as D4 irreps and r^9_n as B4 irreps).
"""

import argparse
import json
import os
import sys
import time
from itertools import product

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from weights import DomTable, dynkin_B4, dynkin_D4, dim_irrep

HERE = os.path.dirname(os.path.abspath(__file__))

V8 = [tuple(2 * s if i == k else 0 for i in range(4)) for k in range(4) for s in (1, -1)]
S8 = [w for w in product((1, -1), repeat=4) if w.count(-1) % 2 == 0]
E1 = (2, 0, 0, 0)
HEIGHT = np.array([10, 3, 2, 1])


def fock_dims(N):
    """Independent integer check: coefficients of prod_m ((1+q^m)/(1-q^m))^8."""
    c = [1] + [0] * N
    for m in range(1, N + 1):
        for _ in range(8):
            for i in range(N, m - 1, -1):      # multiply by (1+q^m)
                c[i] += c[i - m]
            for i in range(m, N + 1):          # divide by (1-q^m)
                c[i] += c[i - m]
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-level", type=int, default=12)
    ap.add_argument("--out", default=os.path.join(HERE, "data", "r9.json"))
    args = ap.parse_args()
    N = args.max_level

    t0 = time.time()
    T = DomTable("D", 4, 4 * N * N)
    print(f"D4 dominant table: {T.size} weights (|x| <= {N})", flush=True)

    # ---------------- 1. Fock space ----------------
    F = [T.zeros()]
    F[0][T.index[(0, 0, 0, 0)]] = 1
    limit = np.iinfo(np.int64).max // 64
    for n in range(1, N + 1):
        acc = T.zeros()
        for j in range(1, n + 1):
            Fm = F[n - j]
            for k in range(1, j + 1):
                if j % k:
                    continue
                c = j // k
                sv = sum(Fm[T.gather(tuple(k * a for a in w))] for w in V8)
                ss = sum(Fm[T.gather(tuple(k * a for a in w))] for w in S8)
                acc += c * (sv - (-1) ** k * ss)
        assert np.all(acc % n == 0), f"exp-log recurrence not divisible at level {n}"
        Fn = acc // n
        Fn[T.ZERO] = 0
        assert np.abs(Fn).max() < limit, "int64 headroom exhausted"
        assert Fn.min() >= 0
        F.append(Fn)
    print(f"Fock space through level {N}: {time.time()-t0:.1f}s", flush=True)

    # orbit sizes for dimension checks
    orbit = np.array([orbit_size_D4(tuple(int(a) for a in w)) for w in T.W] + [0])
    fd = fock_dims(N)
    for n in range(N + 1):
        assert int((F[n] * orbit).sum()) == fd[n], f"Fock dimension mismatch at level {n}"

    # ---------------- 2. divide by 8_v + 8_s ----------------
    D = V8 + S8
    others = [d for d in D if d != E1]
    g_up = T.gather(tuple(-a for a in E1))                     # w + e1
    g_oth = [T.gather(tuple(d[i] - E1[i] for i in range(4))) for d in others]
    g_D = [T.gather(d) for d in D]
    order = np.argsort(-(T.W @ HEIGHT), kind="stable")
    rank = np.empty(T.size + 1, dtype=np.int64)
    rank[order] = np.arange(T.size)
    rank[T.ZERO] = -1
    for g in g_oth:  # strict triangularity of the division
        sel = g[:T.size] != T.ZERO
        assert np.all(rank[g[:T.size]][sel] < rank[:T.size][sel])

    g_up_l = g_up.tolist()
    g_oth_l = [g.tolist() for g in g_oth]

    TB = DomTable("B", 4, 4 * N * N)
    b_to_d = np.array([T.index[tuple(int(a) for a in w)] for w in TB.W] + [T.ZERO])

    levels = {}
    for n in range(1, N + 1):
        Fl = F[n].tolist()
        r = [0] * (T.size + 1)
        for i in order.tolist():
            v = Fl[g_up_l[i]]
            for g in g_oth_l:
                v -= r[g[i]]
            r[i] = v
        r = np.array(r, dtype=np.int64)
        back = sum(r[g] for g in g_D)
        back[T.ZERO] = 0
        assert np.array_equal(back, F[n]), f"(8v+8s) r_n != F_n at level {n}"
        assert r.min() >= 0, f"negative weight multiplicity in r_{n}"
        assert 16 * int((r * orbit).sum()) == fd[n]

        # SO(8) irreps
        d4 = T.decompose(r)
        assert d4.min() >= 0, f"negative D4 multiplicity in r^8_{n}"

        # SO(9): x4 -> -x4 symmetry, then B4 alternation
        flip = T.lookup(T.W * np.array([1, 1, 1, -1]))
        assert np.array_equal(r[:T.size], r[flip]), f"r_{n} not W(B4)-invariant"
        rB = r[b_to_d]
        b4 = TB.decompose(rB)
        assert b4.min() >= 0, f"negative B4 multiplicity in r^9_{n}"

        r8 = irrep_list(T, d4, "D", dynkin_D4)
        r9 = irrep_list(TB, b4, "B", dynkin_B4)
        dimr = int((r * orbit).sum())
        assert sum(e["mult"] * e["dim"] for e in r8) == dimr
        assert sum(e["mult"] * e["dim"] for e in r9) == dimr
        levels[str(n)] = {"dim": dimr, "r8_D4": r8, "r9_B4": r9}
        print(f"level {n:2d}: dim r = {dimr}, #D4 irreps = {len(r8)}, "
              f"#B4 irreps = {len(r9)}   [{time.time()-t0:.1f}s]", flush=True)

    out = {
        "meta": {
            "description": "SO(8) content r^8_n and SO(9) content r^9_n of one-sided "
                           "massive IIA level n, F_n = (8_v+8_s) r_n.",
            "generator": "so10_uplift/derive_r9.py",
            "max_level": N,
            "conventions": "ambient2 = doubled orthogonal highest weight; dynkin = Sage "
                           "WeylCharacterRing(style='coroots') labels; 8_s = D4(0,0,0,1).",
        },
        "levels": levels,
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {args.out}  [{time.time()-t0:.1f}s]")


def irrep_list(T, arr, typ, dynkin):
    out = []
    for i in np.nonzero(arr[:T.size])[0]:
        w = tuple(int(a) for a in T.W[i])
        out.append({"mult": int(arr[i]), "ambient2": list(w), "dynkin": list(dynkin(w)),
                    "dim": dim_irrep(w, typ, T.n)})
    out.sort(key=lambda e: (sum(e["ambient2"]), e["ambient2"]))
    return out


def orbit_size_D4(w):
    from itertools import permutations
    s = set()
    for p in permutations(w):
        for sg in product((1, -1), repeat=4):
            if sg.count(-1) % 2 == 0:
                s.add(tuple(p[i] * sg[i] for i in range(4)))
    return len(s)


if __name__ == "__main__":
    main()
