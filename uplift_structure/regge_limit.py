"""
S15 (model): the Regge limit of the SO(9) spectrum r9_n, from the light-cone Fock
space alone.  No spectrum data are read here; the comparison with the data is
regge_scaling.py.

Coordinates.  An SO(9) weight nu = (x1, mu), mu = (a, b, c) an SO(7) weight, at level
n has top distance D = n - x1.  GMSZ's support cost n0 = l + 3m + c(r,s) splits as
    n0(nu) = x1 + cost7(mu),   cost7(a,b,c) = 2a - 3b + c(b - c, 2c),
so the GMSZ depth is w = n - n0 = D - cost7(mu).

Model.  In the Regge limit x1 -> oo at fixed D, r9_n(x1, mu) -> m_D(mu), the
multiplicity of the SO(7) irrep mu in [q^D] of

    Z_Regge = q * prod_{k>=2} 1/((1-q^k) det_7(1-q^k g)) * prod_{k>=1} det_8(1+q^{k+1/2} g),

i.e. bosons in 1+7 at levels k >= 2 and SO(7)-spinor fermions at levels k+1/2 >= 3/2
(derivation in STATUS.md, S15).  Since 1+7 = Res 8v and 8 = Res 8s (D4 > B3),
Z_Regge is the restriction of the SO(8) Fock space with the same oscillators.

Parts:
  (1) Fock dimensions by D;
  (2) support: least D at which each SO(7) irrep appears, against cost7 (GMSZ
      c(r,s)), and the zeros above it;
  (3) the Regge limit of the all-full lift: a_D(y) = sum_{eps in {0,1}^3} (-1)^|eps|
      m_D(y + eps) (tensor sector), by column family (j,i) at the bottom level:
      the D with a_D(j,i,0) < 0 and the offset d(j,i) = least such D; the excess
      e_D(j,i) = -a_D(j,i,0) - sum_{L=1..i} a_D(j,i,L) and the least D with e > 0
      (level i+1 needed by the same-column count);
  (4) writes the multiplicities m_D(mu) (mu_1 <= LM) to data/regge_limit_D<Dmax>.json;
  (0) checks the one non-standard step of the derivation on weights: (8v + 8c) x
      (level-1 fermions of charge +1/2) = q T(q), T(q) = 256|_{SO(2) x SO(7)} at z = 1/q;
  (5) the SO(8) lift: decomposes Z8 = q prod det_8v(1-q^k g)^-1 det_8s(1+q^{k+1/2} g)
      into D4 irreps for D <= D8, checks Res^{D4}_{B3} [q^D] Z8 = [q^D] Z_Regge irrep by
      irrep (interlacing), and lists its non-full classes by (L, J) = (l4 - floor,
      |l5| - floor) (the repair types of S13).

Usage: regge_limit.py [--dmax 60] [--lm 10] [--d8 20]   Output: results/regge_limit.txt
"""

import argparse
import json
import os
import sys
import time
from itertools import permutations, product

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS

HERE = os.path.dirname(os.path.abspath(__file__))
RHO2 = (5, 3, 1)                       # 2 rho(B3)


def c_rs(r, s):
    """GMSZ (3.6), as in gmsz_compare.py (not imported: that module loads the spectrum data)."""
    return max(6 * r + 5 * s - 6, 6 * r + 4 * s - 2, 6 * r + 3 * s - 1, 2 * s + 1)


def cost7_t(A, B, C):
    """2*cost7 for the SO(7) weight (A, B, C)/2 (doubled coordinates)."""
    return 2 * A - 3 * B + 2 * c_rs((B - C) // 2, C)


def t_identity():
    """(0): weights doubled; a monomial is (2 * q-exponent, SO(7) weight)."""
    from collections import Counter
    v9 = [tuple(s * 2 * (j == i) for j in range(4)) for i in range(4) for s in (1, -1)] + [(0, 0, 0, 0)]
    sp = list(product((1, -1), repeat=4))
    add = lambda *ws: tuple(sum(c) for c in zip(*ws))
    w256 = [add(v9[a], v9[b]) for a in range(9) for b in range(a, 9)]
    w256.remove((0, 0, 0, 0))                                            # 44
    w256 += [add(v9[a], v9[b], v9[c]) for a in range(9) for b in range(a + 1, 9) for c in range(b + 1, 9)]   # 84
    vs = Counter(add(v, w) for v in v9 for w in sp)
    vs.subtract(Counter(sp))
    w256 += list(vs.elements())                                          # 128
    rhs = Counter((2 - W[0], W[1:]) for W in w256)                      # q * T(q): q^(1 - x1)
    v8 = [w for w in v9 if any(w)]
    c8 = [w for w in sp if w.count(-1) % 2 == 1]
    f = [w for w in sp if w.count(-1) % 2 == 0 and w[0] == 1]           # 8s, charge +1/2: depth 1/2
    lhs = Counter()
    for g in v8 + c8:
        for sub in product((0, 1), repeat=4):
            W = add(g, *[f[i] for i in range(4) if sub[i]]) if any(sub) else g
            lhs[(-g[0] + sum(sub), W[1:])] += 1
    return len(w256), lhs == rhs


def fock(Dmax, log):
    """Weight multiplicities of Z_Regge: Z[t] (t = 2D) on p-coordinates, mu = p + (t % 2)/2."""
    T, P = 2 * Dmax, Dmax // 2 + 3
    S = 2 * P + 1
    Z = np.zeros((T + 1, S, S, S), dtype=np.int64)
    Z[2, P, P, P] = 1                                            # the prefactor q (D = 1)

    def sl(d):
        return {1: (slice(1, None), slice(None, -1)), -1: (slice(None, -1), slice(1, None)),
                0: (slice(None), slice(None))}[d]

    def add(td, ts, delta):
        s0, s1, s2 = (sl(d) for d in delta)
        Z[td][s0[0], s1[0], s2[0]] += Z[ts][s0[1], s1[1], s2[1]]

    t0 = time.time()
    bos = [(0, 0, 0), (0, 0, 0)] + [tuple(s if j == i else 0 for j in range(3)) for i in range(3) for s in (1, -1)]
    for k in range(2, Dmax):                                     # 1/(1 - q^k x^w): increasing t
        for d in bos:
            for t in range(2 * k, T + 1):
                add(t, t - 2 * k, d)
    for k in range(1, Dmax):                                     # (1 + q^{k+1/2} x^w): decreasing t
        s = 2 * k + 1
        for eps in product((1, -1), repeat=3):
            for t in range(T, s - 1, -1):
                par = (t - s) % 2                                # parity of the source slice
                add(t, t - s, tuple((e - 1) // 2 if par == 0 else (e + 1) // 2 for e in eps))
    edge = max(int(np.abs(Z[:, [0, -1], :, :]).max()), int(np.abs(Z[:, :, [0, -1], :]).max()),
               int(np.abs(Z[:, :, :, [0, -1]]).max()))
    assert edge == 0, "box too small"
    assert Z.min() >= 0
    log(f"  Fock space built to D = {Dmax} in {time.time() - t0:.0f}s (box |p| <= {P}, boundary empty)")
    return Z, P


def weyl_b3():
    out = []
    for perm in permutations(range(3)):
        inv = sum(1 for i in range(3) for j in range(i + 1, 3) if perm[i] > perm[j])
        for sg in product((1, -1), repeat=3):
            out.append((tuple(sg[i] * RHO2[perm[i]] for i in range(3)), (-1) ** inv * sg[0] * sg[1] * sg[2]))
    return out


def irreducible(Z, P, Dmax, LM):
    """m[t][(A,B,C)] for dominant doubled weights with A <= 2*LM (Racah-Speiser)."""
    W = weyl_b3()
    m = {}
    for t in range(2, 2 * Dmax + 1):
        par = t % 2
        lams = [(A, B, C) for A in range(par, 2 * LM + 1, 2) for B in range(par, A + 1, 2) for C in range(par, B + 1, 2)]
        vals = np.zeros(len(lams), dtype=np.int64)
        S = Z.shape[1]
        for wr, sg in W:
            idx = np.array([[(l[i] + RHO2[i] - wr[i] - par) // 2 + P for l in lams] for i in range(3)])
            ok = ((idx >= 0) & (idx < S)).all(axis=0)                  # outside the box: weight absent
            v = Z[t][tuple(np.clip(idx[i], 0, S - 1) for i in range(3))]
            vals += sg * np.where(ok, v, 0)
        m[t] = {l: int(v) for l, v in zip(lams, vals) if v}
        assert all(v >= 0 for v in m[t].values())
    return m


def so8_lift(D8, m, log):
    """(5): Z8 on the D4 torus (p-coordinates as in fock), D4 irreps, restriction, non-full classes."""
    from collections import Counter
    T, P = 2 * D8, D8 // 2 + 3
    S = 2 * P + 1
    Z = np.zeros((T + 1, S, S, S, S), dtype=np.int64)
    Z[2, P, P, P, P] = 1

    def sl(d):
        return {1: (slice(1, None), slice(None, -1)), -1: (slice(None, -1), slice(1, None)),
                0: (slice(None), slice(None))}[d]

    def add(td, ts, delta):
        ss = [sl(d) for d in delta]
        Z[td][tuple(x[0] for x in ss)] += Z[ts][tuple(x[1] for x in ss)]

    for k in range(2, D8):
        for i in range(4):
            for sg in (1, -1):
                d = tuple(sg * (j == i) for j in range(4))
                for t in range(2 * k, T + 1):
                    add(t, t - 2 * k, d)
    for k in range(1, D8):
        s = 2 * k + 1
        for eps in product((1, -1), repeat=4):
            if eps.count(-1) % 2:
                continue                                                 # 8s
            for t in range(T, s - 1, -1):
                par = (t - s) % 2
                add(t, t - s, tuple((e - 1) // 2 if par == 0 else (e + 1) // 2 for e in eps))
    assert all(int(np.abs(np.take(Z, [0, -1], axis=ax)).max()) == 0 for ax in (1, 2, 3, 4))
    rho = (6, 4, 2, 0)                                                   # 2 rho(D4)
    W = []
    for perm in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if perm[i] > perm[j])
        for sg in product((1, -1), repeat=4):
            if sg.count(-1) % 2 == 0:
                W.append((tuple(sg[i] * rho[perm[i]] for i in range(4)), (-1) ** inv))
    agree, types = True, {}
    for t in range(2, T + 1):
        par, R = t % 2, P - 1
        lams = [(A, B, C, E) for A in range(par, 2 * R + 1, 2) for B in range(par, A + 1, 2)
                for C in range(par, B + 1, 2) for E in range(-C, C + 1, 2)]
        vals = np.zeros(len(lams), dtype=np.int64)
        for wr, sg in W:
            idx = [[(l[i] + rho[i] - wr[i] - par) // 2 + P for l in lams] for i in range(4)]
            ok = [all(0 <= idx[i][j] < S for i in range(4)) for j in range(len(lams))]
            v = Z[t][tuple(np.clip(idx[i], 0, S - 1) for i in range(4))]
            vals += sg * np.where(ok, v, 0)
        m8 = {l: int(v) for l, v in zip(lams, vals) if v}
        assert all(v > 0 for v in m8.values())
        res = Counter()
        for l, c in m8.items():
            for a in range(l[1], l[0] + 1, 2):
                for b in range(l[2], l[1] + 1, 2):
                    for cc in range(abs(l[3]), l[2] + 1, 2):
                        res[(a, b, cc)] += c
        agree &= dict(res) == m[t]
        tc = Counter()
        for l, c in m8.items():
            L, J = (l[2] - par) // 2, (abs(l[3]) - par) // 2
            if J:
                tc[(L, J)] += c
        types[t] = tc
    log(f"(5) SO(8) lift Z8 (8v bosons k >= 2, 8s fermions k+1/2 >= 3/2): decomposed into D4 irreps for D <= {D8}; "
        f"non-negative; Res^D4_B3 [q^D] Z8 = [q^D] Z_Regge irrep by irrep: {agree}")
    log("    its non-full classes, (L, J): count, by D (tensor sector, integer D):")
    for t in range(2, T + 1, 2):
        if types[t]:
            log(f"      D={t // 2}: " + ", ".join(f"{k}: {v}" for k, v in sorted(types[t].items())))
    first = {}
    for t in range(2, T + 1):
        for k in types[t]:
            first.setdefault(k, t / 2)
    log("    first D of each non-full type (L, J), both sectors: " + ", ".join(f"{k}: {v:g}" for k, v in sorted(first.items())))
    assert agree


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dmax", type=int, default=60)
    ap.add_argument("--lm", type=int, default=10)
    ap.add_argument("--d8", type=int, default=20)
    args = ap.parse_args()
    Dmax, LM = args.dmax, args.lm
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    log(f"Regge-limit Fock space Z_Regge to D = {Dmax} (SO(7) irreps with mu_1 <= {LM})")
    nw, ok = t_identity()
    log(f"(0) (8v + 8c) x (level-1 fermions of charge +1/2) = q T(q) on the torus (256 = 44+84+128, {nw} weights): {ok}")
    assert ok
    Z, P = fock(Dmax, log)
    log("(1) Fock dimensions: " + ", ".join(f"D={D}: {int(Z[2 * D].sum()):.3e}" for D in (10, 20, 30, 40, 50, 60) if D <= Dmax))
    m = irreducible(Z, P, Dmax, LM)
    del Z

    # (2) support against GMSZ
    first, zeros, bad = {}, [], []
    for t in sorted(m):
        for l in m[t]:
            first.setdefault(l, t)
    for A in range(0, 2 * LM + 1):
        for B in range(A % 2, A + 1, 2):
            for C in range(A % 2, B + 1, 2):
                l = (A, B, C)
                ct = cost7_t(*l)
                if ct > 2 * Dmax:
                    continue
                if first.get(l) != ct:
                    bad.append((l, first.get(l), ct))
                for t in range(ct, 2 * Dmax + 1, 2):
                    if not m[t].get(l):
                        zeros.append((l, t))
    log(f"(2) least D of each SO(7) irrep mu (mu_1 <= {LM}, cost7 <= {Dmax}) vs cost7(mu) = 2a - 3b + c(b-c, 2c): "
        f"{'all equal' if not bad else f'{len(bad)} differ: {bad[:6]}'}; zero multiplicities at D >= cost7 "
        f"(D <= {Dmax}): " + (", ".join(f"mu={tuple(x / 2 for x in l)} D={t / 2:g}" for l, t in zeros) or "none"))

    # (3) Regge limit of the all-full lift (tensor sector), by family
    def mD(D, y):
        if not (y[0] >= y[1] >= y[2] >= 0):
            return 0
        return m[2 * D].get(tuple(2 * v for v in y), 0)

    def a(D, y):
        return sum((-1) ** sum(e) * mD(D, tuple(p + q for p, q in zip(y, e))) for e in product((0, 1), repeat=3))

    log("(3) a_D(j,i,L) = Regge limit of the all-full lift at columns (n-D, j, i), level L.  Per family (j,i):")
    log("    D with a_D(j,i,0) < 0 [offset d = least such D] | D with excess e_D(j,i) > 0 [-> level i+1]")
    fam_rows = {}
    for j in range(0, LM - 1):
        for i in range(0, j + 1):
            negD = [D for D in range(1, Dmax + 1) if a(D, (j, i, 0)) < 0]
            excD = [D for D in range(1, Dmax + 1) if i >= 1 and -a(D, (j, i, 0)) - sum(a(D, (j, i, L)) for L in range(1, i + 1)) > 0]
            fam_rows[(j, i)] = (negD, excD)
            if negD or excD:
                rng = lambda xs: (f"{xs[0]}..{xs[-1]}" + ("" if xs == list(range(xs[0], xs[-1] + 1)) else " (gaps)")) if xs else "-"
                log(f"    ({j},{i}): neg D = {rng(negD)} [d = {negD[0] if negD else '-'}] | excess D = {rng(excD)}")
    log("    (1,1) family near the edge: D: a_D(1,1,0), a_D(1,1,1), e_D")
    for D in range(24, min(Dmax, 40) + 1):
        a0, a1 = a(D, (1, 1, 0)), a(D, (1, 1, 1))
        log(f"      D={D}: {a0:+d} {a1:+d} e={-a0 - a1:+d}")
    so8_lift(args.d8, m, log)
    os.makedirs(RESULTS, exist_ok=True)
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    with open(os.path.join(HERE, "data", f"regge_limit_D{Dmax}.json"), "w") as f:
        json.dump({"Dmax": Dmax, "LM": LM, "note": "m[t][A,B,C]: multiplicity of the SO(7) irrep (A,B,C)/2 in [q^(t/2)] Z_Regge",
                   "m": {str(t): [[list(l), v] for l, v in sorted(d.items())] for t, d in m.items()}}, f)
    with open(os.path.join(RESULTS, "regge_limit.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
