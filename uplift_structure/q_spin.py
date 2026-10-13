"""
S16 (model): the spin-resolved Regge function Q(q, y), from the light-cone Fock
space alone.  No spectrum data are read here; q_spin_check.py compares with the data.

Definition.  For an SO(9) weight (N, mu) at level n (spin N = x1, top distance D = n - N),
    Q(q, y) = (1 - y) sum_{n>=1} q^n K(r9_n)(z = y/q),   K(chi) = det_7(1 - z^-1 g) chi,
with z the x1-charge.  By the parabolic BGG resolution (S15), for every SO(9)-dominant
(N, mu) and D >= 1,
    r9_{N+D}(N, mu) = sum_{e <= N} c(D, e; mu),   c(D, e; mu) = mult of the SO(7) irrep mu in [q^D y^e] Q.
So for N > mu_1, c(D, N; mu) = r9_n(N, mu) - r9_{n-1}(N-1, mu) counts the trajectories at
(D, mu) that start at spin N; the sum over all e is the plateau, Q(q, 1) = Z_Regge (S15);
and Q to depth Dmax gives the top Dmax units of every level n.

Construction.  With z = y/q, K applied to the light-cone character is H(q, y)/(1 - y):
1/(1 - y) is the alpha^+_{-1} tower, and
    H = det_7(1 - (q/y) g) (8v + 8c)(z = y/q) P,
    P = prod  alpha^+_{-k} (k >= 2): 1/(1 - q^{k-1} y);   alpha^-_{-k}: 1/(1 - q^{k+1}/y);
              alpha^i_{-k} (6v): 1/(1 - q^k x^w);
              S_{-k} (8s): (1 + q^{k-1/2} y^{1/2} x^w) on the 4, (1 + q^{k+1/2} y^{-1/2} x^w) on the 4bar.
Level 0 (not SO(9)-covariant) is removed by P -> P - 1 + y.  Dividing by
256(z = y/q) = (y/q)^2 (1 + R) gives Q = q y^-2 V with V = q (8v + 8c) det_7 (P - 1 + y) / (1 + R).

Arithmetic is exact modulo two primes (int64 storage), combined by CRT.

Parts:
  (0) checks: every (D, e) slice of Q is W(B3)-invariant; sum_e c(D, e; mu) = m_D(mu) of S15;
  (1) Regge monotonicity at top distances D <= Dmax, for all levels: the starts c(D, e; mu), e > mu_1;
  (2) saturation: the last start, l_sat = e_last - mu_1, against the GMSZ depth w = D - cost7(mu);
  (3) start distributions for small shapes, and the top end (last starts) as D grows;
  writes data/q_spin_D<Dmax>.json.

Usage: q_spin.py [--dmax 24] [--lm mu_1 bound, default Dmax]   Output: results/q_spin.txt
"""

import argparse
import json
import os
import sys
import time
from collections import Counter
from itertools import product

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS
from regge_limit import cost7_t, weyl_b3

HERE = os.path.dirname(os.path.abspath(__file__))
RHO2 = (5, 3, 1)


def is_prime(n):
    if n < 2:
        return False
    for a in (2, 3, 5, 7, 11, 13):
        if n % a == 0:
            return n == a
    d, r = n - 1, 0
    while d % 2 == 0:
        d, r = d // 2, r + 1
    for a in (2, 3, 5, 7, 11, 13):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


P1 = 2147483647
P2 = next(m for m in range(P1 - 2, 0, -2) if is_prime(m))


def crt(a1, a2):
    M = P1 * P2
    x = (a1 + P1 * ((a2 - a1) * pow(P1, -1, P2) % P2)) % M
    return x - M if x > M // 2 else x


class Series:
    """Slices Z[t], t = 2*depth, of shape (charge, p1, p2, p3): doubled charge s = smin(t) + 2i,
    SO(7) weight mu = p + (t % 2)/2 with p in [-B(t), B(t)].  Entries mod p."""

    def __init__(self, T, p):
        self.T, self.p = T, p
        self.B = [t // 2 + 3 for t in range(T + 1)]
        self.smin = [-t - 8 for t in range(T + 1)]
        self.Z = [np.zeros((t + 9,) + (2 * self.B[t] + 1,) * 3, dtype=np.int64) for t in range(T + 1)]

    def acc(self, D, td, S, ts, ds, W2, coef=1):
        """D (slice at doubled depth td) += coef * S (slice at ts) shifted by doubled charge ds, weight W2/2."""
        offs = [(td - ts + ds) // 2] + [(W2[k] + ts % 2 - td % 2) // 2 + self.B[td] - self.B[ts] for k in range(3)]
        sd, ss = [], []
        for o, Ls, Ld in zip(offs, S.shape, D.shape):
            a, b = max(0, -o), min(Ls, Ld - o)
            if a >= b:
                return
            ss.append(slice(a, b))
            sd.append(slice(a + o, b + o))
        view = D[tuple(sd)]
        if coef == 1:
            view += S[tuple(ss)]
        elif coef == -1:
            view -= S[tuple(ss)]
        else:
            view += coef * S[tuple(ss)]
        np.remainder(view, self.p, out=view)

    def geom(self, d, ds, W2):                       # 1/(1 - q^{d/2} y^{ds/2} x^{W2/2})
        for t in range(d, self.T + 1):
            self.acc(self.Z[t], t, self.Z[t - d], t - d, ds, W2)

    def binom(self, d, ds, W2, coef):                # (1 + coef q^{d/2} y^{ds/2} x^{W2/2})
        for t in range(self.T, d - 1, -1):
            self.acc(self.Z[t], t, self.Z[t - d], t - d, ds, W2, coef)

    def boundary_empty(self):
        for Zt in self.Z:
            for ax in range(4):
                if np.take(Zt, [0, -1], axis=ax).any():
                    return False
        return True


def w256():
    """Doubled SO(9) weights of 256 = 44 + 84 + 128, as a list with multiplicity."""
    v9 = [tuple(s * 2 * (j == i) for j in range(4)) for i in range(4) for s in (1, -1)] + [(0, 0, 0, 0)]
    sp = list(product((1, -1), repeat=4))
    add = lambda *ws: tuple(sum(c) for c in zip(*ws))
    out = [add(v9[a], v9[b]) for a in range(9) for b in range(a, 9)]
    out.remove((0, 0, 0, 0))
    out += [add(v9[a], v9[b], v9[c]) for a in range(9) for b in range(a + 1, 9) for c in range(b + 1, 9)]
    vs = Counter(add(v, w) for v in v9 for w in sp)
    vs.subtract(Counter(sp))
    out += list(vs.elements())
    assert len(out) == 256
    return out


def build(Dmax, p, log):
    T = 2 * (Dmax - 1)                                           # V to depth Dmax - 1; Q = q y^-2 V
    S = Series(T, p)
    t0 = time.time()
    b0, i0 = S.B[0], (0 - S.smin[0]) // 2
    S.Z[0][i0, b0, b0, b0] = 1
    zero = (0, 0, 0)
    six = [tuple(2 * s * (j == i) for j in range(3)) for i in range(3) for s in (1, -1)]
    four = [w for w in product((1, -1), repeat=3) if w.count(-1) % 2 == 0]       # 8s at charge +1/2
    fourbar = [w for w in product((1, -1), repeat=3) if w.count(-1) % 2 == 1]    # 8s at charge -1/2
    for k in range(2, Dmax + 1):                                 # alpha^+_{-k}
        if 2 * (k - 1) <= T:
            S.geom(2 * (k - 1), 2, zero)
    for k in range(1, Dmax + 1):                                 # alpha^-_{-k}
        if 2 * (k + 1) <= T:
            S.geom(2 * (k + 1), -2, zero)
    for k in range(1, Dmax + 1):                                 # transverse alpha^i_{-k}
        if 2 * k <= T:
            for w in six:
                S.geom(2 * k, 0, w)
    for k in range(1, Dmax + 1):                                 # fermions S_{-k}
        if 2 * k - 1 <= T:
            for w in four:
                S.binom(2 * k - 1, 1, w, 1)
        if 2 * k + 1 <= T:
            for w in fourbar:
                S.binom(2 * k + 1, -1, w, 1)
    log(f"    P built ({time.time() - t0:.0f}s)")
    S.Z[0][i0, b0, b0, b0] = (S.Z[0][i0, b0, b0, b0] - 1) % p    # P -> P - 1 + y
    S.Z[0][i0 + 1, b0, b0, b0] = (S.Z[0][i0 + 1, b0, b0, b0] + 1) % p
    for w in [zero] + six:                                       # det_7(1 - (q/y) g)
        S.binom(2, -2, w, -1)
    ground = [(0, 2, zero), (4, -2, zero)] + [(2, 0, w) for w in six] \
        + [(1, 1, w) for w in fourbar] + [(3, -1, w) for w in four]                # q (8v + 8c) at z = y/q
    for t in range(T, -1, -1):
        new = np.zeros_like(S.Z[t])
        for d, ds, w in ground:
            if d <= t:
                S.acc(new, t, S.Z[t - d], t - d, ds, w)
        S.Z[t] = new
    R = Counter()
    for W in w256():
        if W != (4, 0, 0, 0):
            R[(4 - W[0], W[0] - 4, W[1:])] += 1                  # (q/y)^{2 - x1} x^mu
    R = sorted(R.items())
    for t in range(T + 1):                                       # V = H' / (1 + R), increasing depth
        for (d, ds, W2), c in R:
            if d <= t:
                S.acc(S.Z[t], t, S.Z[t - d], t - d, ds, W2, -c)
    assert S.boundary_empty(), "box too small"
    log(f"    Q built mod {p} ({time.time() - t0:.0f}s, boundary empty)")
    return S


def b3_violations(S):
    bad = 0
    for t in range(S.T + 1):
        par, L = t % 2, 2 * S.B[t] + 1
        for A in S.Z[t]:
            if not A.any():
                continue
            ok = np.array_equal(A, A.transpose(1, 0, 2)) and np.array_equal(A, A.transpose(0, 2, 1))
            if par == 0:
                ok = ok and np.array_equal(A, A[:, :, ::-1])
            else:
                ok = ok and np.array_equal(A[:, :, :L - 1], A[:, :, L - 2::-1]) and not A[:, :, L - 1].any()
            bad += not ok
    return bad


def irreps(S, LM):
    """Racah-Speiser on every slice of V.  Returns, mod p, the lumped edge values
    edge[(tQ, lam2)] = sum_{s <= lam_1} c and the starts st[(tQ, sQ, lam2)] = c for s > lam_1."""
    W = weyl_b3()
    sg = np.array([g for _, g in W])[:, None]
    edge, st = Counter(), {}
    for t in range(S.T + 1):
        par, Bt = t % 2, S.B[t]
        L = 2 * Bt + 1
        lams = [(A, B, C) for A in range(par, 2 * LM + 1, 2) for B in range(par, A + 1, 2) for C in range(par, B + 1, 2)]
        idx = np.array([[[(l[i] + RHO2[i] - wr[i] - par) // 2 + Bt for i in range(3)] for l in lams] for wr, _ in W])
        ok = ((idx >= 0) & (idx < L)).all(axis=2)
        cl = np.clip(idx, 0, L - 1)
        for i, A in enumerate(S.Z[t]):
            if not A.any():
                continue
            v = np.where(ok, A[cl[..., 0], cl[..., 1], cl[..., 2]], 0)
            vals = (sg * v).sum(axis=0) % S.p
            s = S.smin[t] + 2 * i - 4                            # charge of Q = charge of V - 2
            for l, x in zip(lams, vals):
                if x:
                    if s <= l[0]:
                        edge[(t + 2, l)] = (edge[(t + 2, l)] + int(x)) % S.p
                    else:
                        st[(t + 2, s, l)] = int(x)
    return edge, st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dmax", type=int, default=24)
    ap.add_argument("--lm", type=int, default=0, help="largest mu_1 (default Dmax: every shape that occurs)")
    args = ap.parse_args()
    Dmax = args.dmax
    LM = args.lm or Dmax
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    T = 2 * (Dmax - 1)
    size = sum((t + 9) * (2 * (t // 2 + 3) + 1) ** 3 for t in range(T + 1))
    log(f"Q(q, y) to depth D = {Dmax} (SO(7) irreps mu_1 <= {LM}); primes {P1}, {P2}; {size * 8 / 1e6:.0f} MB per prime")
    res, bad = {}, {}
    for p in (P1, P2):
        S = build(Dmax, p, log)
        bad[p] = b3_violations(S)
        res[p] = irreps(S, LM)
        del S
    edge, st = {}, {}
    for k in set(res[P1][0]) | set(res[P2][0]):
        v = crt(res[P1][0].get(k, 0), res[P2][0].get(k, 0))
        if v:
            edge[k] = v
    for k in set(res[P1][1]) | set(res[P2][1]):
        v = crt(res[P1][1].get(k, 0), res[P2][1].get(k, 0))
        if v:
            st.setdefault(k[:1] + k[2:], {})[k[1]] = v
    tot = Counter(edge)
    for k, row in st.items():
        tot[k] += sum(row.values())
    log(f"(0) slices of Q that are not W(B3)-invariant: {bad[P1]} (mod {P1}), {bad[P2]} (mod {P2})")
    import regge_limit
    Zr, Pr = regge_limit.fock(Dmax, log)
    mR = regge_limit.irreducible(Zr, Pr, Dmax, LM)
    del Zr
    keys = {(t, l) for t in mR for l in mR[t]} | {k for k, v in tot.items() if v}
    mism = [k for k in keys if tot.get(k, 0) != mR.get(k[0], {}).get(k[1], 0)]
    log(f"    sum over spin of Q = Z_Regge (S15) for all D <= {Dmax} and all shapes (mu_1 <= {LM}): "
        f"{'yes' if not mism else f'NO, {len(mism)} mismatches, e.g. {sorted(mism)[:4]}'}")

    # (1) starts above the dominance edge; lumped edge values
    allst = [v for row in st.values() for v in row.values()]
    neg = sorted((k, s, v) for k, row in st.items() for s, v in row.items() if v < 0)
    log(f"(1) starts c(D, e; mu), e > mu_1, at D <= {Dmax} (every level n): {len(allst)} non-zero, {len(neg)} negative"
        + (f", e.g. {neg[:4]}" if neg else " -- Regge monotonicity holds at these depths for every n"))
    log(f"    lumped edge values (= r9 at l = 0): {sum(1 for v in edge.values() if v < 0)} negative of {len(edge)}; "
        f"shapes with zero plateau but a non-zero edge value or start: "
        f"{sum(1 for k in set(edge) | set(st) if not tot.get(k))}")

    # (2) saturation
    sat = Counter()
    for (t, l), v in tot.items():
        if v:
            lsat = (max(st.get((t, l), {l[0]: 0})) - l[0]) // 2
            sat[lsat - (t - cost7_t(*l)) // 2] += 1
    log(f"(2) saturation l_sat - w over all (D, mu), D <= {Dmax}, with non-zero plateau: {dict(sorted(sat.items()))}")

    # (3) start distributions
    def row(t, l):
        if not tot.get((t, l)):
            return None
        r = st.get((t, l), {})
        return edge.get((t, l), 0), [r.get(x, 0) for x in range(l[0] + 2, max(r, default=l[0]) + 1, 2)]

    for l, name in (((0, 0, 0), "singlet"), ((2, 2, 0), "(1,1,0)"), ((1, 1, 1), "spinor")):
        log(f"(3) shape {name}: r9 at the edge (l=0) | starts at l = 1, 2, ...  [plateau]")
        for t in range(l[0] % 2 + 2, 2 * Dmax + 1, 2):
            r = row(t, l)
            if r:
                log(f"      D={t / 2:g}: {r[0]} | {r[1]}  [{r[0] + sum(r[1])}]")
    log("    singlet top end: the last starts, at l = l_sat, l_sat - 1, ... as D grows")
    for t in range(4, 2 * Dmax + 1, 2):
        r = row(t, (0, 0, 0))
        if r and len(r[1]) >= 2:
            log(f"      D={t // 2}: {r[1][::-1][:12]}")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(HERE, "data", f"q_spin_D{Dmax}.json"), "w") as f:
        json.dump({"Dmax": Dmax, "LM": LM,
                   "note": "edge [2D, [A,B,C], v]: sum_{e <= mu_1} c(D, e; mu) = r9 at l = 0; "
                           "starts [2D, [A,B,C], [[2e, c], ...]]: c(D, e; mu) for e > mu_1 (mu = (A,B,C)/2)",
                   "edge": [[t, list(l), v] for (t, l), v in sorted(edge.items())],
                   "starts": [[t, list(l), sorted(r.items())] for (t, l), r in sorted(st.items())]}, f)
    with open(os.path.join(RESULTS, "q_spin.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
