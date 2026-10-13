"""
S17 (model): the alpha^+_{-3} residue, and the hierarchy of Regge directions.

Write T_k for the shift (spin +1, level +k): in top-distance coordinates (D, N) -> (D + k - 1, N + 1),
i.e. multiplication by q^{k-1} y, the monomial of the light-cone oscillator alpha^+_{-k}.
For a shape mu (the SO(7) part of the weight) define, at dominant spins N >= mu_1,
    r9(D, N)  = r9_{N+D}(N, mu),
    s(D, N)   = r9(D, N) - r9(D, N-1)          (T_1-differences: the starts of S16),
    O(D, N)   = s(D, N) - s(D-1, N-1)          (T_2-differences of the starts),
    L3(D, N)  = O(D, N) - O(D-2, N-1)          (T_3-differences),
each term counted as zero below the dominance edge N = mu_1.  The stable totals along each
direction are residues of the light-cone generating function at the alpha^+_{-k} poles:
    sum_N s(D, N)          = Z_Regge_D        (alpha^+_{-1}; S15),
    sum along T_2 of O     = Qhat_{D-N}       (alpha^+_{-2}; S16, exact for 2N - D >= -1),
    sum along T_3 of L3    = Ohat_{D-2N}      (alpha^+_{-3}; here),
with Ohat the residue at w = q^2 y = 1 of O = (1 - q y) Q, graded by h = D - 2e:
    Ohat(q) = q^5 (q ground) det_7(1 - q^3 g) P' / (1 + R)     at y = q^-2,
P' the oscillators without alpha^+_{-1,-2,-3}.  This script computes Ohat from the light
cone (no spectrum data) and compares with the starts of data/q_spin_D<Dmax>.json (q_spin.py
output, also first-principles):
  (1) Ohat multiplicities for small shapes; its non-negativity; its dimension series and the
      obvious positive Fock form with the same dimensions;
  (2) exactness: where O(D, N) = mult_mu Ohat_{D-2N} (threshold in 3N - D), and the two-pole
      formula s(D, N) = Qhat_{D-N} - sum_{j>=1} Ohat_{D-2N-j};
  (3) higher Regge monotonicity: signs of the T_2- and T_3-differences O and L3 (D <= Dmax);
  (4) every residue R_k (k = 1..5, grade g = D - (k-1)N) against its predicted positive Fock form,
      and the k = 2 form against Qhat at full depth.

Usage: q_multi_regge.py [--q path] [--hmax 22]   Output: results/q_multi_regge.txt
"""

import argparse
import json
import os
import sys
from collections import Counter
from itertools import product

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import q_top_end
from common import RESULTS
from q_spin import w256
from regge_limit import weyl_b3

HERE = os.path.dirname(os.path.abspath(__file__))
RHO2 = (5, 3, 1)
ZERO, SIX, FOUR, FOURBAR, SPIN8 = q_top_end.ZERO, q_top_end.SIX, q_top_end.FOUR, q_top_end.FOURBAR, q_top_end.SPIN8


class TSeries:
    """Slices Z[t], t = 2h, on SO(7) p-coordinates: mu = p + (t % 2)/2, p in [-B, B]."""

    def __init__(self, tlo, thi, B):
        self.tlo, self.thi, self.B = tlo, thi, B
        S = 2 * B + 1
        self.Z = {t: np.zeros((S, S, S), dtype=np.int64) for t in range(tlo, thi + 1)}

    def shifted(self, A, ts, td, W2):
        S = 2 * self.B + 1
        d = [(W2[k] + ts % 2 - td % 2) // 2 for k in range(3)]
        out = np.zeros_like(A)
        sl = [(slice(max(x, 0), S + min(x, 0)), slice(max(-x, 0), S + min(-x, 0))) for x in d]
        out[sl[0][0], sl[1][0], sl[2][0]] = A[sl[0][1], sl[1][1], sl[2][1]]
        return out

    def geom(self, dt, W2):                                   # 1/(1 - q^{dt/2} x^{W2/2}), dt > 0
        for t in range(self.tlo + dt, self.thi + 1):
            self.Z[t] += self.shifted(self.Z[t - dt], t - dt, t, W2)

    def binom(self, dt, W2, c):                               # (1 + c q^{dt/2} x^{W2/2}), any dt
        if dt == 0:
            for t in self.Z:
                self.Z[t] += c * self.shifted(self.Z[t], t, t, W2)
            return
        rng = range(self.thi, self.tlo + dt - 1, -1) if dt > 0 else range(self.tlo, self.thi + dt + 1)
        for t in rng:
            self.Z[t] += c * self.shifted(self.Z[t - dt], t - dt, t, W2)

    def msum(self, terms):                                    # multiply by sum c q^{dt/2} x^{W2/2}
        new = {t: np.zeros_like(A) for t, A in self.Z.items()}
        for t in self.Z:
            for dt, W2, c in terms:
                if self.tlo <= t - dt <= self.thi:
                    new[t] += c * self.shifted(self.Z[t - dt], t - dt, t, W2)
        self.Z = new

    def divide(self, terms):                                  # divide by 1 + sum c q^{dt/2} x^{W2/2}, dt > 0
        for t in range(self.tlo, self.thi + 1):
            for dt, W2, c in terms:
                if t - dt >= self.tlo:
                    self.Z[t] -= c * self.shifted(self.Z[t - dt], t - dt, t, W2)

    def mult(self, t, lam):
        if t not in self.Z or t % 2 != lam[0] % 2:
            return 0
        S, A, par, tot = 2 * self.B + 1, self.Z[t], t % 2, 0
        for wr, sg in weyl_b3():
            idx = [(lam[k] + RHO2[k] - wr[k] - par) // 2 + self.B for k in range(3)]
            if all(0 <= x < S for x in idx):
                tot += sg * int(A[idx[0], idx[1], idx[2]])
        return tot


def ohat(hmax):
    """Ohat for 2h <= 2*hmax: W = Ohat / q^5 on t in [-8, 2*hmax - 10 + 10] (margin for the
    negative shifts), then t -> t + 10."""
    tlo, thi, B = -8, 2 * hmax, hmax + 8
    S = TSeries(tlo, thi, B)
    S.Z[0][B, B, B] = 1
    K = thi // 2 + 4
    for k in range(4, K):
        S.geom(2 * (k - 3), ZERO)                             # alpha^+_{-k}, k >= 4: h = k - 3
    for k in range(1, K):
        S.geom(2 * (k + 3), ZERO)                             # alpha^-_{-k}: h = k + 3
    for k in range(1, K):
        for w in SIX:
            S.geom(2 * k, w)                                  # transverse: h = k
    for k in range(1, K):
        for w in FOUR:
            S.binom(2 * k - 3, w, 1)                          # charge +1/2: h = k - 3/2
        for w in FOURBAR:
            S.binom(2 * k + 3, w, 1)                          # charge -1/2: h = k + 3/2
    for w in [ZERO] + SIX:
        S.binom(6, w, -1)                                     # det_7(1 - q^3 g)
    S.msum([(-4, ZERO, 1), (8, ZERO, 1)] + [(2, w, 1) for w in SIX]
           + [(-1, w, 1) for w in FOURBAR] + [(5, w, 1) for w in FOUR])     # q ground at y = q^-2
    R = Counter()
    for W in w256():
        if W != (4, 0, 0, 0):
            R[(12 - 3 * W[0], W[1:])] += 1                    # (q/y)^{2-x1}: h = 3(2 - x1)
    S.divide([(dt, w, c) for (dt, w), c in sorted(R.items())])
    for t in list(S.Z):
        for ax in range(3):
            assert not np.take(S.Z[t], [0, -1], axis=ax).any(), "box too small"
    out = TSeries(tlo + 10, thi + 10 - 9, B)                  # top 9 slices of W are incomplete
    out.Z = {t + 10: A for t, A in S.Z.items() if t + 10 <= out.thi}
    return out


class DSeries(TSeries):
    """Like TSeries but on doubled weight coordinates, both parities in every slice (|2 mu| <= B)."""

    def shifted(self, A, ts, td, W2):
        S = 2 * self.B + 1
        out = np.zeros_like(A)
        sl = [(slice(max(x, 0), S + min(x, 0)), slice(max(-x, 0), S + min(-x, 0))) for x in W2]
        out[sl[0][0], sl[1][0], sl[2][0]] = A[sl[0][1], sl[1][1], sl[2][1]]
        return out

    def mult(self, t, lam):
        if t not in self.Z:
            return 0
        S, A, tot = 2 * self.B + 1, self.Z[t], 0
        for wr, sg in weyl_b3():
            idx = [lam[k] + RHO2[k] - wr[k] + self.B for k in range(3)]
            if all(0 <= x < S for x in idx):
                tot += sg * int(A[idx[0], idx[1], idx[2]])
        return tot


def residue(k, gmax):
    """R_k: residue of Q prod_{j=2}^{k-1} (1 - q^{j-1} y) at the alpha^+_{-k} pole y = q^{-(k-1)},
    graded by g = D - (k-1) e.  Returns (DSeries W, shift): R_k at t = 2g is W at t - shift.
    Factors with negative grade are flipped, (1 + q^d x^w) = q^d x^w (1 + q^-d x^-w) (a 4-plet's
    weights sum to zero), and the ground states are multiplied by q^{2a}, so every shift is >= 0."""
    a = k - 1
    shift = 2 + 4 * a - 2 * a                                 # q y^-2, and q^{2a} taken out of the ground
    for m in range(1, k + 1):
        d = 2 * m - 1 - a                                     # the 4 at level m: t = 2m - 1 - a
        if d < 0:
            shift += 4 * d
    thi = 2 * gmax - shift
    S = DSeries(0, thi, thi + 8)
    S.Z[0][S.B, S.B, S.B] = 1
    M = thi // 2 + k + 3
    for m in range(k + 1, M):
        S.geom(2 * (m - k), ZERO)                             # alpha^+_{-m}, m > k
    for m in range(1, M):
        S.geom(2 * (m + k), ZERO)                             # alpha^-_{-m}
        for w in SIX:
            S.geom(2 * m, w)                                  # transverse
    for m in range(1, M):
        d = 2 * m - 1 - a
        for w in FOUR:
            if d < 0:
                S.binom(-d, tuple(-x for x in w), 1)          # flipped: a 4bar at grade -d
            else:
                S.binom(d, w, 1)
        for w in FOURBAR:
            S.binom(2 * m + 1 + a, w, 1)
    for w in [ZERO] + SIX:
        S.binom(2 * k, w, -1)                                 # det_7(1 - (q/y) g)
    S.msum([(0, ZERO, 1), (4 + 4 * a, ZERO, 1)] + [(2 + 2 * a, w, 1) for w in SIX]
           + [(1 + a, w, 1) for w in FOURBAR] + [(3 + 3 * a, w, 1) for w in FOUR])
    R = Counter()
    for W in w256():
        if W != (4, 0, 0, 0):
            R[(k * (4 - W[0]), W[1:])] += 1                   # (q/y)^{2 - x1}
    S.divide([(dt, w, c) for (dt, w), c in sorted(R.items())])
    for t in S.Z:
        for ax in range(3):
            assert not np.take(S.Z[t], [0, -1], axis=ax).any(), "box too small"
    return S, shift


def fock_form(k, gmax):
    """Predicted R_k = q^{p_k} [prod_{m < k/2} det_8(1 + q^{k/2-m} g)] [(1+7+8) if k even]
    prod_{j=1}^{k-1} det_7(1 - q^j g)^-1  prod_{j >= k+1} ((1-q^j) det_7(1-q^j g))^-1  prod_{j >= k/2+1} det_8(1 + q^j g),
    p_k = k + sum_{m < k/2} 4 (m - k/2).  Returns (DSeries, shift) as residue() does."""
    shift = 2 * k + sum(4 * (2 * m - k) for m in range(1, k) if 2 * m < k)       # 2 p_k
    thi = 2 * gmax - shift
    S = DSeries(0, thi, thi + 8)
    S.Z[0][S.B, S.B, S.B] = 1
    for m in range(1, k):
        if 2 * m < k:
            for w in SPIN8:
                S.binom(k - 2 * m, w, 1)                      # spinor at grade k/2 - m
    if k % 2 == 0:
        for w in FOUR:
            S.binom(0, w, 1)                                  # prod_{w in 4} (1 + x^w) = 1 + 7 + 8
    for j in range(1, k):
        for w in [ZERO] + SIX:
            S.geom(2 * j, w)
    for j in range(k + 1, thi // 2 + 2):
        for w in [ZERO, ZERO] + SIX:
            S.geom(2 * j, w)
    for t2 in range(k + 2, thi + 1, 2):
        for w in SPIN8:
            S.binom(t2, w, 1)                                 # spinor at grade t2/2 >= k/2 + 1
    return S, shift


def fock_guess(hmax, B):
    """q det_8(1 + q^{1/2} g) / (det_7(1 - q g) det_7(1 - q^2 g)) prod_{j>=4} 1/((1-q^j) det_7(1-q^j g))
    prod_{m >= 5/2} det_8(1 + q^m g): positive, with the dimension series of Ohat."""
    S = TSeries(2, 2 * hmax, B)
    S.Z[2][B, B, B] = 1
    for w in SPIN8:
        S.binom(1, w, 1)
    for j in (1, 2):
        for w in [ZERO] + SIX:
            S.geom(2 * j, w)
    for j in range(4, hmax + 1):
        for w in [ZERO, ZERO] + SIX:
            S.geom(2 * j, w)
    for m2 in range(5, 2 * hmax + 1, 2):
        for w in SPIN8:
            S.binom(m2, w, 1)
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", default=os.path.join(HERE, "data", "q_spin_D22.json"))
    ap.add_argument("--hmax", type=int, default=22)
    args = ap.parse_args()
    Qd = json.load(open(args.q))
    Dmax, hmax = Qd["Dmax"], args.hmax
    edge = {(t, tuple(l)): v for t, l, v in Qd["edge"]}
    st = {(t, tuple(l)): {s: v for s, v in r} for t, l, r in Qd["starts"]}
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    H2 = q_top_end.qhat(Dmax, 2 * Dmax + 8)                   # Qhat, graded by g = D - N
    H3 = ohat(hmax)                                           # Ohat, graded by t = 2h, h = D - 2N
    log(f"Ohat (residue of (1 - qy) Q at the alpha^+_(-3) pole) for h = D - 2N <= {hmax}; compared with {os.path.basename(args.q)}")

    # (1) Ohat itself
    log("(1) Ohat multiplicities from the first non-zero h:")
    for lam, name in (((0, 0, 0), "singlet"), ((2, 0, 0), "vector"), ((2, 2, 0), "(1,1,0)"), ((1, 1, 1), "spinor")):
        seq = [(t / 2, H3.mult(t, lam)) for t in range(H3.tlo, H3.thi + 1) if t % 2 == lam[0] % 2]
        nz = [i for i, (_, m) in enumerate(seq) if m]
        seq = seq[nz[0]:] if nz else []
        log(f"    {name}: from h = {seq[0][0] if seq else '-'}: {[m for _, m in seq]}")
    lams = [(A, B, C) for par in (0, 1) for A in range(par, 2 * (hmax // 2 + 4) + 1, 2)
            for B in range(par, A + 1, 2) for C in range(par, B + 1, 2)]
    negs = [(t / 2, l) for t in range(H3.tlo, H3.thi + 1) for l in lams if H3.mult(t, l) < 0]
    log(f"    Ohat non-negative for h <= {H3.thi / 2:g}, all shapes: {not negs}" + (f"; negative at {negs[:5]}" if negs else ""))
    dims = {t / 2: int(A.sum()) for t, A in sorted(H3.Z.items()) if A.any()}
    log(f"    dimensions by h: {dict(list(dims.items())[:12])}")
    F = fock_guess(hmax, H3.B)
    diff = next((t / 2 for t in range(2, H3.thi + 1) if t in F.Z and any(H3.mult(t, l) != F.mult(t, l) for l in lams)), None)
    dims_ok = all(int(F.Z[t].sum()) == int(H3.Z[t].sum()) for t in range(2, H3.thi + 1) if t in F.Z)
    log(f"    positive Fock guess q det_8(1+q^(1/2)g)/(det_7(1-qg)det_7(1-q^2g)) prod_(j>=4) [(1-q^j)det_7(1-q^jg)]^-1 "
        f"prod_(m>=5/2) det_8(1+q^m g): dimensions agree {dims_ok}; "
        + ("agrees for all h" if diff is None else f"SO(7) content first differs at h = {diff:g}"))

    # helpers on Q's data (doubled units: t = 2D, X = 2N, l = 2 mu)
    shapes = {l for (_, l) in set(edge) | set(st)}

    def s(t, X, l):
        if t < 2 or X < l[0]:
            return 0
        if X == l[0]:
            return edge.get((t, l), 0)
        return st.get((t, l), {}).get(X, 0)

    def O(t, X, l):
        return s(t, X, l) - s(t - 2, X - 2, l)

    def L3(t, X, l):
        return O(t, X, l) - O(t - 4, X - 2, l)

    pts = []
    for l in shapes:
        for t in range(2 + l[0] % 2, 2 * Dmax + 1, 2):
            last = max([l[0]] + list(st.get((t, l), {})))
            for X in range(l[0], last + 3, 2):
                pts.append((t, X, l))

    # (2) exactness of Ohat and the two-pole formula
    thr3, bad3, thr2, bad2 = Counter(), [], Counter(), []
    for t, X, l in pts:
        o, h2 = O(t, X, l), t - 2 * X                                 # 2h = 2(D - 2N)
        if o != H3.mult(h2, l):
            bad3.append(((3 * X - t) // 2, t / 2, X / 2, tuple(x / 2 for x in l)))
        if X > l[0]:
            pred = H2.mult((t - X) // 2, l) - sum(H3.mult(h2 - 2 * j, l) for j in range(1, (h2 - H3.tlo) // 2 + 1))
            if s(t, X, l) != pred:
                bad2.append(((3 * X - t) // 2, t / 2, X / 2, tuple(x / 2 for x in l)))
    bad3i = [b for b in bad3 if b[2] > b[3][0]]                       # interior points N > mu_1
    m3i = max(b[0] for b in bad3i) if bad3i else None
    m3 = max(b[0] for b in bad3) if bad3 else None
    m2 = max(b[0] for b in bad2) if bad2 else None
    log(f"(2) {len(pts)} points (D <= {Dmax}, all shapes, spins from the edge to two past the last start)")
    log(f"    O(D, N) = Ohat_(D-2N): fails at {len(bad3)} points, the largest 3N - D among them {m3} "
        f"(e.g. D, N, mu = {[b[1:] for b in sorted(bad3)[-3:]]}); so it holds whenever 3N - D >= {m3 + 1 if m3 is not None else '-'}")
    log(f"      at interior points (N > mu_1) only: fails at {len(bad3i)}, largest 3N - D {m3i} "
        f"(e.g. {[b[1:] for b in sorted(bad3i)[-3:]]}); holds whenever 3N - D >= {m3i + 1 if m3i is not None else '-'}")
    log(f"    two-pole formula s(D, N) = Qhat_(D-N) - sum_(j>=1) Ohat_(D-2N-j), N > mu_1: fails at {len(bad2)} points, "
        f"the largest 3N - D among them {m2} (e.g. {[b[1:] for b in sorted(bad2)[-3:]]}); "
        f"so it holds whenever 3N - D >= {m2 + 1 if m2 is not None else '-'}")

    # (3) higher Regge monotonicity
    negO = [(t / 2, X / 2, tuple(x / 2 for x in l), O(t, X, l)) for t, X, l in pts if X > l[0] and O(t, X, l) < 0]
    negL = [(t / 2, X / 2, tuple(x / 2 for x in l), L3(t, X, l)) for t, X, l in pts if X > l[0] and L3(t, X, l) < 0]
    nO = sum(1 for t, X, l in pts if X > l[0] and O(t, X, l))
    nL = sum(1 for t, X, l in pts if X > l[0] and L3(t, X, l))
    log(f"(3) T_2-differences of the starts O(D, N) = s(D, N) - s(D-1, N-1), N > mu_1, D <= {Dmax}: "
        f"{nO} non-zero, {len(negO)} negative" + (f", e.g. {sorted(negO, key=lambda x: x[3])[:4]}" if negO else ""))
    log(f"    T_3-differences L3(D, N) = O(D, N) - O(D-2, N-1): {nL} non-zero, {len(negL)} negative"
        + (f", e.g. {sorted(negL, key=lambda x: x[3])[:4]}" if negL else ""))
    # (4) every residue R_k is a positive Fock space
    log("(4) residue R_k at the alpha^+_(-k) pole (grade g = D - (k-1)N) against its predicted Fock form, g <= 12:")
    lamsD = [(A, B, C) for par in (0, 1) for A in range(par, 29, 2) for B in range(par, A + 1, 2) for C in range(par, B + 1, 2)]
    for k in range(1, 6):
        Rk, sh = residue(k, 12)
        Fk, shf = fock_form(k, 12)
        assert sh == shf
        ts = [t for t in Rk.Z if Rk.Z[t].any()]
        first = None
        for t in sorted(Rk.Z):
            if any(Rk.mult(t, l) != Fk.mult(t, l) for l in lamsD):
                first = (t + sh) / 2
                break
        neg = any(Rk.mult(t, l) < 0 for t in Rk.Z for l in lamsD)
        log(f"    k={k}: lowest grade g = {(min(ts) + sh) / 2 if ts else '-'}; "
            + ("equals the Fock form for all g <= 12" if first is None else f"differs from the Fock form at g = {first}")
            + f"; non-negative: {not neg}")
    G = q_top_end.GSeries(2, Dmax, H2.Bd)                     # k = 2 at full depth, against Qhat
    G.Z[0][H2.Bd, H2.Bd, H2.Bd] = 1
    for w in FOUR:
        G.binom(0, w, 1)
    for w in [ZERO] + SIX:
        G.geom(1, w)
    for j in range(3, Dmax + 1):
        for w in [ZERO, ZERO] + SIX:
            G.geom(j, w)
    for j in range(2, Dmax + 1):
        for w in SPIN8:
            G.binom(j, w, 1)
    lamsG = [(A, B, C) for par in (0, 1) for A in range(par, 2 * (Dmax // 2 + 3) + 1, 2)
             for B in range(par, A + 1, 2) for C in range(par, B + 1, 2)]
    diffG = next((g for g in range(2, min(H2.g0 + H2.n, G.g0 + G.n)) if any(H2.mult(g, l) != G.mult(g, l) for l in lamsG)), None)
    log(f"    k=2 at full depth: Qhat = q^2 (1+7+8) det_7(1-qg)^-1 prod_(j>=3) [(1-q^j)det_7(1-q^jg)]^-1 prod_(j>=2) det_8(1+q^jg) "
        f"for all g <= {min(H2.g0 + H2.n, G.g0 + G.n) - 1}: " + ("yes" if diffG is None else f"NO, first differs at g = {diffG}"))
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "q_multi_regge.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
