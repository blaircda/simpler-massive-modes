"""
S16 (model): the top end of the ramp, as the residue of Q at the alpha^+_{-2} pole.

Near saturation the last trajectory starts are carried by alpha^+_{-2} (weight e_1, level 2;
the shift of S8's T_2), which raises depth and spin together.  Put y = v/q in
Q(q, y) (q_spin.py) and grade by g = D - e: the factor 1/(1 - q^{k-1} y) of alpha^+_{-k}
at k = 2 becomes 1/(1 - v), so at fixed g the coefficient of v^e is eventually constant,
equal to the residue

    Qhat(q) = q^3 (q ground) det_7(1 - q^2 g) Phat / (1 + R)        (everything at y = 1/q),

with Phat = P without alpha^+_{-2}, each factor graded by g = depth - charge.  Prediction:
c(D, e; mu) = mult of mu in [q^(D-e)] Qhat once e is close enough to its top.

This script computes Qhat from the light-cone oscillators (no spectrum data) and compares
it with the starts of Q in data/q_spin_D<Dmax>.json (q_spin.py output, also first-principles):
  (1) Qhat multiplicities for small shapes;
  (2) for every shape and depth, the spin from which the starts of Q equal Qhat (threshold in 2e - D);
  (3) the dimension series of Qhat, and a test of the obvious positive Fock form with the same
      dimensions, q^2 (1+7+8) det_7(1+qg) prod_{j>=2} det_8(1+q^j g)/det_7(1-q^j g) prod_{j>=3} 1/(1-q^j);
  (4) whether Qhat is non-negative.

Usage: q_top_end.py [--q path] [--gmax 22]   Output: results/q_top_end.txt
"""

import argparse
import json
import os
import sys
from collections import Counter
from itertools import product

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS
from q_spin import w256
from regge_limit import weyl_b3

HERE = os.path.dirname(os.path.abspath(__file__))
RHO2 = (5, 3, 1)
ZERO = (0, 0, 0)
SIX = [tuple(2 * s * (j == i) for j in range(3)) for i in range(3) for s in (1, -1)]
FOUR = [w for w in product((1, -1), repeat=3) if w.count(-1) % 2 == 0]       # 8s at charge +1/2
FOURBAR = [w for w in product((1, -1), repeat=3) if w.count(-1) % 2 == 1]    # 8s at charge -1/2
SPIN8 = list(product((1, -1), repeat=3))


class GSeries:
    """Slices Z[i], i = g - g0, on the doubled SO(7) torus (both parities), box |2 mu| <= Bd."""

    def __init__(self, g0, gmax, Bd):
        self.g0, self.n, self.Bd = g0, gmax - g0 + 1, Bd
        S = 2 * Bd + 1
        self.Z = np.zeros((self.n, S, S, S), dtype=np.int64)

    def sh(self, A, w):
        S = 2 * self.Bd + 1
        out = np.zeros_like(A)
        sl = [(slice(max(d, 0), S + min(d, 0)), slice(max(-d, 0), S + min(-d, 0))) for d in w]
        out[sl[0][0], sl[1][0], sl[2][0]] = A[sl[0][1], sl[1][1], sl[2][1]]
        return out

    def geom(self, d, w):                                     # 1/(1 - q^d x^w), d >= 1
        for i in range(d, self.n):
            self.Z[i] += self.sh(self.Z[i - d], w)

    def binom(self, d, w, c):                                 # (1 + c q^d x^w), d >= 0
        if d == 0:
            for i in range(self.n):
                self.Z[i] += c * self.sh(self.Z[i], w)
            return
        for i in range(self.n - 1, d - 1, -1):
            self.Z[i] += c * self.sh(self.Z[i - d], w)

    def mult(self, g, lam):
        i, S, tot = g - self.g0, 2 * self.Bd + 1, 0
        if not 0 <= i < self.n:
            return 0
        for wr, sg in weyl_b3():
            idx = [lam[k] + RHO2[k] - wr[k] + self.Bd for k in range(3)]
            if all(0 <= x < S for x in idx):
                tot += sg * int(self.Z[i][idx[0], idx[1], idx[2]])
        return tot


def qhat(gmax, Bd):
    """Qhat for g <= gmax (the slice g = gmax is complete because ground's g = -1 term is applied to a longer series)."""
    G = gmax - 3 + 1                                          # W = Qhat / q^3 to g = G (one extra slice)
    S = GSeries(-1, G + 1, Bd)
    S.Z[1, Bd, Bd, Bd] = 1                                    # Phat starts at 1 (g = 0)
    for k in range(3, G + 4):
        S.geom(k - 2, ZERO)                                   # alpha^+_{-k}, k >= 3
    for k in range(1, G + 2):
        S.geom(k + 2, ZERO)                                   # alpha^-_{-k}
    for k in range(1, G + 2):
        for w in SIX:
            S.geom(k, w)                                      # transverse alpha^i_{-k}
    for k in range(1, G + 3):
        for w in FOUR:
            S.binom(k - 1, w, 1)                              # charge +1/2 fermions: g = k - 1
    for k in range(1, G + 2):
        for w in FOURBAR:
            S.binom(k + 1, w, 1)                              # charge -1/2 fermions: g = k + 1
    for w in [ZERO] + SIX:
        S.binom(2, w, -1)                                     # det_7(1 - q^2 g)
    ground = [(-1, ZERO), (3, ZERO)] + [(1, w) for w in SIX] + [(0, w) for w in FOURBAR] + [(2, w) for w in FOUR]
    new = np.zeros_like(S.Z)
    for i in range(S.n):
        for d, w in ground:
            if 0 <= i - d < S.n:
                new[i] += S.sh(S.Z[i - d], w)
    S.Z = new
    R = Counter()
    for W in w256():
        if W != (4, 0, 0, 0):
            R[(4 - W[0], W[1:])] += 1                         # (q/y)^{2 - x1} x^mu at y = 1/q: g = 4 - 2 x1
    for i in range(S.n):
        for (d, w), c in sorted(R.items()):
            if i - d >= 0:
                S.Z[i] -= c * S.sh(S.Z[i - d], w)
    assert not any(np.take(S.Z, [0, -1], axis=ax).any() for ax in (1, 2, 3)), "box too small"
    assert np.abs(S.Z).max() < 2 ** 60
    S.g0 += 3                                                 # Qhat = q^3 W
    S.n -= 1                                                  # drop the incomplete top slice
    S.Z = S.Z[:-1]
    return S


def fock_guess(gmax, Bd):
    S = GSeries(2, gmax, Bd)
    for w in [ZERO, ZERO] + SIX + SPIN8:
        S.Z[0] += S.sh(np.pad(np.ones((1, 1, 1), dtype=np.int64), Bd), w)        # q^2 (1 + 7 + 8)
    for w in [ZERO] + SIX:
        S.binom(1, w, 1)                                      # det_7(1 + q g)
    for j in range(2, gmax + 1):
        for w in SPIN8:
            S.binom(j, w, 1)
        for w in [ZERO] + SIX:
            S.geom(j, w)
    for j in range(3, gmax + 1):
        S.geom(j, ZERO)
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", default=os.path.join(HERE, "data", "q_spin_D22.json"))
    ap.add_argument("--gmax", type=int, default=22)
    args = ap.parse_args()
    Qd = json.load(open(args.q))
    st = {(t, tuple(l)): {s: v for s, v in r} for t, l, r in Qd["starts"]}
    Dmax, gmax = Qd["Dmax"], args.gmax
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    Bd = 2 * gmax + 8
    H = qhat(gmax, Bd)
    log(f"Qhat (residue of Q at the alpha^+_(-2) pole) for g = D - e <= {gmax}; compared with {os.path.basename(args.q)}")
    log("(1) Qhat multiplicities from the first non-zero g:")
    for lam, name in (((0, 0, 0), "singlet"), ((2, 0, 0), "vector"), ((2, 2, 0), "(1,1,0)"), ((1, 1, 1), "spinor")):
        seq = [(g, H.mult(g, lam)) for g in range(H.g0, H.g0 + H.n)]
        seq = [(g, m) for g, m in seq if m] if any(m for _, m in seq) else []
        log(f"    {name}: from g = {seq[0][0] if seq else '-'}: {[m for _, m in seq]}")
    # (2) where the starts of Q equal Qhat: the threshold in 2e - D
    thr, diff = Counter(), []
    for (t, l), row in st.items():
        if not row:
            continue
        s_ag = max(row) + 2
        for s in sorted(row, reverse=True):                   # last starts first; g = D - e = (t - s)/2
            if row[s] != H.mult((t - s) // 2, l):
                break
            s_ag = s
        thr[(2 * s_ag - t) // 2] += 1
        diff += [((2 * s - t) // 2, t / 2, tuple(x / 2 for x in l)) for s in row if row[s] != H.mult((t - s) // 2, l)]
    log(f"(2) {len(st)} rows of starts (all shapes, D <= {Dmax}).  First 2e - D from which every later start equals Qhat: "
        f"{dict(sorted(thr.items()))}")
    log(f"    largest 2e - D at which a start differs from Qhat: {max(diff)[0] if diff else '-'} "
        f"(e.g. D, mu = {[(d, m) for _, d, m in sorted(diff)[-3:]]}).  So c(D, e; mu) = mult_mu [q^(D-e)] Qhat "
        f"whenever 2e - D >= {max(diff)[0] + 1 if diff else '-'}.")
    # (3) dimensions and the Fock guess
    dims = [int(H.Z[i].sum()) for i in range(H.n)]
    log(f"(3) dimensions of Qhat, g = {H.g0}..: {dims[:12]}")
    F = fock_guess(gmax, Bd)
    lams = [(A, B, C) for par in (0, 1) for A in range(par, 19, 2) for B in range(par, A + 1, 2) for C in range(par, B + 1, 2)]
    first = None
    for g in range(H.g0, min(H.g0 + H.n, F.g0 + F.n)):
        if any(H.mult(g, l) != F.mult(g, l) for l in lams):
            first = g
            break
    log(f"    Fock guess with the same dimensions: {'agrees for all g' if first is None else f'first differs at g = {first}'}"
        f" (dimensions agree: {all(int(F.Z[i].sum()) == dims[i] for i in range(min(F.n, H.n)))})")
    # (4) non-negativity of Qhat
    lams = [(A, B, C) for par in (0, 1) for A in range(par, 2 * (gmax // 2 + 3) + 1, 2)
            for B in range(par, A + 1, 2) for C in range(par, B + 1, 2)]
    negs = [(g, l) for g in range(H.g0, H.g0 + H.n) for l in lams if H.mult(g, l) < 0]
    log(f"(4) Qhat non-negative for g <= {H.g0 + H.n - 1}, all shapes: {not negs}" + (f"; negative at {negs[:5]}" if negs else ""))
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "q_top_end.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
