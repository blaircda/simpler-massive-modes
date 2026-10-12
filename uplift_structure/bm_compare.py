"""
S10: checks of Basile-Markou, "On the deep superstring spectrum", arXiv:2405.18467
(BM), on our data, and what their framework says about Regge monotonicity.

BM describe the full (not supermultiplet-reduced) open-superstring spectrum via
super-Howe duality so(9,1) x osp(*|*).  A Young diagram Y with diagonal length
K is written in modified Frobenius form (l_1..l_K | h_1..h_K) (arms l_i, legs
h_i including the diagonal box).  Their principal embedding (depth 0) is at
    NS (bosons):   N = sum_i ( i l_i + (i - 1/2) h_i ),  multiplicity 1,
    R (fermions):  N = sum_i i (l_i + h_i),               multiplicity 2^(K-1),
with our level n = N - 1/2 (NS; kept by GSO iff integral) and n = N (R).

The full spectrum is F_n = (44+84+128) (x) r9_n; it is computed here from our
r9 data by the Brauer-Klimyk formula.  The script reports

  (1) for every SO(9) irrep first appearing in F_n (n <= 32): does the first
      level and multiplicity agree with BM's principal embedding (GSO-even), or
      appear only later (GSO-odd)?  Also the K = 2, 3, 4 "square" fermions;
  (2) single-box trajectory moves at fixed K (BM's level costs): a box on arm i
      (cost i); a box on leg i (R: cost i; NS: two leg boxes, cost 2i-1 or
      i+j-1): is F monotone, F_(n+c)(Y') >= F_n(Y)?
  (3) column moves (w1 +1, w2 +3, w3 +6, 2w4 +10) on F, for comparison with the
      super-representation result of S9.

Output: results/bm_compare.txt
"""

import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, SO10
from regge_stress import DEFAULT_EXT, load_R9

sys.path.insert(0, SO10)

N, r9, R9 = load_R9(DEFAULT_EXT)
RHO = (7, 5, 3, 1)


def weights_256():
    from sage.all import WeylCharacterRing
    B4 = WeylCharacterRing("B4", style="coroots")
    out = {}
    for hw in ((2, 0, 0, 0), (0, 0, 1, 0), (1, 0, 0, 1)):
        for w, m in B4(*hw).weight_multiplicities().items():
            k = tuple(int(2 * x) for x in w.to_vector())
            out[k] = out.get(k, 0) + int(m)
    assert sum(out.values()) == 256
    return out


def dominate(lam):
    """Signed-permutation (W(B4)) reflection to the dominant chamber: (sign, weight) or None."""
    a = [abs(x) for x in lam]
    if 0 in a or len(set(a)) < 4:
        return None
    order = sorted(range(4), key=lambda i: -a[i])
    perm, sgn = order[:], 1
    for i in range(4):
        while perm[i] != i:
            j = perm[i]
            perm[i], perm[j] = perm[j], perm[i]
            sgn = -sgn
    return sgn * (-1) ** sum(1 for x in lam if x < 0), tuple(a[i] for i in order)


def full_spectrum(W):
    F = {}
    for n in range(1, N + 1):
        out = {}
        for nu, c in r9[n].items():
            for g, m in W.items():
                r = dominate(tuple(nu[i] + g[i] + RHO[i] for i in range(4)))
                if r is not None:
                    mu = tuple(r[1][i] - RHO[i] for i in range(4))
                    out[mu] = out.get(mu, 0) + r[0] * m * c
        F[n] = {k: v for k, v in out.items() if v}
        assert min(F[n].values()) > 0
    return F


def frobenius(mu):
    """(spinorial?, K, arms, legs) of the tensor part of a doubled B4 weight."""
    spin = mu[0] % 2
    rows = [(x - spin) // 2 for x in mu if (x - spin) // 2 > 0]
    cols = [sum(1 for r in rows if r > j) for j in range(rows[0])] if rows else []
    K = sum(1 for i, r in enumerate(rows) if r >= i + 1)
    return spin, K, tuple(rows[i] - (i + 1) for i in range(K)), tuple(cols[i] - i for i in range(K))


def from_frobenius(spin, K, arms, legs):
    cells = set()
    for i in range(K):
        cells.update((i, j) for j in range(i, i + arms[i] + 1))
        cells.update((r, i) for r in range(i, i + legs[i]))
    nrows = max((r for r, _ in cells), default=-1) + 1
    rows = [sum(1 for (r, _) in cells if r == a) for a in range(nrows)]
    if nrows > 4 or any(rows[a] < rows[a + 1] for a in range(nrows - 1)):
        return None
    if {(r, c) for r in range(nrows) for c in range(rows[r])} != cells:
        return None
    return tuple(2 * x + spin for x in rows + [0] * (4 - nrows))


def principal(mu):
    spin, K, arms, legs = frobenius(mu)
    if spin:
        return "R", K, Fraction(sum((i + 1) * (arms[i] + legs[i]) for i in range(K))), 2 ** (K - 1) if K else 1
    N_ns = sum((i + 1) * arms[i] + Fraction(2 * i + 1, 2) * legs[i] for i in range(K))
    return "NS", K, N_ns - Fraction(1, 2), 1


def main():
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    F = full_spectrum(weights_256())
    log(f"full spectrum F_n = (44+84+128) x r9_n computed for n <= {N} [{time.time() - t0:.0f}s]")

    # (1)
    first = {}
    for n in range(1, N + 1):
        for mu, v in F[n].items():
            first.setdefault(mu, (n, v))
    stats, odd = {}, []
    for mu, (n0, v0) in first.items():
        sec, K, npred, mpred = principal(mu)
        even = npred.denominator == 1
        key = (sec, K, "GSO-even" if even else "GSO-odd")
        st = stats.setdefault(key, [0, 0])
        st[0] += 1
        ok = (n0 == npred and v0 == mpred) if even else n0 > npred
        st[1] += ok
        if not ok:
            odd.append((sec, K, mu, n0, v0, str(npred), mpred))
    log(f"(1) {len(first)} irreps first appear in F_n, n <= {N}; agreement with BM's principal embedding "
        "(GSO-even: first level and multiplicity; GSO-odd: appears only later)")
    for key in sorted(stats):
        log(f"  {key}: {stats[key][1]}/{stats[key][0]}")
    log(f"  disagreements (sector, K, weight, first n, mult, BM n, BM mult): {odd}")
    for K in (2, 3, 4):
        mu = tuple(2 * K + 1 if i < K else 1 for i in range(4))
        seq = [(n, F[n][mu]) for n in range(1, N + 1) if F[n].get(mu)][:3]
        log(f"  square fermion K={K} {mu}: BM first level {principal(mu)[2]}, mult {principal(mu)[3]}; data {seq}")

    # (2)
    log("(2) single-box trajectory moves at fixed K, BM costs: F_(n+c)(Y') >= F_n(Y)?")
    res = {}
    for n in range(1, N + 1):
        for mu, v in F[n].items():
            spin, K, arms, legs = frobenius(mu)
            moves = []
            for i in range(K):
                a = list(arms)
                a[i] += 1
                moves.append((f"arm {i + 1}", (spin, K, tuple(a), legs), i + 1))
                if spin:
                    g = list(legs)
                    g[i] += 1
                    moves.append((f"leg {i + 1}", (spin, K, arms, tuple(g)), i + 1))
                else:
                    g = list(legs)
                    g[i] += 2
                    moves.append((f"leg {i + 1} twice", (spin, K, arms, tuple(g)), 2 * i + 1))
                    for j in range(i + 1, K):
                        g = list(legs)
                        g[i] += 1
                        g[j] += 1
                        moves.append((f"legs {i + 1}+{j + 1}", (spin, K, arms, tuple(g)), i + j + 1))
            for name, tgt, c in moves:
                mu2 = from_frobenius(*tgt)
                if mu2 is None or n + c > N:
                    continue
                r = res.setdefault(("R" if spin else "NS", name), [0, 0])
                r[0] += 1
                r[1] += F[n + c].get(mu2, 0) < v
    for key in sorted(res):
        log(f"  {key[0]:2s} {key[1]:12s}: {res[key][1]:5d} violations of {res[key][0]}")

    # (3)
    log("(3) column moves on F (cf. S9 for r9)")
    for name, w, d in (("w1", (2, 0, 0, 0), 1), ("w2", (2, 2, 0, 0), 3), ("w3", (2, 2, 2, 0), 6),
                       ("2w4", (2, 2, 2, 2), 10), ("w2", (2, 2, 0, 0), 2)):
        viol = tot = 0
        for n in range(d + 1, N + 1):
            for nu, v in F[n - d].items():
                tot += 1
                viol += F[n].get(tuple(x + y for x, y in zip(nu, w)), 0) < v
        log(f"  {name:3s} at +{d:2d}: {viol} violations of {tot}")
    log(f"[{time.time() - t0:.0f}s]")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "bm_compare.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
