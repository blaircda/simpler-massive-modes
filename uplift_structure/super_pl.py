"""
S7 (step A): super-plethystic logarithm of the SO(9) spectrum.

Z_m(q) = sum_n r9_n q^{n-1} (per level) and R(q) = Z_m/(1-q) (cumulative).
Spinorial (half-integer) weights are fermionic, tensorial ones bosonic, so
    PE_super[P] = exp( sum_k (1/k) eps_k psi^k(P) ),   eps_k = 1 (tensor), (-1)^{k+1} (spinor),
and P = PL_super(Z) is solved degree by degree from the exp-log recurrence
    d Z_d = sum_j G_j Z_{d-j},   G_d = sum_{k|d} (d/k) eps_k psi^k(P_{d/k}).
PL_super(R) = PL_super(Z_m) + q.  Each P_d is decomposed into SO(9) irreps
(virtual), and also compared with SO(10)->SO(9) restrictions of small SO(10)
irreps (10, 16/16bar, 45, 54, 120, 126, 144, 210).

Torus characters are W(D4)-invariant arrays on dominant weights (doubled
coordinates); r9_n are recomputed here from the light-cone Fock space
(F_n = (8v+8s) r_n), independently of stored data, and cross-checked
against so10_uplift/data/r9.json.

Usage: super_pl.py [DMAX]      Output: results/super_pl.txt
"""

import json
import os
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RESULTS, SO10
from regge_monotone import adams, divide, exp_log, mult

sys.path.insert(0, SO10)
from derive_r9 import S8, V8
from weights import DomTable, dynkin_B4


def main():
    DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    N = DMAX + 2
    T4 = DomTable("D", 4, 4 * N * N)
    TB = DomTable("B", 4, 4 * N * N)
    b_to_d = np.array([T4.index[tuple(int(a) for a in w)] for w in TB.W] + [T4.ZERO])
    flip = T4.lookup(T4.W * np.array([1, 1, 1, -1]))
    spin = T4.W[:, 0] % 2 == 1

    def so9(f):
        assert np.array_equal(f[:T4.size], f[flip]), "not W(B4)-invariant"
        b = TB.decompose(f[b_to_d])
        return {tuple(int(a) for a in TB.W[i]): b[i] for i in np.nonzero(b[:TB.size])[0]}

    # light-cone Fock F = PE[sum q^m 8v] PE_F[sum q^m 8s]; r_n = F_n/(8v+8s)
    psiV = {k: adams(T4, V8, k) for k in range(1, N + 1)}
    psiS = {k: adams(T4, S8, k) for k in range(1, N + 1)}
    G = [None]
    for j in range(1, N + 1):
        g = T4.zeros()
        for k in range(1, j + 1):
            if j % k == 0:
                g += (j // k) * (psiV[k] + (-1) ** (k + 1) * psiS[k])
        G.append(g)
    F = exp_log(T4, N, G)
    Z = [divide(T4, F[d + 1]) for d in range(DMAX + 1)]       # Z_d = r_{d+1}
    stored = json.load(open(os.path.join(DATA, "r9.json")))["levels"]
    for d in range(DMAX + 1):
        mine = {k: int(v) for k, v in so9(Z[d]).items()}
        ref = {tuple(e["ambient2"]): e["mult"] for e in stored[str(d + 1)]["r9_B4"]}
        assert mine == ref, f"r9_{d + 1} mismatch with stored data"
    log(f"r9_1..r9_{DMAX + 1} recomputed from the Fock space and matched to stored data [{time.time() - t0:.0f}s]")

    def psi(f, k):
        """psi^k on a dominant-weight array (exact, integer)."""
        out = T4.zeros()
        for i in np.nonzero(f[:T4.size])[0]:
            w = tuple(int(k * a) for a in T4.W[i])
            out[T4.index[w]] += f[i]
        return out

    # G_d from Z (exp-log), exact integers
    Gz = [None]
    for d in range(1, DMAX + 1):
        acc = d * Z[d]
        for j in range(1, d):
            acc = acc - mult(T4, Gz[j], Z[d - j]) if Gz[j].any() else acc
        acc[T4.ZERO] = 0
        Gz.append(acc)
    # super-PL inversion: d P_d = G_d - sum_{k|d,k>=2} (d/k) eps_k psi^k(P_{d/k})
    P = [None]
    for d in range(1, DMAX + 1):
        acc = Gz[d].astype(object)
        for k in range(2, d + 1):
            if d % k:
                continue
            Pk = P[d // k]
            ten = np.where(spin, 0, Pk[:T4.size]).astype(object)
            spn = np.where(spin, Pk[:T4.size], 0).astype(object)
            tenv = np.append(ten, 0)
            spnv = np.append(spn, 0)
            sign = (-1) ** (k + 1)
            contrib = psi_obj(T4, tenv, k) + sign * psi_obj(T4, spnv, k)
            acc = acc - (d // k) * contrib
        assert all(int(v) % d == 0 for v in acc), f"PL not integral at degree {d}"
        P.append(np.array([int(v) // d for v in acc], dtype=object))
    log(f"super-PL computed to degree {DMAX} [{time.time() - t0:.0f}s]")

    names = {(2, 0, 0, 0): "9", (0, 0, 0, 0): "1", (1, 1, 1, 1): "16", (2, 2, 0, 0): "36",
             (4, 0, 0, 0): "44", (2, 2, 2, 0): "84", (3, 1, 1, 1): "128", (2, 2, 2, 2): "126",
             (6, 0, 0, 0): "156", (4, 2, 0, 0): "231", (3, 3, 1, 1): "432", (3, 3, 3, 1): "576'"}

    def fmt(dct):
        items = sorted(dct.items(), key=lambda t: (sum(t[0]), t[0]))
        return " ".join(f"{int(v):+d}*{names.get(k, str(dynkin_B4(k)))}" for k, v in items if v)

    for d in range(1, DMAX + 1):
        pz = so9(np.array([int(v) for v in P[d]], dtype=np.int64))
        log(f"  PL(Z_m)_{d}: {fmt(pz)}")
    log("  PL(R) = PL(Z_m) + 1*q  (cumulative factor)")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "super_pl.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


def psi_obj(T4, f, k):
    out = np.zeros(T4.size + 1, dtype=object)
    for i in range(T4.size):
        if f[i]:
            w = tuple(int(k * a) for a in T4.W[i])
            out[T4.index[w]] += f[i]
    return out


if __name__ == "__main__":
    main()
