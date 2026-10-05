"""
Library for the BBMS-graded lift (see derive_graded.py for the method).

    levels_data(N)            -> r[n][nu] (D4 dominant arrays), R9[n] (B4 dicts), T
    build_system(n, r, R9, T, slack, lams=None)
                              -> (A sparse QQ matrix, b, cols, rows)
    phi(lam, D2)              -> image of V_lam x^Delta as {(nu, w4): mult}
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "so10_uplift"))

from derive_r9 import E1, HEIGHT, S8, V8, fock_dims, orbit_size_D4
from uplift import candidates
from weights import DomTable, dynkin_D5

from sage.all import GF, QQ, WeylCharacterRing, matrix, vector

D5 = WeylCharacterRing("D5", style="coroots")
_WM = {}


def d5_weights(lam):
    """Doubled-weight multiplicity dict of the D5 irrep with doubled h.w. lam."""
    if lam not in _WM:
        V = D5(*dynkin_D5(lam))
        _WM[lam] = {tuple(int(2 * a) for a in w.to_vector()): int(m)
                    for w, m in V.weight_multiplicities().items()}
    return _WM[lam]


def graded_fock(T, N):
    """F[n] = {nu: dominant-weight array} for n = 0..N."""
    F = [{0: T.zeros()}]
    F[0][0][T.index[(0, 0, 0, 0)]] = 1
    for n in range(1, N + 1):
        acc = {}
        for j in range(1, n + 1):
            for k in range(1, j + 1):
                if j % k:
                    continue
                c = j // k
                for nu0, Fm in F[n - j].items():
                    sv = sum(Fm[T.gather(tuple(k * a for a in w))] for w in V8)
                    ss = sum(Fm[T.gather(tuple(k * a for a in w))] for w in S8)
                    term = c * (sv - (-1) ** k * ss)
                    acc[nu0 + k] = acc.get(nu0 + k, T.zeros()) + term
        Fn = {}
        for nu, a in acc.items():
            assert np.all(a % n == 0), f"graded exp-log not divisible at level {n}, nu {nu}"
            a = a // n
            a[T.ZERO] = 0
            if a.any():
                assert a.min() >= 0
                Fn[nu] = a
        F.append(Fn)
    return F


def divide(T, f):
    """Exact division of a W(D4)-invariant array by 8_v + 8_s."""
    D = V8 + S8
    others = [d for d in D if d != E1]
    g_up = T.gather(tuple(-a for a in E1)).tolist()
    g_oth = [T.gather(tuple(d[i] - E1[i] for i in range(4))).tolist() for d in others]
    order = np.argsort(-(T.W @ HEIGHT), kind="stable").tolist()
    fl = f.tolist()
    r = [0] * (T.size + 1)
    for i in order:
        v = fl[g_up[i]]
        for g in g_oth:
            v -= r[g[i]]
        r[i] = v
    r = np.array(r, dtype=np.int64)
    back = sum(r[T.gather(d)] for d in D)
    back[T.ZERO] = 0
    assert np.array_equal(back, f), "(8v+8s) does not divide a nu-slice"
    return r


def levels_data(N):
    """Graded r_n(x) and cumulative R^9_n for n <= N (checked against Fock dims)."""
    T = DomTable("D", 4, 4 * N * N)
    TB = DomTable("B", 4, 4 * N * N)
    b_to_d = np.array([T.index[tuple(int(a) for a in w)] for w in TB.W] + [T.ZERO])
    orbit = np.array([orbit_size_D4(tuple(int(a) for a in w)) for w in T.W] + [0])
    fd = fock_dims(N)
    F = graded_fock(T, N)
    for n in range(N + 1):
        assert sum(int((a * orbit).sum()) for a in F[n].values()) == fd[n]
    r, R9 = {}, {0: {}}
    flip = T.lookup(T.W * np.array([1, 1, 1, -1]))
    for n in range(1, N + 1):
        r[n] = {nu: divide(T, a) for nu, a in F[n].items()}
        r[n] = {nu: a for nu, a in r[n].items() if a.any()}
        for a in r[n].values():
            assert a.min() >= 0
        tot = sum(r[n].values())
        assert np.array_equal(tot[:T.size], tot[flip]), "r_n not W(B4)-invariant"
        b4 = TB.decompose(tot[b_to_d])
        assert b4.min() >= 0
        R9[n] = dict(R9[n - 1])
        for i in np.nonzero(b4[:TB.size])[0]:
            mu = tuple(int(a) for a in TB.W[i])
            R9[n][mu] = R9[n].get(mu, 0) + int(b4[i])
    return r, R9, T


def phi(lam, D2):
    """V_lam x^Delta  ->  {(nu, D4-dominant w): mult}, nu = Delta - J, J = x5."""
    out = {}
    for wt, m in d5_weights(lam).items():
        w4, J2 = wt[:4], wt[4]
        if not (w4[0] >= w4[1] >= w4[2] >= abs(w4[3])):
            continue  # W(D4)-invariant: compare at D4-dominant weights only
        nu2 = D2 - J2
        if nu2 % 2:
            continue
        key = (nu2 // 2, tuple(w4))
        out[key] = out.get(key, 0) + m
    return out


def default_lams(n, R9):
    classes = set(candidates(R9[n])) | (set(candidates(R9[n - 1])) if n > 1 else set())
    return sorted(classes | {c[:4] + (-c[4],) for c in classes if c[4]})


def bound2(lam):
    """2 x BBMS SO(10) bound (hep-th/0310292 eq. 2.15): 1 + k + 2l + 3m + 5p/2 + 3q/2."""
    k, l, m, p, q = dynkin_D5(lam)
    return 2 + 2 * k + 4 * l + 6 * m + 5 * p + 3 * q


def build_system(n, r, R9, T, slack, lams=None, gmax=None):
    """
    Domain: if gmax is None, Delta with [Delta - l1, Delta + l1] meeting the
    nu-support of r_n widened by slack; otherwise the 'unitarity window'
    0 <= Delta - bound(lambda) <= gmax.
    """
    nus = sorted(r[n])
    nu_lo, nu_hi = nus[0], nus[-1]
    if lams is None:
        lams = default_lams(n, R9)
    cols = []
    for lam in lams:
        l1, par = lam[0], lam[0] % 2
        if gmax is None:
            rng = range(2 * nu_lo - l1 - 2 * slack, 2 * nu_hi + l1 + 2 * slack + 1)
        else:
            rng = range(bound2(lam), bound2(lam) + 2 * gmax + 1)
        for D2 in rng:
            if D2 % 2 == par:
                cols.append((lam, D2))
    rows = {}
    for nu in nus:
        for i in np.nonzero(r[n][nu][:T.size])[0]:
            rows[(nu, tuple(int(a) for a in T.W[i]))] = len(rows)
    entries = {}
    for j, col in enumerate(cols):
        for key, m in phi(*col).items():
            if key not in rows:
                rows[key] = len(rows)
            entries[(rows[key], j)] = m
    A = matrix(QQ, len(rows), len(cols), entries, sparse=True)
    b = [0] * len(rows)
    for (nu, w), i in rows.items():
        if nu in r[n] and w in T.index:
            b[i] = int(r[n][nu][T.index[w]])
    return A, vector(QQ, b), cols, rows
