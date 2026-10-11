"""
Library for the BBMS-graded lift (see derive_graded.py for the method).

    levels_data(N)            -> r[n][nu] (D4 dominant arrays), R9[n] (B4 dicts), T
    build_system(n, r, R9, T, slack, lams=None, gmax=None)
                              -> (sparse integer entries, b, cols, rows)
    solve_exact(entries, b, nrows, ncols) -> (nullity, integer solution or None)
    phi(lam, D2)              -> image of V_lam x^Delta as {(nu, w4): mult}
    d5_slices(lam)            -> {(J2, D4-dominant w4): mult}  (fast Freudenthal)
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


# ------------------------------------------------------------------
# Fast D5 weight data: Freudenthal on dominant weights only
# ------------------------------------------------------------------

_POS_ROOTS = []
for _i in range(5):
    for _j in range(_i + 1, 5):
        for _sg in (1, -1):
            _r = [0] * 5
            _r[_i], _r[_j] = 2, 2 * _sg
            _POS_ROOTS.append(tuple(_r))
_RHO2 = (8, 6, 4, 2, 0)


def _dom5(v):
    a = sorted((abs(t) for t in v), reverse=True)
    if a[4] and sum(1 for t in v if t < 0) % 2:
        a[4] = -a[4]
    return tuple(a)


def _simple_coeffs2(v):
    """Doubled simple-root coefficients x2 of a doubled D5 vector (D5 simple roots)."""
    c1 = v[0]
    c2 = v[0] + v[1]
    c3 = c2 + v[2]
    c4 = (c3 + v[3] - v[4]) // 2
    c5 = (c3 + v[3] + v[4]) // 2
    return (c1, c2, c3, c4, c5)


_DOMM = {}


def d5_dominant_mults(lam):
    """{dominant doubled weight: multiplicity} for the D5 irrep lam (Freudenthal)."""
    if lam in _DOMM:
        return _DOMM[lam]
    n2 = sum(a * a for a in lam)
    # all dominant mu <= lam: lam - mu a non-negative integer combination of simple roots
    doms = []
    B = lam[0]
    def rec(prefix):
        k = len(prefix)
        if k == 5:
            mu = tuple(prefix)
            diff = tuple(lam[i] - mu[i] for i in range(5))
            c = _simple_coeffs2(diff)
            if all(t >= 0 and t % 2 == 0 for t in c):
                doms.append(mu)
            return
        hi = prefix[-1] if prefix else B
        lo = -prefix[-1] if k == 4 else 0
        for a in range(hi, lo - 1, -1):
            if (a - lam[0]) % 2 == 0 and sum(t * t for t in prefix) + a * a <= n2:
                rec(prefix + [a])
    rec([])
    doms.sort(key=lambda mu: sum(_simple_coeffs2(tuple(lam[i] - mu[i] for i in range(5)))))
    LP = [lam[i] + _RHO2[i] for i in range(5)]
    LP2 = sum(t * t for t in LP)
    m = {lam: 1}
    for mu in doms[1:]:
        MP = [mu[i] + _RHO2[i] for i in range(5)]
        den = LP2 - sum(t * t for t in MP)
        tot = 0
        for al in _POS_ROOTS:
            k = 1
            while True:
                v = tuple(mu[i] + k * al[i] for i in range(5))
                mv = m.get(_dom5(v), 0)
                if not mv:
                    break
                tot += mv * sum(v[i] * al[i] for i in range(5))
                k += 1
        val, rem = divmod(2 * tot, den)
        assert rem == 0, "Freudenthal not integral"
        if val:
            m[mu] = val
    _DOMM[lam] = m
    return m


_SLICES = {}


def d5_slices(lam):
    """{(J2, D4-dominant w4): mult} over all weights (w4, J) of V_lam (J = x5)."""
    if lam in _SLICES:
        return _SLICES[lam]
    out = {}
    for nu, mult in d5_dominant_mults(lam).items():
        absn = [abs(t) for t in nu]
        haszero = 0 in absn
        par = 1 if nu[4] < 0 else 0
        seen = set()
        for i in range(5):
            for s in ((1, -1) if absn[i] else (1,)):
                J = s * absn[i]
                rest = sorted(absn[:i] + absn[i + 1:], reverse=True)
                if haszero:  # a zero coordinate makes every sign pattern reachable
                    forms = {tuple(rest), tuple(rest[:3] + [-rest[3]])}
                else:
                    neg = (par - (1 if s < 0 else 0)) % 2
                    forms = {tuple(rest[:3] + [-rest[3]]) if neg else tuple(rest)}
                for w4 in forms:
                    key = (J, w4)
                    if key not in seen:
                        seen.add(key)
                        out[key] = out.get(key, 0) + mult
    _SLICES[lam] = out
    return out


def phi(lam, D2):
    """V_lam x^Delta  ->  {(nu, D4-dominant w): mult}, nu = Delta - J, J = x5."""
    out = {}
    for (J2, w4), m in d5_slices(lam).items():
        nu2 = D2 - J2
        if nu2 % 2:
            continue
        key = (nu2 // 2, w4)
        out[key] = out.get(key, 0) + m
    return out


def default_lams(n, R9):
    classes = set(candidates(R9[n])) | (set(candidates(R9[n - 1])) if n > 1 else set())
    return sorted(classes | {c[:4] + (-c[4],) for c in classes if c[4]})


def bound2(lam):
    """2 x BBMS SO(10) bound (hep-th/0310292 eq. 2.15): 1 + k + 2l + 3m + 5p/2 + 3q/2."""
    k, l, m, p, q = dynkin_D5(lam)
    return 2 + 2 * k + 4 * l + 6 * m + 5 * p + 3 * q


def bound2_min(lam):
    """2 x min(Delta_+, Delta_-), Delta_+- = 1 + k + 2l + 3m + 2(p+q) +- (p-q)/2
    (Bianchi hep-th/0409304 after eq. 3.17): 2 + L1 + L2 + L3 + L4 - |L5| (doubled)."""
    return 2 + lam[0] + lam[1] + lam[2] + lam[3] - abs(lam[4])


BOUNDS = {"plus": bound2, "min": bound2_min}


def domain_cols(n, r, R9, slack, lams=None, gmax=None, bound="plus"):
    """
    Unknowns (lambda, 2 Delta).  If gmax is None: Delta with [Delta - l1,
    Delta + l1] meeting the nu-support of r_n widened by slack; otherwise the
    'unitarity window' 0 <= Delta - bound(lambda) <= gmax, with bound = BBMS's
    Delta_+ ("plus", eq. 2.15 of hep-th/0310292) or min(Delta_+, Delta_-) ("min").
    """
    bfun = BOUNDS[bound]
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
            rng = range(bfun(lam), bfun(lam) + 2 * gmax + 1)
        for D2 in rng:
            if D2 % 2 == par:
                cols.append((lam, D2))
    return cols


def build_system(n, r, R9, T, slack, lams=None, gmax=None):
    """Global sparse system (small n / diagnostics); see solve_level for large n."""
    nus = sorted(r[n])
    cols = domain_cols(n, r, R9, slack, lams=lams, gmax=gmax)
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
    b = [0] * len(rows)
    for (nu, w), i in rows.items():
        if nu in r[n] and w in T.index:
            b[i] = int(r[n][nu][T.index[w]])
    return entries, b, cols, rows


PRIMES = [8388593, 8388587, 8388581, 8388571, 8388547, 8388539]  # < 2^23: fast LinBox path


def solve_exact(entries, b, nrows, ncols):
    """
    Integer solution of A x = b (A given by sparse integer entries), via dense
    solves over GF(p) for word-size primes, CRT to symmetric residues, and an
    exact integer check of A x = b.  Returns (nullity, x or None); nullity is
    computed mod PRIMES[0]; it equals the rational nullity unless p divides
    every maximal non-zero minor (then it would be an over-estimate).
    """
    from sage.all import CRT_list
    sols, mods, nullity = [], [], None
    for p in PRIMES:
        F = GF(p)
        A = matrix(F, nrows, ncols, entries, sparse=False)
        if nullity is None:
            nullity = int(ncols - A.rank())
        try:
            xp = A.solve_right(vector(F, b))
        except ValueError:
            return nullity, None
        sols.append([int(t) for t in xp])
        mods.append(p)
        M = 1
        for q in mods:
            M *= q
        x = []
        for j in range(ncols):
            v = CRT_list([s[j] for s in sols], mods) if len(mods) > 1 else sols[0][j]
            v = int(v) % M
            x.append(v - M if v > M // 2 else v)
        Ax = [0] * nrows
        for (i, j), a in entries.items():
            if x[j]:
                Ax[i] += a * x[j]
        if Ax == list(b):
            return nullity, x
        if nullity:
            # with a kernel, independent mod-p solves need not be compatible
            return nullity, None
    raise RuntimeError("no exact solution reconstructed; add primes")


# ------------------------------------------------------------------
# Memory-lean block solver (large n)
# ------------------------------------------------------------------
#
# The system is block diagonal with 4 blocks, keyed by
#   column (lambda, 2Delta):  (class bit l1 mod 2,  (sum(lambda) - Delta) mod 2)
#   row    (nu, w):           (class bit w1 mod 2,  (sum(w) - nu) mod 2)
# (central character of Spin(10), and Sum w_i + J = Sum lambda_i mod 2 on all
# weights of V_lambda, with J = Delta - nu).  In the unitarity window the
# second bit is the gap parity.  Each block is built directly as a dense
# GF(p) matrix from compact per-irrep weight arrays, so no global entry list
# is ever held in memory.

_W4ID, _W4 = {}, []


def _wid(w):
    i = _W4ID.get(w)
    if i is None:
        i = _W4ID[w] = len(_W4)
        _W4.append(w)
    return i


_SLA = {}


def d5_slice_arrays(lam):
    """(J2, wid, mult) int64 arrays of the (J, D4-dominant w4) weights of V_lam."""
    if lam not in _SLA:
        sl = d5_slices(lam)
        J2 = np.fromiter((k[0] for k in sl), np.int64, len(sl))
        wid = np.fromiter((_wid(k[1]) for k in sl), np.int64, len(sl))
        m = np.fromiter(sl.values(), np.int64, len(sl))
        _SLA[lam] = (J2, wid, m)
        _SLICES.pop(lam, None)   # keep only the compact form
        _DOMM.pop(lam, None)
    return _SLA[lam]


def col_block(lam, D2):
    return (lam[0] % 2, ((sum(lam) - D2) // 2) % 2)


def row_block(nu, w):
    return (w[0] % 2, ((sum(w) - 2 * nu) // 2) % 2)


_OFF, _BIG = 1 << 20, 1 << 32   # row key = (nu + _OFF) * _BIG + wid


def _col_keys(lam, D2):
    J2, wid, m = d5_slice_arrays(lam)
    nu2 = D2 - J2
    ok = nu2 % 2 == 0
    return (nu2[ok] // 2 + _OFF) * _BIG + wid[ok], m[ok]


def _solve_block_mod(p, cols, cj, tgt_keys, tgt_vals):
    """Echelonize [A_block | b_block] over GF(p).  Returns (residues {j: x_j}, nullity, consistent)."""
    ks, cs, vs = [], [], []
    for a, j in enumerate(cj):
        k, m = _col_keys(*cols[j])
        ks.append(k)
        cs.append(np.full(k.size, a, dtype=np.int64))
        vs.append(m)
    keys = np.concatenate(ks + [tgt_keys])
    uk, inv = np.unique(keys, return_inverse=True)
    nA = keys.size - tgt_keys.size
    ri, ti = inv[:nA], inv[nA:]
    cc, vv = np.concatenate(cs), np.concatenate(vs)
    nr, nc = uk.size, len(cj)
    A = matrix(GF(p), nr, nc + 1, sparse=False)
    for i, j, v in zip(ri.tolist(), cc.tolist(), (vv % p).tolist()):
        A[i, j] = v
    for i, v in zip(ti.tolist(), [t % p for t in tgt_vals]):
        A[i, nc] = v
    del ks, cs, vs, keys, cc, vv
    A.echelonize()
    piv = A.pivots()
    if piv and piv[-1] == nc:
        return None, nc - (len(piv) - 1), False
    res = {cj[c]: int(A[k, nc]) for k, c in enumerate(piv)}
    return res, nc - len(piv), True


def solve_level(n, r, R9, T, gmax, slack=0, primes=PRIMES, log=None, bound="plus"):
    """
    Solve phi(X_n) = r_n(x) on the domain_cols domain, block by block.
    Returns (nullity, {col: integer coefficient} or None, block stats).
    Exactness: the CRT-lifted solution is checked against r_n over Z
    (every (nu, w) key, including those where r_n vanishes).
    """
    from sage.all import CRT_list
    cols = domain_cols(n, r, R9, slack, gmax=gmax, bound=bound)
    tgt = {}
    for nu, a in r[n].items():
        for i in np.nonzero(a[:T.size])[0]:
            w = tuple(int(t) for t in T.W[i])
            tgt[(nu + _OFF) * _BIG + _wid(w)] = int(a[i])
    blocks = {}
    for j, c in enumerate(cols):
        blocks.setdefault(col_block(*c), []).append(j)
    tgt_by_block = {}
    for k, v in tgt.items():
        nu, wid = k // _BIG - _OFF, k % _BIG
        tgt_by_block.setdefault(row_block(nu, _W4[wid]), []).append((k, v))
    for b in tgt_by_block:
        if b not in blocks:
            return None, None, {"error": f"target in block {b} with no unknowns"}
    residues, mods, stats, nullity = [], [], {}, 0
    for p in primes:
        resp = {}
        for b, cj in sorted(blocks.items()):
            tk = tgt_by_block.get(b, [])
            tkeys = np.array([k for k, _ in tk], dtype=np.int64)
            res, nl, ok = _solve_block_mod(p, cols, cj, tkeys, [v for _, v in tk])
            if p == primes[0]:
                stats[str(b)] = {"cols": len(cj), "nullity": nl, "rhs_rows": len(tk)}
                nullity += nl
                if log:
                    log(f"   block {b}: {len(cj)} unknowns, nullity {nl}, consistent {ok}")
            if not ok:
                return nullity, None, stats
            resp.update(res)
        residues.append(resp)
        mods.append(p)
        M = 1
        for q in mods:
            M *= q
        x = {}
        for j in set().union(*residues):
            v = CRT_list([rs.get(j, 0) for rs in residues], mods) if len(mods) > 1 else residues[0][j]
            v = int(v) % M
            v = v - M if v > M // 2 else v
            if v:
                x[cols[j]] = v
        # exact check over Z: sum_j x_j phi(col_j) == r_n on every key
        acc = {}
        for (lam, D2), c in x.items():
            k, m = _col_keys(lam, D2)
            for kk, mm in zip(k.tolist(), m.tolist()):
                acc[kk] = acc.get(kk, 0) + c * mm
        acc = {k: v for k, v in acc.items() if v}
        if acc == tgt:
            return nullity, x, stats
        if nullity:
            return nullity, None, stats
    raise RuntimeError("no exact solution reconstructed; add primes")
