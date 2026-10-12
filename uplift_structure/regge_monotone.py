"""
S6: Regge monotonicity  r9_n(nu + e1) >= r9_{n-1}(nu)  via the null-cone arc space.

Theory (STATUS S6):
  * B = C[J_inf(Q)], Q = {y in C^9 : y.y = 0}, deg y_k = k+1, is a graded domain
    with Hilbert series PE[sum_m q^m 9] prod_{m>=2} (1 - q^m) (Mustata: jet
    schemes of the canonical lci cone Q are irreducible lci of expected
    dimension).  Multiplication by the highest-weight vector of y_0 is
    injective and preserves SO(9) highest-weight vectors, so B is
    Regge-monotone:  m_{B,d+1}(nu + e1) >= m_{B,d}(nu).
  * Crystal lemma: if A is Regge-monotone and M >= 0, then A (x) M is.
  * Hence per-level monotonicity of Z_m(q) = sum_n r9_n q^{n-1} follows from
        M'(q) := Z_m(q) / Hilb(B)(q)  having non-negative SO(9) multiplicities.
    In light-cone (SO(8)) terms Hilb(B) = F_B / (1-q) and
        M' = (F_F - F_B^{-1}) (1 - q) / (q (8v + 8s)),
    with F_B = PE[sum q^m 8v], F_F = prod_m prod_{8s}(1 + q^m z).
  * Lemma (also checked here): for the Regge induction A_n = T(A_{n-1}) + E_n,
    top(A_{n-1}) <= r9_{n-1}, so per-level monotonicity gives half (a)
    of the Regge continuation at every level, for any rule for E_n.

This script computes, for degrees d <= NMAX-1:
  (1) B_d as SO(9) characters: non-negative? Regge-monotone? (theory predicts yes)
  (2) M'_d as SO(9) characters: non-negative?   (the open input)
  (3) sanity: sum_i B_{d-i} (x) M'_i == r9_{d+1} for d <= 6 (Sage tensor products)
  (4) top(A_{n-1}) <= r9_{n-1} along the max-height Regge induction, n <= NMAX.

Output: results/regge_monotone.txt
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, SO10
from regge_stress import DEFAULT_EXT, T, load_R9, residual

sys.path.insert(0, SO10)
from derive_r9 import E1, HEIGHT, S8, V8
from uplift import height, interlacing, solve_ilp
from weights import DomTable, dynkin_B4


def exp_log(T4, N, G):
    """F with log F = sum_j G_j q^j / j  (G_j arrays);  n F_n = sum_j G_j F_{n-j}."""
    F = [T4.zeros()]
    F[0][T4.index[(0, 0, 0, 0)]] = 1
    for n in range(1, N + 1):
        acc = T4.zeros()
        for j in range(1, n + 1):
            Gj = G[j]
            if not Gj.any():
                continue
            # product of W-invariant functions: (G_j * F_{n-j}) via weight shifts of G_j's support
            acc += mult(T4, Gj, F[n - j])
        assert np.all(acc % n == 0)
        Fn = acc // n
        Fn[T4.ZERO] = 0
        F.append(Fn)
    return F


_ORBIT_CACHE = {}


def orbit(w):
    """All weights in the W(D4)-orbit of the dominant doubled weight w."""
    if w in _ORBIT_CACHE:
        return _ORBIT_CACHE[w]
    from itertools import permutations, product
    s = set()
    for p in permutations(w):
        for sg in product((1, -1), repeat=4):
            if sg.count(-1) % 2 == 0:
                s.add(tuple(p[i] * sg[i] for i in range(4)))
    _ORBIT_CACHE[w] = sorted(s)
    return _ORBIT_CACHE[w]


def mult(T4, a, b):
    """Product of two W(D4)-invariant functions given on dominant weights (a small support)."""
    out = T4.zeros()
    for i in np.nonzero(a[:T4.size])[0]:
        w = tuple(int(t) for t in T4.W[i])
        for x in orbit(w):
            out += a[i] * b[T4.gather(x)]
    out[T4.ZERO] = 0
    return out


def adams(T4, weights, k):
    """psi^k of a single W(D4)-orbit character (8v or 8s), as a dominant-weight array:
    multiplicity 1 at the dominant representative of k*w."""
    f = T4.zeros()
    doms = {tuple(int(t) for t in T4.dominate(np.array([[k * c for c in w]]))[0]) for w in weights}
    for d in doms:
        f[T4.index[d]] = 1
    return f


def divide(T4, f):
    """Exact division of a W(D4)-invariant array by 8_v + 8_s (as in so10_uplift/derive_r9.py)."""
    D = V8 + S8
    others = [d for d in D if d != E1]
    g_up = T4.gather(tuple(-a for a in E1)).tolist()
    g_oth = [T4.gather(tuple(d[i] - E1[i] for i in range(4))).tolist() for d in others]
    order = np.argsort(-(T4.W @ HEIGHT), kind="stable").tolist()
    fl = f.tolist()
    r = [0] * (T4.size + 1)
    for i in order:
        v = fl[g_up[i]]
        for g in g_oth:
            v -= r[g[i]]
        r[i] = v
    r = np.array(r, dtype=np.int64)
    back = sum(r[T4.gather(d)] for d in D)
    back[T4.ZERO] = 0
    assert np.array_equal(back, f), "(8v+8s) does not divide"
    return r


def main():
    N, r9, R9 = load_R9(DEFAULT_EXT)
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else N
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    Dmax = NMAX + 1
    T4 = DomTable("D", 4, 4 * Dmax * Dmax)
    TB = DomTable("B", 4, 4 * Dmax * Dmax)
    b_to_d = np.array([T4.index[tuple(int(a) for a in w)] for w in TB.W] + [T4.ZERO])
    flip = T4.lookup(T4.W * np.array([1, 1, 1, -1]))

    def so9(f):
        assert np.array_equal(f[:T4.size], f[flip]), "not W(B4)-invariant"
        b = TB.decompose(f[b_to_d])
        return {tuple(int(a) for a in TB.W[i]): int(b[i]) for i in np.nonzero(b[:TB.size])[0]}

    # G_j for log F_B = sum_m sum_k q^{mk} psi^k(8v)/k  and log F_F, log F_B^{-1}
    psiV = {k: adams(T4, V8, k) for k in range(1, Dmax + 1)}
    psiS = {k: adams(T4, S8, k) for k in range(1, Dmax + 1)}
    GB, GBinv, GF = [None], [None], [None]
    for j in range(1, Dmax + 1):
        gb, gf = T4.zeros(), T4.zeros()
        for k in range(1, j + 1):
            if j % k == 0:
                gb += (j // k) * psiV[k]
                gf += (j // k) * (-1) ** (k + 1) * psiS[k]
        GB.append(gb)
        GBinv.append(-gb)
        GF.append(gf)
    FB = exp_log(T4, Dmax, GB)
    FBinv = exp_log(T4, Dmax, GBinv)
    FF = exp_log(T4, Dmax, GF)
    log(f"Fock series to degree {Dmax} [{time.time() - t0:.0f}s]")

    # (1) B_d = sum_{i<=d} FB_i
    log("(1) B = C[J_inf(null cone of C^9)]: SO(9) content (theory: genuine rep, Regge-monotone)")
    Bso9, cum = [], T4.zeros()
    for d in range(NMAX):
        cum = cum + FB[d]
        Bso9.append(so9(cum))
    negB = [d for d, b in enumerate(Bso9) if min(b.values()) < 0]
    monoB = sum(1 for d in range(1, NMAX) for nu, v in Bso9[d - 1].items()
                if Bso9[d].get((nu[0] + 2,) + nu[1:], 0) < v)
    log(f"    degrees with negative multiplicities: {negB}; Regge-monotonicity violations: {monoB}")

    # (2) M'_d = [(FF - FBinv)(1 - q)]_{d+1} / (8v + 8s)
    log("(2) M' = Z_m / Hilb(B): SO(9) content (non-negativity => per-level Regge monotonicity)")
    diff = [FF[d] - FBinv[d] for d in range(Dmax + 1)]
    Mp, allpos = [], True
    for d in range(NMAX):
        num = diff[d + 1] - diff[d]
        m = so9(divide(T4, num))
        Mp.append(m)
        neg = {dynkin_B4(k): v for k, v in m.items() if v < 0}
        allpos &= not neg
        if neg or d < 8 or d % 4 == 3:
            log(f"    M'_{d:2d}: {sum(m.values()):10d} irreps, {len(m):5d} distinct; negative: "
                + (f"{list(neg.items())[:4]}" if neg else "none"))
    log(f"    M' >= 0 for all degrees < {NMAX}: {allpos}")

    # (3) sanity: sum_i B_{d-i} (x) M'_i == r9_{d+1}, d <= 6
    from sage.all import WeylCharacterRing
    B4 = WeylCharacterRing("B4", style="coroots")

    def ch(dct):
        return sum((v * B4(*dynkin_B4(k)) for k, v in dct.items()), B4.zero())
    ok = True
    for d in range(7):
        tot = sum((ch(Bso9[d - i]) * ch(Mp[i]) for i in range(d + 1)), B4.zero())
        ok &= tot == ch(r9[d + 1])
    log(f"(3) sanity  sum_i B_(d-i) (x) M'_i == r9_(d+1) for d <= 6: {ok}")

    # (4) lemma along the max-height induction
    A, bad = {}, 0
    for n in range(1, NMAX + 1):
        if n >= 2:
            top = {}
            for lam, c in A.items():
                for mu in interlacing(lam):
                    if mu[0] == lam[0]:
                        top[mu] = top.get(mu, 0) + c
            bad += sum(1 for mu, v in top.items() if v > r9[n - 1].get(mu, 0))
        rho = residual(r9[n], A)
        E = solve_ilp(rho, lambda l: height(l[:4]), maximize=True) if rho else {}
        A = T(A)
        for l, c in E.items():
            A[l] = A.get(l, 0) + c
    log(f"(4) top(A_(n-1)) <= r9_(n-1) along the max-height induction, n <= {NMAX}: violations {bad}")
    log(f"[{time.time() - t0:.0f}s]")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "regge_monotone.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
