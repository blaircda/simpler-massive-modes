"""
S8: multi-mode Regge induction.  Heads E_n = F_n + P_n, where the forced part
F_n = sum_{m<n} sum_{j>=1} p_K(n-m, j) T^j(P_m) collects Cartan shifts of earlier
primaries by modes k in K (level +k, spin +e1 each), p_K(N, j) = #partitions of N
into exactly j parts from K; vacuum exempt.  P_n >= 0 (new primaries) must
uplift rho_n - Res(F_n).  K = {2} is the "double Regge induction" (shifts T and
T_2); K = {3} tests a level-3 direction on its own; K = {2,3,...} would make
the uplift a Cartan module over all modes.

Output: results/regge_modes.txt
"""
import sys, time
from functools import lru_cache
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, SO10
sys.path.insert(0, SO10)
from regge_stress import DEFAULT_EXT, load_R9, residual, T
from uplift import restrict, solve_ilp, height
N, r9, R9 = load_R9(DEFAULT_EXT)
def shift(d, j): return {(l[0] + 2 * j,) + l[1:]: c for l, c in d.items()}
def run(K, rule):
    obj, mx = {"max-height": (lambda l: height(l[:4]), True), "min-height": (lambda l: height(l[:4]), False),
               "feasible": (None, False)}[rule]
    @lru_cache(None)
    def pK(Nn, j, mn=0):
        if j == 0: return 1 if Nn == 0 else 0
        return sum(pK(Nn - k, j - 1, k) for k in K if k >= mn and k <= Nn)
    A, Pm = {}, {}
    for n in range(1, N + 1):
        rho = residual(r9[n], A)
        F = {}
        for m, Pd in Pm.items():
            for j in range(1, n - m + 1):
                c = pK(n - m, j)
                if c:
                    for l, v in shift(Pd, j).items(): F[l] = F.get(l, 0) + c * v
        rest = dict(rho)
        for mu, v in restrict(F).items(): rest[mu] = rest.get(mu, 0) - v
        if min(rest.values(), default=0) < 0:
            return f"STUCK at n={n} (forced copies overshoot)"
        rest = {k: v for k, v in rest.items() if v}
        P = solve_ilp(rest, obj, maximize=mx) if rest else {}
        if P is None:
            return f"STUCK at n={n} (remainder not liftable)"
        if n > 1: Pm[n] = P                      # vacuum (n=1) is not a primary
        E = dict(F)
        for l, c in P.items(): E[l] = E.get(l, 0) + c
        A = T(A)
        for l, c in E.items(): A[l] = A.get(l, 0) + c
        assert restrict(A) == R9[n]
    tot = {n: sum(P.values()) for n, P in Pm.items()}
    return f"ok to n={N}; primaries per level: {[tot.get(n, 0) for n in range(2, min(N, 16) + 1)]}"
lines = []
for K in ((2,), (3,), (2, 3), (2, 3, 4), tuple(range(2, 40))):
    for rule in ("max-height", "min-height", "feasible"):
        t0 = time.time()
        s = f"K={K if len(K) < 5 else 'all k>=2'} [{rule}]: {run(K, rule)} [{time.time()-t0:.0f}s]"
        print(s, flush=True)
        lines.append(s)
os.makedirs(RESULTS, exist_ok=True)
with open(os.path.join(RESULTS, "regge_modes.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")
