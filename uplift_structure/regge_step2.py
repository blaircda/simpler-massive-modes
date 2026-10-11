"""
S5 step 2 diagnostics (reads r9 data only; see STATUS S5).

  duals     : optimal LP dual certificates for max_A top(A)(nu) at the low
              weights where the universal certificate y_nu is not optimal
              (singlet, vector, spinor), n in {13, 16, 20}.
  singlet   : r9_n(9) versus R9_{n-1}(1) (the near-maximal SO(10)-singlet
              content of an uplift), and the singlet content of the
              max-height Regge induction, n <= 32.
  monotone  : Regge monotonicity r9_n(nu+e1) >= r9_{n-1}(nu) and its cumulative
              version, all nu, n <= 32.
  wfactor   : W(q) = Z_m(q) PE[-q 9] (equivalently R9(q) Res lambda_{-q}(10),
              so that R9(q) = Res PE[q 10] . W(q)); per level: W_j >= 0?
              W_j SO(10)-liftable (non-negatively)?  (Klimyk on B4.)

Usage: regge_step2.py MODE      Output: results/regge_step2_<mode>.txt
"""

import os
import sys
import time
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS
from regge_stress import DEFAULT_EXT, T, load_R9, residual

from uplift import candidates, height, interlacing, solve_ilp
from weights import dynkin_B4


def mode_duals(N, r9, R9, log):
    from sage.all import MixedIntegerLinearProgram
    for n in (13, 16, 20):
        R = R9[n - 1]
        cands = candidates(R)
        mus = sorted(R)
        cover, topof = {mu: [] for mu in mus}, {}
        for j, lam in enumerate(cands):
            for mu in interlacing(lam):
                cover[mu].append(j)
                if mu[0] == lam[0]:
                    topof.setdefault(mu, []).append(j)
        for nu in [(0, 0, 0, 0), (2, 0, 0, 0), (1, 1, 1, 1)]:
            p = MixedIntegerLinearProgram(maximization=True, solver="GLPK")
            x = p.new_variable(nonnegative=True)
            for mu in mus:
                p.add_constraint(p.sum(x[j] for j in cover[mu]) == R[mu])
            p.solver_parameter("simplex_or_intopt", "simplex_only")
            p.set_objective(p.sum(x[j] for j in topof.get(nu, [])))
            val = p.solve()
            b = p.get_backend()
            y = {mus[i]: b.get_row_dual(i) for i in range(len(mus))}
            y = {k: v for k, v in y.items() if abs(v) > 1e-9}
            # validate: y.R equals the LP value (strong duality)
            assert abs(sum(v * R[k] for k, v in y.items()) - val) < 1e-6
            log(f"  n={n} nu={dynkin_B4(nu)} (Dynkin): LP max top = {val:.0f}; optimal dual certificate "
                + ", ".join(f"{v:+.0f}*{dynkin_B4(k)}" for k, v in sorted(y.items(), key=lambda t: (sum(t[0]), t[0]))))


def mode_singlet(N, r9, R9, log):
    Z, V = (0, 0, 0, 0), (2, 0, 0, 0)
    log("   n    r9_n(9)  R9_{n-1}(1)  ratio | max-height induction: A_{n-1}(SO(10) singlet)")
    A = {}
    for n in range(1, N + 1):
        rho = residual(r9[n], A)
        E = solve_ilp(rho, lambda l: height(l[:4]), maximize=True) if rho else {}
        if n >= 8:
            rv, Rs = r9[n].get(V, 0), R9[n - 1].get(Z, 0)
            log(f"  {n:2d} {rv:10d} {Rs:12d}  {rv / Rs:5.2f} | {A.get((0, 0, 0, 0, 0), 0)}")
        A = T(A)
        for l, c in E.items():
            A[l] = A.get(l, 0) + c


def mode_monotone(N, r9, R9, log):
    bad_cum = bad_lvl = 0
    for n in range(2, N + 1):
        for nu, v in R9[n - 1].items():
            if R9[n].get((nu[0] + 2,) + nu[1:], 0) < v:
                bad_cum += 1
        for nu, v in r9[n - 1].items():
            if r9[n].get((nu[0] + 2,) + nu[1:], 0) < v:
                bad_lvl += 1
    log(f"  per level  r9_n(nu+e1) >= r9_(n-1)(nu): violations {bad_lvl} (all nu, n <= {N})")
    log(f"  cumulative R9_n(nu+e1) >= R9_(n-1)(nu): violations {bad_cum} (all nu, n <= {N})")


RHO_B4 = (7, 5, 3, 1)
VEC9 = [tuple(2 * s if i == k else 0 for i in range(4)) for k in range(4) for s in (1, -1)] + [(0, 0, 0, 0)]
LAMBDA9 = {k: [tuple(sum(w[i] for w in S) for i in range(4)) for S in combinations(VEC9, k)]
           for k in range(10)}


def _dom_sign_B4(v):
    a = [abs(t) for t in v]
    if 0 in a or len(set(a)) < 4:
        return None, 0
    order = sorted(range(4), key=lambda i: -a[i])
    sgn, seen = 1, [False] * 4
    for i in range(4):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = order[j]
                L += 1
            if L % 2 == 0:
                sgn = -sgn
    sgn *= (-1) ** sum(1 for t in v if t < 0)
    return tuple(a[i] for i in order), sgn


def tensor_B4(char, weights, coeff):
    """char (x) U for U given by its weight list (Klimyk), times coeff."""
    out = {}
    for mu, m in char.items():
        for w in weights:
            d, s = _dom_sign_B4(tuple(mu[i] + RHO_B4[i] + w[i] for i in range(4)))
            if d is None:
                continue
            lam = tuple(d[i] - RHO_B4[i] for i in range(4))
            out[lam] = out.get(lam, 0) + coeff * s * m
    return out


def mode_wfactor(N, r9, R9, log):
    t = {k: v for k, v in tensor_B4({(2, 0, 0, 0): 1}, LAMBDA9[1], 1).items() if v}
    assert t == {(0, 0, 0, 0): 1, (2, 2, 0, 0): 1, (4, 0, 0, 0): 1}, "Klimyk check 9x9 failed"
    negs, unl = [], []
    for j in range(N):
        W = {}
        for k in range(0, min(9, j) + 1):
            for lam, v in tensor_B4(r9[j - k + 1], LAMBDA9[k], (-1) ** k).items():
                W[lam] = W.get(lam, 0) + v
        W = {k: v for k, v in W.items() if v}
        neg = {dynkin_B4(k): v for k, v in W.items() if v < 0}
        lift = (solve_ilp(W) is not None) if (W and not neg) else None
        if neg:
            negs.append(j)
        if lift is False:
            unl.append(j)
        log(f"  W_{j:2d}: {sum(W.values()):9d} irreps; negative entries {len(neg)}"
            + (f" {list(neg.items())[:3]}" if neg else f"; non-negatively SO(10)-liftable: {lift}"))
    log(f"  W_j has negative entries at j = {negs}; W_j >= 0 but not liftable at j = {unl}")


def main():
    mode = sys.argv[1]
    N, r9, R9 = load_R9(DEFAULT_EXT)
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    {"duals": mode_duals, "singlet": mode_singlet, "monotone": mode_monotone,
     "wfactor": mode_wfactor}[mode](N, r9, R9, log)
    log(f"[{time.time() - t0:.0f}s]")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, f"regge_step2_{mode}.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
