"""
B6 diagnostic: does ANY graded lift (window of width >= 10) give
chirality-resolved non-negative cumulative sums?

Graded lifts differ from the BBMS window solution by multiples of
K = lambda_{-x}(10) (B2).  For window width G >= 10 the single-term
multiples K V_mu x^Delta (gap 0) lie in the window, and at x = 1 every
K V_mu x^Delta becomes e V_mu, independently of Delta.  Hence the
cumulative sums reachable by graded lifts are exactly

    R_n = R_n^BBMS + sum_mu t_mu e V_mu,      t_mu in Z,

(chirality-resolved; e = sum_k (-1)^k Lambda^k(10)).  Increments at
successive levels add up, so per level the question is independent.
This script decides, by an integer linear program per level, whether some
t makes R_n >= 0 with support inside the (resolved) uplift candidates of
R^9_n (positivity forces that support).  It then bounds the minimal such
|t|_1 (the distance of the BBMS normal form from positivity): LP lower
bound, and the best integer value found within a time limit.

Generators: mu over resolved candidates of R^9_n with supp(e V_mu) inside
the resolved candidates ("inner" generators, the S3b-type set); with
--all-mu, any resolved candidate mu, with c_lambda = 0 imposed outside.

e V_mu is computed by Klimyk's formula over the 243 weights of
e = prod_i (2 - z_i - 1/z_i) and checked against Sage on small mu.

Usage: positive_graded.py [n1 n2 ...] [--all-mu]
Output: results/positive_graded.txt
"""

import json
import os
import sys
import time
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "so10_uplift"))
from uplift import candidates, restrict
from weights import dynkin_D5

RHO2 = (8, 6, 4, 2, 0)
E_WEIGHTS = []   # (doubled weight, coefficient) of e
for eps in product((-1, 0, 1), repeat=5):
    c = 1
    for t in eps:
        c *= 2 if t == 0 else -1
    E_WEIGHTS.append((tuple(2 * t for t in eps), c))


def _dominate_signed(v):
    """(dominant v', det w) with v' = w v strictly dominant for D5, or (None, 0) on a wall."""
    a = [abs(t) for t in v]
    if len(set(a)) < 5:
        return None, 0
    order = sorted(range(5), key=lambda i: -a[i])
    # sign of the sorting permutation
    perm = order[:]
    sgn = 1
    seen = [False] * 5
    for i in range(5):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                L += 1
            if L % 2 == 0:
                sgn = -sgn
    out = [a[i] for i in order]
    if out[4] and sum(1 for t in v if t < 0) % 2:
        out[4] = -out[4]
    return tuple(out), sgn


_EV = {}


def e_times(mu):
    """Resolved decomposition of e (x) V_mu: {lambda (doubled): coeff} (Klimyk)."""
    if mu in _EV:
        return _EV[mu]
    out = {}
    for w, c in E_WEIGHTS:
        v = tuple(mu[i] + RHO2[i] + w[i] for i in range(5))
        d, s = _dominate_signed(v)
        if d is None:
            continue
        lam = tuple(d[i] - RHO2[i] for i in range(5))
        out[lam] = out.get(lam, 0) + s * c
    out = {k: v for k, v in out.items() if v}
    _EV[mu] = out
    return out


def check_klimyk():
    from sage.all import WeylCharacterRing
    D5 = WeylCharacterRing("D5", style="coroots")
    V = D5(1, 0, 0, 0, 0)
    e = sum((-1) ** k * (V.exterior_power(k) if k else D5.one()) for k in range(11))
    fw = D5.fundamental_weights()
    for dyn in [(0, 0, 0, 0, 0), (1, 0, 0, 0, 0), (0, 0, 0, 0, 1), (0, 0, 0, 1, 0), (2, 1, 0, 0, 1),
                (0, 1, 1, 0, 0), (1, 0, 0, 2, 1)]:
        mu = tuple(int(2 * a) for a in sum(c * fw[i + 1] for i, c in enumerate(dyn)).to_vector())
        sage_res = {}
        for w, m in (e * D5(*dyn)).monomial_coefficients().items():
            sage_res[tuple(int(2 * a) for a in w.to_vector())] = int(m)
        assert e_times(mu) == sage_res, f"Klimyk mismatch at {dyn}"


def resolved(classes):
    return sorted(set(classes) | {c[:4] + (-c[4],) for c in classes if c[4]})


def _build(R, gens, cand, integer, objective):
    from sage.all import MixedIntegerLinearProgram
    support = set(R)
    for mu in gens:
        support |= set(e_times(mu))
    p = MixedIntegerLinearProgram(maximization=False, solver="GLPK")
    t = p.new_variable(integer=integer)
    rows = {lam: [] for lam in support}
    for mu in gens:
        for lam, c in e_times(mu).items():
            rows[lam].append((mu, c))
    for lam, terms in rows.items():
        expr = R.get(lam, 0) + p.sum(c * t[mu] for mu, c in terms)
        if lam in cand:
            p.add_constraint(expr >= 0)
        else:
            p.add_constraint(expr == 0)
    if objective:
        u = p.new_variable(nonnegative=True)
        for mu in gens:
            p.add_constraint(u[mu] - t[mu] >= 0)
            p.add_constraint(u[mu] + t[mu] >= 0)
        p.set_objective(p.sum(u[mu] for mu in gens))
    else:
        p.set_objective(None)
    return p, t


def _apply(R, tv):
    c = dict(R)
    for mu, k in tv.items():
        for lam, v in e_times(mu).items():
            c[lam] = c.get(lam, 0) + k * v
    return {k: v for k, v in c.items() if v}


def solve_level(n, R, all_mu, timelimit=120):
    """
    Returns dict: generators, feasible, l1_feasible (|t|_1 of the feasibility
    solution), l1_lp (LP lower bound on min |t|_1), l1_best (best integer
    |t|_1 within the time limit), l1_optimal (bool).
    """
    from sage.numerical.mip import MIPSolverException
    R9 = restrict(R)
    cand = set(resolved(candidates(R9)))
    gens = sorted(cand) if all_mu else [mu for mu in sorted(cand) if set(e_times(mu)) <= cand]
    out = {"generators": len(gens)}
    p, t = _build(R, gens, cand, True, False)
    try:
        p.solve()
    except MIPSolverException:
        out["feasible"] = False
        return out
    tv = {mu: int(round(v)) for mu, v in p.get_values(t).items() if round(v)}
    c = _apply(R, tv)
    assert min(c.values()) >= 0 and restrict(c) == R9 and set(c) <= cand
    out.update(feasible=True, l1_feasible=sum(abs(v) for v in tv.values()))
    # LP lower bound on the minimal correction
    p, t = _build(R, gens, cand, False, True)
    out["l1_lp"] = round(float(p.solve()), 3)
    # integer minimisation, time-limited
    p, t = _build(R, gens, cand, True, True)
    p.solver_parameter("timelimit", timelimit)
    try:
        val = p.solve()
        tv = {mu: int(round(v)) for mu, v in p.get_values(t).items() if round(v)}
        c = _apply(R, tv)
        if min(c.values()) >= 0 and restrict(c) == R9:
            import math
            out["l1_best"] = sum(abs(v) for v in tv.values())
            out["min_mult_best"] = min(c.values())
            # |t|_1 is an integer, so best == ceil(LP bound) proves optimality
            out["l1_optimal"] = out["l1_best"] == math.ceil(out["l1_lp"] - 1e-6)
            sp = sum(abs(v) for mu, v in tv.items() if mu[0] % 2)
            out["support"] = {"spinor_l1": sp, "tensor_l1": out["l1_best"] - sp,
                              "largest": sorted(((abs(v), v, dynkin_D5(mu)) for mu, v in tv.items()),
                                                reverse=True)[:6]}
    except MIPSolverException as ex:
        out["l1_best"] = None
        out["l1_note"] = str(ex)[:80]
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    all_mu = "--all-mu" in sys.argv
    levels = [int(a) for a in args] or list(range(17, 25))
    check_klimyk()
    data = json.load(open(os.path.join(HERE, "data", "bbms_unitarity_g8.json")))["levels"]
    lines = [f"generators: {'all resolved candidates' if all_mu else 'inner (supp e V_mu in candidates)'}"]
    print(lines[0], flush=True)
    for n in levels:
        t0 = time.time()
        R = {tuple(e["ambient2"]): e["mult"] for e in data[str(n)]["R_cumulative_D5"]}
        o = solve_level(n, R, all_mu)
        s = f"level {n:2d}: {o['generators']} generators; resolved-positive graded lift exists: {o['feasible']}"
        if o["feasible"]:
            import math
            lo = math.ceil(o["l1_lp"] - 1e-6)
            if o.get("l1_optimal"):
                s += f"; min |t|_1 = {o['l1_best']} (= ceil LP bound, optimal)"
            else:
                s += f"; min |t|_1 in [{lo}, {o.get('l1_best')}] (time-limited)"
            if "support" in o:
                sup = o["support"]
                s += (f"\n      correction: spinor-mu |t|_1 {sup['spinor_l1']}, tensor-mu |t|_1 "
                      f"{sup['tensor_l1']}; largest t_mu (Dynkin): "
                      + ", ".join(f"{v:+d} e.{d}" for _, v, d in sup["largest"]))
        s += f" [{time.time() - t0:.0f}s]"
        print(s, flush=True)
        lines.append(s)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    name = "positive_graded_allmu.txt" if all_mu else "positive_graded.txt"
    with open(os.path.join(HERE, "results", name), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
