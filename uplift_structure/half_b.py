"""
S11: half (b) of the Regge step -- liftability of rho_n = r9_n - S(top A_(n-1)).

A class-level lift E of an SO(9) character rho is a non-negative combination of
interlacing boxes.  Processing weights top-down, the height-maximal remaining
weight mu must be the top (l1..l4) of whatever lifts it, so a lift is fixed by
choosing l5 in [floor, mu4] for every top (floor = 0 for tensors, 1/2 for
spinors).  This script reports (rho_n along the max-height and min-height Regge
inductions, n <= NMAX; cumulative R9_n, n <= 24):

  (1) top-down greedy lifts with a fixed l5 rule ("max": l5 = mu4, "min":
      l5 = floor) and with an adaptive rule (as many boxes as possible reach
      down the mu4 column): where they fail;
  (2) restricted box families (ILP): l5 extremal (l5 in {floor, l4});
      l3 extremal (l3 in {l2, l4}); and one box that is not a non-negative
      combination of l5-extremal boxes (so the l5-extremal cone is a proper
      subcone);
  (3) the level recursion for l5-extremal lifts.  At level t = mu4 the tops
      (mu1, mu2, mu3) are the UNIQUE finite-difference lift T_t of the level
      remainder; the only choice is the degenerate part D_t <= T_t (l5 = t;
      the rest have l5 = floor), and
          T_(t-1) = Lift(rho_(t-1) - rho_t) + M(D_t),
      where one degenerate top (m1, m2, m3) contributes +[m1,m2,m3] and, at the
      boundary tops mu3 = t-1, -[m1,m2] + [m1,m3-1] - [m2-1,m3-1].
      Rules: no degenerate tops; all degenerate; degenerate only where the
      non-boundary entries force it; and a per-level ILP (min or max total
      degenerate) choosing D_t so that T_(t-1) >= 0, with no lookahead;
  (4) where the per-level choices matter (max-height induction, n <= NMAX;
      cumulative R9_n, n <= 24): start from
      the forced degenerate tops (D >= -Lift(rho_(t-1) - rho_t) at non-boundary
      tops) and record the boundary tops (mu3 = t-1) that are still negative,
      by sector, level transition and k2 - (t-1); and whether the per-level
      ILP's extra degenerate tops repair them directly through +[m1, m3-1].

Output: results/half_b.txt
"""

import os
import sys
import time
from itertools import product

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, SO10
from regge_stress import DEFAULT_EXT, T, load_R9, residual

sys.path.insert(0, SO10)
from uplift import candidates, greedy, height, interlacing, restrict, solve_ilp

N, r9, R9 = load_R9(DEFAULT_EXT)


def adaptive(m):
    """Top-down lift; for each forced top take the pointwise-largest feasible profile down the mu4 column."""
    rem, out = dict(m), {}
    while rem:
        mu = max(rem, key=height)
        r = rem[mu]
        b3 = list(product(range(mu[1], mu[0] + 1, 2), range(mu[2], mu[1] + 1, 2), range(mu[3], mu[2] + 1, 2)))
        levels = list(range(mu[3] % 2, mu[3] + 1, 2))
        cap = {t: min(rem.get(x + (t,), 0) for x in b3) for t in levels}
        if cap[mu[3]] < r:
            return None
        C, run = {}, r
        for t in reversed(levels):
            run = min(run, cap[t])
            C[t] = run
        prev = 0
        for t in levels:
            c, prev = C[t] - prev, C[t]
            if c:
                out[mu + (t,)] = out.get(mu + (t,), 0) + c
                for x in b3:
                    for s in range(t, mu[3] + 1, 2):
                        v = rem.get(x + (s,), 0) - c
                        if v < 0:
                            return None
                        if v:
                            rem[x + (s,)] = v
                        else:
                            rem.pop(x + (s,), None)
    return out


def ilp_family(m, keep):
    from sage.all import MixedIntegerLinearProgram
    from sage.numerical.mip import MIPSolverException
    cands = [l for l in candidates(m) if keep(l)]
    cover = {mu: [] for mu in m}
    for j, lam in enumerate(cands):
        for mu in interlacing(lam):
            cover[mu].append(j)
    if any(not v for v in cover.values()):
        return False
    p = MixedIntegerLinearProgram(solver="GLPK")
    x = p.new_variable(integer=True, nonnegative=True)
    for mu, js in cover.items():
        p.add_constraint(p.sum(x[j] for j in js) == m[mu])
    try:
        p.solve()
        return True
    except MIPSolverException:
        return False


def lift3(g, t):
    """Unique T on tops (m1,m2,m3), m3 >= t, with sum_T 1[m1>=n1>=m2>=n2>=m3>=n3>=t] = g (doubled, step 2)."""
    G = lambda a, b, c: g.get((a, b, c), 0)
    keys = set()
    for (a, b, c) in g:
        keys.update({(a, b, c), (a - 2, b, c), (a, b - 2, c), (a - 2, b - 2, c)})
    out = {}
    for (m1, m2, m3) in keys:
        if not m1 >= m2 >= m3 >= t:
            continue
        h = lambda a, b: G(a, b, m3) - (G(a, b, m3 + 2) if m3 + 2 <= b else 0)
        U = lambda a, b: h(a, b) - (h(a, b + 2) if b + 2 <= a else 0)
        v = U(m1, m2) - U(m1 + 2, m2)
        if v:
            out[(m1, m2, m3)] = v
    return out


def res3(Tt, t):
    out = {}
    for (m1, m2, m3), c in Tt.items():
        for nu in product(range(m2, m1 + 1, 2), range(m3, m2 + 1, 2), range(t, m3 + 1, 2)):
            out[nu] = out.get(nu, 0) + c
    return {k: v for k, v in out.items() if v}


def m_pattern(m, t):
    m1, m2, m3 = m
    s, out = t - 2, {}
    for k, c in ((m, 1), ((m1, m2, s), -1), ((m1, m3 - 2, s), 1), ((m2 - 2, m3 - 2, s), -1)):
        if k[0] >= k[1] >= k[2] >= s:
            out[k] = out.get(k, 0) + c
    return out


def level_lift(rho, rule):
    """l5-extremal lift by the level recursion; returns the lift or None."""
    from sage.all import MixedIntegerLinearProgram
    from sage.numerical.mip import MIPSolverException
    E = {}
    for par in (0, 1):
        part = {k: v for k, v in rho.items() if k[0] % 2 == par}
        if not part:
            continue
        tmax, tmin = max(k[3] for k in part), par

        def slab(t):
            return {k[:3]: v for k, v in part.items() if k[3] == t}

        Tt = lift3(slab(tmax), tmax)
        assert res3(Tt, tmax) == slab(tmax)
        for t in range(tmax, tmin - 2, -2):
            if min(Tt.values(), default=0) < 0:
                return None
            D, Tn = {}, None
            if t > tmin:
                cur, nxt = slab(t), slab(t - 2)
                diff = {}
                for k in set(cur) | set(nxt):
                    if k[2] >= t:
                        diff[k] = nxt.get(k, 0) - cur.get(k, 0)
                    elif k[2] == t - 2:
                        diff[k] = nxt.get(k, 0)
                B = lift3({k: v for k, v in diff.items() if v}, t - 2)
                tops = [mm for mm in Tt if Tt[mm] > 0]
                if rule == "all":
                    D = dict(Tt)
                elif rule == "forced":
                    D = {mm: min(c, max(0, -B.get(mm, 0))) for mm, c in Tt.items()}
                elif rule in ("ilp-min", "ilp-max"):
                    p = MixedIntegerLinearProgram(maximization=(rule == "ilp-max"), solver="GLPK")
                    d = p.new_variable(integer=True, nonnegative=True)
                    contrib = {}
                    for i, mm in enumerate(tops):
                        p.add_constraint(d[i] <= Tt[mm])
                        for k, c in m_pattern(mm, t).items():
                            contrib.setdefault(k, []).append((i, c))
                    for k in set(B) | set(contrib):
                        p.add_constraint(B.get(k, 0) + p.sum(c * d[i] for i, c in contrib.get(k, [])) >= 0)
                    p.set_objective(p.sum(d[i] for i in range(len(tops))) if tops else 0)
                    try:
                        p.solve()
                    except MIPSolverException:
                        return None
                    dv = p.get_values(d)
                    D = {mm: int(round(dv[i])) for i, mm in enumerate(tops) if round(dv[i])}
                D = {mm: c for mm, c in D.items() if c}
                Tn = dict(B)
                for mm, c in D.items():
                    for k, v in m_pattern(mm, t).items():
                        Tn[k] = Tn.get(k, 0) + c * v
                Tn = {k: v for k, v in Tn.items() if v}
            for mm, c in Tt.items():
                dd = D.get(mm, 0)
                if dd:
                    E[mm + (t, t)] = E.get(mm + (t, t), 0) + dd
                if c - dd:
                    E[mm + (t, tmin)] = E.get(mm + (t, tmin), 0) + c - dd
            Tt = Tn if Tn is not None else {}
    assert restrict(E) == rho, "level lift does not restrict correctly"
    return E


def diagnostics(rho, stats):
    """Per-level bookkeeping for the 'ilp-min' recursion (see (4) in the docstring)."""
    from sage.all import MixedIntegerLinearProgram
    for par in (0, 1):
        part = {k: v for k, v in rho.items() if k[0] % 2 == par}
        if not part:
            continue
        tmax, tmin = max(k[3] for k in part), par

        def slab(t):
            return {k[:3]: v for k, v in part.items() if k[3] == t}

        Tt = lift3(slab(tmax), tmax)
        for t in range(tmax, tmin, -2):
            cur, nxt = slab(t), slab(t - 2)
            diff = {}
            for k in set(cur) | set(nxt):
                if k[2] >= t:
                    diff[k] = nxt.get(k, 0) - cur.get(k, 0)
                elif k[2] == t - 2:
                    diff[k] = nxt.get(k, 0)
            B = lift3({k: v for k, v in diff.items() if v}, t - 2)
            forced = {m: max(0, -B.get(m, 0)) for m in Tt if B.get(m, 0) < 0}
            stats["level transitions"] = stats.get("level transitions", 0) + 1
            if any(c > Tt.get(m, 0) for m, c in forced.items()):
                stats["forced degenerate exceeds available tops"] = stats.get("forced degenerate exceeds available tops", 0) + 1
            if any(v < 0 and k[2] >= t and not Tt.get(k) for k, v in B.items()):
                stats["non-boundary negative with no top"] = stats.get("non-boundary negative with no top", 0) + 1
            Tf = dict(B)
            for m, c in forced.items():
                for k, v in m_pattern(m, t).items():
                    Tf[k] = Tf.get(k, 0) + c * v
            negb = {k: v for k, v in Tf.items() if v < 0}
            for k, v in negb.items():
                key = (f"boundary negative after forced: {'spinor' if par else 'tensor'}, "
                       f"mu4 {t / 2:g} -> {(t - 2) / 2:g}, k2 - (t-1) = {(k[1] - t + 2) / 2:g}")
                stats[key] = stats.get(key, 0) - v
            tops = [m for m in Tt if Tt[m] > 0]
            p = MixedIntegerLinearProgram(maximization=False, solver="GLPK")
            d = p.new_variable(integer=True, nonnegative=True)
            contrib = {}
            for i, m in enumerate(tops):
                p.add_constraint(d[i] <= Tt[m])
                p.add_constraint(d[i] >= forced.get(m, 0))
                for k, c in m_pattern(m, t).items():
                    contrib.setdefault(k, []).append((i, c))
            for k in set(B) | set(contrib):
                p.add_constraint(B.get(k, 0) + p.sum(c * d[i] for i, c in contrib.get(k, [])) >= 0)
            p.set_objective(p.sum(d[i] for i in range(len(tops))) if tops else 0)
            p.solve()
            dv = p.get_values(d)
            D = {m: int(round(dv[i])) for i, m in enumerate(tops) if round(dv[i])}
            for m, c in D.items():
                extra = c - forced.get(m, 0)
                if extra:
                    direct = (m[0], m[2] - 2, t - 2) in negb
                    key = "extra degenerate repairing directly via +[m1,m3-1]" if direct else "extra degenerate, indirect"
                    stats[key] = stats.get(key, 0) + extra
            Tn = dict(B)
            for m, c in D.items():
                for k, v in m_pattern(m, t).items():
                    Tn[k] = Tn.get(k, 0) + c * v
            Tt = {k: v for k, v in Tn.items() if v}


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else N
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    ext5 = lambda l: l[4] in (l[3] % 2, l[3])
    ext3 = lambda l: l[2] in (l[1], l[3])
    tests = {
        "(1) greedy l5 = mu4": lambda m: greedy(m, "max")[0] is not None,
        "(1) greedy l5 = floor": lambda m: greedy(m, "min")[0] is not None,
        "(1) adaptive (deepest profile)": lambda m: adaptive(m) is not None,
        "(2) any lift (ILP)": lambda m: solve_ilp(m) is not None,
        "(2) l5 extremal (ILP)": lambda m: ilp_family(m, ext5),
        "(2) l3 extremal (ILP)": lambda m: ilp_family(m, ext3),
        "(3) levels: no degenerate tops": lambda m: level_lift(m, "none") is not None,
        "(3) levels: all degenerate": lambda m: level_lift(m, "all") is not None,
        "(3) levels: degenerate where forced": lambda m: level_lift(m, "forced") is not None,
        "(3) levels: per-level ILP, min": lambda m: level_lift(m, "ilp-min") is not None,
        "(3) levels: per-level ILP, max": lambda m: level_lift(m, "ilp-max") is not None,
    }
    for target in ("rho_n, max-height induction", "rho_n, min-height induction", "cumulative R9_n"):
        cum = target.startswith("cum")
        top = 24 if cum else NMAX
        fails = {k: [] for k in tests}
        A = {}
        for n in range(1, top + 1):
            tgt = R9[n] if cum else residual(r9[n], A)
            if tgt:
                for k, f in tests.items():
                    if not f(tgt):
                        fails[k].append(n)
            if not cum:
                E = solve_ilp(tgt, lambda l: height(l[:4]), maximize=("max" in target)) if tgt else {}
                A = T(A)
                for l, c in E.items():
                    A[l] = A.get(l, 0) + c
        log(f"{target} (n <= {top}):")
        for k, f in fails.items():
            log(f"  {k:40s}: fails at {f if len(f) < 9 else f[:8] + ['...']} ({len(f)} levels)")
        log(f"  [{time.time() - t0:.0f}s]")
    log(f"(4) per-level choices along the max-height induction (n <= {NMAX}):")
    stats, A = {}, {}
    for n in range(1, NMAX + 1):
        rho = residual(r9[n], A)
        if rho:
            diagnostics(rho, stats)
        E = solve_ilp(rho, lambda l: height(l[:4]), maximize=True) if rho else {}
        A = T(A)
        for l, c in E.items():
            A[l] = A.get(l, 0) + c
    for k in sorted(stats):
        log(f"  {k}: {stats[k]}")
    log(f"  [{time.time() - t0:.0f}s]")
    log("(4) the same for the cumulative R9_n (n <= 24):")
    stats = {}
    for n in range(1, 25):
        diagnostics(R9[n], stats)
    for k in sorted(stats):
        log(f"  {k}: {stats[k]}")
    log(f"  [{time.time() - t0:.0f}s]")
    box = {(4, 4, 4, 2): 1, (4, 4, 4, 4): 1}
    log(f"(2) the box of (2,2,2,2,1) = {{(2,2,2,1),(2,2,2,2)}}: lift exists {solve_ilp(box) is not None}, "
        f"l5-extremal lift exists {ilp_family(box, ext5)}")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "half_b.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
