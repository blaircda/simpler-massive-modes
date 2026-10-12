"""
S8 (v): which Cartan shifts are "Regge directions" of the trajectory heads?

A direction d = (s, k) is a Cartan shift lambda -> lambda + s of SO(10) highest
weights together with a level offset k.  T itself is (e1, 1) on the cumulative
A_n; here we ask about directions acting on the heads E_n of the Regge
induction A_n = T(A_{n-1}) + E_n.  No relation is guessed; the script reports

  (1) uniqueness: all non-negative lifts of rho_n along the induction, n <= 8
      (if unique, the low-level heads are rule-free facts);
  (2) a scan: for every dominant shift s with l1 <= 2 (38 shifts; l1 <= 4 in doubled
      coordinates) and k = 1..6, does E_n contain s(E_{n-k}) for all
      k+2 <= n <= 8 (the unique heads)?  Reported with the number of irreps
      actually checked, whether the vacuum step (n = k+1) passes, and how far
      the containment persists under the plain max-/min-height rules (n <= 14);
  (3) forced single-direction inductions to n = NMAX for the scan survivors:
      E_n = F_n + P_n, F_n = monomials in d applied to earlier primaries P_m,
      P_n >= 0 a lift of rho_n - Res F_n (three rules);
  (4) a null test: force a direction only on primaries at levels m >= m0
      (lower ones exempt), m0 = 2, 5, 7, 9.  A spurious direction should fail
      shortly after forcing starts, wherever it starts;
  (5) joint tests for T2 = (e1, 2) and Tw = (e1+e2, 3): a free commutative
      action (overshoot located), and the weaker joint containment
      E_n >= max(T2 E_{n-2}, Tw E_{n-3}) under three greedy rules.

Directions are (s, k, fermionic, vacuum_allowed); fermionic ones act at most
once per monomial.  Class level throughout (l5 >= 0 representatives).

Output: results/regge_directions.txt
"""

import itertools
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, SO10
from regge_stress import DEFAULT_EXT, T, load_R9, residual

sys.path.insert(0, SO10)
from uplift import enumerate_uplifts, height, restrict, solve_ilp
from weights import dim_irrep, dynkin_B4, dynkin_D5

N, r9, R9 = load_R9(DEFAULT_EXT)
RULES = ("max-height", "min-height", "feasible")
E1, W = (2, 0, 0, 0, 0), (2, 2, 0, 0, 0)


def objective(rule):
    return {"max-height": (lambda l: height(l[:4]), True),
            "min-height": (lambda l: height(l[:4]), False),
            "feasible": (None, False)}[rule]


def fold(l):
    return l[:4] + (abs(l[4]),)


def shifted(d, s):
    out = {}
    for l, c in d.items():
        t = fold(tuple(a + b for a, b in zip(l, s)))
        out[t] = out.get(t, 0) + c
    return out


def fmt(d):
    return ", ".join(f"{c}x{dim_irrep(l, 'D', 5)}{list(dynkin_D5(l))}" for l, c in sorted(d.items()) if c)


def plain_heads(rule, NM):
    obj, mx = objective(rule)
    A, heads = {}, {}
    for n in range(1, NM + 1):
        rho = residual(r9[n], A)
        E = solve_ilp(rho, obj, maximize=mx) if rho else {}
        heads[n] = E
        A = T(A)
        for l, c in E.items():
            A[l] = A.get(l, 0) + c
    return heads


def monomials(dirs, Lmax, vacuum):
    """M[L] = {total shift: count} over monomials of level L in the directions."""
    M = [dict() for _ in range(Lmax + 1)]
    M[0] = {(0,) * 5: 1}
    for (s, k, ferm, vok) in dirs:
        if vacuum and not vok:
            continue
        new = [dict(m) for m in M]
        for L in range(Lmax + 1):
            for sv, c in M[L].items():
                cur, LL = sv, L
                while LL + k <= Lmax:
                    LL += k
                    cur = tuple(a + b for a, b in zip(cur, s))
                    new[LL][cur] = new[LL].get(cur, 0) + c
                    if ferm:
                        break
        M = new
    return M


def forced_run(dirs, rule, m0=2, NM=None, keep=False):
    """Free (commutative) forced induction; primaries at 1 < m < m0 are exempt."""
    NM = NM or N
    obj, mx = objective(rule)
    Mall, Mvac = monomials(dirs, NM, False), monomials(dirs, NM, True)
    A, Pm, Eh = {}, {}, {}
    for n in range(1, NM + 1):
        rho = residual(r9[n], A)
        F = {}
        for m, Pd in Pm.items():
            if m == 1:
                Mm = Mvac
            elif m < m0:
                continue
            else:
                Mm = Mall
            for sv, c in Mm[n - m].items():
                for l, v in shifted(Pd, sv).items():
                    F[l] = F.get(l, 0) + c * v
        rest = dict(rho)
        for mu, v in restrict(F).items():
            rest[mu] = rest.get(mu, 0) - v
        neg = {mu: v for mu, v in rest.items() if v < 0}
        if neg:
            return f"STUCK at n={n} (forced copies overshoot)", (n, F, neg, Pm, Eh)
        rest = {mu: v for mu, v in rest.items() if v}
        P = solve_ilp(rest, obj, maximize=mx) if rest else {}
        if P is None:
            return f"STUCK at n={n} (remainder not liftable)", (n, F, None, Pm, Eh)
        Pm[n] = P
        E = dict(F)
        for l, c in P.items():
            E[l] = E.get(l, 0) + c
        Eh[n] = E
        A = T(A)
        for l, c in E.items():
            A[l] = A.get(l, 0) + c
        assert restrict(A) == R9[n]
    return f"ok to n={NM}", (NM, None, None, Pm, Eh)


def max_run(dirs, rule, NM=None):
    """Joint containment E_n >= max_d shift_d(E_{n-k_d}) (images may overlap)."""
    NM = NM or N
    obj, mx = objective(rule)
    A, Eh = {}, {}
    for n in range(1, NM + 1):
        rho = residual(r9[n], A)
        F = {}
        for (s, k, ferm, vok) in dirs:
            src = n - k
            if src < 1 or (src == 1 and not vok):
                continue
            for t, v in shifted(Eh[src], s).items():
                F[t] = max(F.get(t, 0), v)
        rest = dict(rho)
        for mu, v in restrict(F).items():
            rest[mu] = rest.get(mu, 0) - v
        if min(rest.values(), default=0) < 0:
            return f"STUCK at n={n} (forced copies overshoot)"
        rest = {mu: v for mu, v in rest.items() if v}
        P = solve_ilp(rest, obj, maximize=mx) if rest else {}
        if P is None:
            return f"STUCK at n={n} (remainder not liftable)"
        E = dict(F)
        for l, c in P.items():
            E[l] = E.get(l, 0) + c
        Eh[n] = E
        A = T(A)
        for l, c in E.items():
            A[l] = A.get(l, 0) + c
        assert restrict(A) == R9[n]
    return f"ok to n={NM}"


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else N
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    # (1) uniqueness of the heads along the induction
    log("(1) all non-negative lifts of rho_n along the Regge induction")
    A, U = {}, {}
    for n in range(1, 9):
        rho = residual(r9[n], A)
        Es = list(enumerate_uplifts(rho)) if rho else [{}]
        log(f"  n={n}: {len(Es)} head(s)")
        assert len(Es) == 1, "heads not unique; the scan below assumes uniqueness"
        U[n] = Es[0]
        A = T(A)
        for l, c in U[n].items():
            A[l] = A.get(l, 0) + c
    log(f"  [{time.time() - t0:.0f}s]")

    # (2) scan of Cartan shifts against the unique heads
    H = {rule: plain_heads(rule, 14) for rule in ("max-height", "min-height")}
    for rule in H:
        assert all(H[rule][n] == U[n] for n in range(1, 9))
    shifts = set()
    for l in itertools.product(range(0, 5), repeat=4):
        if not l[0] >= l[1] >= l[2] >= l[3]:
            continue
        for l5 in range(-l[3], l[3] + 1):
            w = l + (l5,)
            if len({x % 2 for x in w}) == 1 and any(w):
                shifts.add(w)
    shifts = sorted(shifts, key=lambda w: (sum(abs(x) for x in w), w))
    log(f"(2) scan: {len(shifts)} shifts x k=1..6; listed if E_n >= s(E_(n-k)) for all k+2 <= n <= 8 "
        "with at least one irrep checked")
    survivors = []
    for s in shifts:
        for k in range(1, 7):
            checked, ok = 0, True
            for n in range(k + 2, 9):
                f = shifted(U[n - k], s)
                checked += sum(f.values())
                if any(U[n].get(l, 0) < c for l, c in f.items()):
                    ok = False
                    break
            if not ok or checked == 0:
                continue
            vac = all(U[k + 1].get(l, 0) >= c for l, c in shifted(U[1], s).items()) if k + 1 <= 8 else None
            persist = {}
            for rule in H:
                persist[rule] = next((n for n in range(k + 2, 15)
                                      if any(H[rule][n].get(l, 0) < c
                                             for l, c in shifted(H[rule][n - k], s).items())), None)
            survivors.append((s, k, vac))
            log(f"  s={list(dynkin_D5(s))} k={k}: irreps checked {checked}; vacuum step passes: {vac}; "
                f"plain rules first fail at n = {persist['max-height']} (max), {persist['min-height']} (min)")

    # (3) forced single-direction inductions
    log(f"(3) forced single-direction inductions to n={NMAX} (vacuum treated as in the scan)")
    for s, k, vac in survivors:
        variants = [False, True] if s[0] % 2 else [False]
        for ferm in variants:
            out = [forced_run([(s, k, ferm, bool(vac))], rule, NM=NMAX)[0] for rule in RULES]
            log(f"  s={list(dynkin_D5(s))} k={k}{' fermionic' if ferm else ''}: "
                + " | ".join(f"{r}: {o}" for r, o in zip(RULES, out)))
    log(f"  [{time.time() - t0:.0f}s]")

    # (4) null test: delayed start
    log("(4) null test: direction forced only on primaries at levels m >= m0")
    panel = [(E1, 2), (W, 3), (E1, 3), (E1, 1), (W, 2), ((1, 1, 1, 1, 1), 2), ((2, 2, 2, 0, 0), 3),
             ((4, 0, 0, 0, 0), 3), ((0, 0, 0, 0, 0), 2)]
    for s, k in panel:
        res = []
        for m0 in (2, 5, 7, 9):
            o = [forced_run([(s, k, False, False)], rule, m0=m0, NM=NMAX)[0] for rule in ("max-height", "feasible")]
            res.append(f"m0={m0}: {o[0]}" + ("" if o[0] == o[1] else f" / feasible: {o[1]}"))
        log(f"  s={list(dynkin_D5(s))} k={k}: " + "; ".join(res))
    log(f"  [{time.time() - t0:.0f}s]")

    # (5) joint T2 + Tw
    log("(5) joint T2 = (e1, 2) [vacuum exempt] and Tw = (e1+e2, 3)")
    for vw in (True, False):
        dirs = [(E1, 2, False, False), (W, 3, False, vw)]
        for rule in RULES:
            msg, (n, F, neg, Pm, Eh) = forced_run(dirs, rule, NM=NMAX)
            log(f"  free action, Tw on vacuum {vw} [{rule}]: {msg}")
        if neg:
            log(f"    at n={n}: forced part F = {fmt(F)}")
            log(f"    overshoot at SO(9) weights {[(list(dynkin_B4(mu)), v) for mu, v in neg.items()]}")
            for m in sorted(Pm):
                log(f"    primaries P_{m} = {fmt(Pm[m]) or '-'}")
    for rule in RULES:
        log(f"  joint containment max(T2, Tw) [{rule}]: {max_run([(E1, 2, False, False), (W, 3, False, True)], rule, NM=NMAX)}")
    log(f"[{time.time() - t0:.0f}s]")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "regge_directions.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
