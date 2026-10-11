"""
S5: Regge induction.  T(V_lambda) = V_{lambda + e1} (Regge shift along l1).

Res T(A) = Res A + S(top A), where top(A) = sum_lambda A(lambda) Res V*_lambda,
V*_lambda = V_lambda - V_{lambda - e1}, is the mu1 = l1 slice of Res A, and S
shifts mu1 -> mu1 + 1.  The Regge induction is

    A_n = T(A_{n-1}) + E_n,   E_n >= 0 a non-negative uplift of
    rho_n = R^9_n - Res T(A_{n-1}) = r^9_n - S(top A_{n-1}).

Checks (class level, reading so10_uplift/data only):
 (1) canonical uplifts: is A_n - T(A_{n-1}) >= 0?  is rho_n >= 0?   (n <= 24)
 (2) run the induction from A_1 = 1 with three rules for E_n
     (max-height, min-height, any feasible), n <= 24;
 (3) path independence: for EVERY enumerated uplift A_{n-1} (n-1 <= 10),
     is rho_n >= 0 and liftable?

Output: results/regge_induction.txt
"""

import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RESULTS

from uplift import height, restrict, solve_ilp
from weights import dynkin_D5

L = json.load(open(os.path.join(DATA, "R10.json")))["levels"]
NMAX = 24


def R9(n):
    return {tuple(e["ambient2"]): e["mult"] for e in L[str(n)]["R9_B4"]}


def canonical(n):
    return {tuple(e["ambient2"]): e["mult"] for e in L[str(n)]["canonical_R10_D5"]}


def T(a):
    return {(l[0] + 2,) + l[1:]: c for l, c in a.items()}


def residual(n, A_prev):
    rho = R9(n)
    for mu, v in restrict(T(A_prev)).items():
        rho[mu] = rho.get(mu, 0) - v
    return rho


def main():
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    log("(1) canonical uplifts")
    prev = {}
    for n in range(1, NMAX + 1):
        A = canonical(n)
        TA = T(prev)
        E = {l: A.get(l, 0) - TA.get(l, 0) for l in set(A) | set(TA)}
        nneg = sum(1 for v in E.values() if v < 0)
        rho_ok = min(residual(n, prev).values()) >= 0
        log(f"  n={n:2d}: A_n - T(A_n-1) >= 0: {nneg == 0}"
            + (f" ({nneg} negative entries)" if nneg else "") + f"; rho_n >= 0: {rho_ok}")
        prev = A

    log("(2) Regge induction from A_1 = 1")
    rules = {"max-height": (lambda l: height(l[:4]), True),
             "min-height": (lambda l: height(l[:4]), False),
             "feasible": (None, False)}
    for name, (obj, mx) in rules.items():
        A, t0, status = {}, time.time(), "ok to n=%d" % NMAX
        for n in range(1, NMAX + 1):
            rho = residual(n, A)
            if min(rho.values(), default=0) < 0:
                status = f"stuck at n={n}: rho_n negative"
                break
            rho = {k: v for k, v in rho.items() if v}
            E = solve_ilp(rho, obj, maximize=mx) if rho else {}
            if E is None:
                status = f"stuck at n={n}: rho_n not liftable"
                break
            A = T(A)
            for l, c in E.items():
                A[l] = A.get(l, 0) + c
            assert restrict(A) == R9(n)
        log(f"  rule {name}: {status} [{time.time() - t0:.0f}s]")

    log("(3) path independence: every enumerated uplift A_{n-1}")
    AL = json.load(open(os.path.join(DATA, "R10_all_uplifts.json")))["levels"]
    for m in sorted(map(int, AL)):
        n = m + 1
        cl = [tuple(c) for c in AL[str(m)]["classes_ambient2"]]
        neg = unl = 0
        cache = {}
        for u in AL[str(m)]["uplifts"]:
            A = {cl[i]: k for i, k in enumerate(u) if k}
            rho = residual(n, A)
            if min(rho.values()) < 0:
                neg += 1
                continue
            key = tuple(sorted((k, v) for k, v in rho.items() if v))
            if key not in cache:
                cache[key] = solve_ilp(dict(key)) is not None
            unl += not cache[key]
        log(f"  A_{m} -> level {n}: {len(AL[str(m)]['uplifts'])} uplifts; rho negative for {neg}; "
            f"not liftable for {unl}; distinct residuals {len(cache)}")

    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "regge_induction.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
