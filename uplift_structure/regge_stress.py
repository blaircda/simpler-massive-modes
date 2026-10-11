"""
S5 step 1: stress tests of the Regge continuation conjecture.

Notation (class level, doubled ambient weights): T(V_lambda) = V_{lambda+e1};
top(A) = mu1 = l1 slices of Res A;  rho_n(A) = r9_n - S(top A), where S shifts
mu1 -> mu1 + 1.  Identity used in (B):
    top(A) = Res A - Res T^{-1}(A),
so that, for an uplift A of R9_{n-1} and nu = mu - e1,
    (a) at mu  <=>  Res T^{-1}(A)(nu) >= R9_{n-1}(nu) - r9_n(mu) =: deficit(nu).

Modes
  induction  : run A_n = T(A_{n-1}) + E_n from A_1 = 1 with three rules for
               E_n (max-height / min-height / feasible) up to the data's max level.
  universal-a: rigorous check of (a) for ALL uplifts A of R9_{n-1}: for every nu
               with deficit(nu) > 0, an LP lower bound on min_A Res T^{-1}(A)(nu)
               (LP min <= integer min, so LP >= deficit proves it).
  sample-b   : (b) on sampled uplifts A of R9_{n-1}: random-objective ILP vertices
               plus random walks along the Euler moves h_nu of S3b/S4; for each
               sample check rho_n >= 0 and that rho_n has a non-negative uplift.

Data: r9 from --r9 (default uplift_structure/data/r9_to32.json if present, else
so10_uplift/data/r9.json); levels <= 24 are cross-checked against
so10_uplift/data/r9.json when both are available.

Output: results/regge_stress_<mode>.txt
"""

import argparse
import json
import os
import random
import sys
import time
from itertools import product

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, HERE, RESULTS

from uplift import candidates, height, interlacing, restrict, solve_ilp

DEFAULT_EXT = os.path.join(HERE, "data", "r9_to32.json")


# ------------------------------------------------------------------ data

def load_R9(path):
    levels = json.load(open(path))["levels"]
    N = max(int(k) for k in levels)
    r9, R9 = {}, {0: {}}
    for n in range(1, N + 1):
        r9[n] = {tuple(e["ambient2"]): e["mult"] for e in levels[str(n)]["r9_B4"]}
        R = dict(R9[n - 1])
        for mu, k in r9[n].items():
            R[mu] = R.get(mu, 0) + k
        R9[n] = R
    return N, r9, R9


def T(a):
    return {(l[0] + 2,) + l[1:]: c for l, c in a.items()}


def residual(r9n, A):
    """rho = r9_n - S(top A) = R9_n - Res T(A), computed as r9_n + R9_{n-1} - Res T(A)."""
    out = dict(r9n)
    ResA = restrict(A)
    for mu, v in ResA.items():
        out[mu] = out.get(mu, 0) + v
    for mu, v in restrict(T(A)).items():
        out[mu] = out.get(mu, 0) - v
    return {k: v for k, v in out.items() if v}


# ------------------------------------------------------------------ Euler moves (S3b/S4)

RHO2 = (8, 6, 4, 2, 0)
E_WEIGHTS = []
for eps in product((-1, 0, 1), repeat=5):
    c = 1
    for t in eps:
        c *= 2 if t == 0 else -1
    E_WEIGHTS.append((tuple(2 * t for t in eps), c))


def _dominate_signed(v):
    a = [abs(t) for t in v]
    if len(set(a)) < 5:
        return None, 0
    order = sorted(range(5), key=lambda i: -a[i])
    sgn, seen = 1, [False] * 5
    for i in range(5):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = order[j]
                L += 1
            if L % 2 == 0:
                sgn = -sgn
    out = [a[i] for i in order]
    if out[4] and sum(1 for t in v if t < 0) % 2:
        out[4] = -out[4]
    return tuple(out), sgn


def e_times(mu):
    """Resolved e (x) V_mu by Klimyk (same as bbms_lift/positive_graded.py, validated there)."""
    out = {}
    for w, c in E_WEIGHTS:
        d, s = _dominate_signed(tuple(mu[i] + RHO2[i] + w[i] for i in range(5)))
        if d is None:
            continue
        lam = tuple(d[i] - RHO2[i] for i in range(5))
        out[lam] = out.get(lam, 0) + s * c
    return {k: v for k, v in out.items() if v}


def h_move(nu):
    """Class-level Euler move h_nu (S3b): pi(e V_nu)/2 if l5 = 0, else pi(e V_nu)."""
    merged = {}
    for lam, c in e_times(nu).items():
        key = lam[:4] + (abs(lam[4]),)
        merged[key] = merged.get(key, 0) + c
    merged = {k: v for k, v in merged.items() if v}
    if nu[4] == 0:
        assert all(v % 2 == 0 for v in merged.values())
        merged = {k: v // 2 for k, v in merged.items()}
    return merged


# ------------------------------------------------------------------ modes

def mode_induction(N, r9, R9, log):
    rules = {"max-height": (lambda l: height(l[:4]), True),
             "min-height": (lambda l: height(l[:4]), False),
             "feasible": (None, False)}
    for name, (obj, mx) in rules.items():
        A, t0 = {}, time.time()
        for n in range(1, N + 1):
            rho = residual(r9[n], A)
            if min(rho.values(), default=0) < 0:
                log(f"  [{name}] n={n}: rho_n NEGATIVE -> conjecture (a) fails on this path")
                return False
            E = solve_ilp(rho, obj, maximize=mx) if rho else {}
            if E is None:
                log(f"  [{name}] n={n}: rho_n has NO non-negative uplift -> conjecture (b) fails on this path")
                return False
            A = T(A)
            for l, c in E.items():
                A[l] = A.get(l, 0) + c
            assert restrict(A) == R9[n]
            if n >= 20 or n == N:
                log(f"  [{name}] n={n:2d}: ok, |E_n| = {sum(E.values())} [{time.time() - t0:.0f}s]")
    return True


def mode_universal_a(N, r9, R9, log, nmax):
    from sage.all import MixedIntegerLinearProgram
    allok = True
    for n in range(2, min(N, nmax) + 1):
        t0 = time.time()
        m = R9[n - 1]
        deficit = {}
        for nu, v in m.items():
            mu = (nu[0] + 2,) + nu[1:]
            d = v - r9[n].get(mu, 0)
            if d > 0:
                deficit[nu] = d
        if not deficit:
            log(f"  n={n:2d}: no weight with positive deficit -> (a) holds for every uplift (trivially)")
            continue
        cands = candidates(m)
        cover = {mu: [] for mu in m}
        for j, lam in enumerate(cands):
            for mu in interlacing(lam):
                cover[mu].append(j)
        p = MixedIntegerLinearProgram(maximization=False, solver="GLPK")
        x = p.new_variable(nonnegative=True)
        for mu, js in cover.items():
            p.add_constraint(p.sum(x[j] for j in js) == m[mu])
        p.set_objective(p.sum(x[j] for j in range(len(cands))))  # placeholder
        b = p.get_backend()
        cols = [next(iter(x[j].dict())) for j in range(len(cands))]
        ncols = b.ncols()
        # Res T^{-1}(A)(nu) = sum of A(lambda) over lambda with nu in box(lambda), l1 > nu1
        below = {}
        for j, lam in enumerate(cands):
            for mu in interlacing(lam):
                if lam[0] > mu[0] and mu in deficit:
                    below.setdefault(mu, []).append(j)
        worst = None
        for nu, d in deficit.items():
            coeffs = [0] * ncols
            for j in below.get(nu, []):
                coeffs[cols[j]] = 1
            b.set_objective(coeffs)
            b.solve()
            val = b.get_objective_value()
            slack = val - d
            if worst is None or slack < worst[0]:
                worst = (slack, nu, d, val)
            if val < d - 1e-7:
                allok = False
                log(f"  n={n:2d}: (a) NOT certified at mu = nu + e1, nu = {nu}: "
                    f"LP min Res T^-1 A = {val:.3f} < deficit {d}")
        log(f"  n={n:2d}: {len(deficit)} weights with positive deficit; "
            f"(a) certified for ALL uplifts: {worst[0] >= -1e-7}; tightest slack {worst[0]:.3f} "
            f"at nu={worst[1]} (deficit {worst[2]}) [{time.time() - t0:.0f}s]")
    return allok


def sample_uplifts(m, k_vertices, walk_steps, rng):
    """Diverse non-negative uplifts of m: random-objective ILP vertices + Euler-move walks."""
    cands = candidates(m)
    cset = set(cands)
    moves = []
    for nu in cands:
        h = h_move(nu)
        if h and set(h) <= cset:
            moves.append(h)
    samples = []
    for _ in range(k_vertices):
        w = {lam: rng.uniform(-1, 1) for lam in cands}
        A = solve_ilp(m, lambda l: w[l], maximize=True)
        assert A is not None
        samples.append(A)
        cur = dict(A)
        for step in range(walk_steps):
            h = rng.choice(moves)
            s = rng.choice((1, -1))
            new = dict(cur)
            for lam, v in h.items():
                new[lam] = new.get(lam, 0) + s * v
            if min(new.values()) >= 0:
                cur = {k: v for k, v in new.items() if v}
                if step % 5 == 4:
                    samples.append(dict(cur))
    for A in samples:
        assert restrict(A) == m
    return samples, len(moves)


def mode_sample_b(N, r9, R9, log, levels, k_vertices, walk_steps, seed):
    rng = random.Random(seed)
    allok = True
    for n in levels:
        if n > N:
            break
        t0 = time.time()
        samples, nmoves = sample_uplifts(R9[n - 1], k_vertices, walk_steps, rng)
        distinct = {tuple(sorted(A.items())) for A in samples}
        neg = unl = 0
        cache = {}
        for key in distinct:
            rho = residual(r9[n], dict(key))
            if min(rho.values()) < 0:
                neg += 1
                continue
            rk = tuple(sorted(rho.items()))
            if rk not in cache:
                cache[rk] = solve_ilp(dict(rk)) is not None
            unl += not cache[rk]
        ok = neg == 0 and unl == 0
        allok &= ok
        log(f"  n={n:2d}: {len(distinct)} distinct sampled uplifts of R9_{n - 1} ({nmoves} Euler moves); "
            f"rho negative {neg}; not liftable {unl}; distinct residuals {len(cache)} "
            f"[{time.time() - t0:.0f}s]")
    return allok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["induction", "universal-a", "sample-b"])
    ap.add_argument("--r9", default=DEFAULT_EXT if os.path.exists(DEFAULT_EXT)
                    else os.path.join(DATA, "r9.json"))
    ap.add_argument("--nmax", type=int, default=24, help="universal-a: max level")
    ap.add_argument("--levels", default="12-18", help="sample-b: level range a-b")
    ap.add_argument("--vertices", type=int, default=8)
    ap.add_argument("--walk", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    N, r9, R9 = load_R9(args.r9)
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    base = os.path.join(DATA, "r9.json")
    if os.path.abspath(args.r9) != os.path.abspath(base) and os.path.exists(base):
        Nb, r9b, _ = load_R9(base)
        same = all(r9b[n] == r9[n] for n in range(1, min(N, Nb) + 1))
        log(f"data: {args.r9} (max level {N}); identical to so10_uplift r9 for n <= {min(N, Nb)}: {same}")
        assert same
    else:
        log(f"data: {args.r9} (max level {N})")

    if args.mode == "induction":
        ok = mode_induction(N, r9, R9, log)
    elif args.mode == "universal-a":
        ok = mode_universal_a(N, r9, R9, log, args.nmax)
    else:
        a, b = (int(t) for t in args.levels.split("-"))
        ok = mode_sample_b(N, r9, R9, log, range(a, b + 1), args.vertices, args.walk, args.seed)
    log(f"RESULT: {'no counterexample' if ok else 'COUNTEREXAMPLE / not certified'}")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, f"regge_stress_{args.mode}.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
