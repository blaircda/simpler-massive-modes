"""
SO(9) -> SO(10) uplift machinery for the cumulative Bars III.B scheme.

Branching Spin(10) -> Spin(9) is multiplicity free:
    V_lambda|_{B4} = sum_{mu interlacing lambda} V_mu,
    l1 >= mu1 >= l2 >= mu2 >= l3 >= mu3 >= l4 >= mu4 >= |l5|,
with mu_i = l_i mod 1.  All weights here are doubled ambient tuples.

lambda and its conjugate (l5 -> -l5) restrict identically, so the uplift
problem is posed on conjugation classes, represented by l5 >= 0; a class
with l5 > 0 and total multiplicity c may be split c = c_+ + c_- freely
between the two chiralities.
"""

from itertools import product


def height(w):
    return 10 * w[0] + 3 * w[1] + 2 * w[2] + w[3]


def interlacing(lam):
    """All B4 highest weights mu (doubled) in the restriction of D5 lam."""
    l1, l2, l3, l4, l5 = lam
    a5 = abs(l5)
    rngs = [range(l2, l1 + 1, 2), range(l3, l2 + 1, 2), range(l4, l3 + 1, 2),
            range(a5, l4 + 1, 2)]
    return [tuple(m) for m in product(*rngs)]


def restrict(c):
    """Restriction of a D5 character {lam: mult} to B4 {mu: mult}."""
    out = {}
    for lam, k in c.items():
        for mu in interlacing(lam):
            out[mu] = out.get(mu, 0) + k
    return {mu: k for mu, k in out.items() if k}


def candidates(m):
    """
    Conjugation classes lambda = (L, l5), l5 >= 0, whose whole restriction
    lies inside supp(m).  L ranges over supp(m) since the restriction of
    V_lambda contains V_L.
    """
    supp = set(m)
    out = []
    for L in supp:
        for l5 in range(L[3], -1, -2):
            lam = L + (l5,)
            if all(mu in supp for mu in interlacing(lam)):
                out.append(lam)
    return out


def greedy(m, rule):
    """
    Top-down greedy lift.  Repeatedly take the height-maximal mu in the
    remainder (it is componentwise maximal, so it must be the top component
    (l1..l4) of whatever lifts it) and lift it with
        rule = "max": l5 = mu4   (smallest restriction)
        rule = "min": l5 = 0 or 1/2 (largest restriction).
    Returns (dict or None, failure info).
    """
    rem = dict(m)
    out = {}
    while rem:
        mu = max(rem, key=height)
        k = rem[mu]
        l5 = mu[3] if rule == "max" else mu[3] % 2
        lam = mu + (l5,)
        for nu in interlacing(lam):
            v = rem.get(nu, 0) - k
            if v < 0:
                return None, {"lambda": lam, "mult": k, "missing": nu, "have": rem.get(nu, 0)}
            if v:
                rem[nu] = v
            else:
                rem.pop(nu, None)
        out[lam] = out.get(lam, 0) + k
    return out, None


def solve_ilp(m, objective=None, maximize=False, solver=None):
    """
    Exact non-negative integer uplift via Sage MILP.
    objective: optional function lam -> weight (default: none, feasibility).
    Returns dict {lam: mult} or None if infeasible.
    """
    from sage.all import MixedIntegerLinearProgram
    from sage.numerical.mip import MIPSolverException

    cands = candidates(m)
    cover = {mu: [] for mu in m}
    for j, lam in enumerate(cands):
        for mu in interlacing(lam):
            cover[mu].append(j)
    if any(not v for v in cover.values()):
        return None
    p = MixedIntegerLinearProgram(maximization=maximize, solver=solver)
    x = p.new_variable(integer=True, nonnegative=True)
    for mu, js in cover.items():
        p.add_constraint(p.sum(x[j] for j in js) == m[mu])
    if objective is not None:
        p.set_objective(p.sum(objective(lam) * x[j] for j, lam in enumerate(cands)))
    try:
        p.solve()
    except MIPSolverException:
        return None
    vals = p.get_values(x)
    sol = {cands[j]: int(round(v)) for j, v in vals.items() if round(v)}
    assert restrict(sol) == m, "ILP solution does not restrict correctly"
    return sol


def kernel_dimension(m):
    """dim of {c in Q^cands : restrict(c) = 0}; the affine dimension of the solution set."""
    from sage.all import matrix, QQ
    cands = candidates(m)
    mus = sorted(m)
    idx = {mu: i for i, mu in enumerate(mus)}
    entries = {}
    for j, lam in enumerate(cands):
        for mu in interlacing(lam):
            entries[(idx[mu], j)] = 1
    A = matrix(QQ, len(mus), len(cands), entries, sparse=True)
    return len(cands) - A.rank(), len(cands)


# ------------------------------------------------------------------
# Exact enumeration of all non-negative uplifts
# ------------------------------------------------------------------
#
# Every uplift is produced by the following top-down search: the
# height-maximal mu in the remainder is componentwise maximal, hence equal
# to the top component (l1..l4) of every lambda that covers it; its
# multiplicity must therefore be split among lambda = (mu, l5),
# l5 in {mu4, mu4-1, ...} (>= 0).  The search is complete and exact.

def _compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for k in range(total + 1):
        for rest in _compositions(total - k, parts - 1):
            yield (k,) + rest


def _subtract(rem, lam, k):
    out = dict(rem)
    for nu in interlacing(lam):
        v = out.get(nu, 0) - k
        if v < 0:
            return None
        if v:
            out[nu] = v
        else:
            del out[nu]
    return out


def _branches(rem):
    mu = max(rem, key=height)
    k = rem[mu]
    l5s = list(range(mu[3], -1, -2))
    for comp in _compositions(k, len(l5s)):
        cur = rem
        chosen = []
        for l5, c in zip(l5s, comp):
            if c:
                cur = _subtract(cur, mu + (l5,), c)
                if cur is None:
                    break
                chosen.append((mu + (l5,), c))
        if cur is not None:
            yield chosen, cur


def count_uplifts(m, cap=None):
    """Number of non-negative uplifts (conjugation classes).  Memoised."""
    memo = {}

    def rec(rem):
        if not rem:
            return 1
        key = frozenset(rem.items())
        if key in memo:
            return memo[key]
        tot = 0
        for _, cur in _branches(rem):
            tot += rec(cur)
            if cap is not None and tot > cap:
                break
        memo[key] = tot
        return tot

    return rec(dict(m))


def enumerate_uplifts(m, limit=None):
    """Yield all non-negative uplifts {lambda: mult} (conjugation classes)."""
    found = [0]

    def rec(rem, acc):
        if limit is not None and found[0] >= limit:
            return
        if not rem:
            found[0] += 1
            yield dict(acc)
            return
        for chosen, cur in _branches(rem):
            for lam, c in chosen:
                acc[lam] = acc.get(lam, 0) + c
            yield from rec(cur, acc)
            for lam, c in chosen:
                acc[lam] -= c
                if not acc[lam]:
                    del acc[lam]

    yield from rec(dict(m), {})


def is_unique_optimum(m, objective, sol, maximize=True, solver=None):
    """
    Certify that `sol` is the unique optimum of `objective` over all
    non-negative integer uplifts of m.  On the optimal face
        {x in Z_{>=0}^cands : A x = m, objective(x) = best}
    each coordinate is minimised and maximised; sol is unique iff every
    coordinate is pinned.  LP-relaxation bounds are used first and exact
    ILPs only for coordinates the LP does not pin (LP bounds contain the
    integer bounds, so this is rigorous).  Returns True iff unique.
    """
    from sage.all import MixedIntegerLinearProgram

    cands = candidates(m)
    cover = {mu: [] for mu in m}
    for j, lam in enumerate(cands):
        for mu in interlacing(lam):
            cover[mu].append(j)
    best = sum(objective(l) * c for l, c in sol.items())

    def extreme(j, integer, sense_max):
        p = MixedIntegerLinearProgram(maximization=sense_max, solver=solver)
        x = p.new_variable(integer=integer, nonnegative=True)
        for mu, js in cover.items():
            p.add_constraint(p.sum(x[i] for i in js) == m[mu])
        p.add_constraint(p.sum(objective(lam) * x[i] for i, lam in enumerate(cands)) == best)
        p.set_objective(x[j])
        return p.solve()

    eps = 1e-6
    for j, lam in enumerate(cands):
        v = sol.get(lam, 0)
        lo = extreme(j, False, False)
        hi = extreme(j, False, True)
        if hi - lo < eps:
            continue
        if extreme(j, True, False) < v - 0.5 or extreme(j, True, True) > v + 0.5:
            return False
    return True


def class_order_key(lam):
    """Total order on conjugation classes used for lexicographic tie-breaking."""
    return (height(lam[:4]), lam[4], lam)


def canonical_uplift(m, solver=None):
    """
    Canonical uplift: (1) maximise H(c) = sum_lambda c_lambda height(l1..l4);
    (2) on that optimal face take the lexicographic maximum of c, classes
    ordered by class_order_key descending.  Lex-max over a finite set is
    unique, so this is always well defined.  Sequential exact ILPs.
    Returns (sol, best_height), or (None, None) if no uplift exists.
    """
    from sage.all import MixedIntegerLinearProgram

    first = solve_ilp(m, lambda l: height(l[:4]), maximize=True, solver=solver)
    if first is None:
        return None, None
    best = sum(height(l[:4]) * c for l, c in first.items())
    cands = sorted(candidates(m), key=class_order_key, reverse=True)
    cover = {mu: [] for mu in m}
    for j, lam in enumerate(cands):
        for mu in interlacing(lam):
            cover[mu].append(j)

    # One model; variables are fixed through bounds, only the objective changes.
    p = MixedIntegerLinearProgram(maximization=True, solver=solver)
    x = p.new_variable(integer=True, nonnegative=True)
    for mu, js in cover.items():
        p.add_constraint(p.sum(x[i] for i in js) == m[mu])
    p.add_constraint(p.sum(height(l[:4]) * x[i] for i, l in enumerate(cands)) == best)
    p.set_objective(x[0])
    b = p.get_backend()
    ncols = b.ncols()
    assert ncols == len(cands)
    fixed = {}
    for j in range(len(cands)):
        c = _col_of(p, x[j])
        coeffs = [0] * ncols
        coeffs[c] = 1
        b.set_objective(coeffs)
        b.solve()
        v = int(round(b.get_variable_value(c)))
        fixed[j] = v
        b.variable_lower_bound(c, v)
        b.variable_upper_bound(c, v)
    sol = {cands[j]: v for j, v in fixed.items() if v}
    assert restrict(sol) == m
    return sol, best


def _col_of(p, var):
    """Backend column index of a single MILP variable."""
    d = var.dict()
    assert len(d) == 1
    return next(iter(d))
