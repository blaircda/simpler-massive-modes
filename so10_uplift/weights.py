"""
Weight-lattice helpers for D4, B4, D5 in ambient (orthogonal) coordinates.

Convention: every weight is stored as a tuple of *doubled* ambient
coordinates, so that vector weights (1,0,0,0) -> (2,0,0,0) and spinor
weights (1/2,...,1/2) -> (1,1,1,1) are both integer tuples.  A tuple is a
valid weight iff all entries have the same parity.

Characters are stored by their W-invariant weight multiplicities, restricted
to dominant weights:  dict {dominant doubled weight: multiplicity}.

Irreducible decompositions use the Weyl alternation (Brauer/Klimyk) formula

    mult(V_lambda in chi) = sum_{w in W} det(w) m_chi(lambda + rho - w rho),

evaluated only at dominant lambda.
"""

from itertools import permutations, product


# ------------------------------------------------------------------
# Dominant representatives
# ------------------------------------------------------------------

def dom_B(x):
    """Dominant representative under W(B_n): sorted absolute values."""
    return tuple(sorted((abs(a) for a in x), reverse=True))


def dom_D(x):
    """Dominant representative under W(D_n): even numbers of sign flips."""
    neg = sum(1 for a in x if a < 0) % 2
    y = sorted((abs(a) for a in x), reverse=True)
    if neg and y[-1] != 0:
        y[-1] = -y[-1]
    return tuple(y)


def is_dom_B(x):
    return all(x[i] >= x[i + 1] for i in range(len(x) - 1)) and x[-1] >= 0


def is_dom_D(x):
    n = len(x)
    return all(x[i] >= x[i + 1] for i in range(n - 2)) and x[n - 2] >= abs(x[n - 1])


# ------------------------------------------------------------------
# Weyl groups as signed permutations, with determinant
# ------------------------------------------------------------------

def _perm_sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


def weyl_group(typ, n):
    """List of (perm, signs, det) for W(B_n) or W(D_n)."""
    out = []
    for p in permutations(range(n)):
        ps = _perm_sign(p)
        for signs in product((1, -1), repeat=n):
            nneg = signs.count(-1)
            if typ == "D" and nneg % 2:
                continue
            det = ps * (-1) ** nneg
            out.append((p, signs, det))
    return out


def act(w, x):
    p, signs, _ = w
    return tuple(signs[i] * x[p[i]] for i in range(len(x)))


RHO2 = {  # doubled Weyl vectors
    ("D", 4): (6, 4, 2, 0),
    ("B", 4): (7, 5, 3, 1),
    ("D", 5): (8, 6, 4, 2, 0),
}


def alternation_shifts(typ, n):
    """
    Aggregated {rho - w rho : sum of det(w)} over W, doubled coordinates.
    """
    rho = RHO2[(typ, n)]
    acc = {}
    for w in weyl_group(typ, n):
        wr = act(w, rho)
        s = tuple(rho[i] - wr[i] for i in range(n))
        acc[s] = acc.get(s, 0) + w[2]
    return {s: c for s, c in acc.items() if c != 0}


def decompose(mults, typ, n):
    """
    Irreducible decomposition of a W-invariant character given by dominant
    weight multiplicities (doubled coordinates).  Returns {hw: multiplicity}
    including possibly negative entries (caller checks positivity).
    """
    dom = dom_B if typ == "B" else dom_D
    isdom = is_dom_B if typ == "B" else is_dom_D
    shifts = alternation_shifts(typ, n)
    out = {}
    for lam in mults:
        if not isdom(lam):
            continue
        tot = 0
        for s, c in shifts.items():
            v = mults.get(dom(tuple(lam[i] + s[i] for i in range(n))), 0)
            if v:
                tot += c * v
        if tot:
            out[lam] = tot
    return out


# ------------------------------------------------------------------
# Conversions to Sage "coroots" Dynkin labels
# ------------------------------------------------------------------

def dynkin_D4(x):
    return ((x[0] - x[1]) // 2, (x[1] - x[2]) // 2, (x[2] - x[3]) // 2, (x[2] + x[3]) // 2)


def dynkin_B4(x):
    return ((x[0] - x[1]) // 2, (x[1] - x[2]) // 2, (x[2] - x[3]) // 2, x[3])


def dynkin_D5(x):
    return ((x[0] - x[1]) // 2, (x[1] - x[2]) // 2, (x[2] - x[3]) // 2,
            (x[3] - x[4]) // 2, (x[3] + x[4]) // 2)


# ------------------------------------------------------------------
# Weyl dimension formula (exact, rational-free with doubled coordinates)
# ------------------------------------------------------------------

def dim_irrep(lam, typ, n):
    from fractions import Fraction
    rho = RHO2[(typ, n)]
    l = [lam[i] + rho[i] for i in range(n)]
    r = list(rho)
    num = Fraction(1)
    for i in range(n):
        for j in range(i + 1, n):
            num *= Fraction((l[i] - l[j]) * (l[i] + l[j]), (r[i] - r[j]) * (r[i] + r[j]))
        if typ == "B":
            num *= Fraction(l[i], r[i])
    assert num.denominator == 1
    return int(num)


# ------------------------------------------------------------------
# Vectorised table of dominant weights, with gather maps
# ------------------------------------------------------------------

class DomTable:
    """
    All dominant weights (doubled coords) of type (typ, n) with squared
    doubled norm <= r2max, stored as rows of an int array.  A W-invariant
    function supported inside the ball is an int64 array of length
    len(table)+1 (the last slot is a permanent zero, used as the target of
    lookups falling outside the ball).

    gather(shift)[i] = index of dom(weights[i] - shift)  (or the zero slot).
    Hence  (f * e^{orbit})  style sums become  sum_s c_s f[gather(s)].
    """

    def __init__(self, typ, n, r2max):
        import numpy as np
        self.np = np
        self.typ, self.n, self.r2max = typ, n, r2max
        isdom = is_dom_B if typ == "B" else is_dom_D
        rows = []
        bound = int(r2max ** 0.5) + 1
        for par in (0, 1):
            def rec(prefix):
                if len(prefix) == n:
                    if isdom(prefix) and sum(a * a for a in prefix) <= r2max:
                        rows.append(tuple(prefix))
                    return
                lo = -bound if len(prefix) == n - 1 else 0
                hi = prefix[-1] if prefix else bound
                for a in range(lo, hi + 1):
                    if (a - par) % 2 == 0 and sum(b * b for b in prefix) + a * a <= r2max:
                        rec(prefix + [a])
            rec([])
        self.W = np.array(sorted(set(rows)), dtype=np.int64)
        self.size = len(self.W)
        self.ZERO = self.size
        self.off = bound + 2
        self.base = 2 * self.off + 1
        keys = self._key(self.W)
        order = np.argsort(keys)
        self._skeys = keys[order]
        self._sidx = order
        self.index = {tuple(int(a) for a in w): i for i, w in enumerate(self.W)}
        self._cache = {}

    def _key(self, X):
        k = self.np.zeros(len(X), dtype=self.np.int64)
        for c in range(self.n):
            k = k * self.base + (X[:, c] + self.off)
        return k

    def dominate(self, X):
        np = self.np
        A = np.sort(np.abs(X), axis=1)[:, ::-1].copy()
        if self.typ == "D":
            neg = (X < 0).sum(axis=1) % 2 == 1
            flip = neg & (A[:, -1] != 0)
            A[flip, -1] *= -1
        return A

    def lookup(self, X):
        """Indices of the (already dominant) rows X, ZERO if not in table."""
        np = self.np
        inball = (X * X).sum(axis=1) <= self.r2max
        Xc = np.clip(X, -self.off, self.off)
        k = self._key(Xc)
        pos = np.searchsorted(self._skeys, k)
        pos = np.minimum(pos, self.size - 1)
        hit = (self._skeys[pos] == k) & inball
        out = np.full(len(X), self.ZERO, dtype=np.int64)
        out[hit] = self._sidx[pos[hit]]
        return out

    def gather(self, shift):
        shift = tuple(shift)
        g = self._cache.get(shift)
        if g is None:
            X = self.W - self.np.array(shift, dtype=self.np.int64)
            g = self.lookup(self.dominate(X))
            g = self.np.append(g, self.ZERO)  # zero slot maps to zero slot
            self._cache[shift] = g
        return g

    def zeros(self):
        return self.np.zeros(self.size + 1, dtype=self.np.int64)

    def to_dict(self, f):
        nz = self.np.nonzero(f[:self.size])[0]
        return {tuple(int(a) for a in self.W[i]): int(f[i]) for i in nz}

    def from_dict(self, d):
        f = self.zeros()
        for w, v in d.items():
            f[self.index[w]] = v
        return f

    def decompose(self, f):
        """Weyl-alternation irreducible multiplicities, as an array on the table."""
        out = self.zeros()
        for s, c in alternation_shifts(self.typ, self.n).items():
            out += c * f[self.gather(tuple(-a for a in s))]
        out[self.ZERO] = 0
        return out
