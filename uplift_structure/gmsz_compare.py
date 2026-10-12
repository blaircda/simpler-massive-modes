"""
S9: comparison with Georgoudis, Minahan, Strom, Zoumis, "The growth of SO(9)
super-representations in Type II string theory", arXiv:2609.04107 (GMSZ).

GMSZ's a_{[l,m,r,s],n} is our r9_n (SO(9) content of level n with the 256 =
44+84+128 divided out).  This script tests their empirical statements on our
independently derived data (data/r9_to32.json) and relates their shift
relations to the Regge structures of S5-S8:

  (1) our r9_7 vs their worked example (E.19) (labelled "level 6" there);
  (2) their support rule, Table 1 / eq. (3.6), incl. the gap [n-2,0,0,0],
      for 7 <= n <= 32 (level offsets -1, 0, +1);
  (3) their shift relations (3.32) l -> l+k at level +k, (3.33) m -> m+k at
      +3k, (3.34) r -> r+k at +6k, in their stated regions;
  (4) sharpness: for k = 1, equality as a function of the distance ("margin")
      from their region, and the inequality a_n >= a_(n-d) everywhere;
  (5) monotonicity r9_n(nu + w) >= r9_(n-d)(nu) for all nu, for column moves
      w = w1, w2, w3, 2w4 (a column of height 1..4), single boxes e_i, and
      neighbouring level offsets as controls;
  (6) where the SO(10) heads E_n sit relative to the l-stable range of (3.32),
      and whether top(A_(n-1)) = r9_(n-1) on that range (both rules, n <= 20);
  (7) SO(10) analogues of the r-shift and of a 4-form shift as head
      directions (containment under plain rules; forced inductions; null test),
      with a self-dual 5-form shift as a control for how discriminating the
      forced tests are at large level offsets.

The literature data used here (E.19, Table 1, (3.32)-(3.35)) are typed in from
the paper; nothing is imported from baseline files.

Output: results/gmsz_compare.txt
"""

import itertools
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import regge_directions as RD
from common import RESULTS
from regge_stress import T

from uplift import interlacing, restrict
from weights import dim_irrep, dynkin_B4

N, r9 = RD.N, RD.r9

# GMSZ (E.19): the level with top [6,0,0,0]
E19 = {(6, 0, 0, 0): 1, (4, 0, 0, 1): 1, (3, 1, 0, 0): 1, (4, 0, 0, 0): 1, (3, 0, 0, 1): 1, (2, 0, 1, 0): 1,
       (2, 1, 0, 0): 2, (3, 0, 0, 0): 1, (1, 1, 0, 1): 1, (2, 0, 0, 1): 3, (0, 2, 0, 0): 1, (1, 0, 0, 2): 1,
       (1, 0, 1, 0): 1, (1, 1, 0, 0): 2, (2, 0, 0, 0): 3, (0, 1, 0, 1): 2, (1, 0, 0, 1): 4, (0, 0, 0, 2): 1,
       (0, 0, 1, 0): 3, (0, 1, 0, 0): 2, (1, 0, 0, 0): 2, (0, 0, 0, 1): 2, (0, 0, 0, 0): 2}


def c_rs(r, s):
    """GMSZ (3.6): [l,m,r,s] is allowed at level n iff n - l - 3m >= c(r, s)."""
    return max(6 * r + 5 * s - 6, 6 * r + 4 * s - 2, 6 * r + 3 * s - 1, 2 * s + 1)


def lmax(n, m, r, s):
    return n - 3 * m - c_rs(r, s)


def allowed(n, k):
    return min(k) >= 0 and k[0] <= lmax(n, *k[1:])


def eps(s):
    return 0 if s > 3 else (1 if s >= 1 else 2)


def labels(n):
    L = n + 1
    for k in itertools.product(range(L + 1), range(L // 3 + 1), range(L // 5 + 2), range(L + 1)):
        if allowed(n, k):
            yield k


def main():
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    a = {n: {tuple(dynkin_B4(k)): v for k, v in r9[n].items()} for n in range(1, N + 1)}

    def A(n, k):
        return a[n].get(k, 0) if 1 <= n <= N else 0

    # (1)
    dim7 = sum(v * dim_irrep(k, 'B', 4) for k, v in r9[7].items())
    log(f"(1) our r9_7 == GMSZ (E.19): {a[7] == E19}  (23 irreps; total dimension {dim7}, GMSZ (E.20): 26136)")

    # (2)
    log("(2) support rule (3.6) incl. the gap [n-2,0,0,0], levels 7..%d" % N)
    for off in (0, -1, 1):
        bad = 0
        for n in range(7, N + 1):
            sup = {k for k, v in a[n].items() if v}
            pred = {k for k in labels(n + off) if k != (n + off - 2, 0, 0, 0)}
            bad += sup != pred
        log(f"  level offset {off:+d}: support differs at {bad} of {N - 6} levels")

    # (3) and (4)
    log("(3) shift relations as stated; (4) k=1: equality by margin from the stated region "
        "(region = margin > 0), and the inequality a_n >= a_(n-d) on all labels")
    for name in ("l", "m", "r"):
        tot = viol = 0
        bymargin, ge, ge_tot = {}, 0, 0
        for n in range(2, N + 1):
            for k in labels(n):
                l, m, r, s = k
                LM = lmax(n, m, r, s)
                for kk in range(1, n):
                    if name == "l":
                        k2, n2, cond, marg = (l - kk, m, r, s), n - kk, 2 * l > LM + kk, 2 * l - (LM + 1)
                    elif name == "m":
                        k2, n2, cond, marg = (l, m - kk, r, s), n - 3 * kk, l > LM - m + kk - 1, l - (LM - m)
                    else:
                        k2, n2 = (l, m, r - kk, s), n - 6 * kk
                        cond, marg = l > LM - r + (kk - 1) + eps(s), l - (LM - r + eps(s))
                    if n2 < 1 or not allowed(n2, k2):
                        continue
                    if cond:
                        tot += 1
                        viol += A(n, k) != A(n2, k2)
                    if kk == 1:
                        t = bymargin.setdefault(marg, [0, 0])
                        t[0] += 1
                        t[1] += A(n, k) == A(n2, k2)
                        ge_tot += 1
                        ge += A(n, k) >= A(n2, k2)
        row = ", ".join(f"{mg}: {bymargin[mg][1]}/{bymargin[mg][0]}" for mg in sorted(bymargin) if -3 <= mg <= 2)
        log(f"  {name}-shift: {tot} instances in the stated regions, {viol} violations")
        log(f"     k=1 equality by margin  {row}")
        log(f"     k=1 inequality holds in {ge}/{ge_tot}")

    # (5)
    log("(5) monotonicity r9_n(nu + w) >= r9_(n-d)(nu), all nu in r9_(n-d), n <= %d" % N)
    moves = [("column 1 = w1", (2, 0, 0, 0), (1, 2)), ("column 2 = w2", (2, 2, 0, 0), (2, 3)),
             ("column 3 = w3", (2, 2, 2, 0), (5, 6)), ("column 4 = 2w4", (2, 2, 2, 2), (9, 10)),
             ("box row 2 = e2", (0, 2, 0, 0), (2, 3)), ("box row 3 = e3", (0, 0, 2, 0), (3, 4)),
             ("box row 4 = e4", (0, 0, 0, 2), (3, 4)), ("spinor w4", (1, 1, 1, 1), (3, 4)), ("level only", (0, 0, 0, 0), (1,))]
    for name, w, ds in moves:
        for d in ds:
            viol, tot, lead = 0, 0, 0
            for n in range(d + 1, N + 1):
                for nu, v in r9[n - d].items():
                    t = tuple(x + y for x, y in zip(nu, w))
                    if not t[0] >= t[1] >= t[2] >= t[3]:
                        continue
                    tot += 1
                    if r9[n].get(t, 0) < v:
                        viol += 1
                        lead += tuple(dynkin_B4(nu))[1:] == (0, 0, 0) and nu[0] == 2 * (n - d - 1)
            log(f"  {name:16s} at +{d:2d}: {viol:5d} violations of {tot:5d} "
                f"({lead} of them at the top [n-d-1,0,0,0] of level n-d)")

    # (6)
    log("(6) SO(10) heads vs the l-stable range of (3.32) (margin = 2l - (lmax+1); equality for margin >= 0)")

    def margin_l(n, mu):
        l, m, r, s = dynkin_B4(mu)
        return 2 * l - (lmax(n, m, r, s) + 1)

    for rule in ("max-height", "min-height"):
        H = RD.plain_heads(rule, 20)
        Acum, worst, sat = {}, [], 0
        for n in range(1, 21):
            if n >= 2:
                top = {}
                for lam, cc in Acum.items():
                    for mu in interlacing(lam):
                        if mu[0] == lam[0]:
                            top[mu] = top.get(mu, 0) + cc
                for mu, v in r9[n - 1].items():
                    if margin_l(n, (mu[0] + 2,) + mu[1:]) >= 0 and top.get(mu, 0) != v:
                        sat += 1
            resE = restrict(H[n]) if H[n] else {}
            worst.append(max((margin_l(n, mu) for mu in resE), default=None))
            Acum = T(Acum)
            for l, cc in H[n].items():
                Acum[l] = Acum.get(l, 0) + cc
        log(f"  [{rule}] max margin of Res E_n, n=1..20: {worst}")
        log(f"  [{rule}] on margin >= 0: top(A_(n-1)) != r9_(n-1) at {sat} weights")
    log(f"  [{time.time() - t0:.0f}s]")

    # (7)
    log("(7) SO(10) head directions with the r-shift and 4-form costs (control: self-dual 5-form at +10)")
    R3, F4, SD5 = (2, 2, 2, 0, 0), (2, 2, 2, 2, 0), (2, 2, 2, 2, 2)
    H = {rule: RD.plain_heads(rule, 16) for rule in ("max-height", "min-height")}
    for name, s, k in (("3-form [0,0,1,0,0]", R3, 6), ("4-form [0,0,0,1,1]", F4, 10), ("sd 5-form [0,0,0,0,2]", SD5, 10)):
        for rule in H:
            E = H[rule]
            fails = [n for n in range(k + 1, 17) if any(E[n].get(l, 0) < c for l, c in RD.shifted(E[n - k], s).items())]
            log(f"  {name} at +{k} [{rule}]: plain containment fails at n = {fails or 'none'} (n <= 16)")
        out = [RD.forced_run([(s, k, False, True)], rule)[0] for rule in RD.RULES]
        log(f"  {name} at +{k}: forced " + " | ".join(out))
        nt = [RD.forced_run([(s, k, False, False)], "max-height", m0=m0)[0] for m0 in (2, 5, 7, 9)]
        log(f"  {name} at +{k}: null test m0=2,5,7,9: " + " | ".join(nt))
    log(f"[{time.time() - t0:.0f}s]")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "gmsz_compare.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
