"""
S13 (part 2): why and where repair levels are needed, from data and from
asymptotics.

All-full lift a = Delta^4 R of the cumulative R9_n (cumulative_lift.py).  A
level-L class skipping mu4 = 0, (m1,m2,m3,L,1) instead of (m1,m2,m3,L,0),
changes the bottom-level coefficients by
    +[m1,m2,m3]  and  -[m1,m2,L-1] + [m1,m3-1,L-1] - [m2-1,m3-1,L-1],
so a bottom deficit at a column (k,c,c) is covered by same-column repairs at
levels 1..c or by pattern feed from level c+1.  Hence the indicator

  e_c(n, k) = -a(k,c,c,0) - sum_{L=1..c} a(k,c,c,L)
  (> 0: level c+1 needed by the same-column count)

reported (1) as max over k, and (2) as the relative excess e/capacity at bulk
columns k = c..c+7, for c = 0..3 and every n of the data.

(3) Asymptotics at fixed weights: multiplicities ~ C(n) dim(nu), so
a(l) ~ C(n) D4(l), D4(l) = sum_eps (-1)^|eps| dim(l + eps) (dominant terms).
The script checks the sign of D4 on all dominant l with l1 <= 20 and prints the
asymptotic relative excess for columns (k,c,c), c = 0..5.

Usage: repair_asymptotics.py [--data path]     Output: results/repair_asymptotics[_<N>].txt
"""

import argparse
import os
import sys
from fractions import Fraction
from itertools import product

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cumulative_lift as CL
from common import RESULTS
from regge_stress import load_R9

HERE = os.path.dirname(os.path.abspath(__file__))


def dim_b4(x):
    rho = (Fraction(7, 2), Fraction(5, 2), Fraction(3, 2), Fraction(1, 2))
    y = [x[i] + rho[i] for i in range(4)]
    num = den = Fraction(1)
    for i in range(4):
        num *= y[i]
        den *= rho[i]
        for j in range(i + 1, 4):
            num *= (y[i] - y[j]) * (y[i] + y[j])
            den *= (rho[i] - rho[j]) * (rho[i] + rho[j])
    return num / den


def d4_dim(l):
    def R(x):
        return dim_b4(x) if x[0] >= x[1] >= x[2] >= x[3] >= 0 else 0
    return sum((-1) ** sum(e) * R(tuple(p + q for p, q in zip(l, e))) for e in product((0, 1), repeat=4))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "data", "r9_to40.json"))
    args = ap.parse_args()
    N, r9, R9 = load_R9(args.data)
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    log(f"data {os.path.basename(args.data)}: n <= {N} (tensor sector)")
    log("(1) max_k e_c(n,k) (argmax k):  n | c=0 | c=1 | c=2 | c=3")
    rel = {c: {} for c in range(4)}
    for n in range(6, N + 1):
        a, _ = CL.all_full(R9[n], 0)
        row = []
        for c in range(4):
            best = None
            for k in sorted({y[0] for y in a if y[1:3] == (c, c)}):
                if k < c:
                    continue
                d = -a.get((k, c, c, 0), 0)
                cap = sum(a.get((k, c, c, L), 0) for L in range(1, c + 1))
                v = d - cap
                if best is None or v > best[0]:
                    best = (v, k)
                if c <= k <= c + 7 and cap:
                    rel[c].setdefault(n, {})[k] = (d - cap) / cap
            row.append(f"{best[0]} ({best[1]})" if best else "-")
        log(f"  {n:2d} | " + " | ".join(row))
    log("(2) relative excess (deficit - capacity)/capacity at bulk columns (k,c,c), k = c..c+7")
    for c in (1, 2, 3):
        log(f"  c={c}:")
        for n in sorted(rel[c]):
            if n % 4 == 0 or n == N:
                log(f"    n={n:2d}: " + " ".join(f"{rel[c][n].get(k, float('nan')):+.3f}" for k in range(c, c + 8)))
    neg = sum(1 for l in product(range(21), repeat=4) if l[0] >= l[1] >= l[2] >= l[3] and d4_dim(l) < 0)
    mn = min(d4_dim(l) for l in product(range(21), repeat=4) if l[0] >= l[1] >= l[2] >= l[3])
    log(f"(3) D4(dim) on dominant l with l1 <= 20: {neg} negative values; minimum {mn}")
    for c in range(6):
        row = []
        for k in range(max(c, 1), max(c, 1) + 8):
            d = -d4_dim((k, c, c, 0))
            cap = sum(d4_dim((k, c, c, L)) for L in range(1, c + 1))
            row.append(f"{float((d - cap) / cap):+.3f}" if cap else f"d={float(d):+.0f}")
        log(f"  c={c}: asymptotic relative excess at (k,{c},{c}), k={max(c, 1)}..: " + " ".join(row))
    os.makedirs(RESULTS, exist_ok=True)
    name = "repair_asymptotics.txt" if N == 40 else f"repair_asymptotics_{N}.txt"
    with open(os.path.join(RESULTS, name), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
