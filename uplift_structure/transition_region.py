"""
S14: locating the transition region of the repairs (cumulative R9_n, tensor sector).

  (1) shape of the spectrum: ln(R9_n(nu)/dim nu) - ln R9_n(0) along the symmetric
      tensors (k,0,0,0) and along (k,k,0,0), with the mean slope per unit k over
      k = 4..12 and 16..40 (a dimension-proportional regime would give slope ~ 0);
      (1b) a heat-kernel fit ln(R/dim) = A - beta*C2(nu) over the bulk
      (nu1 <= 0.9 sqrt(n)), with R^2 and beta*n;
  (2) the negative bottom-level coefficients of the all-full lift a = Delta^4 R,
      by column family (j,i) (columns (k,j,i), level 0): their k-ranges, and the
      excess band e = -a(k,j,i,0) - sum_{L=1..i} a(k,j,i,L) > 0 (needs level i+1);
  (3) the edges over n: the offset d(j,i) = n - k_max of each family's negative
      range from the top, and the (1,1) excess band [k_lo, k_hi] with n - k_hi.

Usage: transition_region.py [--data path]   Output: results/transition_region[_<N>].txt
"""

import argparse
import math
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cumulative_lift as CL
from common import RESULTS
from regge_stress import load_R9

HERE = os.path.dirname(os.path.abspath(__file__))
RHO = (3.5, 2.5, 1.5, 0.5)
FAMS = [(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (3, 2), (3, 3), (4, 4), (5, 5), (6, 6)]


def dimf(x):
    y = [x[i] + RHO[i] for i in range(4)]
    num = den = 1.0
    for i in range(4):
        num *= y[i]
        den *= RHO[i]
        for j in range(i + 1, 4):
            num *= (y[i] - y[j]) * (y[i] + y[j])
            den *= (RHO[i] - RHO[j]) * (RHO[i] + RHO[j])
    return num / den


def ranges(ks):
    ks, out, s, p = sorted(ks), [], None, None
    for k in ks:
        if s is None:
            s = p = k
        elif k == p + 1:
            p = k
        else:
            out.append((s, p))
            s = p = k
    if s is not None:
        out.append((s, p))
    return ",".join(f"{a}" if a == b else f"{a}-{b}" for a, b in out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "data", "r9_to40.json"))
    args = ap.parse_args()
    N, r9, R9 = load_R9(args.data)
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    log(f"data {os.path.basename(args.data)}: n <= {N}")
    log("(1) ln(R/dim) - ln R(0): mean slope per unit k over k=4..12 and k=16..40 (sym = (k,0,0,0); 2-row = (k,k,0,0))")
    for n in range(24, N + 1, 4):
        y0 = math.log(R9[n][(0, 0, 0, 0)])

        def f(x):
            v = R9[n].get(tuple(2 * c for c in x), 0)
            return math.log(v / dimf(x)) - y0 if v else None
        s1 = [f((k, 0, 0, 0)) for k in (4, 12, 16, min(40, n - 4))]
        s2 = [f((k, k, 0, 0)) for k in (4, 12)]
        sl = lambda a, b, ka, kb: (b - a) / (kb - ka) if a is not None and b is not None else float("nan")
        log(f"  n={n}: sym slope {sl(s1[0], s1[1], 4, 12):+.3f} (k=4..12), {sl(s1[2], s1[3], 16, min(40, n - 4)):+.3f} "
            f"(k=16..{min(40, n - 4)}); 2-row slope {sl(s2[0], s2[1], 4, 12):+.3f}; values at k=1,2,4: "
            + ", ".join(f"{f((k, 0, 0, 0)):+.2f}" for k in (1, 2, 4)))
    import numpy as np
    log("(1b) heat-kernel fit ln(R/dim) = A - beta*C2 over the bulk nu1 <= 0.9 sqrt(n): beta*n, R^2")
    C2 = lambda x: sum(x[i] * (x[i] + 2 * RHO[i]) for i in range(4))
    for n in range(16, N + 1, 4):
        X, Y = [], []
        for k, v in R9[n].items():
            x = tuple(w / 2 for w in k)
            if k[0] % 2 or v <= 0 or x[0] > 0.9 * math.sqrt(n):
                continue
            X.append(C2(x))
            Y.append(math.log(v / dimf(x)))
        X, Y = np.array(X), np.array(Y)
        b, a0 = np.polyfit(X, Y, 1)
        res = Y - (a0 + b * X)
        log(f"  n={n}: beta*n = {-b * n:.3f}, beta*n/sqrt(n) = {-b * math.sqrt(n):.3f}, "
            f"R^2 = {1 - (res ** 2).sum() / ((Y - Y.mean()) ** 2).sum():.4f}")
    edge_rows = []
    for n in range(20, N + 1):
        a, _ = CL.all_full(R9[n], 0)
        neg, exc = defaultdict(list), defaultdict(list)
        for y, v in a.items():
            if y[3] != 0 or v >= 0:
                continue
            k, j, i = y[:3]
            neg[(j, i)].append(k)
            if i >= 1 and -v - sum(a.get((k, j, i, L), 0) for L in range(1, i + 1)) > 0:
                exc[(j, i)].append(k)
        if n == N or n == N - 4:
            log(f"(2) n={n}: negative bottom columns (k,j,i), by family (j,i): k-ranges [excess k-ranges -> needs level i+1]")
            for ji in sorted(neg):
                ex = f"   EXCESS k={ranges(exc[ji])} -> level {ji[1] + 1}" if ji in exc else ""
                log(f"    {ji}: k={ranges(neg[ji])}{ex}")
        ex11 = exc.get((1, 1), [])
        edge_rows.append(f"  {n:2d} | " + " ".join(f"{n - max(neg[f]):6d}" if neg.get(f) else "     -" for f in FAMS)
                         + (f" | [{min(ex11)},{max(ex11)}] n-k_hi={n - max(ex11)}" if ex11 else " | -")
                         + (f" | other excess families: {sorted(set(exc) - {(1, 1)})}" if set(exc) - {(1, 1)} else ""))
    log("(3) offsets d(j,i) = n - k_max of the negative range, families " + " ".join(str(f) for f in FAMS)
        + "; (1,1) excess band")
    for r in edge_rows:
        log(r)
    os.makedirs(RESULTS, exist_ok=True)
    name = "transition_region.txt" if N == 40 else f"transition_region_{N}.txt"
    with open(os.path.join(RESULTS, name), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
