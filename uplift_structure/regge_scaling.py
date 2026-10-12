"""
S15 (verification): the spectrum data against the Regge-limit model of regge_limit.py
(which reads no data; its output data/regge_limit_D<Dmax>.json is loaded here).

Weight nu = (x1, mu) at level n, top distance D = n - x1, l = x1 - mu_1, GMSZ depth
w = D - cost7(mu).  Model: r9_n(nu) -> m_D(mu) as l -> oo at fixed (D, mu).

  (1) per level: every (n, x1, mu) with mu_1 <= LM, D <= Dmax, n <= N: r9_n against
      m_D(mu): counts equal / below / above; monotonicity in l at fixed (D, mu);
      the saturation point l_sat(D, mu) (r9 = m_D for all l >= l_sat in the data)
      against the GMSZ depth w;
  (2) cumulative: R9_n(nu) against G_D(mu) = sum_{D' <= D} m_D'(mu), same counts,
      and the relative deficit 1 - R9/G against l - w;
  (3) the all-full lift a = Delta^4 R9_n (tensor sector) at columns (n-D, j, i):
      the offsets d_n(j,i) = n - k_max of the bottom negatives against the model
      offsets, and the (1,1) excess e(n, n-D) against e_D for D = 24..36;
  (4) where the data equal the model exactly: for each n, the least D at which some
      a(n-D, j, i, L) (j <= 6) differs from a_D(j,i,L).

Usage: regge_scaling.py [--data path] [--model path]   Output: results/regge_scaling[_<N>].txt
"""

import argparse
import json
import os
import sys
from collections import Counter
from itertools import product

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS
from regge_limit import cost7_t
from regge_stress import load_R9

HERE = os.path.dirname(os.path.abspath(__file__))
FAMS = [(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (3, 2), (3, 3), (4, 4), (5, 5), (6, 6)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "data", "r9_to40.json"))
    ap.add_argument("--model", default=os.path.join(HERE, "data", "regge_limit_D60.json"))
    args = ap.parse_args()
    N, r9, R9 = load_R9(args.data)
    M = json.load(open(args.model))
    Dmax, LM = M["Dmax"], M["LM"]
    m = {int(t): {tuple(l): v for l, v in rows} for t, rows in M["m"].items()}
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    log(f"data {os.path.basename(args.data)} (n <= {N}) against model {os.path.basename(args.model)} "
        f"(D <= {Dmax}, mu_1 <= {LM})")
    G = {}
    for par in (0, 1):
        acc = {}
        for t in range(2 + par, 2 * Dmax + 1, 2):
            for l, v in m.get(t, {}).items():
                acc[l] = acc.get(l, 0) + v
            G[t] = dict(acc)

    # (1), (2): every (n, x1, mu) in range
    for name, data, model in (("(1) per level r9_n vs m_D", r9, m), ("(2) cumulative R9_n vs G_D", R9, G)):
        cnt, mono, sat, unsat, viol, ramp = Counter(), 0, Counter(), 0, 0, {}
        mus = sorted({l for t in model for l in model[t]})
        for t in range(2, min(2 * Dmax, 2 * N) + 1):
            for mu in mus:
                if mu[0] > 2 * LM or mu[0] % 2 != t % 2:
                    continue
                mv = model[t].get(mu, 0)
                w = (t - cost7_t(*mu)) // 2              # t - cost7_t is even
                seq = []
                for X in range(mu[0], 2 * N - t + 1, 2):            # X = 2 x1 >= 2 mu_1; n = (X + t)/2 <= N
                    n = (X + t) // 2
                    v = data[n].get((X,) + mu, 0)
                    seq.append(v)
                    cnt["equal" if v == mv else ("below" if v < mv else "above")] += 1
                    viol += v != mv and (X - mu[0]) // 2 >= w
                if t % 2 == 0 and mu in ((0, 0, 0), (2, 2, 0)) and w >= 8 and w % 4 == 0 and mv:
                    ramp[(mu, w)] = [v / mv for v in seq]
                if not seq:
                    continue
                mono += sum(1 for a, b in zip(seq, seq[1:]) if b < a)
                if mv == 0:
                    continue
                if seq[-1] == mv:
                    ls = len(seq) - 1
                    while ls > 0 and seq[ls - 1] == mv:
                        ls -= 1
                    sat[ls - w] += 1
                else:
                    unsat += 1
        log(f"{name}: {sum(cnt.values())} comparisons: {dict(cnt)}; decreases along l at fixed (D,mu): {mono}")
        log(f"    saturation point l_sat - w over (D,mu) saturated in the data: {dict(sorted(sat.items()))}; "
            f"not yet saturated at the largest l: {unsat}; (D,mu,l) with l >= w and data != model: {viol}")
        log("    ramp at fixed (mu, w): least l with data/model >= 1/2, >= 0.99, = 1 (- : not reached at the largest l)")
        for (mu, w), q in sorted(ramp.items()):
            f = lambda c: next((l for l, x in enumerate(q) if x >= c), "-")
            log(f"      mu={tuple(x // 2 for x in mu)} w={w:2d}: {f(0.5)}, {f(0.99)}, {f(1.0)}  (l <= {len(q) - 1})")

    # (3), (4): the all-full lift
    def a_model(D, y):
        t = 2 * D
        tot = 0
        for e in product((0, 1), repeat=3):
            z = tuple(p + q for p, q in zip(y, e))
            if z[0] >= z[1] >= z[2] >= 0:
                tot += (-1) ** sum(e) * m.get(t, {}).get(tuple(2 * v for v in z), 0)
        return tot

    def a_data(n, y):
        tot = 0
        for e in product((0, 1), repeat=4):
            z = tuple(p + q for p, q in zip(y, e))
            if z[0] >= z[1] >= z[2] >= z[3] >= 0:
                tot += (-1) ** sum(e) * R9[n].get(tuple(2 * v for v in z), 0)
        return tot

    off_model = {}
    for f in FAMS:
        neg = [D for D in range(1, Dmax + 1) if a_model(D, f + (0,)) < 0]
        off_model[f] = neg[0] if neg else None
    log("(3) offsets d(j,i): model (Regge limit) | data d_n(j,i) for n = N-4..N, families " + " ".join(map(str, FAMS)))
    log("    model : " + " ".join(f"{off_model[f] if off_model[f] else '-':>6}" for f in FAMS))
    for n in range(max(20, N - 8), N + 1):
        row = []
        for f in FAMS:
            ks = [k for k in range(f[0], n + 1) if a_data(n, (k,) + f + (0,)) < 0]
            row.append(f"{n - max(ks):6d}" if ks else "     -")
        log(f"    n={n:2d}: " + " ".join(row))
    log("    (1,1) excess e = -a(.,1,1,0) - a(.,1,1,1) at top distance D: data e(n, n-D) for n = N-8, N-6, .., N | model e_D")
    for D in range(24, 37):
        ed = []
        for n in range(N - 8, N + 1, 2):
            k = n - D
            ed.append(f"{-a_data(n, (k, 1, 1, 0)) - a_data(n, (k, 1, 1, 1)):+10d}" if k >= 1 else "         -")
        em = -a_model(D, (1, 1, 0)) - a_model(D, (1, 1, 1))
        log(f"      D={D}: " + " ".join(ed) + f" | {em:+10d}")
    log("(4) least D at which the data's all-full lift a(n-D, j, i, L) (j <= 6, all L) differs from the model a_D(j,i,L)")
    out = []
    for n in range(10, N + 1):
        Dd = None
        for D in range(1, n + 1):
            if any(a_data(n, (n - D, j, i, L)) != a_model(D, (j, i, L))
                   for j in range(0, 7) for i in range(0, j + 1) for L in range(0, i + 1) if n - D >= j):
                Dd = D
                break
        out.append(f"n={n}: {Dd}")
    log("    " + ", ".join(out))
    os.makedirs(RESULTS, exist_ok=True)
    name = "regge_scaling.txt" if N == 40 else f"regge_scaling_{N}.txt"
    with open(os.path.join(RESULTS, name), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
