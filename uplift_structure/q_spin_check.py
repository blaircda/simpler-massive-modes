"""
S16 (verification): the spectrum data against Q (q_spin.py, which reads no data;
its output data/q_spin_D<Dmax>.json is loaded here).

  (1) every r9_n(N, mu) with mu_1 <= LM and 1 <= D = n - N <= Dmax, n <= N_data (all dominant
      (N, mu), zeros included) against sum_{e <= N} c(D, e; mu);
  (2) every start r9_n(N, mu) - r9_{n-1}(N-1, mu) (N > mu_1) against c(D, N; mu);
  (3) the upper-half theorem (S18) on the raw data, independently of Q: for every level n and every
      weight (N, mu) with 3N >= n - 1 and N > mu_1, the start equals the multiplicity of mu in
      [q^(n-2N)] of the positive Fock space Qhat (q_top_end.qhat_fock).

Usage: q_spin_check.py [--data path] [--q path]   Output: results/q_spin_check[_<N>].txt
"""

import argparse
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS
import q_top_end
from regge_stress import load_R9

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "data", "r9_to40.json"))
    ap.add_argument("--q", default=os.path.join(HERE, "data", "q_spin_D22.json"))
    args = ap.parse_args()
    N, r9, _ = load_R9(args.data)
    Qd = json.load(open(args.q))
    Dmax, LM = Qd["Dmax"], Qd["LM"]
    edge = {(t, tuple(l)): v for t, l, v in Qd["edge"]}
    st = {(t, tuple(l)): dict(r) for t, l, r in Qd["starts"]}
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    log(f"data {os.path.basename(args.data)} (n <= {N}) against {os.path.basename(args.q)} (D <= {Dmax}, mu_1 <= {LM})")
    cnt, scnt, bad, sbad, nz = Counter(), Counter(), [], [], 0
    for t in range(2, 2 * Dmax + 1):
        par = t % 2
        for A in range(par, 2 * LM + 1, 2):
            for B in range(par, A + 1, 2):
                for C in range(par, B + 1, 2):
                    l = (A, B, C)
                    row, e0 = st.get((t, l), {}), edge.get((t, l), 0)
                    for X in range(A, 2 * N - t + 1, 2):         # doubled spin; level n = (X + t)/2 <= N
                        n = (X + t) // 2
                        v = r9[n].get((X,) + l, 0)
                        nz += v > 0
                        pred = e0 + sum(x for s, x in row.items() if s <= X)
                        cnt[pred == v] += 1
                        if pred != v:
                            bad.append((n, X, l, v, pred))
                        if X > A:
                            prev = r9[n - 1].get((X - 2,) + l, 0) if n > 1 else 0
                            scnt[row.get(X, 0) == v - prev] += 1
                            if row.get(X, 0) != v - prev:
                                sbad.append((n, X, l, v - prev, row.get(X, 0)))
    log(f"(1) r9_n(N, mu) = sum_(e <= N) c(D, e; mu): {cnt[True]} equal, {cnt[False]} different "
        f"({nz} of the compared multiplicities are non-zero)" + (f"; e.g. {bad[:4]}" if bad else ""))
    log(f"(2) starts r9_n(N, mu) - r9_(n-1)(N-1, mu) = c(D, N; mu): {scnt[True]} equal, {scnt[False]} different"
        + (f"; e.g. {sbad[:4]}" if sbad else ""))
    # (3) the upper-half theorem on the raw data
    G = (N + 2) // 3 + 1
    Hf = q_top_end.qhat_fock(G, 2 * G + 8)
    pred = {}
    for g in range(2, G + 1):
        for par in (0, 1):
            for A in range(par, 2 * g + 5, 2):
                for B in range(par, A + 1, 2):
                    for C in range(par, B + 1, 2):
                        m = Hf.mult(g, (A, B, C))
                        if m:
                            pred[(g, (A, B, C))] = m
    pts = set()
    for n in range(1, N + 1):
        for nu in r9[n]:
            pts.add((n, nu))
            if n < N:
                pts.add((n + 1, (nu[0] + 2,) + nu[1:]))
        for (g, l) in pred:
            X = n - g                                            # doubled spin X = 2N, g = n - 2N
            pts.add((n, (X,) + l))
    cnt3, bad3, nz3 = Counter(), [], 0
    for n, nu in pts:
        X, l = nu[0], nu[1:]
        if 3 * X < 2 * (n - 1) or X <= l[0] or X % 2 != l[0] % 2 or n > N:
            continue
        start = r9[n].get(nu, 0) - (r9[n - 1].get((X - 2,) + l, 0) if n > 1 else 0)
        g2 = 2 * n - 2 * X                                       # 2g = 2(n - 2N)
        pr = pred.get((g2 // 2, l), 0) if g2 % 2 == 0 else None
        nz3 += start != 0
        cnt3[start == pr] += 1
        if start != pr:
            bad3.append((n, X / 2, tuple(x / 2 for x in l), start, pr))
    log(f"(3) upper-half theorem on the raw data (3N >= n - 1, N > mu_1, n <= {N}): start = mult of mu in "
        f"[q^(n-2N)] Qhat (Fock form): {cnt3[True]} equal, {cnt3[False]} different ({nz3} non-zero starts)"
        + (f"; e.g. {bad3[:4]}" if bad3 else ""))
    os.makedirs(RESULTS, exist_ok=True)
    name = "q_spin_check.txt" if N == 40 else f"q_spin_check_{N}.txt"
    with open(os.path.join(RESULTS, name), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
