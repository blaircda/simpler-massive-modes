"""
SO(9) -> SO(10) uplift of the cumulative Bars III.B sums
    R^9_n = sum_{i<=n} r^9_i
reading data/r9.json (produced by derive_r9.py).  Does not use any
baseline file.

For each level n it records
  * canonical uplift: maximise H = sum_lambda c_lambda * height(l1..l4)
    (height = 10,3,2,1 pairing); ties broken by lexicographic maximum over
    classes ordered by (height(l1..l4), l5) descending.  Whether the
    max-H face is a single point is certified separately (per-coordinate
    exact bounds);
  * further extremal uplifts: min height, max / min  sum c_lambda * l5;
  * greedy (max-l5 / min-l5) outcomes;
  * kernel dimension of the restriction map on admissible candidates;
  * exact number of uplifts (conjugation classes) for n <= --count-max;
  * all uplifts for n <= --enumerate-max (separate file).

Uplifts are recorded on conjugation classes (l5 >= 0 representative); a
class with l5 > 0 and multiplicity c may be split between lambda and its
conjugate in any way (c+1 choices), all with identical SO(9) restriction.
"""

import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from uplift import (canonical_uplift, count_uplifts, enumerate_uplifts, greedy, height,
                    is_unique_optimum, kernel_dimension, restrict, solve_ilp)
from weights import dim_irrep, dynkin_B4, dynkin_D5


def d5_list(sol):
    out = []
    for lam, c in sol.items():
        out.append({"mult": c, "ambient2": list(lam), "dynkin": list(dynkin_D5(lam)),
                    "dim": dim_irrep(lam, "D", 5), "conjugate_pair": lam[4] > 0})
    out.sort(key=lambda e: (sum(e["ambient2"]), e["ambient2"]))
    return out


def b4_list(m):
    out = [{"mult": c, "ambient2": list(mu), "dynkin": list(dynkin_B4(mu)),
            "dim": dim_irrep(mu, "B", 4)} for mu, c in m.items()]
    out.sort(key=lambda e: (sum(e["ambient2"]), e["ambient2"]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "data", "r9.json"))
    ap.add_argument("--max-level", type=int, default=None)
    ap.add_argument("--count-max", type=int, default=11)
    ap.add_argument("--enumerate-max", type=int, default=10)
    ap.add_argument("--out", default=os.path.join(HERE, "data", "R10.json"))
    ap.add_argument("--out-all", default=os.path.join(HERE, "data", "R10_all_uplifts.json"))
    args = ap.parse_args()

    levels_in = json.load(open(args.data))["levels"]
    N = max(int(k) for k in levels_in)
    if args.max_level:
        N = min(N, args.max_level)

    t0 = time.time()
    R = {}
    levels, all_out = {}, {}
    hobj = lambda lam: height(lam[:4])
    l5obj = lambda lam: lam[4]
    for n in range(1, N + 1):
        for e in levels_in[str(n)]["r9_B4"]:
            mu = tuple(e["ambient2"])
            R[mu] = R.get(mu, 0) + e["mult"]

        canon, best_h = canonical_uplift(R)
        if canon is None:
            print(f"level {n}: NO NON-NEGATIVE UPLIFT (obstruction)", flush=True)
            levels[str(n)] = {"R9_B4": b4_list(R), "feasible": False}
            continue
        unique = is_unique_optimum(R, hobj, canon, maximize=True)
        extremal = {
            "min_height": solve_ilp(R, hobj, maximize=False),
            "max_l5": solve_ilp(R, l5obj, maximize=True),
            "min_l5": solve_ilp(R, l5obj, maximize=False),
        }
        for name, s in list(extremal.items()) + [("canonical", canon)]:
            assert restrict(s) == R, f"{name} uplift at level {n} fails restriction"
        gmax, _ = greedy(R, "max")
        gmin, _ = greedy(R, "min")
        kdim, ncand = kernel_dimension(R)
        entry = {
            "feasible": True,
            "R9_B4": b4_list(R),
            "canonical_R10_D5": d5_list(canon),
            "canonical_height": best_h,
            "max_height_face_is_single_point": unique,
            "extremal_R10_D5": {k: d5_list(v) for k, v in extremal.items()},
            "num_classes_sum_c": sum(canon.values()),
            "greedy_max_l5_succeeds": gmax is not None,
            "greedy_min_l5_succeeds": gmin is not None,
            "num_candidates": ncand,
            "kernel_dim": kdim,
        }
        if n <= args.count_max:
            entry["num_uplifts_mod_chirality"] = count_uplifts(R)
        if n <= args.enumerate_max:
            sols = list(enumerate_uplifts(R))
            assert len(sols) == entry["num_uplifts_mod_chirality"]
            assert all(restrict(s) == R for s in sols)
            assert canon in sols
            assert all(sum(s.values()) == entry["num_classes_sum_c"] for s in sols)
            cl = sorted({lam for s in sols for lam in s})
            all_out[str(n)] = {"classes_ambient2": [list(l) for l in cl],
                               "classes_dynkin": [list(dynkin_D5(l)) for l in cl],
                               "uplifts": [[s.get(l, 0) for l in cl] for s in sols]}
        levels[str(n)] = entry
        print(f"level {n:2d}: #B4 irreps {len(R)}, canonical #classes {len(canon)}, "
              f"sum c = {entry['num_classes_sum_c']}, max-H unique = {unique}, "
              f"kernel dim {kdim}/{ncand}, "
              f"#uplifts = {entry.get('num_uplifts_mod_chirality', '-')}  "
              f"[{time.time()-t0:.1f}s]", flush=True)

    meta = {
        "description": "SO(10) uplifts R^10_n of cumulative SO(9) sums R^9_n = sum_{i<=n} r^9_i "
                       "(Bars hep-th/9503228 Sec. III.B).",
        "generator": "so10_uplift/uplift_so10.py",
        "input": "so10_uplift/data/r9.json",
        "max_level": N,
        "conventions": "ambient2 = doubled orthogonal highest weight with l5 >= 0 "
                       "(conjugation-class representative); dynkin = Sage "
                       "WeylCharacterRing('D5', style='coroots'). For conjugate_pair "
                       "classes the multiplicity may be split freely between lambda and "
                       "its conjugate. Canonical = maximiser of H = sum c * (10 l1 + 3 l2 + 2 l3 + l4), "
                       "ties broken by lex-max over classes ordered by (height, l5) "
                       "descending. The canonical uplift equals the hardcoded baseline "
                       "for n <= 10 (checked in verify_R10.py, not used here).",
    }
    with open(args.out, "w") as f:
        json.dump({"meta": meta, "levels": levels}, f, indent=1)
    if all_out:
        meta_all = dict(meta, description="All non-negative SO(10) uplifts (conjugation "
                        "classes) for low levels; each uplift is a multiplicity vector over "
                        "classes_ambient2.")
        with open(args.out_all, "w") as f:
            json.dump({"meta": meta_all, "levels": all_out}, f)
    print(f"wrote {args.out}  [{time.time()-t0:.1f}s]")


if __name__ == "__main__":
    main()
