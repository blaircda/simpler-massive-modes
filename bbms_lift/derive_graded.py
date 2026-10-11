"""
BBMS-graded SO(10) x SO(2)_Delta lift of the chiral string levels.

Reference: Beisert, Bianchi, Morales, Samtleben, hep-th/0310292, sec. 2
(and Bianchi, hep-th/0409304, eqs. 3.13-3.17).

1. Light-cone Fock character refined by the oscillator number nu (every
   alpha_{-n}^I and S_{-n}^a counts 1), via the graded exp-log recurrence
        n F_n(x) = sum_j G_j(x) F_{n-j}(x),
        G_j(x)  = sum_{k|j} (j/k) x^k [psi^k 8_v - (-1)^k psi^k 8_s].
2. r_n(x) = F_n(x) / (8_v + 8_s), slice by slice in nu (no shift of nu),
   so that e.g. r_2 = 1.x + 8_v.x^2, matching BBMS's vac_2.
3. Solve exactly over Q for an integer combination
        X_n = sum_{lambda, Delta} c(lambda, Delta) V_lambda x^Delta
   of D5 irreps (both chiralities) such that, on branching
   SO(10) -> SO(8) x SO(2)_J and setting nu = Delta - J,
        phi(X_n) = r_n(x)        (as SO(8) x nu characters).
   SO(2)_J rotates the last ambient coordinate x5; SO(8) acts on x1..x4,
   the SO(8) inside the interlacing SO(9) of so10_uplift/.
   Unknown domain: lambda among the class candidates of R^9_n or R^9_{n-1}
   (both signs of l5); Delta with [Delta - l1, Delta + l1] meeting the
   nu-support of r_n, widened by --slack on each side.
4. Report the nullity of the graded system.  If it is zero, the graded
   solution is unique in the domain; otherwise the level is flagged and the
   solution returned by the solver is an arbitrary representative.
   Linear algebra (graded.solve_level): the system is block diagonal in 4
   blocks (central character x gap parity); each block is echelonized as a
   dense GF(p) matrix (word-size primes), residues are CRT-lifted to
   symmetric integers, and the solution is checked exactly against r_n.
   Forget Delta, form cumulative sums R_n = X_1 + ... + X_n, and check
   (i) non-negativity and (ii) R_n|_{B4} = R^9_n (computed here).

Positive (not super-) counting: fermionic states carry + signs.  The map
preserves the spinor/tensor class, so this is BBMS's supercharacter with
spinorial coefficients negated.

Output: data/bbms_graded.json (rewritten after every level; --resume continues
from the last completed level of an existing output file).
"""

import argparse
import json
import os
import resource
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from graded import HERE, dynkin_D5, levels_data, solve_level
from uplift import restrict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-level", type=int, default=7)
    ap.add_argument("--slack", type=int, default=1, help="extra Delta range on each side")
    ap.add_argument("--gmax", type=int, default=None,
                    help="use the unitarity window 0 <= Delta - bound(lambda) <= gmax instead")
    ap.add_argument("--bound", choices=["plus", "min"], default="plus",
                    help="window lower edge: BBMS Delta_+ (default) or min(Delta_+, Delta_-)")
    ap.add_argument("--out", default=os.path.join(HERE, "data", "bbms_graded.json"))
    ap.add_argument("--resume", action="store_true",
                    help="continue after the last completed level stored in --out")
    args = ap.parse_args()
    N = args.max_level
    t0 = time.time()

    r, R9, T = levels_data(N)
    print(f"graded Fock space and r_n(x) through level {N} [{time.time()-t0:.1f}s]", flush=True)

    meta = {"generator": "bbms_lift/derive_graded.py", "max_level": N,
            "slack": args.slack, "gmax": args.gmax, "bound": args.bound,
            "conventions": "ambient2 doubled D5 weights; J = x5; SO(8) on x1..x4; "
                           "nu = Delta - J; Delta2 = 2*Delta; positive counting. "
                           "If unique_in_domain is false the stored solution is an "
                           "arbitrary representative."}

    def save():
        out = {"meta": meta, "levels": levels}
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        tmp = args.out + ".tmp"
        with open(tmp, "w") as f:
            json.dump(out, f, indent=1)
        os.replace(tmp, args.out)

    levels, cum, start = {}, {}, 1
    if args.resume and os.path.exists(args.out):
        old = json.load(open(args.out))
        assert (old["meta"]["gmax"] == args.gmax and old["meta"]["slack"] == args.slack
                and old["meta"].get("bound", "plus") == args.bound), \
            "resume: domain parameters differ from the stored run"
        done = [int(k) for k, v in old["levels"].items() if "R_cumulative_D5" in v]
        last = 0
        while last + 1 in done:
            last += 1
        levels = {str(k): old["levels"][str(k)] for k in range(1, last + 1)}
        if last:
            cum = {tuple(e["ambient2"]): e["mult"]
                   for e in old["levels"][str(last)]["R_cumulative_D5"]}
        start = last + 1
        print(f"resuming after level {last}", flush=True)

    for n in range(start, N + 1):
        tn = time.time()
        nullity, x, stats = solve_level(n, r, R9, T, args.gmax, slack=args.slack,
                                        log=lambda s: print(s, flush=True), bound=args.bound)
        solvable = x is not None
        ncols = sum(v.get("cols", 0) for v in stats.values() if isinstance(v, dict))
        rec = {"n_unknowns": ncols, "blocks": stats, "nullity_modp": nullity,
               "solvable": solvable, "unique_in_domain": nullity == 0,
               "seconds": round(time.time() - tn, 1)}
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
        print(f"level {n}: unknowns {ncols}, nullity {nullity}, solvable {solvable} "
              f"[{time.time()-tn:.1f}s, total {time.time()-t0:.1f}s, peak RSS {rss:.2f} GB]",
              flush=True)
        if not solvable:
            levels[str(n)] = rec
            save()
            print("   no graded solution in this domain; stopping cumulative sums")
            break
        graded = sorted((col, int(v)) for col, v in x.items())
        X = {}
        for (lam, D2), k in graded:
            X[lam] = X.get(lam, 0) + k
        X = {lam: k for lam, k in X.items() if k}
        for lam, k in X.items():
            cum[lam] = cum.get(lam, 0) + k
        cum = {lam: k for lam, k in cum.items() if k}
        res_ok = restrict(cum) == {mu: k for mu, k in R9[n].items() if k}
        nonneg = min(cum.values()) >= 0

        def entries(d):
            return [{"ambient2": list(lam), "dynkin": list(dynkin_D5(lam)), "mult": k}
                    for lam, k in sorted(d.items())]

        rec.update({
            "graded": [{"ambient2": list(lam), "dynkin": list(dynkin_D5(lam)),
                        "Delta2": D2, "mult": k} for (lam, D2), k in graded],
            "X_D5": entries(X), "R_cumulative_D5": entries(cum),
            "restriction_ok": res_ok, "cumulative_nonneg": nonneg,
        })
        print(f"   cumulative: restriction to B4 ok {res_ok}, non-negative {nonneg}, "
              f"#D5 irreps {len(cum)}", flush=True)
        levels[str(n)] = rec
        save()
    print(f"wrote {args.out} [{time.time()-t0:.1f}s]")


if __name__ == "__main__":
    main()
