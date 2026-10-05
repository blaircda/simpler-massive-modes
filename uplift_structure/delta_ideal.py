"""
S4: the restriction kernel R(Spin10) -> R(Spin9) is the principal ideal (delta),
delta = S+ - S-, and the Euler class is e = -delta^2.

The proof is algebraic (see STATUS.md).  This script checks the concrete
identities used there:
 * delta and e restrict to zero on B4;
 * e = sum_k (-1)^k Lambda^k(10) equals -delta^2;
 * Lambda^k(10) restricts to Lambda^k(9) + Lambda^(k-1)(9) (triangularity);
 * S+ and S- both restrict to the B4 spinor;
 * the level-6 BBMS-vs-baseline chirality difference reported in
   previous_efforts.md equals delta * (10 + 2), and vanishes on
   conjugation classes;
 * on sample D5 irreps, V - Vbar is outer-anti-invariant and restricts
   to zero, and delta * (V - Vbar) is outer-invariant;
 * the key step for the conjugation-class lattice: with pi the merge map
   to conjugation classes and z sigma-invariant,
       pi(delta * S+ * z) = -1/2 pi(e * z),
   so that pi(e V_nu)/2 (l5 = 0) and pi(e V_nu) (l5 > 0) lie in
   pi((delta)) = ker of the class-level restriction (the S3b generators).
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS

from sage.all import WeylCharacterRing

D5 = WeylCharacterRing("D5", style="coroots")
B4 = WeylCharacterRing("B4", style="coroots")
ONE = D5(0, 0, 0, 0, 0)
V = D5(1, 0, 0, 0, 0)
Sp = D5(0, 0, 0, 1, 0)
Sm = D5(0, 0, 0, 0, 1)
delta = Sp - Sm


def res(x):
    return x.branch(B4, rule="symmetric")


def ext(rep, k):
    return rep.exterior_power(k) if k > 0 else rep.parent().one()


def sigma(x):
    """Outer automorphism of D5: flip the sign of the last ambient coordinate."""
    out = D5.zero()
    for w, m in x.monomial_coefficients().items():
        c = list(w.to_vector())
        out += m * D5(D5.space()(c[:4] + [-c[4]]))
    return out


def merge(x):
    """pi: D5 character -> {conjugation class (ambient weight, last coord >= 0): coeff}."""
    out = {}
    for w, m in x.monomial_coefficients().items():
        c = list(w.to_vector())
        c[4] = abs(c[4])
        k = tuple(c)
        out[k] = out.get(k, 0) + m
    return {k: v for k, v in out.items() if v != 0}


def main():
    lines = []

    def log(s):
        print(s)
        lines.append(s)

    V9 = B4(1, 0, 0, 0)
    log("delta restricts to 0: %s" % (res(delta) == 0))
    e = sum((-1) ** k * ext(V, k) for k in range(11))
    log("e = sum (-1)^k Lambda^k(10) = %s" % e)
    log("e == -delta^2: %s" % (e == -delta * delta))
    log("e restricts to 0: %s" % (res(e) == 0))
    for k in range(1, 5):
        ok = res(ext(V, k)) == ext(V9, k) + ext(V9, k - 1)
        log("Lambda^%d(10)|B4 == Lambda^%d(9) + Lambda^%d(9): %s" % (k, k, k - 1, ok))
    log("S+|B4 == S-|B4 == spinor 16: %s"
        % (res(Sp) == B4(0, 0, 0, 1) and res(Sm) == B4(0, 0, 0, 1)))

    diff = Sp + D5(1, 0, 0, 1, 0) - Sm - D5(1, 0, 0, 0, 1)
    log("level-6 BBMS - baseline difference == delta*(10+2): %s" % (diff == delta * (V + 2 * ONE)))
    log("  ... and sigma-merged it is zero: %s" % (diff + sigma(diff) == 0))

    samples = [(0, 0, 0, 1, 0), (1, 0, 0, 1, 0), (0, 0, 0, 2, 0), (0, 1, 0, 1, 0), (1, 0, 0, 3, 0)]
    for lab in samples:
        x = D5(*lab)
        anti = x - sigma(x)
        log("V%s - sigma: anti-invariant %s; delta*(V - Vbar) sigma-invariant %s; restricts to 0 %s"
            % (lab, sigma(anti) == -anti, sigma(delta * anti) == delta * anti, res(anti) == 0))

    e_ = -delta * delta
    zs = [(0, 0, 0, 0, 0), (1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (2, 0, 1, 0, 0),
          (0, 0, 0, 1, 0), (1, 0, 0, 1, 0), (0, 0, 1, 2, 0), (1, 1, 0, 0, 3)]
    for lab in zs:
        v = D5(*lab)
        z = v if lab[3] == lab[4] else v + sigma(v)
        lhs = merge(delta * Sp * z)
        rhs = {k: -c / 2 for k, c in merge(e_ * z).items()}
        even = all(c % 2 == 0 for c in merge(e_ * v).values())
        log("z from V%s: pi(delta S+ z) == -pi(e z)/2: %s; pi(e V) even: %s"
            % (lab, lhs == rhs, even))

    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "delta_ideal.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
