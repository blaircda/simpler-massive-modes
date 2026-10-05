"""
Exact identities for the grading map
    phi : R(Spin10)[x^{1/2}] -> R(Spin8)[x^{1/2}],   V_lam x^Delta |-> x^Delta Res(V_lam)(y = 1/x),
(SO(2)_J fugacity y, nu = Delta - J), using graded.phi.

 (a) phi(K) = 0 for K = lambda_{-x}(10) = sum_k (-1)^k Lambda^k(10) x^k
     (graded Euler class; K(x=1) = e).
 (b) phi(10 x - 1 - x^2) = 8_v x                  (bosonic letter)
 (c) phi(N) = x (x^2 - 1) (8_v - 8_s)            (supersymmetric letter)
     with N = 1 - 10 x + 16 x^{3/2} - 16bar x^{5/2} + 10 x^3 - x^4,
     16 = D5(0,0,0,0,1), 16bar = D5(0,0,0,1,0).
     N(t^2) is the numerator of the D = 10 pure-spinor partition function.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from graded import HERE, phi

from sage.all import WeylCharacterRing

D5 = WeylCharacterRing("D5", style="coroots")


def amb(dyn):
    fw = D5.fundamental_weights()
    return tuple(int(2 * a) for a in sum(c * fw[i + 1] for i, c in enumerate(dyn)).to_vector())


def image(terms):
    """terms: list of (Dynkin label, Delta2, coeff) -> phi image {(nu, w4): m}."""
    out = {}
    for dyn, D2, c in terms:
        for key, m in phi(amb(dyn), D2).items():
            out[key] = out.get(key, 0) + c * m
    return {k: v for k, v in out.items() if v}


def target(terms):
    """terms: list of (D4 dominant doubled weight, nu, coeff)."""
    out = {}
    for w, nu, c in terms:
        out[(nu, w)] = out.get((nu, w), 0) + c
    return {k: v for k, v in out.items() if v}


V8, S8, C8, ONE = (2, 0, 0, 0), (1, 1, 1, 1), (1, 1, 1, -1), (0, 0, 0, 0)
# dominant-weight multiplicity functions: 8v = {(2,0,0,0):1}; 8s = {(1,1,1,1):1}


def main():
    LAM = {0: [(0, 0, 0, 0, 0)], 1: [(1, 0, 0, 0, 0)], 2: [(0, 1, 0, 0, 0)], 3: [(0, 0, 1, 0, 0)],
           4: [(0, 0, 0, 1, 1)], 5: [(0, 0, 0, 2, 0), (0, 0, 0, 0, 2)]}
    K = []
    for k in range(11):
        for dyn in LAM[min(k, 10 - k)]:
            K.append((dyn, 2 * k, (-1) ** k))
    lines = []
    ok_a = image(K) == {}
    ok_b = image([((1, 0, 0, 0, 0), 2, 1), ((0, 0, 0, 0, 0), 0, -1), ((0, 0, 0, 0, 0), 4, -1)]) \
        == target([(V8, 1, 1)])
    N = [((0, 0, 0, 0, 0), 0, 1), ((1, 0, 0, 0, 0), 2, -1), ((0, 0, 0, 0, 1), 3, 1),
         ((0, 0, 0, 1, 0), 5, -1), ((1, 0, 0, 0, 0), 6, 1), ((0, 0, 0, 0, 0), 8, -1)]
    ok_c = image(N) == target([(V8, 3, 1), (S8, 3, -1), (V8, 1, -1), (S8, 1, 1)])
    for s in (f"(a) phi(lambda_-x(10)) = 0: {ok_a}",
              f"(b) phi(10x - 1 - x^2) = 8v x: {ok_b}",
              f"(c) phi(1 - 10x + 16 x^3/2 - 16bar x^5/2 + 10x^3 - x^4) = x(x^2-1)(8v - 8s): {ok_c}"):
        print(s)
        lines.append(s)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", "letter_identities.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    assert ok_a and ok_b and ok_c


if __name__ == "__main__":
    main()
