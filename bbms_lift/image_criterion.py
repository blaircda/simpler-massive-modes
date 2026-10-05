"""
B5: which graded SO(8) characters lift to SO(10) x SO(2)_Delta.

Theorem (proof in STATUS.md): f in R(Spin8)[x^{+-1}] lies in phi(R(Spin10)[x^{+-1/2}])
iff
  (x = +1)  f(1) is W(B4)-invariant  (an SO(9) character), and
  (x = -1)  f(-1) has W(B4)-invariant tensor part and W(B4)-anti-invariant
            spinor part (the twisted character of the reflection of e9).

This script checks:
 (1) the generators used in the sufficiency proof:
        phi(10 + 1 - x^{-1} - x)             = 8v + 1 = Res 9
        phi(16 x^{1/2} + 16bar x^{1/2})      = (1 + x)(8s + 8c)
        phi(16 x^{1/2} - 16bar x^{1/2})      = (1 - x)(8s - 8c)
 (2) both conditions for the graded string levels r_n(x), n <= NMAX,
     i.e. that every BBMS graded lift X_n exists.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from graded import HERE, levels_data, phi
from letter_identities import amb, image, target

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
V8, S8, C8, ONE = (2, 0, 0, 0), (1, 1, 1, 1), (1, 1, 1, -1), (0, 0, 0, 0)


def main():
    lines = []

    def log(s):
        print(s)
        lines.append(s)

    g9 = image([((1, 0, 0, 0, 0), 0, 1), ((0, 0, 0, 0, 0), 0, 1),
                ((0, 0, 0, 0, 0), -2, -1), ((0, 0, 0, 0, 0), 2, -1)])
    log(f"phi(10 + 1 - 1/x - x) = 8v + 1: {g9 == target([(V8, 0, 1), (ONE, 0, 1)])}")
    sp = image([((0, 0, 0, 0, 1), 1, 1), ((0, 0, 0, 1, 0), 1, 1)])
    log(f"phi(16 x^1/2 + 16bar x^1/2) = (1+x)(8s+8c): "
        f"{sp == target([(S8, 0, 1), (C8, 0, 1), (S8, 1, 1), (C8, 1, 1)])}")
    sm = image([((0, 0, 0, 0, 1), 1, 1), ((0, 0, 0, 1, 0), 1, -1)])
    log(f"phi(16 x^1/2 - 16bar x^1/2) = (1-x)(8s-8c): "
        f"{sm == target([(S8, 0, 1), (C8, 0, -1), (S8, 1, -1), (C8, 1, 1)])}")

    r, R9, T = levels_data(NMAX)
    flip = T.lookup(T.W * np.array([1, 1, 1, -1]))
    spin = T.W[:, 0] % 2 == 1
    allok = True
    for n in range(1, NMAX + 1):
        res = []
        for sgn in (1, -1):
            f = sum((sgn ** nu) * a for nu, a in r[n].items())[:T.size]
            ft, fs = np.where(spin, 0, f), np.where(spin, f, 0)
            ok = np.array_equal(ft, ft[flip]) and np.array_equal(fs, sgn * fs[flip])
            res.append(ok)
            allok &= ok
        log(f"level {n:2d}: r_n(+1) SO(9)-covariant {res[0]}; r_n(-1) twisted-covariant {res[1]}")
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", "image_criterion.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    assert allok


if __name__ == "__main__":
    main()
