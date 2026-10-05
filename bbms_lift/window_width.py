"""
B3: nullity of the graded system versus the unitarity-window width G
(domain 0 <= Delta - b(lambda) <= G), at a few levels.  G = 8 (equivalently
9, gaps being even) is the widest kernel-free window; from G = 10 multiples
of K = lambda_{-x}(10) enter.

Output: results/window_width.txt
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from graded import HERE, build_system, levels_data, solve_exact

LEVELS = (4, 6, 10, 13, 14)
WIDTHS = (8, 10, 12, 16, 20)


def main():
    r, R9, T = levels_data(max(LEVELS))
    lines = []
    for n in LEVELS:
        out = []
        for G in WIDTHS:
            e, b, cols, rows = build_system(n, r, R9, T, 0, gmax=G)
            nl, x = solve_exact(e, b, len(rows), len(cols))
            out.append(f"G={G}: nullity {nl}, solvable {x is not None}")
        s = f"level {n:2d}: " + "; ".join(out)
        print(s, flush=True)
        lines.append(s)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", "window_width.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
