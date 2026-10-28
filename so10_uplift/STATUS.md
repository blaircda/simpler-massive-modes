# so10_uplift — independent SO(9) → SO(10) uplift (Sage-based)

## Pipeline

| step | script | output | runtime |
|---|---|---|---|
| Fock space → $r^8_n$ → $r^9_n$ | `derive_r9.py --max-level 24` | `data/r9.json` | ~30 s |
| $R^9_n=\sum_{i\le n} r^9_i$ → SO(10) uplifts | `uplift_so10.py` | `data/R10.json`, `data/R10_all_uplifts.json` | ~5.5 min |
| check vs baseline + Sage branching | `verify_r9_against_baseline.py`, `verify_R10.py` | — | ~30 s, ~3 min |
| loader for later scripts | `repdata.py` (`RepData().r9(n)`, `.R9(n)`, `.R10(n)`, `.R10_all(n)`, …) | — | — |

Library code: `weights.py` (dominant-weight tables, Weyl alternation, dims), `uplift.py` (interlacing branching, candidates, ILP, greedy, exact enumeration, canonical uplift).

Derivation scripts never import baseline files; comparisons live only in the `verify_*` scripts.

## Method

- **Fock space:** exp–log recurrence $nF_n=\sum_j G_jF_{n-j}$, where $G_j=\sum_{k|j}(j/k)[\psi^k 8_v-(-1)^k\psi^k 8_s]$ and $\psi^k$ are Adams operations. It runs on W(D4)-dominant weight multiplicities. The dimension is checked against $\prod_m((1+q^m)/(1-q^m))^8$.
- **Division by $8_v+8_s$:** exact triangular division, ordered by the height $h=(10,3,2,1)$ under which $e_1$ is the unique top weight of $8_v+8_s$. Checked by re-multiplying.
- **SO(9) step:** Spin(8) and Spin(9) share a maximal torus, so $r^9_n$ has the same weights as $r^8_n$. SO(9) covariance is the weight symmetry $x_4\to-x_4$, which is checked. The decomposition is then unique, via Weyl alternation. No greedy search is needed.
- **SO(10) step:** branching is multiplicity-free interlacing, $l_1\ge\mu_1\ge l_2\ge\mu_2\ge l_3\ge\mu_3\ge l_4\ge\mu_4\ge|l_5|$. Uplifts are posed on conjugation classes ($l_5\ge 0$). Existence of a non-negative uplift is decided by an exact ILP (Sage/GLPK), and every solution is re-verified by restriction.

## Results (2026-09-30)

1. $r^8_n$ and $r^9_n$ agree with the baseline for $n\le10$.
2. **A non-negative SO(10) uplift of $R^9_n$ exists for every $n\le 24$.** No obstruction was found.

   > **Update (2026-10-12): existence now holds for every $n\le40$.** This is shown by the Regge induction in [`../uplift_structure/`](../uplift_structure/STATUS.md) (S12, `beyond32.py`), under both the max-height and min-height rules, with every $A_n$ verified by restriction. Its input is `derive_r9.py --max-level 40`, stored as `uplift_structure/data/r9_to40.json`. The uplifts for $n>24$ are recomputed (~2.5 min), not stored here. The $n\le24$ figure above is the range of this directory's stored uplifts.
3. **Non-uniqueness.** The number of non-negative uplifts (conjugation classes, before chirality splits) is exactly:

   | n | 1–6 | 7 | 8 | 9 | 10 | 11 |
   |---|---|---|---|---|---|---|
   | #uplifts | 1 | 2 | 8 | 72 | 6057 | 6 722 076 |

   The affine dimension of the solution set, taken over candidate classes, is:

   | n | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 | 22 | 24 |
   |---|---|---|---|---|---|---|---|---|---|---|
   | dim | 1 | 3 | 8 | 16 | 30 | 52 | 83 | 129 | 190 | 271 |

   Kernel sources include $\lambda\leftrightarrow\bar\lambda$ and $\sum_k(-1)^k\Lambda^k(10)\otimes V$. For example, $2(1+45+210)$ and $2(10+120)+126+\overline{126}$ restrict identically.
4. **Canonical uplift.** The hardcoded baseline $R^{10}_n$ ($n\le10$) is the *unique* maximiser of $H=\sum_\lambda c_\lambda(10l_1+3l_2+2l_3+l_4)$ ("most top-heavy"). This holds for $n\le 11$. From $n=12$ onward the max-$H$ face is not a point. The canonical choice stored is therefore max-$H$ first, then the lex-max over classes ordered by $(\mathrm{height}, l_5)$ descending. This reproduces the baseline exactly for $n\le10$.
5. **Greedy lifts fail.** Top-down greedy with $l_5$ maximal fails for every $n\ge 7$. With $l_5$ minimal it fails at $n=6$ and for every $n\ge8$. A constructive proof therefore cannot be this naive greedy.

Structural explorations of this output (invariants, fibres, kernel moves) live separately in [`../uplift_structure/`](../uplift_structure/STATUS.md).

## Open / next (computational)

- Push the stored uplifts beyond level 24. (Existence alone is done to $n=40$; see the update under item 2.) Fock + SO(9) is cheap. The canonical-uplift sequential ILP is the bottleneck (~110 s at n=24, growing ×1.5/level). Existence alone (single ILP) is much cheaper.
- A smarter greedy with provable success (item 5 rules out the naive ones).
