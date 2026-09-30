# STATUS

Goal: all-orders existence/admissibility proof for the minimal cumulative Bars III.B scheme (hep-th/9503228): for every level $n$ there is a genuine (non-negative) SO(10) character $R^{10}_n$ with

$$R^{10}_n \downarrow_{SO(9)} = R^9_n := \sum_{i=1}^{n} r^9_i ,$$

where the physical IIA level-$n$ (one-sided) states are $r^9_n \otimes (44+84+128)$.

## Established (as of 2026-09-30)

Pipeline `baseline_SO8_to_SO10.py` (runtime ~55 s, passes):

1. `baseline_SO8_construct.py`: SO(8) light-cone Fock space $F_n$ through level 10, then division by the SUSY factor $8_v+8_s$ to get $r^8_n$ (exact reconstruction asserted).
2. `baseline_SO8_to_SO9.py`: greedy SO(8)→SO(9) lift $r^8_n \to r^9_n$ (branching asserted equal).
3. `baseline_hardcoded_level10_SO10.py`: previously found $R^{10}_n$, $n=1..10$ (**baseline data**; derivation scripts must not import it). Verified: $R^{10}_n\downarrow B_4 = R^9_n$ for $n\le 10$.

## Research directions

- [`so10_uplift/`](so10_uplift/STATUS.md): independent Sage-based pipeline Fock → $r^9_n$ → SO(10) uplifts, through level 24. This is the hard computational evidence. Stored data: `so10_uplift/data/{r9,R10,R10_all_uplifts}.json`; loader `so10_uplift/repdata.py`.
  - A non-negative uplift exists for all $n\le24$.
  - The uplift is unique for $n\le6$, then has rapidly growing multiplicity (6 722 076 at $n=11$).
  - The baseline $R^{10}_n$ ($n\le10$) is the unique "max-height" uplift. That uplift stays unique through $n=11$ but not from $n=12$; the stored canonical choice adds a lex tie-break.
- [`uplift_structure/`](uplift_structure/STATUS.md): exploratory analysis of that output (reads it, never modifies it).
  - The SO(10) irrep count is uplift-invariant; the $(l_2,l_4)$ fibre sizes are fixed by $R^9_n$.
  - Fibres are coupled, but restriction is injective within each fibre.
  - The canonical uplift maximises $\sum c\,l_1$.
  - The uplift kernel is generated over $\mathbb Z$ by Euler-class moves $e\otimes V_\nu$ (halved when $l_5(\nu)=0$, checked $n\le16$), with $e=\sum(-1)^k\Lambda^k(10)$.

## Open / next

- No all-orders argument yet (see `uplift_structure/STATUS.md` for ideas).

## Minor code notes

- Assertion messages in `baseline_SO8_to_SO10.py` reference an undefined `i` (would raise NameError instead of the intended message if an assertion ever failed).
