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
  - **Proved at all levels (S4):** the kernel of restriction is the principal ideal $(\delta)$, with $\delta=S^+-S^-$ and $e=-\delta^2$. On conjugation classes this gives exactly the halved Euler moves above. So existence is a class-level question, and chirality is purely a selection question.
  - **Regge induction (S5, main line from 2026-10-11).** $A_n=T(A_{n-1})+E_n$, with $T:\lambda\mapsto\lambda+e_1$ and $E_n\ge0$ the new trajectory heads. It runs to $n=24$ for every tested rule, and continues from *every* enumerated uplift ($n\le11$). Conjecture: one step always works, which would prove existence at all $n$ by induction. Stress tests found no counterexample: the induction runs to $n=32$, and half (a) is certified for every uplift at $n\le24$. Step 2 (S5) proves a universal certificate lemma that reduces half (a) to an SO(9)-only inequality (a*); (a*) fails only at a few weights near the origin. There is a heuristic warning that the every-uplift form may fail at very large $n$, so the target is the existential induction (max-height rule). **S6:** half (a) of the induction is *proved* to follow from per-level Regge monotonicity, for any rule. Monotonicity is proved for the bosonic analogue: the cumulative bosonic spectrum is the coordinate ring of the arc space of the SO(9) null cone, a domain (via Mustaţă). The superstring needs a fermionic (super arc space) extension, since the naive factorisation through the bosonic algebra fails. **S7–S8:**
    - Natural super-null-cone arc-space models are ruled out at degree 2.
    - The trajectory heads are unique through $n=8$. They multiply ($E_5\supseteq\Lambda^216$, $E_6\supseteq16\otimes45$, $E_7\supseteq\mathrm{Sym}^245-1$), but form no finite algebra.
    - Two **secondary Regge directions** pass a null test: $T_2$ (spin $+e_1$, level $+2$) and $T_w$ (spin $+e_1+e_2$, level $+3$). Each forced induction runs to $n=32$ for every rule.
    - But $T_2$ and $T_w$ do not act freely together (strict obstruction at $n=9$), and a level-3 $e_1$ direction fails at $n=7$. So the naive two-vector-oscillator (Howe-dual) model is ruled out.
    - **S9, comparison with GMSZ (2609.04107):**
      - Our data reproduce their worked example, their support rule (all $7\le n\le32$) and their three shift relations, with zero violations.
      - Their shift relations are the stable-range equalities of monotone sequences that hold at *all* depths on our data. The monotone moves are "column moves" at triangular costs: $(\omega_1,+1)$ [S6], $(\omega_2,+3)$, $(\omega_3,+6)$, and a 4-form at $+10$.
      - $(\omega_2,+3)$ is the SO(9) shadow of $T_w$.
      - The vacuum exemption of $T_2$ is exactly their gap $[n-2,0,0,0]$.
      - The SO(10) heads sit just outside their $\ell$-stable range, at its edge.
      - Lead: Howe-duality "deep spectrum" constructions (Markou–Skvortsov; Basile–Markou) as a mechanism for proving the monotonicities.
    - **S10, Basile–Markou (2405.18467) read and checked:**
      - Their principal embedding (lowest level and multiplicity; R multiplicity $2^{K-1}$) holds for all 2 798 irreps first appearing in the full spectrum for $n\le32$. This includes $K=3,4$, beyond their own check.
      - Super-representations have a principal-embedding analogue: at GMSZ's first level, the multiplicity depends only on $(r,s)$.
      - But the paper gives no injective map for half (a). It is RNS, with no osp description of dividing out $44+84+128$, and multiplying by column minors does not preserve the constraints on dressed states.
    - **S11, half (b)** (liftability of $\rho_n$):
      - Local greedy lifts fail.
      - Lifts can always be taken with $l_5$ extremal ($\{0\text{ or }\tfrac12,l_4\}$). This is a property of the data, not of lifting in general.
      - An extremal lift is then an exact level-by-level recursion in $\mu_4$. At each level the tops are unique (finite differences); the only choice is which tops are degenerate.
      - Choosing them by a small per-level ILP with no lookahead succeeds for the max-height induction ($n\le32$), and for the cumulative $R^9_n$ directly ($n\le24$).
      - Above the bottom level only forced choices occur. All real repairs sit at the bottom transition, next to the symmetric tensors $[k,0,0,0,0]$ and the spinors $[k,0,0,0,1]$, i.e. next to the leading trajectory, as in half (a).
      - For the cumulative $R^9_n$ ($n\le24$) the only repairs are at the symmetric tensors $[k,0,0,0,0]$. Existence then reduces (empirically) to explicit finite-difference inequalities on $R^9_n$, with no Regge induction.
- [`bbms_lift/`](bbms_lift/STATUS.md): the BBMS SO(10)×SO(2)$_\Delta$-graded lift (hep-th/0310292) as a uniform-in-$n$ existence vehicle.
  - Reproduces BBMS's published levels 1–3, chirality included.
  - The grading map's kernel is the ideal generated by $K=\sum_k(-1)^k\Lambda^k(10)x^k$ (proved), and $K(1)=e$. So the grading fixes chirality, and leaves only Euler moves, which are fixed by the Δ-window.
  - In BBMS's unitarity window the graded lift is unique for all $n\le24$ and restricts correctly. **But it is non-negative only for $n\le18$** (B6): chirality-resolved negatives appear from $n=19$, and class-level negatives (singlet, 210) from $n=23$. So the window normal form is *not* an all-orders existence vehicle. Existence at $n\le24$ still holds via `so10_uplift/`. But chirality-resolved positive *graded* lifts do exist for all $n\le24$, at window width 10: the minimal fix is a few $K$-moves, mainly $e\cdot\overline{16}$, then $e\cdot[k,0,0,1,0]$ and $e\cdot1$. So the issue is the choice of normal form, not the grading. The literature-motivated alternative window, using Bianchi's $\min(\Delta_\pm)$, fails at $n=15$; no further window tuning is pursued.
  - The letter has a finite SO(10) preimage, the D=10 pure-spinor numerator.
  - **Theorem (B5).** A graded SO(8) character lifts if and only if it is SO(9)-covariant at $x=1$ and twisted-covariant at $x=-1$. The string levels satisfy both: at $x=-1$, $(-1)^\nu$ acts as the $e_9$ reflection. So the graded lift exists at all levels.

## Open / next

- No all-orders argument yet (see `uplift_structure/STATUS.md` for ideas).
- First-pass literature search: [`LITERATURE_SEARCH.md`](LITERATURE_SEARCH.md). Key lead: Beisert–Bianchi–Morales–Samtleben (hep-th/0310292) independently observed the III.B partial-sum SO(10) structure. They give an SO(10)×SO(2)$_\Delta$ grading, and an SO(10,2) match to free N=4 SYM (Pólya counting) up to $\Delta=10$.

## Minor code notes

- Assertion messages in `baseline_SO8_to_SO10.py` reference an undefined `i` (would raise NameError instead of the intended message if an assertion ever failed).
