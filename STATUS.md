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

- [`so10_uplift/`](so10_uplift/STATUS.md): independent Sage-based pipeline Fock → $r^9_n$ → SO(10) uplifts, with stored uplifts through level 24 ($r^9_n$ derived to $n=40$; existence to $n=40$ is in `uplift_structure/` S12). This is the hard computational evidence. Stored data: `so10_uplift/data/{r9,R10,R10_all_uplifts}.json`; loader `so10_uplift/repdata.py`.
  - A non-negative uplift exists for all $n\le24$ (extended to $n\le40$ by the Regge induction, `uplift_structure/` S12).
  - The uplift is unique for $n\le6$, then has rapidly growing multiplicity (6 722 076 at $n=11$).
  - The baseline $R^{10}_n$ ($n\le10$) is the unique "max-height" uplift. That uplift stays unique through $n=11$ but not from $n=12$; the stored canonical choice adds a lex tie-break.
- [`uplift_structure/`](uplift_structure/STATUS.md): exploratory analysis of that output (reads it, never modifies it).
  - The SO(10) irrep count is uplift-invariant; the $(l_2,l_4)$ fibre sizes are fixed by $R^9_n$.
  - Fibres are coupled, but restriction is injective within each fibre.
  - The canonical uplift maximises $\sum c\,l_1$.
  - The uplift kernel is generated over $\mathbb Z$ by Euler-class moves $e\otimes V_\nu$ (halved when $l_5(\nu)=0$, checked $n\le16$), with $e=\sum(-1)^k\Lambda^k(10)$.
  - **Proved at all levels (S4):** the kernel of restriction is the principal ideal $(\delta)$, with $\delta=S^+-S^-$ and $e=-\delta^2$. On conjugation classes this gives exactly the halved Euler moves above. So existence is a class-level question, and chirality is purely a selection question.
  - **Regge induction (S5, main line from 2026-10-11).** $A_n=T(A_{n-1})+E_n$, with $T:\lambda\mapsto\lambda+e_1$ and $E_n\ge0$ the new trajectory heads. It runs to $n=24$ for every tested rule, and continues from *every* enumerated uplift ($n\le11$). Conjecture: one step always works, which would prove existence at all $n$ by induction. Stress tests found no counterexample: the induction runs to $n=32$ (to $n=40$ with the new data, under both rules; S12), and half (a) is certified for every uplift at $n\le24$. Step 2 (S5) proves a universal certificate lemma that reduces half (a) to an SO(9)-only inequality (a*); (a*) fails only at a few weights near the origin. There is a heuristic warning that the every-uplift form may fail at very large $n$, so the target is the existential induction (max-height rule). **S6:** half (a) of the induction is *proved* to follow from per-level Regge monotonicity, for any rule. Monotonicity is proved for the bosonic analogue: the cumulative bosonic spectrum is the coordinate ring of the arc space of the SO(9) null cone, a domain (via Mustaţă). The superstring needs a fermionic (super arc space) extension, since the naive factorisation through the bosonic algebra fails. **S7–S8:**
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
    - **S11, half (b)** (liftability of $\rho_n$; **corrected**: the $l_5$-extremal structure holds only for $n\le33$, see S12):
      - Local greedy lifts fail.
      - Lifts can always be taken with $l_5$ extremal ($\{0\text{ or }\tfrac12,l_4\}$). This is a property of the data, not of lifting in general.
      - An extremal lift is then an exact level-by-level recursion in $\mu_4$. At each level the tops are unique (finite differences); the only choice is which tops are degenerate.
      - Choosing them by a small per-level ILP with no lookahead succeeds for the max-height induction ($n\le32$), and for the cumulative $R^9_n$ directly ($n\le24$).
      - Above the bottom level only forced choices occur. All real repairs sit at the bottom transition, next to the symmetric tensors $[k,0,0,0,0]$ and the spinors $[k,0,0,0,1]$, i.e. next to the leading trajectory, as in half (a).
      - For the cumulative $R^9_n$ ($n\le24$) the only repairs are at the symmetric tensors $[k,0,0,0,0]$. Existence then reduces (empirically) to explicit finite-difference inequalities on $R^9_n$, with no Regge induction.
      - **Correction (S12, data to $n=40$).** From $n=34$ no $l_5$-extremal lift of $\rho_n$ exists. Lifts need classes $(l_1,l_2,2,2,1)$, so the repairs climb to the next level. The pattern held for $n\le33$ and then changed.
    - **S12, the cumulative route and existence to $n=40$:**
      - Lifting $R^9_n$ directly (all-full lift $\Delta^4R$ plus bottom repairs) works for $n\le33$, with stable-range minimum slack and singlet-sized chain margins.
      - It breaks at $n=34$: the singlet node is drained, then forced coverage fails at column $(5,1,1)$.
      - With new data `r9_to40.json`, **the Regge induction runs to $n=40$ under both rules**, so non-negative uplifts exist for all $n\le40$.
    - **S13, repair levels of half (b):**
      - Onsets $n_1=6$, $n_2=34$ (31 under the min-height rule). No level 3 is needed by $n=44$, and repairs never skip more than one bottom level.
      - A repair level $c+1$ becomes necessary when a bottom column $(k,c,c)$ has more deficit than its same-column capacity. This is local; aggregate capacity is not the trigger.
      - At fixed weights the all-full lift is asymptotically positive ($D_4(\dim)\ge70$). So the repairs are pre-asymptotic, in a transition region between the stable range and fixed weights.
      - Level 3 extrapolates to $n\sim100$. A proof needs uniform asymptotics in that region plus the exact stable range.
    - **S14, the transition region:**
      - At $n\le44$ the multiplicities fall roughly exponentially in spin; there is no dimension-proportional regime yet.
      - The negative all-full coefficients fill the bulk of each column family, up to a fixed distance below the top.
      - The level-2 band is $k\in[4,n-30]$, in the family $(1,1)$ only.
      - So the region is the whole Regge-scaling middle, between fixed edges. A proof needs uniform asymptotics there, where GMSZ-type cost and shift structure is the natural tool.
    - **S15, the Regge limit in closed form** (`regge_limit.py` reads no data; `regge_scaling.py` checks it against the data to $n=44$):
      - In GMSZ's cost coordinates ($D=n-x_1$, depth $w=D-{\rm cost}_7(\mu)$), the large-spin limit of $r^9_n$ is a positive SO(7) Fock space. It is derived from the light cone; the 256 cancels exactly. It has bosons $1+7$ at levels $\ge2$ and spinor fermions $8$ at levels $\ge\tfrac32$.
      - The data never exceed it, and equal it exactly whenever $\ell\ge w$ (321 101 comparisons, zero exceptions). GMSZ's support costs $c(r,s)$ and their gap are its lowest energies.
      - It reproduces all of S14's offsets and the band's upper edge $n-30$ as $n$-independent quantities. It predicts a second level-2 family $(2,1)$ for $D\ge50$, and no forced level 3 near the top up to $D=60$.
      - The model is the restriction of the SO(8) Fock space ($8_v$ bosons, $8_s$ fermions), so in the Regge limit the uplift exists in closed form. At finite $n$ those classes overshoot in the ramp $\ell<w$. What remains of half (b) is when each class is created, which the ramp governs (corrected 2026-10-13: the first version said the uplift near the top was settled).
    - **S16, the ramp** (`q_spin.py`, `q_top_end.py` read no data; `q_spin_check.py`):
      - The ramp is the cumulative count of Regge trajectories by starting spin. Their starts are the coefficients of one generating function $Q(q,y)$ (y grades spin), built exactly from the light cone. Its sum over spin is S15's Fock space.
      - $Q$ is computed to depth 22. It reproduces all 1.2 million multiplicities and 1.1 million starts in the data ($n\le40$) with zero exceptions.
      - Because $Q$ covers every level, it **proves Regge monotonicity (half (a)'s input) at every weight within 22 of the top, for every $n$**. It also proves GMSZ's shift-relation saturation (last start at $\ell=w$ or $w-1$) there.
      - The upper half of every ramp (spin $\ge(D-1)/2$) is a second closed form, the residue at the $\alpha^+_{-2}$ pole (S8's $T_2$ shift). It is non-negative to $g\le23$; a positive formula for it is open.
      - Half (a) at all orders is now positivity of one explicit function. A multi-Regge expansion (one residue per $\alpha^+_{-k}$) is the next test.
    - **S17, the Regge directions** (`q_multi_regge.py`, reads no data):
      - The stable total along each Regge direction $T_k$ (spin $+1$, level $+k$, the oscillator $\alpha^+_{-k}$) is a residue $R_k$. Each has an explicit positive Fock-space formula, derived by S15's mechanism and checked for $k\le5$. $R_2$ supplies the positive formula S16 left open.
      - At $D\le22$, for every level, the starts are $T_2$-monotone except on the vacuum line (GMSZ's gap; S8's exemption). $T_3$-monotonicity fails (389 negatives), as S8 found for the SO(10) heads.
      - Successive residues give the starts exactly in successive spin bands: $\hat Q$ for $2N-D\ge-1$, and $\hat Q$ minus a sum of $\hat O$ for $3N-D\ge-2$. With a degree bound, Regge monotonicity on the upper half of every ramp follows at all orders. The lower part needs a reason for $T_2$-monotonicity.
    - **S18, theorem: Regge monotonicity for spin $\ge(n-1)/3$** (`q_degree_bound.py`; checks in `q_spin_check.py` (3)):
      - For every level $n$ and every weight $\nu=(N,\mu)$ with $3N\ge n-1$, $r^9_n(\nu)-r^9_{n-1}(\nu-e_1)$ equals the multiplicity of $\mu$ at grade $n-2N$ in an explicit positive SO(7) Fock space $\hat Q$. So it is $\ge0$.
      - The proof is an excess count (2·spin − depth), factor by factor, in the light-cone formula for $(1-qy)Q$. It needs no cancellations, only standard facts (light cone, $256\otimes r^9_n$, parabolic BGG).
      - Checked weight by weight to depth 18, and on all 16 683 starts with $3N\ge n-1$ in the data ($n\le40$).
      - Half (a)'s input is now proved at every level for spin $\ge(n-1)/3$ and within 22 of the top. The open part is low spin, $\nu_1<(n-1)/3$, beyond $n=40$.
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
