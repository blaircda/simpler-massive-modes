# uplift_structure — exploring the structure of the SO(10) uplifts

This directory holds exploratory analysis of the computational output in [`../so10_uplift/`](../so10_uplift/STATUS.md): invariants, fibre structure, and the kernel of the uplift problem. These directions may or may not lead to the all-orders proof. The hard computational evidence (existence through level 24, counts, canonical choice) stays in `so10_uplift/`.

All scripts here only *read* `so10_uplift/data/*.json` and reuse `so10_uplift/uplift.py`, `weights.py` (via `common.py`). They never import baseline files and never write to `so10_uplift/`.

| script | topic | output | runtime |
|---|---|---|---|
| `count_invariant.py` | S1 | `results/count_invariant.txt` | ~1 min |
| `fibre_analysis.py` | S2 | `results/fibre_analysis.txt` | ~1 min |
| `canonical_in_fibres.py` | S3a | `results/canonical_in_fibres.txt` | ~10 min |
| `kernel_moves.py` | S3b | `results/kernel_moves.txt` | ~15 min |
| `delta_ideal.py` | S4 | `results/delta_ideal.txt` | ~5 s |
| `regge_induction.py` | S5 | `results/regge_induction.txt` | ~20 s |
| `regge_stress.py` | S5 step 1: induction to $n=32$; (a) for all uplifts (LP); (b) on samples | `results/regge_stress_*.txt` | ~2 min / ~2 min / ~10 s |
| `regge_extremal.py` | S5 step 2: per-weight LP maxima of top content vs min/max-height uplifts | `results/regge_extremal.txt` | ~1 min ($n\le18$) |
| `regge_a_star.py` | S5 step 2: universal certificate (lemma), (a*) to $n=32$, certificate optimality | `results/regge_a_star.txt` | ~1 min |
| `regge_monotone.py` | S6: arc-space algebra $B$; $M'=Z_m/\mathrm{Hilb}(B)$; induction lemma | `results/regge_monotone.txt` | ~3 min ($n\le32$) |
| `regge_heads.py` | S8: heads $E_n$ under two rules; product test | `results/regge_heads.txt` | ~5 s |
| `head_pl.py` | S8: SO(10) super-PL of the head series, chiralities fixed by the algebra | `results/head_pl_max-height.txt` | ~1 min (degree 8) |
| `regge_modes.py` | S8: multi-mode Regge induction (shifts $T_k$, $k\in K$) | `results/regge_modes.txt` | ~15 s |
| `regge_directions.py` | S8 (v): head uniqueness; Cartan-direction scan; forced inductions; null test; joint $T_2$, $T_w$ | `results/regge_directions.txt` | ~3 min |
| `gmsz_compare.py` | S9: GMSZ (2609.04107) support rule, shift relations, monotonicities; heads vs stable range | `results/gmsz_compare.txt` | ~2.5 min |
| `bm_compare.py` | S10: Basile–Markou principal embedding and trajectory moves on the full spectrum | `results/bm_compare.txt` | ~10 s |
| `half_b.py` | S11: half (b): greedy rules, restricted box families, the level recursion for $l_5$-extremal lifts | `results/half_b.txt` | ~15 s |
| `cumulative_lift.py` | S12: all-full lift of $R^9_n$, forced recursion, bottom chains, uplift by the level recursion; `--bosonic` for $B$ | `results/cumulative_lift*.txt` | ~10 s (n ≤ 40) |
| `beyond32.py` | S11/S12 beyond 32 (data `r9_to40.json`): Regge induction to 40; extremal lifts of $\rho_n$; intermediate-$l_5$ classes | `results/beyond32.txt` | ~2.5 min |
| `repair_levels.py` | S13: repair levels $K_L(n)$, $K_J(n)$ of $\rho_n$ (both rules) and minimal top-level repairs | `results/repair_levels*.txt` | ~25 min (n ≤ 40, both rules) |
| `repair_asymptotics.py` | S13: column indicator $e_c(n,k)$ for $R^9_n$; asymptotic $D_4(\dim)$ | `results/repair_asymptotics*.txt` | ~1 min |
| `repair_capacity.py` | S13: least vs largest number of repairs per level | `results/repair_capacity.txt` | ~6 min |
| `transition_region.py` | S14: spectrum profile; negative region and excess band of the all-full lift; edges over $n$ | `results/transition_region*.txt` | ~1 min |
| `super_pl.py` | S7: super-plethystic log of the SO(9) spectrum | `results/super_pl.txt` | ~5 s (degree 7) |
| `regge_step2.py` | S5 step 2: optimal duals at low weights; singlet margins; Regge monotonicity; $W$ factorisation | `results/regge_step2_*.txt` | ~30 s total |

Notation: an uplift is $c=(c_\lambda)$ over D5 conjugation classes $\lambda=(l_1,\dots,l_5)$, $l_5\ge0$, with $\sum_\lambda c_\lambda\,V_\lambda|_{B4}=R^9_n=\sum_\mu m_n(\mu)V_\mu$. Branching is interlacing: $l_1\ge\mu_1\ge l_2\ge\mu_2\ge l_3\ge\mu_3\ge l_4\ge\mu_4\ge|l_5|$.

## Findings (2026-09-30)

### S1. Irrep-count invariant and fibre sizes (`count_invariant.py`)

$N_n=\sum_\lambda c_\lambda$ is the total number of SO(10) irreps in $R^{10}_n$. Chirality splits don't change it. It is the same for every uplift:

$$N_n=\sum_{a\ge b\ge 0} m_n(a,a,b,b),$$

i.e. the total multiplicity in $R^9_n$ of B4 irreps with Dynkin labels $(0,p,0,q)$.

*Proof.* If $\mu$ interlaces $\lambda$ and $\mu_1=\mu_2$, then $\mu_1\ge l_2\ge\mu_2$ forces $\mu_1=\mu_2=l_2$. Likewise $\mu_3=\mu_4$ forces $\mu_3=\mu_4=l_4$. So $V_\lambda|_{B4}$ contains exactly one such constituent, $(l_2,l_2,l_4,l_4)$.

**Fibre refinement.** For each $(a,b)$, the number of classes in *any* uplift with $l_2=a$, $l_4=b$ (the fibre $F(a,b)$) equals $m_n(a,a,b,b)$.

Checks:
- The functional $f$ with $\sum_{\mu\in I(\lambda)}f(\mu)=1$ for all D5 classes with $l_1\le5$ is unique and equals $1[\mu_1=\mu_2,\mu_3=\mu_4]$.
- The all-ones vector lies in the restriction row space at every $n\le24$.
- The fibre sizes hold for all stored uplifts ($n\le24$) and for all enumerated uplifts ($n\le10$).

$N_n$ = 1, 1, 2, 3, 6, 11, 21, 39, 75, 141, 266, 495, 921, 1688, 3086, 5587, 10048, 17925, 31779, 55921, 97804, 169904, 293400, 503579.

### S2. Fibres are coupled; the freedom lies only between fibres (`fibre_analysis.py`)

- **Coupling.** Fibre $(a,b)$ restricts into the box $\mu_2\le a\le\mu_1$, $\mu_4\le b\le\mu_3$. B4 irreps reachable from more than one fibre carry 25% of $R^9_n$ at $n=4$, 73% at $n=10$ and 90% at $n=24$. The fibre graph has exactly 2 connected components: tensors and spinors.
- **Not a product.** For $n\le10$, the set of all uplifts is not a product of per-fibre choices. At $n=10$ there are 6057 uplifts, against about $9.6\times10^{10}$ combinations of the per-fibre options that occur.
- **Restriction is injective on each fibre** (the fibre-local kernel is 0 for all $n\le24$). *Proof:* within $F(a,b)$, $\mathrm{Res}\,V_\lambda=1[a\le\mu_1\le l_1]\cdot1[l_3\le\mu_2\le a]\,1[b\le\mu_3\le l_3]\cdot1[l_5\le\mu_4\le b]$. This is a tensor product of three families indexed by $l_1$, $l_3$, $l_5$, each triangular.
- **Consequence.** Once the SO(9) content has been divided among the fibres, each fibre's SO(10) content is unique, recovered by finite differences. The uplift problem is therefore *how to divide the shared $\mu$ among fibres* so that every fibre's inverse is non-negative.

### S3. The canonical uplift, and the kernel (`canonical_in_fibres.py`, `kernel_moves.py`)

**(a) What the canonical uplift is.**
- $\sum c\,l_2$ and $\sum c\,l_4$ are invariant (by S1). So the height objective reduces to $10\sum c\,l_1+2\sum c\,l_3$ plus a constant.
- The canonical uplift maximises $\sum c\,l_1$ at every $n\le24$. That maximiser is unique for $n\le10$.
- Lex on $(\sum c\,l_1,\sum c\,l_3,\sum c\,l_5)$ reproduces the canonical uplift exactly for $n\le11$. It is not unique from $n=12$, where the tie is a move with all linear label moments zero.
- **Where shared $\mu$ go.** For each non-degenerate box direction, measure where the multiplicity lands.
  - In the $l_4$ direction ($\mu_3>\mu_4$), canonical puts **74–82%** at $b=\mu_4$ ($n=8$–$24$). Min-height splits roughly 50/50.
  - In the $l_2$ direction ($\mu_1>\mu_2$), both favour $a=\mu_2$ (canonical 55–69%, min-height 48–64%). That part looks like branching geometry.
  - So the canonical uplift prefers the *lowest* available $l_4$ fibre.
  - *Correction (2026-09-30):* an earlier corner histogram (35–58% at $(\mu_2,\mu_3)$) was an artefact. On degenerate boxes several corners coincide, and the label depended on the order of the checks. It has been replaced by these marginals.

**(b) Kernel = Euler-class moves.**
- Let $e=\sum_k(-1)^k\Lambda^k(10)$, which restricts to zero on Spin(9). Conjugation-merged, $e=2(1-10+45-120+210-[126])$.
- For a D5 class $\nu$, $e\otimes V_\nu$ is even iff $l_5(\nu)=0$. Define $h_\nu=\tfrac12e\otimes V_\nu$ ($l_5=0$) or $e\otimes V_\nu$ ($l_5>0$).
- For every $n\le16$, the $h_\nu$ supported in the admissible candidates are linearly independent, their number equals the kernel dimension, and they generate the full **integer** kernel. So every uplift is uniquely $R^{10}_{\rm can}+\sum_\nu t_\nu h_\nu$ with $t\in\mathbb Z$.
- For $n\le9$ all uplifts have $t_\nu\ge0$, so canonical is the corner of the orthant. At $n=10$ this fails: 792 of 6057 uplifts have $t_{45}=-1$, where $\nu=(1,1,0,0,0)$.

### S4. The restriction kernel is the ideal $(\delta)$, at all levels (`delta_ideal.py`; added 2026-10-05)

This is a proof, not an empirical observation. The script checks the concrete identities it uses.

Let $\delta=S^+-S^-=[0,0,0,1,0]-[0,0,0,0,1]$ and let $\sigma$ be the outer automorphism of D5 (it swaps $S^\pm$).

1. **$\ker(\mathrm{Res}:R(\mathrm{Spin}_{10})\to R(\mathrm{Spin}_9))=(\delta)$.**
   - The two rings are polynomial: $R(\mathrm{Spin}_{10})=\mathbb Z[\Lambda^1,\Lambda^2,\Lambda^3,S^+,S^-]$ and $R(\mathrm{Spin}_9)=\mathbb Z[\Lambda^1,\Lambda^2,\Lambda^3,\Delta]$.
   - Restriction sends $\Lambda^k_{10}\mapsto\Lambda^k_9+\Lambda^{k-1}_9$ and $S^\pm\mapsto\Delta$.
   - Change generators to $\mu_k=\Lambda^k_{10}-\mu_{k-1}$ (with $\mu_0=1$), a triangular change invertible over $\mathbb Z$. Then restriction is $\mu_k\mapsto\Lambda^k_9$, $S^+\mapsto\Delta$, $\delta\mapsto0$.
2. **Euler class.** $e=\sum_k(-1)^k\Lambda^k(10)=-\delta^2$.
3. **The $\sigma$-invariant part of $(\delta)$ is $e\cdot R^\sigma$.**
   - An anti-invariant $y(S^+,S^-)$ vanishes at $S^+=S^-$, so $y=\delta z$ with $z\in R^\sigma$. Here $R^\sigma=A[S^++S^-,S^+S^-]$, with $A$ the polynomial ring in the $\mu_k$.
   - Conversely $\delta z=(1-\sigma)(S^+z)$.
4. **Conjugation-class lattice (the setting of `so10_uplift/`).**
   - Let $\pi$ merge $V_\lambda,V_{\bar\lambda}$. Restriction factors through $\pi$, so the class-level kernel is $\pi((\delta))$.
   - Using $\pi=\pi\circ\sigma$ and item 3, $\pi(\delta y)=-\tfrac12\pi(e z)$ with $z=(y-\sigma y)/\delta$, and every $z\in R^\sigma$ arises (take $y=S^+z$).
   - $R^\sigma$ has $\mathbb Z$-basis $V_\nu$ ($l_5=0$) and $V_\nu+V_{\bar\nu}$ ($l_5>0$). Hence the class-level kernel is generated over $\mathbb Z$ by $\tfrac12\pi(eV_\nu)$ ($l_5=0$) and $\pi(eV_\nu)$ ($l_5>0$).
   - That is exactly the $h_\nu$ of S3b, now valid at every level on the full (unbounded) lattice.
   - The S3b statement that the $h_\nu$ *supported inside the admissible candidates* generate the finite-window kernel is a further support statement. It is still empirical ($n\le16$).

**Consequences.**
- **Existence.** Existence of a non-negative uplift is a purely conjugation-class question. Any non-negative class-level solution can be split into chiralities arbitrarily.
- **Chirality is a selection question.** It is controlled by the $\delta$-odd part of the kernel, which is invisible on classes.
- **The BBMS difference is an example.** The level-6 difference between BBMS and the baseline reported in `../previous_efforts.md` is $[0,0,0,1,0]+[1,0,0,1,0]-[0,0,0,0,1]-[1,0,0,0,1]=\delta\,(10+2)$. It vanishes under $\pi$.

### S5. Regge induction (`regge_induction.py`; added 2026-10-11)

**Set-up (exact).**
- Let $V^*_\lambda=V_\lambda-V_{\lambda-e_1}$ (the "starred" reps of BBMS, hep-th/0310292 eq. 2.13), with $V^*_\lambda=V_\lambda$ when $l_1=l_2$.
- By interlacing, $\mathrm{Res}\,V^*_\lambda$ is exactly the **top slice** $\{\mu:\mu_1=l_1,\ (\mu_2,\mu_3,\mu_4)\prec(l_2,..,l_5)\}$, which is the SO(8)→SO(7) branching of $(l_2,..,l_5)$ placed at $\mu_1=l_1$. In the starred basis, restriction is diagonal in $\mu_1=l_1$.
- Let $T:V_\lambda\mapsto V_{\lambda+e_1}$ (Regge shift), $\mathrm{top}(A)=\sum_\lambda A_\lambda\,\mathrm{Res}\,V^*_\lambda$, and $S$ the shift $\mu_1\mapsto\mu_1+1$. Then
$$\mathrm{Res}\,T(A)=\mathrm{Res}\,A+S(\mathrm{top}\,A).$$
- **Regge induction:** $A_n=T(A_{n-1})+E_n$ with $E_n\ge0$ a non-negative uplift of
$$\rho_n:=R^9_n-\mathrm{Res}\,T(A_{n-1})=r^9_n-S(\mathrm{top}\,A_{n-1}).$$
- Unrolled, this is a decomposition of every level into **SO(10) Regge trajectories**: $r^9_n=\sum_\lambda h_n(\lambda)\,\mathrm{Res}\,V_\lambda+\sum_{j\ge1}\sum_\lambda h_{n-j}(\lambda)\,\mathrm{Res}\,V^*_{\lambda+je_1}$ with heads $h_m=E_m\ge0$. These are Bars' "copies of earlier levels", organised as trajectories.
- $|E_n|=N_n-N_{n-1}$ for any rule, by S1, since $T$ preserves the irrep count.

**Results.**
1. **Canonical uplifts.** $A^{\rm can}_n\supseteq T(A^{\rm can}_{n-1})$ holds only for $n\le8$; the canonical uplift is not inductive. But $\rho_n\ge0$ for all $n\le24$: the shifted canonical uplift never overshoots $R^9_n$.
2. **The induction runs to $n=24$** from $A_1=1$, for each of three rules for $E_n$: max-height, min-height, and an arbitrary feasible ILP solution. Each run takes seconds.
3. **Path independence.** For *every* enumerated uplift $A_{n-1}$ ($n-1\le10$; 6057 uplifts at $n-1=10$), $\rho_n\ge0$ and $\rho_n$ is liftable. The residual depends on $A_{n-1}$ only through its top slices: 90 distinct residuals at $n=11$.

**Step 1 stress tests (2026-10-11; `regge_stress.py`).** Data: $r^9_n$ for $n\le32$ in `data/r9_to32.json`, from `so10_uplift/derive_r9.py`. It is identical to `so10_uplift/data/r9.json` for $n\le24$.
- **Induction to $n=32$** for all three rules for $E_n$, with $|E_n|=N_n-N_{n-1}$ (12,294,174 at $n=32$).
- **(a) for *all* uplifts, rigorously, for $n\le24$.**
  - The identity $\mathrm{top}(A)=\mathrm{Res}\,A-\mathrm{Res}\,T^{-1}(A)$ turns (a) at $\mu=\nu+e_1$ into $\min_A\mathrm{Res}\,T^{-1}(A)(\nu)\ge R^9_{n-1}(\nu)-r^9_n(\mu)$. An LP lower bound on the minimum suffices.
  - It is certified at every weight with positive deficit (484 weights at $n=24$).
  - It is **tight but never violated**: slack exactly 0, always at symmetric tensors $\nu=(k,0,0,0)$, on the leading trajectory. The deficits there come in repeating pairs: 1,1,1,1,2,2,3,3,7,7,12,12,25,25,43,43,83,83,145,145,267,267 for $n=3..24$.
- **(b) on sampled uplifts, $n=12..18$.** Samples are random-objective ILP vertices plus random walks along the S3b Euler moves: 36–111 distinct uplifts per level. Every residual is non-negative and liftable. The samples are modest because the walks accept few moves.

**Step 2(a): the continuation inequality (2026-10-11).**

*Lemma (proved; pointwise check on 1584 classes in `regge_a_star.py`).*
- **Statement.** For a B4 weight ν let $\tilde\nu=(\nu_1+1,\nu_1+1,\nu_3,\nu_4)$ and $y_\nu=\delta_\nu-\delta_{\nu+e_1}+\delta_{\tilde\nu}$. Then $y_\nu\cdot\mathrm{Res}\,V_\lambda\ge[\text{top slice of }\lambda\ni\nu]$ for every class λ.
- **Proof, by cases on $l_2$.**
  - If $l_2\le\nu_1$, then $\tilde\nu\notin$ box and the left side is exactly $[l_1=\nu_1]\cdot[\nu'\prec l']$.
  - If $l_2=\nu_1+1$, then $\nu\notin$ box and $[\tilde\nu\in\text{box}]\ge[\nu+e_1\in\text{box}]$, because the $\tilde\nu$ conditions are weaker ($l_3\le\nu_1+1$ instead of $l_3\le\nu_2$).
  - If $l_2\ge\nu_1+2$, all three terms vanish.
- **Consequence.** For every non-negative uplift $A$ of $R=R^9_{n-1}$: $\mathrm{top}(A)(\nu)\le R(\nu)-R(\nu+e_1)+R(\tilde\nu)$. So half (a) holds for **all** uplifts at every $\mu=\nu+e_1$ where the uplift-free inequality
$$(\mathrm a^*)\qquad R^9_n(\nu+e_1)\ \ge\ R^9_{n-1}(\nu)+R^9_{n-1}(\tilde\nu)$$
holds.

*Data.*
- **(a*) holds** at every weight for $n\le12$.
- **Failures are few and near the origin.** For $13\le n\le32$ it fails only at a few weights near the origin: the singlet, the vector, the spinor, and later $(2,0,0,0)$; at most 5 of 1,426 weights at $n=32$.
- **Equality is common.** (a*) holds with *equality* at about half of all weights (748/1426 at $n=32$), presumably the Regge-edge region of GMSZ's shift relations.
- **The certificate is optimal** (equals the LP maximum of top content) everywhere except those few low weights ($n\le20$).
- **Better certificates at the low weights** (`regge_step2.py duals`):
  - The spinor's optimal certificate is the $n$-independent 5-term chain $+[0,0,0,1]-[1,0,0,1]+[0,1,0,1]-[0,0,1,1]+[0,0,0,3]$ (Dynkin). It extends $y_\nu$ by also raising coordinates 3 and 4.
  - The singlet's is $+[0]$ plus a 4-term chain placed near the top of the support ($\mu_1=n-6$).
- The per-weight maxima are integral ($n\le18$). No single uplift attains them all: the min-height uplift does only for $n\le10$.

*Universal versus existential.*
- **Some uplift hoards nearly every possible singlet.** $\max_A A(\text{singlet})\approx R^9_{n-1}(1)$ (636 vs 639 at $n=20$). So the universal (a) at the vector needs $r^9_n(9)\gtrsim R^9_{n-1}(1)$.
  - The ratio rises from 1.25 to 1.91 over $n=8..32$.
  - But fixed-irrep asymptotics give $r_n(9)/r_n(1)\to9$, while $R_{n-1}(1)/r_n(1)$ grows like $\sqrt n$. That *suggests* (heuristic, unproved) that the every-uplift form fails eventually, at $n$ of order $10^3$ or beyond.
- **A proof should therefore target the existential induction with a specific rule.** The max-height induction never creates an SO(10) singlet ($A_{n-1}(1)=0$ for all $n\le32$), so it has the full margin $r_n(9)$ there.

*Regge monotonicity* (`regge_step2.py monotone`). $r^9_n(\nu+e_1)\ge r^9_{n-1}(\nu)$ for **all** ν and $n\le32$, and hence also cumulatively. This would follow from an injective "multiply by the mode-1 vector's highest-weight vector" map if the spectrum were a free $\mathrm{Sym}(q\cdot9)$-module.

*Plethystic factorisation* (`regge_step2.py wfactor`).
- **Definition.** $R^9(q)=\mathrm{Res}\,\mathrm{PE}[q\cdot10]\cdot W(q)$ with $W=R^9\cdot\mathrm{Res}\,K(q)=Z_m\cdot\mathrm{PE}[-q\cdot9]$, where $K(q)=\sum_k(-1)^k\Lambda^k(10)q^k$ is the graded Euler class of `../bbms_lift/` B2.
- **Positivity.** $W_j\ge0$ for all $j\ge8$ (and $j=5,6$). It is negative only at $j=2,3,4,7$, with one entry each: $-1$, $-16$, $-9$, $-1$.
- **Liftability.** $W_j$ is non-negatively SO(10)-liftable at $j=10$ and for all $12\le j\le31$. It is not liftable at $j=5,6,8,9,11$.
- So the free-module structure holds only eventually. An explicit $\mathrm{PE}[q\cdot10]\otimes W^{(10)}$ construction would need a low-order correction.

*Emerging proof architecture for (a).*
- **Edge:** the universal certificate, plus (a*) as an inequality on the string spectrum itself. It is often an equality, so it calls for an injection or shift-relation argument.
- **Bulk** (low weights): a bound specific to the construction, e.g. the max-height rule keeps low-$l_1$ SO(10) content small.
- Half (b) has not yet been attacked theoretically. Note that (a)+(b) together $\iff$ $T(A)\le$ some uplift of $R^9_n$.

**Conjecture (Regge continuation).** For all $n\ge2$ and every non-negative uplift $A$ of $R^9_{n-1}$:
- (a) $r^9_n\ge S(\mathrm{top}\,A)$, i.e. continuing all trajectories one step never overshoots level $n$;
- (b) $\rho_n=r^9_n-S(\mathrm{top}\,A)$ has a non-negative uplift (the new heads).

By induction from $A_1=1$, this implies existence at all $n$. Path independence means a proof needs only a one-step lemma for an arbitrary uplift.

### S6. Regge monotonicity: reduction, and a proof for the bosonic analogue (2026-10-12)

**(i) Half (a) of the Regge induction follows from per-level Regge monotonicity (proved).**
- *Lemma.* Along any Regge induction $A_n=T(A_{n-1})+E_n$ (any rule for $E_n\ge0$), $\mathrm{top}(A_{n-1})\le r^9_{n-1}$.
- *Proof.*
  - $\mathrm{top}(TA)=S(\mathrm{top}A)$, so $\mathrm{top}(A_{n-1})=S(\mathrm{top}A_{n-2})+\mathrm{top}(E_{n-1})$.
  - The induction step gives $\mathrm{Res}\,E_{n-1}=r^9_{n-1}-S(\mathrm{top}A_{n-2})$, hence
    $$\mathrm{top}(A_{n-1})=r^9_{n-1}-\mathrm{Res}\,T^{-1}(E_{n-1})\le r^9_{n-1}.$$
  - The base case is $A_1=1$.
  - Checked along the max-height induction for $n\le32$.
- *Consequence.* If **per-level Regge monotonicity** $r^9_n(\nu+e_1)\ge r^9_{n-1}(\nu)$ holds, then $\rho_n=r^9_n-S(\mathrm{top}A_{n-1})\ge r^9_n-S(r^9_{n-1})\ge0$. So half (a) holds at every level, whatever rule picks the new irreps. No certificates are needed, and the possible large-$n$ failure of the every-uplift form (S5) is irrelevant.
- Monotonicity itself holds for all ν and $n\le32$ (S5, `regge_step2.py monotone`).

**(ii) Monotonicity for the bosonic analogue (proved, modulo standard results).**
- *The algebra.* Let $Q=\{y\in\mathbb C^9:y\cdot y=0\}$ and $B=\mathbb C[J_\infty(Q)]$, the arc space $\{Y(t)=\sum_k y_kt^k: Y\cdot Y=0\}$, graded by $\deg y_k=k+1$.
- *Regularity.* $Q$ is a canonical lci singularity (a Gorenstein rational quadric cone). By Mustaţă (Invent. Math. 2001), the jet schemes $J_m(Q)$ are irreducible lci of the expected dimension $8(m+1)$. So the coefficients of $Y\cdot Y$ form a regular sequence, and $\mathrm{Hilb}(B)=\mathrm{PE}[\sum_mq^m9]\prod_{m\ge2}(1-q^m)$.
  - This is exactly the cumulative bosonic spectrum $\mathrm{PE}[\sum_mq^m8_v]/(1-q)$.
  - Side remark: it is a genuine SO$(D-1)$ character in *any* dimension, unlike the per-level spectrum.
- *$B$ is a domain.* It is a union of the domains $\mathbb C[J_m(Q)]$, each reduced (CM and generically reduced) and irreducible.
- *Monotonicity of $B$.* Multiplication by the highest-weight vector $y_0^+$ (weight $e_1$, degree 1) is injective. It commutes with the SO(9) raising operators, which act as derivations killing $y_0^+$. So it maps highest-weight vectors of weight ν injectively to ones of weight $\nu+e_1$, giving $m_{B,d+1}(\nu+e_1)\ge m_{B,d}(\nu)$.
- *Checked to degree 31:* $B$ is non-negative and monotone.

**(iii) Crystal lemma (proved).**
- If $A$ is Regge-monotone and $M\ge0$, then $A\otimes M$ is Regge-monotone.
- Reason: $c^{\lambda+e_1}_{\mu+e_1,\kappa}\ge c^{\lambda}_{\mu,\kappa}$. By Kashiwara's tensor rule, $c^\lambda_{\mu\kappa}=\#\{b\in B(\kappa):\mathrm{wt}\,b=\lambda-\mu,\ \varepsilon_i(b)\le\langle\mu,\alpha_i^\vee\rangle\}$, and adding the dominant $\omega_1=e_1$ to μ only relaxes the conditions.

**(iv) Obstruction: the superstring does not factor through $B$.**
- *The test.* Per-level monotonicity would follow from $M'=Z_m/\mathrm{Hilb}(B)\ge0$. In light-cone terms, $M'=(F_F-F_B^{-1})(1-q)/(q(8_v+8_s))$.
- *The result.* $M'$ is negative at every degree $2\le d\le31$, on a small fixed set: $9$, $[0,0,1,0]$, $[1,1,0,0]$, $[1,0,0,2]$, with coefficients growing to about $-5\times10^4$ at $d=31$. The cumulative quotient $M=M'/(1-q)$ is also negative.
- *The cause.* $B$ has a mode-2 vector ($y_1$, a 9 at degree 2), whereas the superstring has a spinor 16 there ($r^9_3=44+16$). The fermions are entangled with the bosonic modes.
- *What is needed.* An SO(9)-covariant model of the superstring spectrum itself: plausibly the arc space of a super null cone (bosonic part $Q$, an 8-dimensional fermionic fibre such as $\ker\slashed y\subset16$, or a κ-symmetry quotient) that is torsion-free over $B$. Then the same highest-weight argument applies. Possibly relevant: the pure-spinor character literature (Berkovits–Nekrasov, hep-th/0503075), and the pure-spinor numerator met in `../bbms_lift/` B4.

### S7. Toward a fermionic model: literature, and the super-plethystic log (2026-10-12)

**Literature** (see `../LITERATURE_SEARCH.md`, follow-up section).
- The arc space of a quadric cone as the γ-field generating function is standard (Chandrasekhar–Ray; Aisaka–Arroyo). That supports S6(ii).
- For pure-spinor and quadric βγ systems the full Hilbert space has extra *fermionic*, non-globally-defined states (Aisaka–Arroyo–Berkovits–Nekrasov). So "functions on a supervariety" is likely too naive.
- No SO(9)-covariant non-negative model of the massive spectrum was found.

**Super-PL** (`super_pl.py`).
- *Method.* Spinorial weights are treated as fermionic and tensorial ones as bosonic. The spectrum is recomputed from the Fock space and matched to stored data.
- *Result.*
  - $\mathrm{PL}(Z_m)$: $9$; $16-1$; $9+36-16$; $-9-36-84+2\cdot16+128$; then growing alternating sums at degrees 5–7.
  - So there is no finite generator/relation structure, and no complete-intersection model.
- *Consequence.* Every super arc space whose bosonic generators are consecutive modes $y_0,y_1,\dots$ of a vector fails already at degree 2. Its $y_1$ gives an extra $9$, but $\mathrm{PL}_2=16-1$. The superstring has no independent mode-2 vector: the light-cone mode-2 oscillators are absorbed by the SUSY factor $8_v+8_s$.

**An SO(10) pattern in the cumulative spectrum.**
- $\mathrm{PL}(R)=\mathrm{PL}(Z_m)+q=\mathrm{Res}\big(10\,q+(16-1)q^2+(45-\overline{16})q^3+\dots\big)$. That is SO(10) generators $10$ (degree 1), $16$ (2), $45$ (3), with relations $1$ (2) and $\overline{16}$ (3).
- These are the shapes of $y\cdot y=0$ and $\slashed y\theta=0$ for a null 10-vector and a spinor, plus an adjoint, like angular momentum $x\wedge p$ for a massless particle.
- The canonical uplift reproduces $A_4=210'+144+45$ exactly from these terms.
- This suggests looking for an **SO(10)-equivariant** algebra (massless 10D superparticle / phase-space type) whose SO(9) restriction is $R$. That would settle existence outright, with positivity manifest, rather than proving SO(9) monotonicity.
- Caveat: `previous_efforts.md` item 1 found no simple structure in the PL of the baseline uplift to level 10. The pattern beyond degree 4 is untested.

### S8. The model test: structure of the trajectory heads (2026-10-12)

The test asks, without guessing relations, whether the new irreps (heads) $E_n$ of the Regge induction have "superparticle" structure (S7).

**(i) Heads multiply** (`regge_heads.py`). The heads are rule-independent through $n=8$; the max- and min-height rules first differ at $n=9$.
- $E_3=16$ and $E_4=45$.
- $E_5=1+120+144$, and it contains $\Lambda^216=120$.
- $E_6=16+126+144+320+560$, and it contains all of $16\otimes45=16+144+560$.
- $E_7$ contains $\mathrm{Sym}^245=1+54+210+770$ *except the singlet*, which acts like a relation "$w\cdot w=0$".

So products of earlier heads reappear as heads, consistent with an algebra $A$ with $H=A/(y)$.

**(ii) But not a small algebra** (`head_pl.py`). This is the SO(10) super-PL of the head series. Chiralities are fixed degree by degree to match predicted products; the naive "all spinors are 16" convention creates spurious $\pm$ pairs.
- Generators: $16$ (level 3); $45$ (level 4); $1+144$ (level 5); $126+320$ (level 6).
- From level 7: mixtures of new generators and relations ($-1$ from $w\cdot w$, $-45$, $-54$, $-210$, $-1050$, …), with no sign of termination through level 9.

There is no finite superparticle-type presentation.

**(iii) A second Regge direction** (`regge_heads.py` part 3, `regge_modes.py`). The second-wave generators $144=[1,0,0,0,1]$ (level 5) and $320=[1,1,0,0,0]$ (level 6) are θ and $w$ shifted by $+e_1$ at level $+2$.
- Define $T_2$ as the same Cartan shift as $T$ ($\lambda\mapsto\lambda+e_1$), but at level $+2$.
- $E_n\supseteq T_2(E_{n-2})$ holds for all $4\le n\le8$, where the heads are unique (see (v)).
- The vacuum is exempt: $T_2(E_1)=10$ is not in $E_3=16$.
- The plain rules do not enforce this containment, and lose it from $n=11$ (max-height) and $n=9$ (min-height).
- **Double Regge induction runs to $n=32$ for all three tested rules (max-height, min-height, feasible).** Here $E_n=T_2(E_{n-2})+P_n$ with $P_n\ge0$, vacuum exempt. So the uplift can be chosen with two commuting Regge directions, $T$ (level $+1$) and $T_2$ (level $+2$).
- The primary counts are 0, 1, 1, 2, 4, 7, 13, 26, 48, 89, 163, 301, 538, 972, 1734 for $n=2..16$. They are *automatically* rule-independent: $|P_n|=|E_n|-|E_{n-2}|$, and $|E_n|=|A_n|-|A_{n-1}|$ is fixed by the uplift-invariant irrep count (S1). They are about 70% of $|E_n|$.

**(iv) No level-3 direction.** Forcing $T_3$ (shift $+e_1$, level $+3$) makes the induction **fail at $n=7$ for every rule**. This holds for $T_3$ on its own, together with $T_2$, and for all $T_k$ with $k\ge2$.
- On its own the forced copy is $T_3(E_4)=320$, and no non-negative lift of $\rho_7$ contains it.
- $A_6$ is the unique uplift, so this is a strict obstruction.
- So the structure is not a Cartan module over all string modes (an arc space).

**(v) Which Cartan shifts are directions? Scan, calibration, null test** (`regge_directions.py`).
- **The heads are unique through $n=8$.** At each $n\le8$, $\rho_n$ has exactly one non-negative lift along the induction, so the low-level containments above are rule-free facts.
- *Scan.* All 38 dominant shifts $s$ with $l_1\le2$ and level offsets $k=1..6$ were tested for $E_n\supseteq s(E_{n-k})$ on the unique heads ($k+2\le n\le8$).
  - 13 pairs pass, but only three are tested on at least 5 irreps:
    - $T_2=(e_1,2)$: 10 irreps; vacuum exempt.
    - $T_w=(e_1+e_2,3)$: 5 irreps; vacuum included, since $E_4=45=T_w(1)$.
    - $(\omega_5,3)$: 5 irreps.
  - The other 10 are tested on only 1–2 irreps.
- *Forced single-direction inductions to $n=32$ (three rules).*
  - $T_2$ and $T_w$ survive.
  - $(\omega_5,3)$ fails at $n=9$, or at $n=13$ if it acts at most once, like a fermion.
  - Of the 10 weakly tested pairs, six survive and three fail. One more, $(\omega_4,4)$, survives only as a fermion (acting at most once). The survivors include the composites $T_2^2=(2e_1,4)$ and $T_2T_w=(2e_1+e_2,5)$.
  - So surviving to $n=32$ alone is weak evidence when the low-level test is weak.
- *Null test (delayed start).* Each direction is forced only on primaries at levels $m\ge m_0$, for $m_0=2,5,7,9$.
  - $T_2$ and $T_w$ survive to $n=32$ from every start.
  - These wrong directions fail no later than two levels after they can first act, for every start: $T_3$, $(e_1,1)$, $(e_1+e_2,2)$, $(\omega_5,2)$, $(e_1+e_2+e_3,3)$ and $(2e_1,3)$.
  - The exception is the singlet shift $(0,2)$: it fails for $m_0\le5$ but survives from $m_0\ge7$.
  - So for non-zero shifts the test discriminates: **$T_2$ and $T_w$ are genuine secondary Regge directions.**
- *But they do not act freely together.*
  - Forcing the free commutative action of $T_2$ and $T_w$ fails at $n=9$ for every rule. This holds with the vacuum exempt from $T_2$, and with $T_w$ either acting on the vacuum or not.
  - The forced part needs three copies of $3696_s=[1,1,0,0,1]$: $T_2$ of the two primaries $560$ at level 7, and $T_w$ of the primary $144$ at level 6. Only two fit.
  - The heads up to $n=8$ are unique, so this is a strict obstruction.
  - The weaker joint containment $E_n\supseteq\max(T_2E_{n-2},T_wE_{n-3})$ also fails under the three greedy rules, at $n=12$–$16$. That failure is not strict, because the rules make choices from $n=9$ on.
  - The excess copy would disappear if $144$ (level 6) $=T_2(x)$ and $560$ (level 7) $=T_w(x)$ for a common $x=16$ at level 4, because then $T_2T_w(x)$ is counted once. But $E_4=45$ has no 16.

**Reading (tentative).**
- $T=(e_1,1)$, $T_2=(e_1,2)$ and $T_w=(e_1+e_2,3)$ are exactly the highest-weight vectors $y_0^+$, $y_1^+$ and $(y_0\wedge y_1)^+$ of two vector oscillators: $y_0$ in degree 1 and $y_1$ in degree 2.
- Freely combined, their monomials would give the SO(10) highest weights of the harmonics of $\mathbb C[y_0,y_1]$, graded by GL(2) characters at $(q,q^2)$. This is Howe duality for $O(10)\times\mathfrak{sp}(4)$.
- The $45=w$ at level 4 fits $y_0\wedge y_1$, an angular-momentum-like $x\wedge p$.
- The data confirm each generator separately, but **rule out the free combination at $n=9$**. They also rule out $y_1$ acting on the vacuum: there is no 10 at level 3.
- So the naive two-oscillator harmonic model fails. A super-extension in which a 16 at level 4 is cancelled (the ghost $x$ above) is conceivable but untested.

### S9. Comparison with GMSZ, arXiv:2609.04107 (2026-10-12) (`gmsz_compare.py`)

GMSZ's $a_{[\ell,m,r,s],n}$ is our $r^9_n$, with the same level convention (top irrep $[n-1,0,0,0]$).

**The data agree.**
- Our $r^9_7$ equals their worked example (E.19) exactly. Their appendix labels it "level 6"; the offset is confined to the appendix.
- Their support rule (Table 1 / (3.6)), including the gap $[n-2,0,0,0]$, holds at every level $7\le n\le32$ of our data. With the levels shifted by $\pm1$ it fails at every level.
- Their shift relations hold with zero violations in the stated regions: 37 757 instances for $\ell$ (3.32), 7 223 for $m$ (3.33) and 620 for $r$ (3.34).

**Their shift relations are the saturated ends of monotone sequences.**
- The three shifts are depth-preserving moves: a column of height 1, 2 or 3 is added to the Young diagram, at a cost of 1, 3 or 6 levels.
- On our data the inequality $r^9_n(\nu+w)\ge r^9_{n-d}(\nu)$ holds at **all** depths, not only in their stable ranges:
  - $(\omega_1,+1)$: 0 violations of 11 746. This is S6's per-level Regge monotonicity, which implies half (a) of the Regge step.
  - $(\omega_2,+3)$: 0 of 9 033.
  - $(\omega_3,+6)$: 0 of 5 906.
  - A fourth move that GMSZ do not list, a column of height 4 ($2\omega_4=[0,0,0,2]$) at $+10$: 0 of 3 121.
- One level less fails every time: 8 021, 2 477 and 72 violations.
- So the costs are the triangular numbers 1, 3, 6, 10. This is the bosonic depth-0 cost $\sum_i i\,s_i$ of a column, i.e. an antisymmetrised $\alpha_{-1}\wedge\dots\wedge\alpha_{-k}$ (Markou–Skvortsov).
- Single boxes behave differently by row:
  - Row 2 is monotone only at $+3$, not at the bosonic cost $+2$ (1 078 violations).
  - Row 3 is monotone only at $+4$, not at the bosonic cost $+3$ (269 violations).
  - Row 4 is monotone at its bosonic cost $+4$ (0 violations; $+3$ gives 123).
- How sharp their regions are, at $k=1$:
  - (3.33) and (3.34) are sharp: equality holds exactly in the stated region.
  - (3.32) also holds one step outside its stated region (748/748 at the boundary), and mostly fails one step further (15/678 equalities).

**Dictionary with the SO(10) structures of S5–S8.**
- *$\ell$-shift and $T$.*
  - For $n\le20$, under both rules, every SO(9) component of $\mathrm{Res}\,E_n$ lies outside the $\ell$-stable range of (3.32), and touches its edge (margin $-1$) at every level.
  - Inside that range $\mathrm{top}(A_{n-1})=r^9_{n-1}$ exactly, so there the uplift is pure Regge continuation, with no heads.
  - Their "halfway point" is the frontier of the new SO(10) irreps. It moves by one unit of $l_1$ every two levels, i.e. along $T_2$.
- *$m$-shift and $T_w$.* Restricting $E_n\supseteq T_w(E_{n-3})$ (and so $A_n\supseteq T_w(A_{n-3})$) to SO(9) gives exactly the $(\omega_2,+3)$ inequality. So (3.33) is the SO(9) shadow of $T_w$: same shift vector $e_1+e_2$, same level offset 3.
- *$T_2$ and their gap.*
  - $T_2$, with the vacuum exempt, restricts to $r^9_n(\nu+e_1)\ge r^9_{n-2}(\nu)$ everywhere except on the leading trajectory.
  - The data show exactly 30 violations, all of the form "top $[n-3,0,0,0]$ of level $n-2$ maps to $[n-2,0,0,0]$ at level $n$". That is GMSZ's gap.
  - So **the vacuum exemption of $T_2$ is the SO(10) form of their empirical gap.**
  - $T_2$ raises the depth by one, which is why it is absent from their depth-preserving list.
- *$r$-shift and the 4-form.*
  - The SO(10) head directions $(\omega_3,+6)$ and $([0,0,0,1,1],+10)$ survive forced inductions to $n=32$ and the null test.
  - So does a control, a self-dual 5-form shift $([0,0,0,0,2],+10)$, even though the plain heads violate it at $n=11$–$13$.
  - At these offsets the forced part is too small to discriminate, so these SO(10) directions remain undecided.

**Open on both sides.**
- GMSZ list proving the support rule and the "repeat" relations as open problems (§7). We need the $\ell$-inequality for half (a) (S6).
- Their multiplicities go to level 226 (all irreps) and 501 (single irreps), with data on Zenodo (their ref. [17]). That would test all four monotonicities far beyond $n=32$.

**Lead for half (a).** GMSZ attribute the support rule to explicit depth-0 highest-weight constructions, which use a Howe dual of the Lorentz algebra.
- For the bosonic string the dual is $sp(2K)$ (Markou–Skvortsov 2309.15988); for the superstring it is $osp$ (Basile–Markou 2405.18467).
- Deeper states come from "trajectory-shifting" (dressing) operators.
- Markou–Skvortsov note that adding a box to row $i$ raises the level by $i$, with dressing coefficients smooth in the row lengths.
- An *injective* version of this at fixed depth would prove the column monotonicities, and hence half (a).
- Caveat: for super-representations, single boxes in rows 2 and 3 are not monotone at the bosonic costs, so the superstring mechanism must differ there.

### S10. Basile–Markou, arXiv:2405.18467: reading and checks (2026-10-12) (`bm_compare.py`)

**What the paper does.**
- It gives a covariant (RNS) construction of whole trajectories of the open superstring, based on super-Howe duality.
  - The Lorentz algebra pairs with $osp(2M|2N)$ in the NS sector and $osp(2M+1|2N)$ in the R sector, acting on oscillator mode labels.
  - The Fock space is multiplicity-free under each pair: each Lorentz irrep $Y$ comes with one lowest-weight $osp$-module.
- In the transverse subspace the super-Virasoro constraints ($G_{1/2},G_{3/2}$ in NS; $G_0,G_1$ in R) are combinations of $osp$ lowering operators, plus a "mass" operator in R.
- Physical states of spin $Y$ at depth $w$ are "dressings" of the principal embedding by polynomials in $osp$ raising operators.
  - The coefficients depend on the arm and leg lengths of $Y$. At $w=\tfrac12$ there is one linear condition in the hook lengths.
  - Trajectories are labelled by $(K,w)$, where $K$ is the diagonal length of the Young diagram, rather than by the number of rows as for the bosonic string.
- The principal embedding (depth 0), with $Y=(\ell_1..\ell_K|h_1..h_K)$ (arms; legs including the diagonal):
  - NS: level $\sum_i\big(i\,\ell_i+(i-\tfrac12)h_i\big)$, multiplicity 1.
  - R: level $\sum_i i(\ell_i+h_i)$, multiplicity $2^{K-1}$. They checked only $K=2$.
- They do not discuss monotonicity or stabilization. They work with the full spectrum, not the supermultiplet-reduced $r^9$.

**Checks on our data.** The full spectrum is $F_n=(44+84+128)\otimes r^9_n$, computed by Brauer–Klimyk for $n\le32$.
- *Principal embedding.* Of the 2 798 irreps first appearing at $n\le32$, all but two agree with their prediction: the first level and multiplicity for GSO-even cases, and appearance only later for GSO-odd ones. The two exceptions are the massless vector and spinor (principal level 0).
- This includes the R multiplicities $2^{K-1}$ beyond their check:
  - $K=3$: 346 irreps; the $K=3$ square first appears at level 14 with multiplicity 4.
  - $K=4$: the square first appears at level 30 with multiplicity 8.
- *Trajectory moves at fixed $K$, at their level costs.*
  - Monotone: a box on arm 1 (cost 1) in both sectors (0 violations of 13 577 NS and 9 626 R), and in NS two boxes on leg 1 (cost 1) or legs 1+2 (cost 2).
  - Not monotone (hundreds of violations): boxes on arms 2 and 3, on legs 2 and 3, and on R leg 1.
- *Column moves on $F$.* Monotone at the triangular costs (0 violations); $\omega_2$ at $+2$ is not. This is the same pattern as for $r^9$ (S9).

**A super-representation analogue of the principal embedding** (`gmsz_compare.py` data, $n\le32$; empirical).
- Every super-representation $[\ell,m,r,s]$ first appears exactly at GMSZ's level $\ell+3m+c(r,s)$.
- Its multiplicity there depends only on $(r,s)$, not on $\ell$ or $m$ (rows $r$, columns $s$):

| $r\backslash s$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 1 | 1 | 2 | 1 | 3 | 4 | 4 |
| 1 | 1 | 1 | 3 | 3 | 1 | 3 | 4 | |
| 2 | 2 | 1 | 3 | 3 | 1 | 3 | | |
| 3 | 2 | 1 | 3 | 3 | 1 | | | |
| 4 | 2 | 1 | 3 | | | | | |
| 5 | 2 | 1 | | | | | | |

**Assessment for half (a): no ready-made injection.**
- *No route to $r^9$.* The framework is RNS without manifest spacetime supersymmetry: the supercharges involve spin fields, not $osp$ bilinears. So dividing out $44+84+128$ has no $osp$ description.
- *Not even for the full spectrum.* The natural candidate is multiplication by the highest-weight column minor $(\alpha_{-1}\wedge\dots\wedge\alpha_{-k})\cdot(\varepsilon_1\wedge\dots\wedge\varepsilon_k)$. It does not commute with the constraints.
  - Already $[G_{1/2},\varepsilon\cdot\alpha_{-1}]=\varepsilon\cdot b_{-1/2}$.
  - This kills principally embedded states but not dressed ones; for example, their depth-$\tfrac12$ clone of the spin 2.
  - So the dressing must change with the row lengths, and single-box moves beyond row and column 1 are genuinely non-monotone.
- *What it does give.*
  - The natural coordinates $(K,w)$.
  - The reason column moves have uniform triangular costs: they lengthen only the arms of the first $k$ hooks.
  - A plausible origin of GMSZ's stabilization: the dressing coefficients depend smoothly on the lengths, and degenerate only for short arms.
- The core problem is unchanged from S6/S7: an SO(9)-covariant model of the super-representations with manifest positivity.

### S11. Half (b): liftability of the residual $\rho_n$ (2026-10-12) (`half_b.py`)

> **Correction (2026-10-12, `beyond32.py`, data to $n=40$).** The $l_5$-extremal property below holds only up to $n=33$.
> - From $n=34$ on, $\rho_n$ (max-height induction) has **no** $l_5$-extremal lift: the global ILP is infeasible for every $n=34..40$.
> - The no-lookahead level recursion fails accordingly.
> - A lift of $\rho_{34}$ needs at least 22 487 classes with intermediate $l_5$, all of shape $(l_1,l_2,2,2,1)$, Dynkin $(a_4,a_5)=(1,3)$. At $n=35$ and $36$ it needs 134 661 and 370 181 of the same family.
> - So the repairs climb to the next level as $n$ grows. The clean two-valued structure below is a low-level artefact. The pattern held for every $n\le33$ and then changed, which is exactly the trap `previous_efforts.md` warns about.
> - The Regge induction itself (unrestricted $l_5$) is unaffected: it runs to $n=40$ (S12).

**Set-up.**
- A class-level lift of an SO(9) character $\rho$ is a non-negative sum of interlacing boxes.
- Processing weights from the top down, the height-maximal remaining weight must be the top $(l_1,\dots,l_4)$ of whatever lifts it.
- So a lift is fixed by one choice per top: $l_5\in[\text{floor},\mu_4]$, where the floor is 0 for tensors and ½ for spinors.

**Local rules fail.** Three greedy lifts fail from $n=6$–$10$, for $\rho_n$ under both inductions and for the cumulative $R^9_n$:
- always $l_5=\mu_4$;
- always $l_5=$ floor;
- an adaptive profile that sends as many boxes down the $\mu_4$ column as the remainders allow.

The choices of $l_5$ need global coordination.

**$l_5$ can always be taken extremal.**
- Every $\rho_n$ ($n\le32$, both inductions) and every $R^9_n$ ($n\le24$) has a lift using only $l_5\in\{\text{floor},l_4\}$. In the $\mu_4$ direction each box is then either "full" or "degenerate".
- This is a property of the string data, not of liftability in general: the extremal cone is a proper subcone. For example, the box of $(2,2,2,2,1)$, $\{(2,2,2,1),(2,2,2,2)\}$, has no extremal lift.
- $l_3$ cannot be restricted the same way: $l_3\in\{l_2,l_4\}$ fails from $n=8$.

**Exact level recursion for extremal lifts.** Slice by $t=\mu_4$, from the top down.
- At level $t$ the tops $(\mu_1,\mu_2,\mu_3)$ are the *unique* lift $T_t$ of the level remainder. It is given by explicit third finite differences, because the level restriction is unitriangular.
- The only choice is the degenerate part $0\le D_t\le T_t$ (with $l_5=t$); the remaining tops get $l_5=$ floor.
- The next level is $T_{t-1}=\mathrm{Lift}(\rho_{t-1}-\rho_t)+M(D_t)$, which depends only on $\rho_{t-1}$, $\rho_t$ and $D_t$. One degenerate top $(m_1,m_2,m_3)$ contributes:
  - $+[m_1,m_2,m_3]$;
  - at the boundary tops ($\mu_3=t-1$): $-[m_1,m_2]+[m_1,m_3-1]-[m_2-1,m_3-1]$.
- So an extremal lift exists iff some sequence $D_t$ keeps every $T_t\ge0$.
- Example: the unique $E_8$ needs both mechanisms. One copy of the top $(1,1,1,1)$ is degenerate (126) and one is full (210), to cancel a $-1$ of $\mathrm{Lift}(\rho_0-\rho_1)$ at $(1,1,1)$. The top $(2,1,1,1)$ is degenerate (1050), to cancel a $-1$ at the boundary top $(2,0,0)$ via $+[m_1,m_3-1]$.

**Results of the recursion.**
- Rules with no degenerate tops, all degenerate, or degenerate only where the non-boundary entries force it all fail. The first two reproduce the restricted-family results.
- A per-level ILP choice of $D_t$ (no lookahead; minimising or maximising the number of degenerate tops) **succeeds for every $\rho_n$ of the max-height induction ($n\le32$) and every cumulative $R^9_n$ ($n\le24$).**
- Under the min-height induction it hits dead ends at $n=29$–$32$, although extremal lifts exist there; those levels need lookahead.
- So, for the max-height rule, half (b) reduces (empirically) to a sequence of small per-level problems whose data are finite differences of $\rho_n$.
- The same reduction works for the cumulative $R^9_n$ directly, bypassing the Regge induction.

**Where the choices matter** (max-height induction, $n\le32$; part (4) of the script). Start from the *forced* degenerate tops, $D_t\ge-\mathrm{Lift}(\rho_{t-1}-\rho_t)$ at non-boundary tops.
- The forced degenerate tops never exceed the available tops. Non-boundary negatives are always fixable (0 of 96 level transitions fail either way).
- The boundary tops that remain negative all sit at the **bottom transition**:
  - tensors, $\mu_4\colon1\to0$: at $(k,0,0,0)$, i.e. the symmetric tensors $[k,0,0,0,0]$ (1 538 192 units), and at $(k,1,0,0)$ (99 776);
  - spinors, $\mu_4\colon\tfrac32\to\tfrac12$: at $(k,\tfrac12,\tfrac12,\tfrac12)$, i.e. $[k-\tfrac12,0,0,0,1]$ (82 987).
- All but 42 of the 1 720 997 extra degenerate units repair these directly, through $+[m_1,m_3-1]$: degenerate tops $(k,m_2,1,1)$ feed the symmetric tensor $[k,0,0,0,0]$.
- So, as with half (a) (whose tight spots were the leading trajectory), everything non-trivial in half (b) sits next to the leading-trajectory irreps $[k,0,0,0,0]$ and $[k,0,0,0,1]$, at the lowest $\mu_4$ level. Higher levels need only the forced choices.
- **For the cumulative $R^9_n$ itself** ($n\le24$) the picture is cleaner still.
  - The only negatives left after the forced choices are at the symmetric tensors $[k,0,0,0,0]$, at the tensor $\mu_4\colon1\to0$ transition: 79 956 units, all repaired directly.
  - There are no spinor repairs, no $(k,1,0,0)$ repairs and no indirect extras.
  - So existence of an uplift of $R^9_n$ reduces (empirically) to: forced choices available at every level, plus symmetric-tensor repairs at the bottom. This route does not use the Regge induction.

**Next.**
- Prove the two per-level statements found here:
  - (i) at every level, the forced degenerate tops are available, i.e. $T_t+\mathrm{Lift}(\rho_{t-1}-\rho_t)\ge0$ at each top;
  - (ii) at the bottom transition, the negatives at $[k,0,0,0,0]$ (and $[k,0,0,0,1]$) can be repaired by degenerate tops $(k,m_2,1,1)$ without creating new negatives.
- Both are explicit inequalities in finite differences of $\rho_n$, and (by the same recursion) of the cumulative $R^9_n$.

### S12. The cumulative route, and the data to $n=40$ (2026-10-12) (`cumulative_lift.py`, `beyond32.py`)

Idea: lift $R^9_n$ directly, through the S11 level recursion, with no Regge induction.

**The all-full lift.**
- The classes with minimal $|l_5|$ (0 for tensors, ½ for spinors) restrict to a unitriangular basis. So $R^9_n$ has a unique lift by them, $a(l)=\Delta^4R(l)=\sum_{\varepsilon\in\{0,1\}^4}(-1)^{|\varepsilon|}R(l+\varepsilon)$.
- Per SO(9) irrep this is $V_\mu\mapsto\sum_\varepsilon(-1)^{|\varepsilon|}V_{(\mu-\varepsilon,0)}$.
- Degenerate classes ($l_5=l_4$) repair its negatives level by level, as in S11.

**Up to $n=33$ the picture is clean.**
- **(A)** $a<0$ only at the bottom level, $l_4=0$ (or ½ for spinors). Spinors have no negatives at all until $n=28$.
- **(B)** The repairs need only degenerate classes at $l_4=l_5=1$ (or $\tfrac32$).
  - Forced coverage holds. The tightest case sits at a fixed distance from the top, column $(n-14,2,2)$, with slack 132 for every $n\ge24$.
  - The bottom deficits are fed along the leading-trajectory chains $[k,0,0,0,0]$, $(k,1,0,0)$ and $[k,0,0,0,1]$.
  - Their tail sums from the top are $n$-independent: $1,1,1,1,2,3,3,2,-4,-21,\dots$ for $[k,0,0,0,0]$. So the prefix condition reads $R^9_n(\text{singlet})\ge3$ plus small corrections.
- Uplifts of $R^9_n$ found this way for every $n\le33$, in seconds, use only $l_5$-minimal classes plus $[a_1,a_2,a_3,0,2]$ (and $[a_1,a_2,a_3,0,3]$ from $n=28$).

**It breaks at $n=34$.**
- The forced degenerate tops in columns $(m_1,1,1)$ drain the singlet node of the symmetric-tensor chain. Its prefix sum goes negative ($-60\,960$ at $n=34$, $-1.0\times10^6$ at $n=40$).
- From $n=35$, forced coverage itself fails, at column $(5,1,1)$: slack $-77\,285$ at $n=35$, $-4.7\times10^6$ at $n=40$.
- So "all-full lift plus bottom repairs" is not an all-orders mechanism.

**Bosonic analogue** (`--bosonic 22`). For the cumulative bosonic spectrum $B$ (S6), the all-full lift is non-negative at degrees 7, 8 and 10–22. At degrees 4, 5, 6 and 9 it has a single $-1$ at the bottom level. So for $B$ no repairs are needed in this range, and the bottom-level negatives of the superstring are a fermionic effect.

**Existence to $n=40$.**
- `so10_uplift/derive_r9.py --max-level 40` gives `data/r9_to40.json` (8 min). It is identical to `r9_to32.json` for $n\le32$.
- With it, **the Regge induction runs to $n=40$ under both the max-height and the min-height rules.** Half (a) holds at every step, half (b) is solved by unrestricted ILP, and every $A_n$ is verified by restriction.
- $|A_{40}|=1\,457\,936\,065$ for both rules, as S1 requires.
- So non-negative uplifts of $R^9_n$ exist for all $n\le40$; previously this was established for $n\le32$.
- The S9 column monotonicities also hold to $n=40$ (`beyond32.py` (4)): $(\omega_1,+1)$, which is Regge monotonicity and the input to half (a), $(\omega_2,+3)$, $(\omega_3,+6)$ and $(2\omega_4,+10)$, with 0 violations. $(\omega_1,+2)$ fails only at the 38 gap entries $[n-2,0,0,0]$.

### S13. Mapping the repair levels (2026-10-12) (`repair_levels.py`, `repair_asymptotics.py`, `repair_capacity.py`)

**Definitions.**
- A class is *full* when $l_5$ is minimal (0 for tensors, ½ for spinors).
- A non-full class skips the bottom $J=l_5-\text{floor}$ levels of $\mu_4$ and reaches up to level $L=l_4-\text{floor}$.
- For each residual $\rho_n$:
  - $K_L(n)$ is the least $K$ such that $\rho_n$ lifts with non-full classes only at levels $L\le K$;
  - $K_J(n)$ is the least $J$ such that every class can have skip depth at most $J$.
- The onset $n_K$ of repair level $K$ is the first $n$ with $K_L(n)\ge K$.

**The map.** Data: both rules for $n\le40$; max-height also for $n\le44$, using data from `so10_uplift/derive_r9.py --max-level 44` (15 min, not stored in the repo).
- **Max-height rule.**
  - No repairs at all for $n\le5$ and $n=7$.
  - Level 1 from $n=6$.
  - Level 2 from $n=34$, through 44. No level 3 by $n=44$.
  - $K_J=1$ throughout: no lift ever needs to skip more than one bottom level.
  - Level-2 repairs are all of type $(L,J)=(2,1)$: 22 487 at $n=34$, 5.9 M at 40, 51.9 M at 44.
- **Min-height rule.** Level 2 is already needed from $n=31$, and its minimal lifts also use the type $(2,2)$. So the onsets depend on the rule.

**Mechanism (exact).** Switch a level-$L$ class to skip $\mu_4=0$, i.e. $(m_1,m_2,m_3,L,1)$ instead of $(\dots,L,0)$. This changes the bottom-level coefficients by:
- $+1$ at its own column $(m_1,m_2,m_3)$;
- $-[m_1,m_2,L-1]+[m_1,m_3-1,L-1]-[m_2-1,m_3-1,L-1]$ on bottom columns with third coordinate $L-1$.

So a bottom deficit at a column $(k,c,c)$ is covered only by same-column repairs at levels $1..c$, or by level-$(c+1)$ feed through $+[m_1,m_3-1,L-1]$. This gives a local indicator for the cumulative $R^9_n$:
$$e_c(n,k)=-a(k,c,c,0)-\sum_{L\le c}a(k,c,c,L).$$
- $c=0$: positive from $n=9$. The symmetric-tensor chain has no same-column capacity.
- $c=1$: the relative excess at bulk columns $(4..8,1,1)$ rises in accelerating steps and crosses 0 at $n\approx34$–35; this is the level-2 onset. For example, at $k=5$ it goes $-0.96$ ($n=20$), $-0.02$ (34), $+1.98$ (44).
- $c=2$: $\approx-0.58$ at $n=44$, rising steadily by about 0.009 per level. A linear extrapolation of this single column puts the level-3 onset near $n\approx100$.
- $c=3$: slower still.

**The onset is local, not aggregate.** When level 2 first became necessary, the minimum number of level-1 repairs over all lifts was only 57% of the maximum possible (0.41 → 0.57 over $n=16..33$). For level 2 the ratio is 0.011 → 0.130 over $n=34..40$.

**Asymptotics at fixed weights.**
- For fixed weights, multiplicities become proportional to SO(9) dimensions (Curtright–Thorn). So $a(l)\sim C(n)\,D_4(l)$ with $D_4(l)=\sum_\varepsilon(-1)^{|\varepsilon|}\dim(l+\varepsilon)$.
- $D_4(l)\ge70$ for every dominant $l$ with $l_1\le20$. The minimum is at the origin: $1-9+36-84+126$.
- So at fixed weights the all-full lift is *eventually* positive and needs no repairs: the deficits are pre-asymptotic.
- The approach is slow. At $n=44$ some columns are still moving away from their asymptotic value; for example $(1,1,1)$ is at $-4.2$ against an asymptotic $-3.3$.

**Reading.**
- The repairs live in a transition region between two well-behaved regions: the stable range near the leading trajectory (exactly structured, $n$-independent) and fixed weights (asymptotically positive).
- On the evidence, the needed repair level grows with $n$: level 2 from $n=31$–34, level 3 plausibly near $n\sim100$. So no finite repair pattern covers all $n$.
- A proof of half (b), or of existence directly, will need uniform asymptotics in the transition region together with an exact treatment of the stable range.

### S14. Locating the transition region (2026-10-12) (`transition_region.py`)

Cumulative $R^9_n$, tensor sector. Data to $n=40$ (repo) and to $n=44$ (regenerated).

**The spectrum is not dimension-proportional at these levels.**
- Along the symmetric tensors $(k,0,0,0)$, $\ln(R^9_n/\dim)$ falls almost linearly in $k$. The slope over $k=4..12$ is $-1.15$ at $n=24$ and $-1.00$ at $n=44$; further out it is about $-0.8$. It shrinks only slowly with $n$.
- Even at $k=1$, $R/\dim$ is $e^{-0.75}$ times its value at the singlet ($n=44$).
- So the regime where the asymptotic positivity $D_4(\dim)\ge70$ applies (S13) is not reached for any weight at these $n$.
- A heat-kernel fit $\ln(R/\dim)=A-\beta C_2$ over the bulk ($\nu_1\le0.9\sqrt n$) is only rough: $R^2=0.92$–$0.95$. The fitted $\beta\approx0.85/\sqrt n$, so the spin scale is about $n^{1/4}$, not $\sqrt n$.

**The negative region (bottom level only).** For each column family $(j,i)$, the coefficients $a(k,j,i,0)<0$ form one contiguous range of $k$:
- it starts near the dominance edge, at $k\approx\max(j,4)$;
- it ends at a **fixed distance $d(j,i)$ below the top**, the same for every $n$ once the family appears.

| $(j,i)$ | (0,0) | (1,1) | (2,2) | (3,2) | (3,3) | (1,0) | (2,1) | (4,4) | (5,5) | (6,6) |
|---|---|---|---|---|---|---|---|---|---|---|
| $d(j,i)$ | 8 | 15 | 14 | 17 | 19 | 22 | 26 | 26 | 31 | 36 |

So the negatives fill the bulk, and the stable range cuts them off.

**The level-2 band.**
- Only the family $(1,1)$, the classes $[k-1,0,1,0,0]$, ever has excess beyond its own levels.
- It first appears at $n=35$ ($k=5,6$). For every $n\ge37$ it is exactly $k\in[4,\,n-30]$: a fixed lower edge and a fixed distance from the top.
- So the band widens by one column per level and fills the middle.
- Up to $n=44$, every family with $i\ge2$ is covered by its own levels.

**Reading.**
- The transition region is the whole middle band between a fixed-weight edge and the stable range: the Regge-scaling regime, where spin grows with the level ($k\approx xn$). It is not a thin shell at a thermal spin scale.
- A uniform analysis has to work in that regime. There, GMSZ's shift relations and support costs control the structure: multiplicities depend on $n$ minus the cost of the weight.
- Asymptotic positivity says the lower edge must eventually recede at fixed weights, but this has not started by $n=44$.

### S15. The Regge limit: a closed-form model of the top of the spectrum (2026-10-12) (`regge_limit.py`, `regge_scaling.py`)

`regge_limit.py` derives the model from the light-cone Fock space and reads no spectrum data. `regge_scaling.py` compares it with the data (to $n=40$ in the repo, to $n=44$ regenerated). Results: `results/regge_limit.txt`, `results/regge_scaling.txt`, `results/regge_scaling_44.txt`.

**Coordinates from GMSZ's costs.** Write an SO(9) weight as $\nu=(x_1,\mu)$ with $\mu=(a,b,c)$ an SO(7) weight, and $\ell=x_1-a$.
- At level $n$ the top distance is $D=n-x_1$.
- GMSZ's cost splits as $n_0(\nu)=x_1+{\rm cost}_7(\mu)$, with ${\rm cost}_7(a,b,c)=2a-3b+c(b-c,2c)$.
- So the GMSZ depth is $w=n-n_0=D-{\rm cost}_7(\mu)$.

**The model.** As $x_1\to\infty$ at fixed $D$, $r^9_n(x_1,\mu)\to m_D(\mu)$, the multiplicity of the SO(7) irrep $\mu$ in $[q^D]Z_{\rm Regge}$, where

$$Z_{\rm Regge}=q\prod_{k\ge2}\frac{1}{(1-q^k)\det_7(1-q^kg)}\prod_{k\ge1}\det_8\big(1+q^{k+1/2}g\big).$$

This is a positive Fock space: bosons in $1+7$ at levels $k\ge2$, and SO(7)-spinor fermions at levels $k+\tfrac12\ge\tfrac32$.

Derivation:
1. On the common torus, the SO(9) character of the level-$n$ states $F_n=256\otimes r^9_n$ is the light-cone character: bosons $8_v$, fermions $8_s$, ground states $8_v+8_c$.
2. Restrict to SO(2)$_{12}\times$SO(7), with charge $z=x_1$. The Koszul transform $K(\chi)=\det_7(1-z^{-1}g)\,\chi$ extracts irreducible multiplicities exactly. By the parabolic BGG resolution, for every dominant $(N,\mu)$, $m(N,\mu)$ is the coefficient of $z^N\chi_\mu$ in $K(\chi)$; all other terms have SO(7) highest weights with first entry $\ge N+2$. Since $K$ is multiplicative, $K(F)=256|\cdot K(r^9)$.
3. Put $z=y/q$. At each order in $q$, $\sum_nq^nK(F_n)$ is a Laurent polynomial in $y$ divided by $1-y$, the pole of $\alpha^+_{-1}$. So the coefficient of $y^N$ is constant for large $N$, equal to the residue at $y=1$: this is the Regge limit.
4. Take the residue. The Koszul factor removes the level-1 SO(6) vector and one scalar. What remains of the bosons is two scalars and the SO(6) vector at each level $k\ge2$, i.e. $1+7$. The fermions of charge $\pm\tfrac12$ pair into $8$ at depths $j+\tfrac12\ge\tfrac32$.
5. The ground states times the level-1 fermions of charge $+\tfrac12$ give exactly $q\,T(q)$, with $T(q)=256|_{z=1/q}$ (`regge_limit.py` (0), checked on weights). So the 256 cancels and leaves $Z_{\rm Regge}$.

**GMSZ's support costs follow from the Fock space** (`regge_limit.py` (2), to $D=60$).
- For every SO(7) irrep with $\mu_1\le10$ and cost $\le60$, the least $D$ with $m_D(\mu)>0$ is exactly ${\rm cost}_7(\mu)$.
- The only zero at $D\ge{\rm cost}_7$ is $\mu=0$ at $D=2$: GMSZ's gap $[n-2,0,0,0]$.
- So $c(r,s)$ is the lowest energy of the irrep $[m,r,s]$ in this Fock space.

**Data against the model** (`regge_scaling.py`; every statement holds on both data sets).
- 321 101 comparisons ($n\le44$, $\mu_1\le10$). In all of them $r^9_n\le m_D$, and $r^9_n$ is non-decreasing in $\ell$ at fixed $(D,\mu)$.
- $r^9_n=m_D$ exactly whenever $\ell\ge w$, with zero exceptions. Saturation sets in at $\ell=w$ (2151 cases) or $\ell=w-1$ (21 cases). So GMSZ's shift relation (3.32), iterated, ends exactly on this model.
- The same holds for the cumulative $R^9_n$ against $G_D=\sum_{D'\le D}m_{D'}$.
- Below saturation (the ramp), data/model reaches $\tfrac12$ at $\ell\approx2$–$4$ and $0.99$ at $\ell\approx6$–$12$, both growing slowly with $w$.
- The all-full lift agrees exactly near the top: $a(n-D,j,i,L)=a_D(j,i,L)$ (all $j\le6$, all $L$) for every $D\le\lfloor n/2\rfloor+1$ and every $n\le44$. The first difference is at $D=\lfloor n/2\rfloor+2$.

**The upper edges of the transition region are Regge-limit quantities.**
- The model's offsets ($d(j,i)$ = least $D$ with $a_D(j,i,0)<0$) are 8, 22, 15, 26, 14, 17, 19, 26, 31, 36 for $(0,0),(1,0),(1,1),(2,1),(2,2),(3,2),(3,3),(4,4),(5,5),(6,6)$. These are exactly S14's measured offsets; the data reach all of them by $n=43$.
- In the model each family stays negative for every $D$ from $d(j,i)$ to 60.
- The model's $(1,1)$ excess is $e_{29}=-42\,890$ and $e_{30}=+89\,517$, and $e_D>0$ for all $30\le D\le60$. So the band's upper edge $n-30$ is the sign change of a closed-form, $n$-independent quantity.
- The data at fixed $D$ converge to $e_D$. At $D=24$: $-43\,588$ at $n=44$, against $-43\,589$.
- The lower edge $k=4$ is not a Regge property. There $\ell=3$, which is about where the ramp of $\mu=(1,1,0)$ reaches half the Regge value.

**Predictions beyond $n=44$, near the top.**
- Family $(2,1)$ has excess for $D\ge50$: a second family that needs level 2. In the data it should appear only from $n\approx55$, depending on the ramp.
- No family with $i\ge2$ has excess for $D\le60$. So near the top, the same-column count never forces level 3 up to $D=60$. S13's level-3 extrapolation concerns bulk columns at small $k$, which the model does not cover.
- More families turn negative, e.g. $(2,0)$ from $D=49$. The offsets grow like costs: $d(j,i)=2j+3i+6$ for $4\le i\le j\le8$.

**The Regge-limit uplift problem is solved.**
- $Z_{\rm Regge}$ is the restriction of the SO(8) Fock space

$$Z_8=q\prod_{k\ge2}\det_{8_v}(1-q^kg)^{-1}\prod_{k\ge1}\det_{8_s}\big(1+q^{k+1/2}g\big),$$

because $8_v\to7+1$ and $8_s\to8$. So every Regge slice $\Phi_D=[q^D]Z_{\rm Regge}$ is the restriction of a non-negative SO(8) character. `regge_limit.py` (5) checks this irrep by irrep for $D\le20$.
- In SO(10) terms, the classes $(n-D,\lambda_8)$, $\lambda_8\in[q^D]Z_8$, restrict to the cumulative $R^9$ in the Regge limit.
- This explicit lift uses non-full classes $(L,J)$ (S13's repair types): $(1,1)$ from $D=5$, $(2,1)$ from $D=12$ and $(2,2)$ from $D=13$, with rapidly growing counts (at $D=20$: 26 903, 1 912 and 280). So it is a valid but costly lift. The same-column count forces level 2 only from $D=30$, and the data's minimal lifts never skip two bottom levels ($K_J=1$, S13).

**What this means for half (b)** (derived from the above, not computed separately; corrected 2026-10-13).
- In the Regge induction, $T$ keeps the top distance fixed. Also $\rho_n$ vanishes wherever $r^9_n$ equals the model and the slices of $A_{n-1}$ already restrict to it.
- So the classes at a top distance $D$ are frozen once $D$ is in the wedge $\ell\ge w$, and their SO(8) parts restrict to $\Phi_D$. That is a necessary condition on the induction, and $Z_8$ shows it can be met for every $D$.
- **Correction.** The first version said that near the top half (b) *is* the problem of lifting each $\Phi_D$, and that $Z_8$ solves it. That overstated it. The slice condition says nothing about *when* each class is created, and that is where the ramp enters.
  - Think of the SO(9) multiplicities at fixed (shape, top distance) as counting trajectories. Each trajectory starts at some spin and continues forever (Regge monotonicity), and the ramp is the distribution of starting spins.
  - A class created at the first level where it fits ($l_1=l_2$) starts all its trajectories at spin $l_2$. One of them has a shape $(l_2,\mu_3,\mu_4)$ and starts at its own dominance edge, and every such shape costs at least $3l_2+1$ levels (cheapest $[0,l_2,0,0]$). So these classes start trajectories only at spin $N\le(D-1)/2$.
  - The data have later starts: up to $N=D-2$ for the singlet (at $w=20$ the starts are $1,1,2,5,\dots$ at $\ell=19,18,17,16,\dots$). They must come from classes created after their first fit, which at their creation level fill a whole column segment $x_1\in[l_2,l_1]$ of $r^9_n$.
  - In particular, $Z_8$'s classes carried to finite $n$ from their first fit saturate by spin about $D/2$, while the data saturate only at $\ell=w$. So they overshoot $R^9_n$ in the ramp, and $Z_8$ is not by itself a finite-$n$ lift.
- The open part of half (b) is therefore the creation schedule. It lives in the ramp, $\ell<w$, roughly $D>n/2$, which is also where the residuals are. There the data stay below the Regge values in a controlled way: monotone in $\ell$, and exactly saturating at $\ell=w$.

### S16. The ramp: the spin-resolved function $Q(q,y)$ (2026-10-13) (`q_spin.py`, `q_spin_check.py`, `q_top_end.py`)

`q_spin.py` and `q_top_end.py` derive everything from the light-cone Fock space and read no spectrum data. `q_spin_check.py` compares with the data. Results: `results/q_spin.txt`, `results/q_spin_check.txt`, `results/q_top_end.txt`; data `data/q_spin_D22.json`.

**Trajectories and starts.** Fix a shape $\mu$ (the SO(7) part) and a top distance $D=n-x_1$. Along a Regge trajectory (spin $N$ and level both $+1$) the multiplicity counts the trajectories that have started by spin $N$. The *start* at spin $N$ is $r^9_n(N,\mu)-r^9_{n-1}(N-1,\mu)$, the number of new trajectories. Regge monotonicity (half (a)) is non-negativity of every start. The plateau of S15 is the total number of trajectories, and the ramp is their cumulative count.

**The function.** Let $z$ be the $x_1$-charge and $K(\chi)=\det_7(1-z^{-1}g)\,\chi$ the Koszul transform of S15, and define

$$Q(q,y)=(1-y)\sum_{n\ge1}q^n\,K(r^9_n)\big|_{z=y/q}.$$

By the parabolic BGG argument of S15 (exact at every dominant weight, not only at large spin), with $c(D,e;\mu)$ the multiplicity of $\mu$ in $[q^Dy^e]Q$:

$$r^9_{N+D}(N,\mu)=\sum_{e\le N}c(D,e;\mu).$$

So $c(D,N;\mu)$ is the start at spin $N$ (for $N>\mu_1$), the sum over $e$ is the plateau ($Q(q,1)=Z_{\rm Regge}$), and $Q$ to depth $D_{\max}$ gives the top $D_{\max}$ units of *every* level.

**Construction.**
- Put $z=y/q$ in the light-cone character. The $\alpha^+_{-1}$ tower becomes $1/(1-y)$.
- Every other oscillator carries its charge as a power of $y$: $\alpha^+_{-k}\to q^{k-1}y$, $\alpha^-_{-k}\to q^{k+1}/y$, transverse $\alpha^i_{-k}\to q^k$, fermions $\to q^{k\mp1/2}y^{\pm1/2}$.
- Level 0 is removed, since it is not SO(9)-covariant.
- Divide by $256(z=y/q)$. Unlike at $y=1$, this division is not a clean cancellation.
- Arithmetic is exact modulo two primes, combined by CRT. It is a 5-index computation (depth, charge, SO(7) weight); depth 22 took about 12 min and 482 MB per prime.

**Checks (depth $D\le22$, every shape).**
- Every $(D,e)$ slice of $Q$ is $W(B_3)$-invariant, as it must be. The intermediate factors are only SO(6)-covariant, so this is a real test.
- $\sum_e c(D,e;\mu)$ equals S15's $Z_{\rm Regge}$ for all $D\le22$ and every shape.
- Against the data ($n\le40$): all 1 224 614 multiplicities $r^9_n(N,\mu)$ with $D\le22$ (17 750 of them non-zero) and all 1 135 476 starts agree, with zero exceptions. The $n\le44$ data were in scratch and were lost in a reboot; they were not re-checked.
- Example: the singlet at $D=21$ ($w=20$) has edge value 544 and starts $2241,4122,4730,4057,2898,1865,1130,659,375,211,116,63,34,18,9,5,2,1,1$. Its plateau is 23 081, exactly the ramp of S15.

**Results valid at every level $n$ (top distance $D\le22$).**
- **Regge monotonicity holds:** all 4 382 non-zero starts are positive, and all 858 edge values are non-negative. So half (a)'s input, $r^9_n(\nu+e_1)\ge r^9_{n-1}(\nu)$, holds for every $n$ at every weight within 22 of the top. This is a finite computation, but the statement covers infinitely many levels.
- **GMSZ's shift relation (3.32) is exactly saturation:** for all 863 (depth, shape) pairs, the last start is at $\ell=w$ (843 pairs) or $\ell=w-1$ (20), with no exception at any level.
- **Shape of the starts:** a hump at low spin (peak at $\ell\approx2$–$4$), then a tail that ends in a sequence independent of $D$.

**The top end of the ramp is the residue at the $\alpha^+_{-2}$ pole** (`q_top_end.py`).
- Put $y=v/q$ and grade by $g=D-e$. The $\alpha^+_{-2}$ factor becomes $1/(1-v)$, so at fixed $g$ the starts converge to a residue $\hat Q(q)$. It is computed like S15 (one variable plus SO(7) weights).
- $\alpha^+_{-2}$ has weight $e_1$ at level 2: the shift of S8's secondary direction $T_2$.
- **Exact range:** $c(D,e;\mu)=\mathrm{mult}_\mu[q^{D-e}]\hat Q$ for every start with $2e-D\ge-1$, i.e. spin at least $(D-1)/2$. This holds in all 749 rows. The last disagreement is at $2e-D=-2$ (singlet, $D=18,20,22$). So the upper half of every ramp is closed-form.
- $\hat Q$ multiplicities from the first non-zero $g$:
  - singlet: $1,1,2,5,9,18,34,63,116,211,376,667,\dots$;
  - vector: $1,1,3,6,15,28,60,\dots$;
  - $(1,1,0)$: $1,1,5,10,24,50,108,\dots$;
  - spinor: $1,1,3,8,16,33,68,\dots$
- $\hat Q$ is non-negative for all shapes and $g\le23$.
- **$\hat Q$ is not the obvious Fock space.** Its dimension series factorises like one, $16\,u^2(1+u)^7\prod_{j\ge2}(1+u^j)^8/\big(\prod_{j\ge2}(1-u^j)^7\prod_{j\ge3}(1-u^j)\big)$, starting with $1+7+8$ (the ground states). But the Fock space $q^2(1+7+8)\det_7(1+qg)\prod_{j\ge2}\det_8(1+q^jg)/\det_7(1-q^jg)\prod_{j\ge3}(1-q^j)^{-1}$, which has those dimensions, differs from $\hat Q$ at $g=4$. **S17 finds the positive formula:** this guess had the right dimensions but the wrong split.

**Reading.**
- The ramp is no longer a black box.
  - Its total is S15's positive Fock space (residue at $y=1$, the $\alpha^+_{-1}$ tower).
  - Its upper half is a second closed form ($\alpha^+_{-2}$).
  - Each $\alpha^+_{-(j+1)}$ has depth $j$ per unit of spin. So the next poles should govern successively lower bands of spin, roughly $e\approx D/j$: a multi-Regge structure. Only the first two are confirmed.
- Half (a) at all orders is now a positivity statement about one explicit generating function, $c(D,e;\mu)\ge0$ for $e>\mu_1$. It holds for $D\le22$ at every level.
- On the upper half of each ramp it reduces further to $\hat Q\ge0$, a one-variable statement. $\hat Q$ is cheap to extend, but only to $g\le23$ so far.

### S17. The Regge directions $T_k$ and their residues (2026-10-13) (`q_multi_regge.py`)

`q_multi_regge.py` derives the residues from the light-cone Fock space and reads no spectrum data. It compares them with $Q$ (S16, also first-principles). Results: `results/q_multi_regge.txt`.

**The hierarchy.** Write $T_k$ for the shift (spin $+1$, level $+k$). In top-distance coordinates it is $(D,N)\mapsto(D+k-1,N+1)$, i.e. multiplication by $q^{k-1}y$, the monomial of the light-cone oscillator $\alpha^+_{-k}$. For each shape $\mu$ and dominant spin $N$, define:
- $s(D,N)=r^9(D,N)-r^9(D,N-1)$: the starts ($T_1$-differences, S16);
- $O(D,N)=s(D,N)-s(D-1,N-1)$: the $T_2$-differences of the starts;
- $L_3(D,N)=O(D,N)-O(D-2,N-1)$: the $T_3$-differences.

The stable total along $T_k$ is the residue $R_k$ of the light-cone generating function at the $\alpha^+_{-k}$ pole $y=q^{-(k-1)}$, graded by $g=D-(k-1)N$. So $R_1=Z_{\rm Regge}$ (S15), $R_2=\hat Q$ (S16), and $R_3=\hat O$ (here).

**Every residue is a positive Fock space.** This is derived below, and checked against the residue computed directly from the light cone: for $k=1,\dots,5$ (all shapes, $g\le12$), for $k=2$ to $g\le22$, and for $k=3$ to $h\le22.5$.

$$R_k=q^{p_k}\prod_{1\le m<k/2}\det_8\big(1+q^{k/2-m}g\big)\,\Big[\prod_{w\in4}(1+x^w)\Big]_{k\ \mathrm{even}}\prod_{j=1}^{k-1}\frac{1}{\det_7(1-q^jg)}\prod_{j\ge k+1}\frac{1}{(1-q^j)\det_7(1-q^jg)}\prod_{j\ge k/2+1}\det_8\big(1+q^jg\big)$$

- The exponent is $p_k=k-4\sum_{1\le m<k/2}(k/2-m)$, so $p_k=1,2,1,0,-3$ for $k=1,\dots,5$, matching the lowest grades found.
- The fermion grades in the last product run in steps of 1 from $k/2+1$.
- $\prod_{w\in4}(1+x^w)=1+7+8$ as an SO(7) character: the ground states.
- *Derivation.*
  - At $y=q^{-(k-1)}$ (so $z=q^{-k}$), the identity of S15 applies to the level-$k$ fermions: ground states $\times$ one level of charge-$\tfrac12$ fermions $=z^{-1}\,256(z)$. So the 256 cancels exactly.
  - The Koszul factor $\det_7(1-q^kg)$ removes the grade-$k$ bosons.
  - The charge-$\tfrac12$ fermions of levels $m<k/2$ have negative grade. Flipped, they become $\bar 4$'s and pair with level $k-m$ into SO(7) spinors.
- $R_2=\hat Q=q^2(1+7+8)\det_7(1-qg)^{-1}\prod_{j\ge3}[(1-q^j)\det_7(1-q^jg)]^{-1}\prod_{j\ge2}\det_8(1+q^jg)$.
  - This corrects S16. The guess there had the right dimensions, but put a fermionic 7 at $g=1$ where there is a bosonic one, and missed the boson gap at $g=2$.
- $R_3=\hat O=q\det_8(1+q^{1/2}g)\,[\det_7(1-qg)\det_7(1-q^2g)]^{-1}\prod_{j\ge4}[(1-q^j)\det_7(1-q^jg)]^{-1}\prod_{m\ge5/2}\det_8(1+q^mg)$.
- So every stable total in the hierarchy is manifestly non-negative. For $k\ge2$, though, $R_k$ contains lone 7's (bosons at grades $1,\dots,k-1$), so it is not obviously the restriction of an SO(8) character, unlike $R_1$.

**Exact ranges** (all shapes, $D\le22$).
- $s(D,N)=\hat Q_{D-N}$ for $2N-D\ge-1$ (S16).
- $O(D,N)=\hat O_{D-2N}$ at every interior point ($N>\mu_1$) with $3N-D\ge2$. Edge points of large shapes need more.
- The two-pole formula $s(D,N)=\hat Q_{D-N}-\sum_{j\ge1}\hat O_{D-2N-j}$ holds for every start with $3N-D\ge-2$.
- So successive poles govern successive bands of spin: $N\gtrsim D/2$ by $\alpha^+_{-2}$ alone, $D/3\lesssim N\lesssim D/2$ by $\alpha^+_{-2}$ and $\alpha^+_{-3}$, and so on.

**Higher Regge monotonicity** (all shapes, $D\le22$, every level).
- **The starts are $T_2$-monotone.** Of 1 446 non-zero $T_2$-differences, exactly one is negative: $O(2,1)=-1$ for the singlet. That is the vacuum line. The leading trajectory starts at $(D,N)=(1,0)$, but GMSZ's gap leaves nothing at $(2,1)$. This is S8's vacuum exemption for $T_2$.
- **$T_3$ is not a monotone direction.** 389 of 1 070 non-zero $T_3$-differences are negative, matching S8's failed level-3 direction.
- So the SO(9) spectrum is $T_1$- and $T_2$-monotone (with the vacuum exemption) but not $T_3$-monotone. S8 found the same for the SO(10) heads.

**Reading for half (a).**
- Regge monotonicity follows from $T_2$-monotonicity of the starts. Along each $T_2$ line the starts are partial sums of $O$, beginning at a non-negative edge value. The vacuum line is fine directly: $s=1,0,0,\dots$
- On the upper half of every ramp ($2N-D\ge-1$) the starts equal $\hat Q$, a positive Fock space. If the exactness range is proved, this gives Regge monotonicity there at every level. The proof would be a bound on the $y$-degree; so far it is checked to $D=22$.
- Below that, the formulas are differences of positive Fock spaces. Since $T_3$-monotonicity fails, positivity there cannot come from iterating the residues. It needs a reason for $O\ge0$.

### S18. Theorem: Regge monotonicity for every weight with spin at least a third of the level (2026-10-13) (`q_degree_bound.py`, `q_spin_check.py` (3))

**Theorem.** For every level $n\ge1$ and every SO(9) weight $\nu=(N,\mu)$ ($N=\nu_1$, $\mu$ the SO(7) part) with $N>\mu_1$ and $3N\ge n-1$:

$$r^9_n(\nu)-r^9_{n-1}(\nu-e_1)=\mathrm{mult}_\mu\,[q^{\,n-2N}]\,\hat Q\ \ge\ 0,$$

$$\hat Q=q^2\prod_{w\in4}(1+x^w)\,\frac{1}{\det_7(1-qg)}\prod_{j\ge3}\frac{1}{(1-q^j)\det_7(1-q^jg)}\prod_{j\ge2}\det_8\big(1+q^jg\big).$$

Here $\prod_{w\in4}(1+x^w)=1+7+8$, so $\hat Q$ is a product of characters of genuine SO(7) representations.
- So Regge monotonicity, the input to half (a) (S6), holds at every level for all weights whose spin is at least $(n-1)/3$.
- Moreover the number of new trajectories there depends only on $n-2N$ and the shape.

**Proof.**
1. *Setup (S15, S16).* Let $Q(q,y)=(1-y)\sum_{n\ge1}q^nK(r^9_n)|_{z=y/q}$, with $K=\det_7(1-z^{-1}g)\,\cdot$ and $c(D,e;\mu)=\mathrm{mult}_\mu[q^Dy^e]Q$. By the parabolic BGG resolution, $r^9_{N+D}(N,\mu)=\sum_{e\le N}c(D,e;\mu)$ at every dominant $(N,\mu)$. So for $N>\mu_1$ the start equals $c(n-N,N;\mu)$.
2. *Factorization.* Since $F_n=256\otimes r^9_n$ and $K$ is multiplicative:

   $$Q=q^2y^{-2}\,(8_v+8_c)(y/q)\,\det_7\big(1-\tfrac qy g\big)\,\frac{P-1+y}{1+R},\qquad (y/q)^2(1+R)=256(y/q).$$

   Here $P$ is the light-cone oscillator Fock space without $\alpha^+_{-1}$, and $-1+y$ removes level 0. Write $P=P'/(1-qy)$, where $1/(1-qy)$ is $\alpha^+_{-2}$.
3. *$T_2$-differences.* $O:=(1-qy)Q$ has coefficients $O(D,e)=c(D,e)-c(D-1,e-1)$, and

   $$O=q^2y^{-2}\,(8_v+8_c)(y/q)\,\det_7\big(1-\tfrac qy g\big)\,\frac{P'-1+qy+y-qy^2}{1+R}.$$

4. *Excess bound.* Give $q^Dy^e$ (any weight) the excess $\varepsilon=2e-D$, which is additive. The largest excess in each factor:
   - $q^2y^{-2}$: $-6$.
   - The ground states $y/q,\ q/y,\ 6,\ (y/q)^{1/2}\bar 4,\ (q/y)^{1/2}4$ have $\varepsilon=3,-3,0,\tfrac32,-\tfrac32$: maximum 3, attained only by $y/q$.
   - $\det_7(1-\tfrac qyg)$: each power of $q/y$ has $\varepsilon=-3$, so the maximum is 0.
   - $1/(1+R)=\sum_j(-R)^j$: $R$ consists of $(q/y)^{2-x_1}x^\mu$ with $x_1\le\tfrac32$, so $\varepsilon\le-\tfrac32$ there and the maximum is 0.
   - The bracket. In $P'$ the oscillator monomials have $\varepsilon=3-k\le0$ ($\alpha^+_{-k}$, $k\ge3$), $-k-3$ ($\alpha^-_{-k}$), $-k$ (transverse), $\tfrac32-k$ (charge-$+\tfrac12$ fermions) and $-k-\tfrac32$ (charge-$-\tfrac12$ fermions). The only positive ones are the four level-1 charge-$+\tfrac12$ fermions, each $\tfrac12$, so $\varepsilon(P')\le2$. The other terms $-1,qy,y,-qy^2$ have $\varepsilon=0,1,2,3$. So the maximum is 3, attained only by $-qy^2$, and everything else has $\varepsilon\le2$.

   The total is $\le0$. Equality holds only for $q^2y^{-2}\cdot(y/q)\cdot(-qy^2)=-q^2y$ at weight 0, the vacuum line's $T_2$-difference (S17). Changing any factor to a non-maximal term lowers $\varepsilon$ by at least 1, so every other monomial has $\varepsilon\le-1$.
5. *Conclusion.* Fix $g=D-e$. Then $c(g+e,e)=\sum_{e'\le e}O(g+e',e')$, and $\hat Q_g:=\sum_{e'}O(g+e',e')$ is a finite sum and the stable value along $T_2$. So

   $$c(D,e)=\hat Q_{D-e}-\sum_{e'>e}O(g+e',e').$$

   If $2e-D\ge-1$, i.e. $e\ge g-1$, every $e'>e$ has $\varepsilon=e'-g\ge0$. So the tail vanishes, except the vacuum term, which enters only at $(D,e)=(1,0)$. With $D=n-N$, the condition $2N-D\ge-1$ is $3N\ge n-1$, and $D-N=n-2N$.
6. *$\hat Q$ is a positive Fock space (S17).* Set $qy=1$, so $z=q^{-2}$. The identity of S15,

   $$(8_v+8_c)(z)\prod_{w\in4}(1+z^{-1/2}x^w)=z^{-1}\,256(z),$$

   holds as Laurent polynomials in $z^{1/2}$ (checked on weights). Applied to the level-2 charge-$\tfrac12$ fermions (weight $q=z^{-1/2}$), it cancels the 256. Then $\det_7(1-q^2g)$ removes the grade-2 bosons, and the level-1 fermions give $\prod_{w\in4}(1+x^w)$. What remains is the product above. $\square$

The proof rests on standard facts: the light-cone character, $F_n=256\otimes r^9_n$, and the parabolic BGG resolution. Beyond those, only explicit algebra is used.

**Checks.**
- The excess lemma holds weight by weight, exactly (mod two primes), to $D=18$ (`q_degree_bound.py`). Apart from the vacuum term, which is exactly $-1$ at weight 0, no slice of $O$ with $\varepsilon\ge0$ is non-zero. The largest excess attained is $-1$ (tensor) and $-\tfrac32$ (spinor), so the bound is sharp.
- The conclusion holds on the raw data, independently of $Q$ (`q_spin_check.py` (3)): for $n\le40$, all 16 683 starts with $3N\ge n-1$ (2 974 of them non-zero) equal the multiplicities of the Fock form of $\hat Q$.
- It also agrees with S16 (2): every start of $Q$ with $2e-D\ge-1$ equals $\hat Q$ ($D\le22$).

**Status of half (a).** Regge monotonicity, $r^9_n(\nu+e_1)\ge r^9_{n-1}(\nu)$, now holds:
- at every level for spin $\nu_1\ge(n-1)/3$ (this theorem);
- at every level within 22 of the top (S16);
- everywhere for $n\le40$ (data).

What remains is low spin, $\nu_1<(n-1)/3$, at $n>40$. In the next band, $(n-2)/4\le\nu_1<(n-1)/3$, the starts are $\hat Q_{n-2N}-\sum_{j\ge1}\hat O_{n-3N-j}$ (S17, exact for $D\le22$). There positivity is an inequality between two explicit positive Fock spaces.

## Open / next ideas

- **The next band (from S18).** For $(n-2)/4\le\nu_1<(n-1)/3$, prove the exactness of the two-pole formula by a second excess bound (on $(1-q^2y)O$, with excess $3e-D$; edge points of large shapes need care). Then prove the inequality $\hat Q_{g}(\mu)\ge\sum_{h<g-N}\hat O_h(\mu)$ between the two explicit Fock spaces, by an injection or by bounds.
- **Why are the starts $T_2$-monotone? (from S17)** $O\ge0$ (vacuum aside) would give half (a) at all orders. $O$ tends to $\hat O\ge0$ along $T_3$, but not monotonically. So a different structure is needed, possibly an algebra on which $\alpha^+_{-1}$ and $\alpha^+_{-2}$ (S8's $T$ and $T_2$) act injectively.
- **The band beyond the data.** The all-full lift at level $n$ needs $Q$ to depth about $n$, while the direct computation reaches $D\approx22$–24 on this machine. Going beyond the data needs the multi-Regge structure above, not more brute force.
- **Minimal repair levels for the Regge problem.** An ILP per $D$ on $\Phi_D$ alone is $n$-independent; compare it with $K_L(n)$ of S13.

- **Prove the Regge continuation conjecture (S5).** This is the current main line; the plan is in the root STATUS.
- **Existence in move coordinates.** The uplift set is $\{t\in\mathbb Z^K: R_{\rm can}+\sum t_\nu h_\nu\ge0\}$. Since the moves are multiples of $e$, an all-orders argument might construct one explicit uplift and control the $h_\nu$ algebraically.
- **Transport formulation (from S2).** Split each $m_n(\mu)$ among the fibres whose box contains $\mu$, so that every fibre's finite-difference inverse is non-negative.
- **Generating-function routes to an explicit uplift.** $\sum_n R^9_n q^n=(1-q)^{-1}\sum_n r^9_n q^n$; an SO(10)-covariant formula for $(1-q)^{-1}\times$(light-cone character). In the Regge limit this is realised (S15): the SO(8) Fock space $Z_8$ gives the uplift there, but not at finite $n$ (see the S15 correction).
- Write the explicit per-fibre finite-difference inverse, and check what non-negativity means in those terms.
