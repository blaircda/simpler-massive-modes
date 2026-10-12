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
- $E_n\supseteq T_2(E_{n-2})$ holds for all $4\le n\le8$ under both rules (the rule-independent range).
- The vacuum is exempt: $T_2(E_1)=10$ is not in $E_3=16$.
- The plain rules do not enforce this containment, and lose it from $n=11$ (max-height) and $n=9$ (min-height).
- **Double Regge induction runs to $n=32$ for all three tested rules (max-height, min-height, feasible).** Here $E_n=T_2(E_{n-2})+P_n$ with $P_n\ge0$, vacuum exempt. So the uplift can be chosen with two commuting Regge directions, $T$ (level $+1$) and $T_2$ (level $+2$).
- The primary counts are 0, 1, 1, 2, 4, 7, 13, 26, 48, 89, 163, 301, 538, 972, 1734 for $n=2..16$. They are *automatically* rule-independent: $|P_n|=|E_n|-|E_{n-2}|$, and $|E_n|=|A_n|-|A_{n-1}|$ is fixed by the uplift-invariant irrep count (S1). They are about 70% of $|E_n|$.

**(iv) No level-3 direction.** Forcing $T_3$ (shift $+e_1$, level $+3$) makes the induction **fail at $n=7$ for every rule**. This holds for $T_3$ on its own, together with $T_2$, and for all $T_k$ with $k\ge2$.
- On its own the forced copy is $T_3(E_4)=320$, and no non-negative lift of $\rho_7$ contains it.
- $A_6$ is the unique uplift, so this is a strict obstruction.
- So the structure is not a Cartan module over all string modes (an arc space).

**Tentative reading (unproved).**
- Two commuting first-row shifts at levels 1 and 2 are what the first-row highest-weight monomials $(y_0^+)^a(y_1^+)^b$ of *two* vector oscillators $y_0$ (degree 1) and $y_1$ (degree 2) would produce, in the spirit of $\mathfrak{so}(10)\times\mathfrak{sp}(4)$ Howe duality.
- The $45$ at level 4 also fits $w=y_0\wedge y_1$ (degree 3), an angular-momentum-like $x\wedge p$.
- But a bare $y_1$ (a $10$ at level 3) is absent, so any such model must suppress $y_1$ on the vacuum.

## Open / next ideas

- **Prove the Regge continuation conjecture (S5).** This is the current main line; the plan is in the root STATUS.
- **Existence in move coordinates.** The uplift set is $\{t\in\mathbb Z^K: R_{\rm can}+\sum t_\nu h_\nu\ge0\}$. Since the moves are multiples of $e$, an all-orders argument might construct one explicit uplift and control the $h_\nu$ algebraically.
- **Transport formulation (from S2).** Split each $m_n(\mu)$ among the fibres whose box contains $\mu$, so that every fibre's finite-difference inverse is non-negative.
- **Generating-function routes to an explicit uplift.** $\sum_n R^9_n q^n=(1-q)^{-1}\sum_n r^9_n q^n$; an SO(10)-covariant formula for $(1-q)^{-1}\times$(light-cone character).
- Write the explicit per-fibre finite-difference inverse, and check what non-negativity means in those terms.
