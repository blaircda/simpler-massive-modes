# uplift_structure — exploring the structure of the SO(10) uplifts

This directory holds exploratory analysis of the computational output in [`../so10_uplift/`](../so10_uplift/STATUS.md): invariants, fibre structure, and the kernel of the uplift problem. These directions may or may not lead to the all-orders proof. The hard computational evidence (existence through level 24, counts, canonical choice) stays in `so10_uplift/`.

All scripts here only *read* `so10_uplift/data/*.json` and reuse `so10_uplift/uplift.py`, `weights.py` (via `common.py`). They never import baseline files and never write to `so10_uplift/`.

| script | topic | output | runtime |
|---|---|---|---|
| `count_invariant.py` | S1 | `results/count_invariant.txt` | ~1 min |
| `fibre_analysis.py` | S2 | `results/fibre_analysis.txt` | ~1 min |
| `canonical_in_fibres.py` | S3a | `results/canonical_in_fibres.txt` | ~10 min |
| `kernel_moves.py` | S3b | `results/kernel_moves.txt` | ~15 min |

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

## Open / next ideas

- **Existence in move coordinates.** The uplift set is $\{t\in\mathbb Z^K: R_{\rm can}+\sum t_\nu h_\nu\ge0\}$. Since the moves are multiples of $e$, an all-orders argument might construct one explicit uplift and control the $h_\nu$ algebraically.
- **Transport formulation (from S2).** Split each $m_n(\mu)$ among the fibres whose box contains $\mu$, so that every fibre's finite-difference inverse is non-negative.
- **Generating-function routes to an explicit uplift.** $\sum_n R^9_n q^n=(1-q)^{-1}\sum_n r^9_n q^n$; an SO(10)-covariant formula for $(1-q)^{-1}\times$(light-cone character).
- Write the explicit per-fibre finite-difference inverse, and check what non-negativity means in those terms.
