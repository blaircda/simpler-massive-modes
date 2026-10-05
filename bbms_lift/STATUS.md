# bbms_lift — the BBMS SO(10)×SO(2)$_\Delta$-graded lift as an existence vehicle

Research line started 2026-10-05, following `../previous_efforts.md` item 2 and `../LITERATURE_SEARCH.md` A1/A2. The BBMS reference is Beisert, Bianchi, Morales, Samtleben, [hep-th/0310292](https://arxiv.org/abs/hep-th/0310292), §2.

**Aim.**
- Use the BBMS grading to define *one* uplift uniformly in $n$, ideally by a closed formula, and then try to prove its cumulative sums are non-negative.
- This is about existence, not canonicity. Whether the BBMS representative is physically preferred for IIA/M-theory is a separate question (see `../previous_efforts.md`). By `../uplift_structure/` S4, chirality is irrelevant to existence anyway.

Derivation scripts here never import baseline files. Comparison with `so10_uplift/` output happens only in `compare_uplifts.py`.

| script | what | output | runtime |
|---|---|---|---|
| `graded.py` | library: ν-graded Fock space, $r_n(x)$, the map $\varphi$, linear system | — | — |
| `derive_graded.py` | solve the graded lift level by level; cumulative checks | `data/*.json` | ~35 s to $n=9$ |
| `compare_uplifts.py` | verification: compare with `so10_uplift/data/R10*.json` | stdout | ~2 s |
| `letter_identities.py` | exact identities for $\varphi$ | `results/letter_identities.txt` | ~3 s |
| `image_criterion.py` | B5: generator identities, and the $x=\pm1$ conditions for $r_n(x)$ | `results/image_criterion.txt` | ~10 s |

Data:
- `data/bbms_graded.json`: `--max-level 7 --slack 1` (ν-window domain). It shows the kernel at $n=7$.
- `data/bbms_unitarity_g8.json`: `--max-level 12 --gmax 8` (unitarity-window domain). This is the main dataset.

## Set-up

- **Graded chiral spectrum.** We use the light-cone chiral Fock space, refined by oscillator number ν: every $\alpha^I_{-m}$ and $S^a_{-m}$ counts 1. Then $r_n(x)=F_n(x)/(8_v+8_s)$, sliced in ν with no shift. This reproduces BBMS's $\mathrm{vac}_2 = 1\cdot x+8_v x^2$.
- **The grading map.**
  $$\varphi:\ R(\mathrm{Spin}_{10})[x^{\pm1/2}]\to R(\mathrm{Spin}_8)[x^{\pm1/2}],\qquad V_\lambda x^\Delta\mapsto x^\Delta\,\mathrm{Res}^{SO(10)}_{SO(8)\times SO(2)_J}(V_\lambda)\big|_{y=1/x},$$
  i.e. $\nu=\Delta-J$.
  - $J$ rotates the 5th ambient coordinate, and SO(8) acts on $x_1..x_4$, which is the SO(8) inside the interlacing SO(9) of `so10_uplift/`.
  - $\varphi$ is a ring homomorphism and commutes with Adams operations ($\psi^k x=x^k$).
- **The graded problem.** Find $X_n$ with $\varphi(X_n)=r_n(x)$. Forgetting Δ ($x\to1$), $X_n$ restricts to $r^9_n$ on SO(9), and $R_n=\sum_{i\le n}X_i$ is a candidate uplift.
- **Counting convention.** Positive counting: fermions carry $+$. This is BBMS's supercharacter with spinor-class coefficients negated.

## Findings (2026-10-05)

### B1. Reproduction of BBMS
Levels 1–3 reproduce hep-th/0310292 eq. (2.14) exactly, including the chirality:
$$X_2=[1,0,0,0,0]_2-[0,0,0,0,0]_3,\qquad X_3=[2,0,0,0,0]_3-[1,0,0,0,0]_4+[0,0,0,0,1]_{5/2}.$$

### B2. Kernel of the grading map is $(K)$, $K=\lambda_{-x}(10)=\sum_k(-1)^k\Lambda^k(10)\,x^k$ (proved)
- **Proof.** $\varphi$ evaluates a W(D5)-invariant Laurent polynomial on the hypersurface $xy=1$.
  - Vanishing there implies vanishing on every Weyl image $x\,y_i^{\pm1}=1$.
  - So the function is divisible, in the UFD $\mathbb Z[y_i^{\pm1/2},x^{\pm1/2}]$, by $\prod_{i,\pm}(1-x y_i^{\pm1})=K$.
  - The quotient is again W-invariant.
- **Numerical check.** `letter_identities.py` (a) checks $\varphi(K)=0$. At $n=7$, the kernel vectors found numerically in wide Δ-windows are $x^{-1}K$ and multiples (`data/bbms_graded.json`, nullity 1 at slack 1, 4 at slack 2).
- **Consequences.**
  - $K(1)=e=-\delta^2$, which is $\sigma$-invariant. So the grading fixes chirality completely: there is no δ-odd freedom.
  - It leaves exactly the Euler-class ($e$-move) freedom of `../uplift_structure/` S3b/S4.
  - "Unique in a finite domain" (as in `previous_efforts.md`) is a statement about the Δ-window. $K$ spans 10 units of Δ.

### B3. The unitarity window gives a unique, non-negative lift (empirical, $n\le12$)
- **The domain.** BBMS's SO(10) bound (their eq. 2.15) is
  $$b(\lambda)=1+k+2l+3m+\tfrac52p+\tfrac32q\qquad(\lambda=[k,l,m,p,q]).$$
  Restrict to $0\le\Delta-b(\lambda)\le G$, with λ among the class candidates of $R^9_n,R^9_{n-1}$ in both chiralities. The bound is chirality-asymmetric, consistent with B1.
- **Results.**
  - With $G=8$ the graded system has nullity 0 and is solvable at every $n\le12$. With $G=4$ the same holds for $n\le9$ (the window is too narrow from $n=10$). Where both apply the solutions coincide, since the $G=4$ domain is inside the $G=8$ one.
  - The cumulative sums $R_n$ are **non-negative** and restrict exactly to $R^9_n$ for all $n\le12$. Level 12 took ~5 min.
  - All gaps $\Delta-b(\lambda)$ are even. The positive terms are mostly on the bound.
- **Risk.** The maximal gap grows slowly with $n$: 0, 0, 2, 2, 2, 4, 4, 4, 4, 6, 6, 6 for $n=1..12$. $K$ spans 10 units of Δ. If gaps keep growing, the window must widen, and it may then admit multiples of $K$ (non-uniqueness) or need a sharper normalisation. Uniqueness is only established for windows narrower than that.
- **Comparison with `so10_uplift/`** (`compare_uplifts.py`):
  - At class level, $R_n$ equals the canonical (= baseline) uplift for $n\le8$.
  - For $n=9..12$ it is a *different* non-negative uplift. At $n=9,10$ it is among the enumerated ones. At $n=9$ the difference is a single Euler-type move that lowers top-heaviness: $+[0]-[3,0,0,0,0]+[2,1,0,0,0]-[2,0,1,0,0]+[2,0,0,1,1]-[2,0,0,0,2]$.
  - So BBMS is not the max-height uplift: it is less top-heavy.
  - The level-6 and level-7 differences from the "protected continuation" in `previous_efforts.md` are pure chirality, invisible on classes; for example $\delta(10+2)$ at $n=6$.

### B4. The BBMS letter has a finite SO(10)-covariant preimage: the pure-spinor numerator (proved identity)
$$\varphi\big(1-10x+16x^{3/2}-\overline{16}x^{5/2}+10x^3-x^4\big)=x(x^2-1)(8_v-8_s),\qquad \varphi(10x-1-x^2)=8_v\,x.$$
Here $16=[0,0,0,0,1]$, the same chirality that BBMS assign at level 3. $N(t^2)$ is the numerator $1-10t^2+16t^3-16t^5+10t^6-t^8$ of the D=10 pure-spinor partition function, i.e. the character of the linearised 10D SYM BV multiplet. Since $\varphi$ commutes with Adams operations, one can write a formal closed form
$$X(q,x)\;\equiv\;\frac{x(x^2-1)}{N(x)}\Big(\prod_{m\ge1}\mathrm{PE}\Big[q^m\frac{N(x)}{x^2-1}\Big]-1\Big)\pmod K$$
(supercharacter conventions) in the fraction field of $R(\mathrm{Spin}_{10})[x^{\pm1/2}]/(K)$.

**Assessment (2026-10-05): this formula is automatically true but not computational.**
- **Automatically true.** By B2, $\varphi$ is injective on the quotient, which is a domain. Both sides have the same $\varphi$-image by construction. So there is nothing to "check".
- **Not computational.** $N$ and $x^2-1$ are not units, and the individual letters have no polynomial preimage. For example $x\,8_s$ fails the $x=1$ condition of B5, because $8_s$ is not an SO(9) character.
- **What it would take.** A closed formula producing $X_n$ directly would have to use building blocks that each satisfy B5. That is the analogue of Hanany–Forcella–Troost using NSR letters $[1,0,0,0]_9-1$ instead of GS ones; see Open/next.

### B5. Exactly which graded characters lift (theorem)

**Theorem.** $f\in R(\mathrm{Spin}_8)[x^{\pm1}]$ lies in $\varphi\big(R(\mathrm{Spin}_{10})[x^{\pm1/2}]\big)$ if and only if both of these hold:
- ($x=+1$) $f(1)$ is $W(B_4)$-invariant, i.e. an SO(9) character;
- ($x=-1$) $f(-1)$ has $W(B_4)$-invariant tensor part and $W(B_4)$-*anti*-invariant spinor part.

Here $W(B_4)/W(D_4)$ is generated by $\sigma_4: w_4\mapsto w_4^{-1}$, which swaps $8_s\leftrightarrow8_c$.

*Necessity.* Write $f(w,s)=F(w,z_5^{1/2}=1/s)$ with $s=x^{1/2}$, and use the D5 Weyl element $(z_4,z_5)\mapsto(z_4^{-1},z_5^{-1})$.
- At $s=\pm1$ it maps the hypersurface to itself, giving $f(w)=f(\sigma_4w)$.
- At $s^2=-1$ it lands on $z_5^{1/2}=-1/s$. Multiplying $z_4^{1/2},z_5^{1/2}$ both by $-1$ is the identity element of Spin(10), but on D4 weights it acts as $(-1)^{2\mu_4}$. That gives the sign on the spinor part.

*Sufficiency (constructive).*
- **Decomposition.** $R(\mathrm{Spin}_8)=R^{\sigma_4}\oplus(8_s-8_c)R^{\sigma_4}$. Here $R^{\sigma_4}=\mathrm{Res}\,R(\mathrm{Spin}_9)=R(SO_9)\oplus 16\cdot R(SO_9)$; the argument is the same as in `../uplift_structure/` S4.
- **Divisibility.** The two conditions give the following divisibility by $1\pm x$, class by class:

  | class | parity | condition | divisible by |
  |---|---|---|---|
  | tensor | invariant | none | — |
  | spinor | invariant | vanishes at $x=-1$ | $1+x$ |
  | spinor | anti-invariant | vanishes at $x=1$ | $1-x$ |
  | tensor | anti-invariant | vanishes at $x=\pm1$ | $1-x^2$ |

- **Generators.** All the needed pieces lie in the image, which is a ring and closed under λ-operations, since $\varphi$ is a λ-ring map:
  - $\varphi(11-x^{-1}-x)=\mathrm{Res}\,9$, and hence all of $R(SO_9)$;
  - $\varphi(16x^{1/2}+\overline{16}x^{1/2})=(1+x)(8_s+8_c)$;
  - $\varphi(16x^{1/2}-\overline{16}x^{1/2})=(1-x)(8_s-8_c)$;
  - their product gives $(1-x^2)(8_s^2-8_c^2)$.

  `image_criterion.py` checks these identities.

**Application to the string levels.**
- **$x=+1$.** This is the usual SO(9) covariance of the light-cone spectrum.
- **$x=-1$.** This says $(-1)^\nu$ acts on each mass level as the reflection of the longitudinal (9th) axis, i.e. $e_9\in\mathrm{Pin}(9)$; its twisted character has exactly this symmetry.
  - Physical argument: in light-cone gauge the rest-frame generators $S^{i9}$ come from the oscillator part of $M^{i-}$, which is cubic ($\alpha\alpha\alpha$, $\alpha SS$). So they flip ν-parity, while $S^{ij}$ preserves it.
  - So $(-1)^\nu$ commutes with so(8) and anticommutes with $S^{i9}$. By Schur it is $\pm e_9$ on each SO(9) isotypic block.
- **Checked.** `image_criterion.py` verifies both conditions for $r_n(x)$, $n\le12$.

**Consequence.**
- **A graded BBMS lift exists at every level.** Modulo the light-cone Lorentz-algebra facts above, $X_n$ exists for all $n$, unique modulo $K$. So the BBMS construction is well defined at all orders.
- **What remains open:**
  - choosing the representative (the unitarity window, B3);
  - positivity of the cumulative sums, which is the actual existence content.

## Open / next

1. Push the unitarity-window solve beyond $n=12$, watching the maximal gap against the width-10 span of $K$. The bottleneck is the exact sparse solve over $\mathbb Q$; a triangular solve by highest weight should be much faster.
2. Prove that the unitarity window contains no non-zero element of $(K)$, which would give uniqueness at all $n$. Also find out why solutions exist there, i.e. why the gaps stay small.
3. ~~Check the closed formula of B4~~ It is automatically true (see the B4 assessment). Replacement: build a *computational* closed formula from letters satisfying B5. For example, use an NSR-type grading in which every letter is $8_v=\varphi(10x-1-x^2)$, following Hanany–Forcella–Troost. This would be a different refinement from BBMS's, but any graded lift serves for existence.
4. **Existence.** Positivity of $R_n=\sum_{i\le n}X_i$ is the real content. Look for a manifestly positive rewriting of $(1-q)^{-1}X(q,x)$ at $x=1$ after reduction. Literature to check: pure-spinor partition functions and characters (Berkovits; Berkovits–Nekrasov; Movshev–Schwarz; Eager–Noja–Senghaas–Walcher [arXiv:2404.07167](https://arxiv.org/abs/2404.07167), which also cites Pengpan–Ramond).
