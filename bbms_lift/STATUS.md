# bbms_lift — the BBMS SO(10)×SO(2)$_\Delta$-graded lift as an existence vehicle

Research line started 2026-10-05, following `../previous_efforts.md` item 2 and `../LITERATURE_SEARCH.md` A1/A2. The BBMS reference is Beisert, Bianchi, Morales, Samtleben, [hep-th/0310292](https://arxiv.org/abs/hep-th/0310292), §2.

**Aim.**
- Use the BBMS grading to define *one* uplift uniformly in $n$, ideally by a closed formula, and then try to prove its cumulative sums are non-negative.
- This is about existence, not canonicity. Whether the BBMS representative is physically preferred for IIA/M-theory is a separate question (see `../previous_efforts.md`). By `../uplift_structure/` S4, chirality is irrelevant to existence anyway.

Derivation scripts here never import baseline files. Comparison with `so10_uplift/` output happens only in `compare_uplifts.py`.

| script | what | output | runtime |
|---|---|---|---|
| `graded.py` | library: ν-graded Fock space, $r_n(x)$, the map $\varphi$, linear system | — | — |
| `derive_graded.py` | solve the graded lift level by level (block solver, B7); cumulative checks; `--resume` | `data/*.json` | ~5 min to $n=24$ |
| `positive_graded.py` | B6 diagnostic: minimal $K$-correction to a resolved-positive lift (ILP) | `results/positive_graded.txt` | ~5 min to $n=24$ |
| `level_stats.py` | per-level gaps, $N_n$ check, chirality-resolved and class-level positivity | `results/level_stats.txt` | ~2 s |
| `compare_uplifts.py` | verification: compare with `so10_uplift/data/R10*.json` | stdout | ~2 s |
| `letter_identities.py` | exact identities for $\varphi$ | `results/letter_identities.txt` | ~3 s |
| `window_width.py` | nullity versus window width $G$ | `results/window_width.txt` | ~10 min |
| `fundamental_domain.py` | window as normal form for $(K)$: single-term tests | `results/fundamental_domain.txt` | ~2 min ($L=4,5$) |
| `image_criterion.py` | B5: generator identities, and the $x=\pm1$ conditions for $r_n(x)$ | `results/image_criterion.txt` | ~10 s |

Data:
- `data/bbms_graded.json`: `--max-level 7 --slack 1` (ν-window domain). It shows the kernel at $n=7$.
- `data/bbms_unitarity_g8.json`: `--max-level 24 --gmax 8` (unitarity-window domain). This is the main dataset. Levels 1–18 come from the earlier solver, 19–24 from the block solver.

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

### B3. The unitarity window: unique lift for $n\le24$, non-negative only for $n\le18$ (see B6)
- **The domain.** BBMS's SO(10) bound (their eq. 2.15) is
  $$b(\lambda)=1+k+2l+3m+\tfrac52p+\tfrac32q\qquad(\lambda=[k,l,m,p,q]).$$
  Restrict to $0\le\Delta-b(\lambda)\le G$, with λ among the class candidates of $R^9_n,R^9_{n-1}$ in both chiralities. The bound is chirality-asymmetric, consistent with B1.
- **In ambient coordinates.** $b(\lambda)=1+l_1+l_2+l_3+l_4-l_5=1+\langle\lambda,\xi\rangle$ with $\xi=(1,1,1,1,-1)=2\omega_{\overline{16}}$.
  - So the bound is 1 + the maximal $U(1)$ charge for $SO(10)\supset U(5)$.
  - $g'=\Delta-\langle\lambda,\xi\rangle$ (= gap + 1) is superadditive under tensor products.
- **Results ($G=8$).**
  - The graded system has nullity 0 and is solvable at every $n\le24$.
  - The cumulative sums $R_n$ restrict exactly to $R^9_n$ for all $n\le24$.
  - They are **non-negative only for $n\le18$**; see B6 for the obstruction.
  - Every class has multiplicity $\ge1$.
  - The total irrep count equals the uplift invariant $N_n$ of `../uplift_structure/` S1 (a consistency check).
  - All gaps are even.
- **Window width (from `gw` test, $n=4,6,10,13,14$).**
  - $G=8$ (equivalently 9) is exactly the widest window with nullity 0.
  - At $G=10$ multiples of $K$ enter. For example the nullity at $n=14$ is 38, 114, 266, 418 for $G=10,12,16,20$.
  - This matches $K$'s profile: its terms have $g'\in\{0,2,4,6,8,10\}$, with the top $g'=10$ part being the single term $1\cdot x^{10}$.
- **Maximal gap used,** $n=1..24$: 0, 2, 2, 2, 2, 4, 4, 4, 4, 6, 6, 6, 8, then 8 at every level 14–24.
  - It reaches the largest kernel-free value at $n=13$ and stays there through 24.
  - The number of terms at the edge (gap 8) keeps growing (`results/level_stats.txt`).
- **Normal-form structure** (`fundamental_domain.py`; probe of high even gaps in the 2026-10-05 session).
  - **Parity grading.** Roots pair with $\xi$ in $\{0,\pm2\}$, so $g'$ is additive mod 2 under tensor products. Every term of $K$ has even $g'$ (0, 2, …, 10). So $(K)$ respects the $g'$-parity grading, and the problem splits into an even-gap sector and an odd-gap sector.
  - **Uniqueness.** No non-zero element of $(K)$ lies in the window: nullity 0 for all window terms with $\lambda_1\le5$, both sectors.
  - **Even-gap sector** (odd $g'$, where every string-level term found so far lies). The window has $g'\in\{1,3,5,7,9\}$: five values, matching $K$'s even-$g'$ span.
    - Single terms $V_\mu x^\Delta$ with gaps 10, 12, 14 (and 16 for small $\mu$) all have window representatives.
    - The failures start at gap 18–22, and the threshold drops as $\mu$ grows. This is consistent with truncation at $\lambda_1\le5$: from gap 20, about six nested reductions are needed, and each can raise $\lambda_1$ by one.
  - **Odd-gap sector.** The window has only four $g'$ values, so gaps 9 and 11 have no representative, as the parity count predicts.
  - **Below the window.** No term with gap $<0$ has a representative.
  - **Conjecture (refined).** In the even-gap sector, every element supported at gaps $\ge0$ is congruent mod $K$ to a unique element supported in $\{0\le\Delta-b(\lambda)\le8\}$.
  - **What would then remain.** Existence and uniqueness of the BBMS normal form at all $n$ would further need the class of $r_n$ to have a preimage supported at even gaps $\ge0$.
    - The natural letter preimages sit at gap $\ge-1$: $10x-1-x^2$ has gaps $\{-1,-1,1\}$, and the numerator $N$ has gaps $\{-1,-1,-1,-1,1,3\}$.
    - So this is a genuine constraint, presumably the PSU(2,2|4) unitarity bound that BBMS's window encodes. It is empirically true through $n=18$.
- **Comparison with `so10_uplift/`** (`compare_uplifts.py`):
  - At class level, $R_n$ equals the canonical (= baseline) uplift for $n\le8$.
  - For $n=9..18$ it is a *different* non-negative uplift. At $n=9,10$ it is among the enumerated ones. At $n=9$ the difference is a single Euler-type move that lowers top-heaviness: $+[0]-[3,0,0,0,0]+[2,1,0,0,0]-[2,0,1,0,0]+[2,0,0,1,1]-[2,0,0,0,2]$.
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

### B6. Obstruction: the unitarity-window lift is not positive beyond $n=18$ (2026-10-11)

Data: `data/bbms_unitarity_g8.json` to $n=24$; `results/level_stats.txt`.
- Levels 1–18 were computed with the earlier solver.
- Levels 19–24 used the block solver (B7). Its output is identical to the earlier solver for $n\le16$.

**Exactness of the statement.** At every level $n\le24$ the window system has full column rank mod $p$ (hence over $\mathbb Q$). The solution is checked exactly against $r_n$, and its restriction equals $R^9_n$. So the failure below is a property of *the* unique window lift, not a numerical artefact. The S1 invariant $N_n$ matches at every level.

**Chirality-resolved failure from $n=19$.** Only the spinor class is involved at first:

| $n$ | 16 | 17 | 18 | 19 | 20 | 21 |
|---|---|---|---|---|---|---|
| $\overline{16}=[0,0,0,1,0]$ | 19 | 17 | 4 | **−43** | **−146** | **−346** |
| $16=[0,0,0,0,1]$ | 18 | 40 | 87 | 190 | 390 | 761 |
| merged class | 37 | 57 | 91 | 147 | 244 | 415 |

- Through $n=16$ the two chiralities are nearly balanced.
- From $n=17$ the balance runs away. The merged count still grows smoothly, so at $n=19$–22 the lift is still a valid uplift at class level (minimum class multiplicity 1), as S4 allows.

**Class-level failure from $n=23$.**
- $n=23$: singlet $-88$.
- $n=24$: singlet $-512$, $210=[0,0,0,1,1]$ $-888$, plus six chirality-resolved negatives.
- These irreps are self-conjugate, so no chirality split helps. **The window normal form is not a valid uplift at $n=23,24$.**
- Existence itself is unaffected: `so10_uplift/` has non-negative uplifts for all $n\le24$.

**Interpretation (tentative).** The lower edge of the window is BBMS's unitarity bound, which is physical. The upper edge (gap $\le8$) is not: we chose it as the widest window with no multiples of $K$.
- From $n=13$ the solution presses against that edge, and the edge terms keep growing.
- From $n=17$ the chirality balance tilts. That is a symptom of the reduction modulo $K$: $K$ is σ-invariant, but $K\cdot Y$ is not when $Y$ is spinorial.
- The natural reading is that the "physical" graded lift wants gaps $>8$ from about $n=17$. Squeezing it into the width-8 window then introduces $K$-multiples that are unphysical and eventually negative.

**Consequences for this research line.**
- At class level, every uplift is the $x=1$ image of some graded lift ($\mathrm{BBMS}+K\cdot Y$, with $K(1)=e$, by S4). So the grading by itself imposes **no** constraint relevant to existence. Its only use was the specific window normal form, and that fails.
- What survives, and is independent of the window: B2 (kernel $=(K)$), B5 (exactly which graded characters lift, and that the string levels do at all $n$), B4 (the pure-spinor letter), and S4.
- To be useful for existence, a graded approach needs a *different, principled* choice of representative mod $K$ that stays positive. One option is BBMS's own procedure: lift to SO(10) first, then assign Δ per SO(8)×U(1)$_J$ component. Their unitarity bound is "satisfied but not saturated" at higher levels, with no upper bound. That would need a uniqueness principle other than a width-8 window.

**Diagnostic: positive graded lifts do exist (`positive_graded.py`, `results/positive_graded.txt`).**

*Set-up.* For window width $G\ge10$, the single-term multiples $K\cdot V_\mu x^{\Delta}$ at gap 0 lie in the window, and each becomes $e\cdot V_\mu$ at $x=1$, whatever $\Delta$ is. So the chirality-resolved cumulative sums reachable by graded lifts are $R_n^{\rm BBMS}+\sum_\mu t_\mu\,e\cdot V_\mu$ with $t_\mu\in\mathbb Z$. This is the set $R^{\rm BBMS}+\delta^2R$: the grading fixes chirality data modulo $\delta^2$. Increments at successive levels add, so each level can be tested independently.

*Method.*
- An integer program per level, with generators $\mu$ having $\mathrm{supp}(e\cdot V_\mu)$ inside the resolved candidates.
- $e\cdot V_\mu$ is computed by Klimyk's formula and validated against Sage.
- The minimal $\|t\|_1$ is certified by the LP bound (an integer at ⌈LP⌉ is optimal).

*Results.*
- **A chirality-resolved non-negative graded lift (window width 10) exists at every $n\le24$.**
- The minimal correction to the BBMS normal form is $\|t\|_1=0,0,8,25,58,217,733$ for $n=17..23$ (all certified optimal), and $\in[1934,1936]$ at $n=24$.
- The correction is concentrated in a few low-lying moves:
  - mostly $e\cdot\overline{16}$, with coefficient $+4,+21,+57,+111,+356,+843$;
  - then $e\cdot[k,0,0,1,0]$ for $k=1,2,3$ (the $\overline{16}$ series);
  - and, from $n=23$, $e\cdot1$.
- These are exactly the directions in which the window solution went negative.

*Reading.*
- The grading is **not** incompatible with positivity. The BBMS window normal form just picks the wrong representative, by $K$-multiples in the $\overline{16}$-series and singlet directions. These need gap up to 10, just beyond the kernel-free width 8.
- A principled positive normal form therefore needs a different choice of fundamental domain for $(K)$, at least in the $\overline{16}$-series directions.
- **Literature-motivated variant: tested, failed (pre-registered single test, 2026-10-11).**
  - *Motivation.* Bianchi's lectures (hep-th/0409304, after eq. 3.17) quote two SO(10) unitarity bounds, $\Delta_\pm=1+k+2l+3m+2(p+q)\pm\frac12(p-q)$. His spinor example $[0,0,0,0,1]_{5/2}$ saturates the smaller one.
  - *The variant.* Use $\min(\Delta_\pm)=1+l_1+l_2+l_3+l_4-|l_5|$ as the lower edge (`derive_graded.py --bound min`, `data/bbms_unitarity_min_g8.json`). It is chirality-symmetric, equals our $\Delta_+$ for $l_5\ge0$, and is lower by $2|l_5|$ for $\overline{16}$-type irreps.
  - *Success criterion (fixed before running).* Uniquely solvable, with non-negative cumulative sums, for all $n\le24$; stop at the first failure.
  - *Result.* For $n\le14$ the solution is identical to the $\Delta_+$ window: the extra low-Δ slots are never used. At $n=15$ the even-gap spinor block is **inconsistent**: there is no graded solution in this window.
  - *Why.* Lowering the $\overline{16}$-type windows also lowers their top edge, and the spinor block needs that top slot. This matches the diagnostic above, which says positivity needs *more* room above (gap up to 10) in the $\overline{16}$ directions. The physically motivated bound moves the window the other way.
  - Following the caveat in `../previous_efforts.md`, no further window variants are pursued. Shaping the window to the data would be fitting, not a principle.

### B7. Block-diagonal solver (2026-10-11)

- **Block structure.** The graded system splits exactly into 4 blocks. A column $(\lambda,\Delta)$ is keyed by ($l_1\bmod2$, $(\sum\lambda_i-\Delta)\bmod2$), and a row $(\nu,w)$ by ($w_1\bmod2$, $(\sum w_i-\nu)\bmod2$). The first bit is the central character; the second is the gap parity. At $n=16$, zero of the 1,012,554 entries cross between blocks.
- **Where the data sits.** The right-hand side lies only in the two even-gap blocks. All four blocks have full column rank through $n=24$.
- **Implementation.** `graded.solve_level` builds each block directly as a dense GF(p) matrix from compact per-irrep weight arrays, with no global entry list. It echelonizes $[A\,|\,b]$ once per block and verifies the lifted solution exactly against $r_n$.
- **Cost.** Levels 19–24 took 4.5 min in total, at 1.9 GB peak (85 s for $n=24$).
- **Robustness.** `derive_graded.py` saves after every level and supports `--resume`.

## Open / next

1. ~~Faster solver~~ Done (2026-10-05).
   - D5 weights now come from a dominant-weight Freudenthal recursion with cached $(J,w_4)$ slices; it is validated against Sage on 13 irreps.
   - The linear algebra is dense GF(p) with a CRT lift, followed by an exact integer check of $Ax=b$.
   - It gives identical output for $n\le12$ and is about 15× faster.
   - Superseded by the block solver of B7, which is another ~8× faster and gives identical output for $n\le16$.
2. Prove the refined normal-form conjecture of B3, for the even-gap sector, gaps $\ge0$. A plausible route is top-down division by $K$, whose top $g'$ part is the single monic term $x^{10}$, controlled by superadditivity of $g'$. Separately, show that $r_n$ has a preimage at even gaps $\ge0$ (unitarity).
3. ~~Check the closed formula of B4~~ It is automatically true (see the B4 assessment). Replacement: build a *computational* closed formula from letters satisfying B5. For example, use an NSR-type grading in which every letter is $8_v=\varphi(10x-1-x^2)$, following Hanany–Forcella–Troost. This would be a different refinement from BBMS's, but any graded lift serves for existence.
4. **Existence.** ~~Positivity of the window lift~~ This fails from $n=23$ (B6). A principled positive representative mod $K$ is still lacking: the $\min(\Delta_\pm)$ window fails at $n=15$ (B6), and further window tuning is ruled out as fitting. Remaining leads:
   - why the $\overline{16}$ direction needs a correction that grows roughly 2–3× per level;
   - returning to the class-level structure (`../uplift_structure/`), using the Euler-move freedom with $e\cdot\overline{16}$-series moves as the main lever. Look for a manifestly positive rewriting of $(1-q)^{-1}X(q,x)$ at $x=1$ after reduction. Literature to check: pure-spinor partition functions and characters (Berkovits; Berkovits–Nekrasov; Movshev–Schwarz; Eager–Noja–Senghaas–Walcher [arXiv:2404.07167](https://arxiv.org/abs/2404.07167), which also cites Pengpan–Ramond).
