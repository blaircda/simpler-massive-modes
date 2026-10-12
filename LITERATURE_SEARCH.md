# Literature search — first pass (2026-09-30)

Seed: I. Bars, *Stringy evidence for D = 11 structure in strongly coupled type IIA superstring*, [hep-th/9503228](https://arxiv.org/abs/hep-th/9503228), [INSPIRE 393910](https://inspirehep.net/literature/393910). It has 22 citations on INSPIRE.

Method: we followed the skill in `.agents/skills/hepth-arxiv-inspire-literature-exploration/`. We listed all 22 citers of the seed, filtered them by topic and author continuity, then followed citers/references of the promising ones:

- Bars–Yankielowicz and Bars 1996;
- Curtright–Thorn 1986 (reached via the 2026 SO(9) growth paper);
- Pengpan–Ramond 1998;
- Hanany–Forcella–Troost 2010;
- Beisert–Bianchi–Morales–Samtleben 2003.

We also ran a few INSPIRE full-text queries. This pass is not exhaustive. The papers marked ★ are the ones to read first.

Notation below follows `STATUS.md`:
- $r^9_n$ is the SO(9) "super-representation" at level $n$;
- $R^9_n=\sum_{i\le n} r^9_i$;
- $R^{10}_n$ is the sought SO(10) uplift.

---

## A. Most relevant leads for the all-orders goal

### ★ A1. Beisert, Bianchi, Morales, Samtleben — *On the spectrum of AdS/CFT beyond supergravity* (2003)
[hep-th/0310292](https://arxiv.org/abs/hep-th/0310292) · [INSPIRE 632137](https://inspirehep.net/literature/632137)

**Not a citer of Bars 1995**, but the closest independent rediscovery of the III.B structure. We found it through the 2026 SO(9) growth paper (A3), which cites it for exactly this point.

- **The III.B statement.** §2.1 says: "the partial sums of chiral string excitations $\sum_{\ell=1}^m \mathrm{vac}_\ell$ naturally assemble into true SO(10) representations". For example $1$, then $10$, then $54-16_s$. This is the III.B statement in different language.
- **The level-$\ell$ content as a virtual SO(10) character.** They write each chiral level as a *virtual* SO(10) character, using $[k,l,m,p,q]^* = [k,\dots]-[k-1,\dots]$. Each term carries an extra SO(2) label $\Delta$ (a conformal dimension). For example
  $$\mathrm{vac}_3=[2,0,0,0,0]_3-[1,0,0,0,0]_4+[0,0,0,0,1]_{5/2}.$$
  So $\mathrm{vac}_\ell = R^{10}_\ell - R^{10}_{\ell-1}$, refined by a grading.
- **The rule that fixes $\Delta$.** Take $\Delta-J=\nu$, where $\nu$ is the oscillator number and $J\in SO(2)\subset SO(10)$. This is a BMN-type rule, and it would give an *extra structure that might pick out a canonical uplift*. The grading should be compared against our "max-height" canonical choice.
- **An independent all-orders formula (the main lead).** They compute the free $\mathcal N=4$ SYM single-trace spectrum by Pólya counting, then strip the BPS part and a $T_{SO(10,2)}$ factor. The result matches $Z_{SO(10,2)}=\sum_\ell V^L_\ell\times V^R_\ell$ up to $\Delta=10$.
  - Free SYM counting is an explicit all-orders formula.
  - If the chiral factor can be isolated, it would give an SO(10)×SO(2)-covariant generating function for the $\mathrm{vac}_\ell$ at all levels.
  - Caveats: they only check the L×R product, they note an SO(10) spinor-chirality ambiguity in the lift, and positivity of the partial sums would still need proof.

### ★ A2. Bianchi — *Higher spins and stringy AdS₅×S⁵* (2004, lectures) — **cites Bars 1995**
[hep-th/0409304](https://arxiv.org/abs/hep-th/0409304) · [INSPIRE 660702](https://inspirehep.net/literature/660702)

- A review of A1. It explicitly credits Bars 1995 (and Bars 2T-physics) for the SO(10) lift: "the partial sums … form SO(10) multiplets".
- It has the cleanest statement of the $\Delta_0=J+\nu$ assignment, eqs. (3.13)–(3.17).
- It gives the $T_{SO(10,2)}$ generating functions (4.9)–(4.15).
- It tabulates the SO(10)×SO(2)$_\Delta$ content up to string level 5.
- It speculates that the SO(10,2) structure is linked to Bars's two-time superstring.

The companion paper is Beisert, Bianchi, Morales, Samtleben, *Higher spin symmetry and N=4 SYM* [hep-th/0405057](https://arxiv.org/abs/hep-th/0405057) · [INSPIRE 649935](https://inspirehep.net/literature/649935). It **cites Bars 1995** as the origin of the hidden so(10), and organises the spectrum into hs(2,2|4) multiplets.

Related:
- Morales–Samtleben, *Higher spin holography for SYM in d dimensions*, [hep-th/0411246](https://arxiv.org/abs/hep-th/0411246) · [INSPIRE 665382](https://inspirehep.net/literature/665382). The same hidden SO(10,2) pattern appears for maximal SYM in any $d$, which suggests a universal (d-independent) structure.
- Bianchi–Dolan–Heslop–Osborn, *N=4 superconformal characters and partition functions*, [hep-th/0609179](https://arxiv.org/abs/hep-th/0609179) · [INSPIRE 727090](https://inspirehep.net/literature/727090). Character technology (reduced Verma modules, Weyl symmetry) that would be needed to make A1 rigorous.

### ★ A3. Georgoudis, Minahan, Ström, Zoumis — *The growth of SO(9) super-representations in Type II string theory* (Sept 2026)
[arXiv:2609.04107](https://arxiv.org/abs/2609.04107) · [INSPIRE 3199641](https://inspirehep.net/literature/3199641)

A very recent large-scale study of exactly our $r^9_n$.

- **Support rule.** They give an empirical rule (Table 1 / eq. 3.6) for which SO(9) irreps $[\ell,m,r,s]$ occur in $r^9_n$ for $n>6$. The top symmetric tensor $[n-2,0,0,0]$ is the only gap.
- **Shift relations.** They give empirical relations (3.32)–(3.35), for example $a_{[\ell,m,r,s],n}=a_{[\ell-k,m,r,s],n-k}$ for large $\ell$.
  - These are "repeat" patterns that look closely tied to the cumulative $(1-q)^{-1}$ structure. They may help control $R^9_n$ in the large-$l_1$ region.
  - They could also help prove that the fibre sizes $m_n(a,a,b,b)$ behave as needed.
- **Data and methods.** A theta-function form of the refined partition function (3.5), and multiplicities computed to level 226 (all reps) and 501 (individual reps).
- **Open problems.** They list proving the support rule and the shift relations from first principles.
- They note (§5) the SO(10) grouping of A1 and ask whether multiplicities organise into SO(10,2).

### ★ A4. Hanany, Forcella, Troost — *The covariant perturbative string spectrum* (2010)
[arXiv:1007.2622](https://arxiv.org/abs/1007.2622) · [INSPIRE 861347](https://inspirehep.net/literature/861347)

- **The formula.** A manifestly SO(9)-covariant generating function for the massive spectrum, eqs. (3.16)–(3.18). It is built from plethystic exponentials of $\frac{q}{1-q}([1,0,0,0]_9-1)$, i.e. it uses $\chi_{8_v}=\chi_{9}-1$ on restriction.
- **Template for our problem.** It is the direct template for the "generating-function route" in `uplift_structure/STATUS.md`. The analogous object would be an SO(10)-covariant expression for $(1-q)^{-1}\times$(light-cone character), for instance using $\chi_{8_v}=\chi_{10}-2$.
- **Same open problem one step down.** They remark that positivity of the SO(9) coefficients is "very plausible" but "a challenge to prove". That is the SO(8)→SO(9) analogue of our SO(9)→SO(10) problem, so any proof technique for one should inform the other.
- They also observe a selection rule (third Dynkin label ≤ 1, for levels ≤ 9).

### A5. Lüst, Mekareeya, Schlotterer, Thomson — *Refined partition functions for open superstrings with 4, 8 and 16 supercharges* (2012)
[arXiv:1211.1018](https://arxiv.org/abs/1211.1018) · [INSPIRE 1198465](https://inspirehep.net/literature/1198465)

- Super-Poincaré covariant refined partition functions.
- It is the source (ref. [7] of A3) of the observation that multiplicities at level $n$ "repeat" at level $n+1$ with $\ell\to\ell+1$.
- It includes a conjecture on large-$s$ behaviour and asymptotics of leading and subleading Regge trajectories.

---

## B. Bars's own follow-ups (direct citers of the seed)

### ★ B1. Bars & Yankielowicz — *U-duality multiplets and nonperturbative superstring states* (1995)
[hep-th/9511098](https://arxiv.org/abs/hep-th/9511098) · [INSPIRE 402246](https://inspirehep.net/literature/402246)

- In toroidal compactifications, U-duality independently demands extra non-perturbative states. The internal $k$-group has to be extended, for example $(4_{\rm int})_L\to(5_{\rm int})_L$.
- For several compactifications these extra states coincide with the III.B minimal completion $9\to10$ at every level they checked.
- **Why it matters.** It gives a second, duality-based derivation of the same cumulative structure. That could supply an alternative all-orders mechanism: compare the $K$-covariance of the compactified perturbative spectrum.

### ★ B2. Bars — *Consistency between 11-D and U duality* (1996)
[hep-th/9601164](https://arxiv.org/abs/hep-th/9601164) · [INSPIRE 404610](https://inspirehep.net/literature/404610)

- A concise restatement of the III.B scheme, eq. (11): $|l\rangle_p\oplus|l-1\rangle_{np}\oplus\cdots\oplus|1\rangle_{np}$ forms complete SO(10)$_{L,R}$ multiplets. It was "verified explicitly up to $l=5$" and "conjectured to be true at all levels".
- For $(d,c)=(10,0),(9,1),(8,2),(6,4)$, U-duality and the 11D minimal completion imply each other.
- For $(7,3),(5,5),(4,6),(3,7)$ they differ, and a larger (non-minimal) set of states is needed. This bears on III.B versus III.C.

### B3. Bars — S-theory series (1996)
- *Supersymmetry, p-brane duality and hidden space-time dimensions*, [hep-th/9604139](https://arxiv.org/abs/hep-th/9604139) · [INSPIRE 417903](https://inspirehep.net/literature/417903)
- *Duality and hidden dimensions*, [hep-th/9604200](https://arxiv.org/abs/hep-th/9604200) · [INSPIRE 418147](https://inspirehep.net/literature/418147)
- *S theory*, [hep-th/9607112](https://arxiv.org/abs/hep-th/9607112) · [INSPIRE 420682](https://inspirehep.net/literature/420682)
- *Algebraic structure of S theory*, [hep-th/9608061](https://arxiv.org/abs/hep-th/9608061) · [INSPIRE 421827](https://inspirehep.net/literature/421827)

Representation theory of the 32-supercharge superalgebra with all central extensions (SO(10,2) and hidden 12th/13th dimensions). This is background for why SO(10) (and SO(10,2)) should act on the index space. It is probably not directly useful for the combinatorial proof, but it is the conceptual link to A1's SO(10,2).

### B4. Bars — 2T-physics / twistor superstring (2002–2005)
- *Hidden twelve-dimensional structures in AdS₅×S⁵ and M⁴×R⁶ supergravities*, [hep-th/0208012](https://arxiv.org/abs/hep-th/0208012) · [INSPIRE 591848](https://inspirehep.net/literature/591848)
- *Twistor superstring in 2T-physics*, [hep-th/0407239](https://arxiv.org/abs/hep-th/0407239) · [INSPIRE 655260](https://inspirehep.net/literature/655260). A 10+2 dimensional generalisation whose extra-boson excitations are proposed as the origin of the SO(10) structure. It is cited by A1/A2 for this reason.

This line may not help the combinatorics, but it is the only proposed *dynamical* origin of the SO(10)/SO(10,2) structure.

---

## C. Representation-theoretic tools (reached by traversal)

### ★ C1. Curtright & Thorn (with Goldstone, Ghandour) — SO(D−1) content of string spectra (1986)
- Curtright–Thorn, *Symmetry patterns in the mass spectra of dual string models*, Nucl. Phys. B274 (1986) · [INSPIRE 18086](https://inspirehep.net/literature/18086)
- Curtright–Thorn–Goldstone, *Spin content of the bosonic string*, Phys. Lett. B175 (1986) · [INSPIRE 18393](https://inspirehep.net/literature/18393)
- Curtright–Ghandour–Thorn, *Spin content of string models*, Phys. Lett. B182 (1986) · [INSPIRE 18991](https://inspirehep.net/literature/18991)

These are the classic generating functions for the multiplicity of each O(D−1) irrep at level $n$, with the derivation valid for general odd D−1. The superstring version uses GL(∞|∞) supercharacters.

- The CTG method turns SO(D−2) light-cone characters into SO(D−1) multiplicities, and so is the natural thing to try to push one rank up: SO(9) cumulative → SO(10).
- Curtright–Thorn also extract "simple empirical rules relating different representations" from data.

(Pre-arXiv; need library or journal access.)

### C2. Pengpan & Ramond — *M(ysterious) patterns in SO(9)* (1998)
[hep-th/9808190](https://arxiv.org/abs/hep-th/9808190) · [INSPIRE 475611](https://inspirehep.net/literature/475611)

The supergravity triplet $44+84+128$ (the factor multiplying $r^9_n$ in our physical states) extends to infinite families of "Euler triplets" via the three embeddings SO(9)⊂F₄.

- **Physical states as F₄/SO(9) data.** The physical states $r^9_n\otimes(44+84+128)$ are Kostant/Euler-triplet objects. That is a possible structural handle on the full (not "super-reduced") spectrum.
- **Possible link to III.C.** This may be more relevant to the III.C higher-multiplicity scheme.

Follow-ups:
- Brink–Ramond–Xiong, *Supersymmetry and Euler multiplets*, [hep-th/0207253](https://arxiv.org/abs/hep-th/0207253) · [INSPIRE 591548](https://inspirehep.net/literature/591548)
- Ramond, *Algebraic dreams*, [hep-th/0112261](https://arxiv.org/abs/hep-th/0112261) · [INSPIRE 569035](https://inspirehep.net/literature/569035)
- Ramond, *The sirens of eleven dimensions*, [hep-th/0411024](https://arxiv.org/abs/hep-th/0411024) · [INSPIRE 663344](https://inspirehep.net/literature/663344)

### C3. Gross, Kostant, Ramond, Sternberg — *The Weyl character formula, the half-spin representations, and equal rank subgroups* (1998)
[math/9808133](https://arxiv.org/abs/math/9808133) · [INSPIRE 478320](https://inspirehep.net/literature/478320)

- **The identity.** $V_\lambda\otimes(S^+-S^-)=\sum_{w\in W_G/W_H}\pm U_{w\bullet\lambda}$ for equal-rank $H\subset G$.
- **Link to our kernel.** Our uplift kernel is generated by Euler-class moves $e\otimes V_\nu$ with $e=\sum(-1)^k\Lambda^k(10)$, which vanishes on Spin(9). For SO(9)⊂SO(10) (unequal rank), the analogous statement is $e|_{SO(9)}=0$.
- **How it could help.** Their methods (cubic Dirac operator, Weyl-group alternating sums) could give closed forms for $e\otimes V_\nu$, and hence a handle on the move lattice. It also underlies C2.

### C4. Cremmer & Gervais — *T duality of perturbative characters for closed bosonic and type II string theories* (1997)
[hep-th/9707107](https://arxiv.org/abs/hep-th/9707107) · [INSPIRE 445577](https://inspirehep.net/literature/445577)

A citer of Bars–Yankielowicz.

- It factorises the long SUSY multiplet out of the transverse-rotation characters for all massive levels at once, which is our $F_n/(8_v+8_s)$ step.
- It also uses O(8) triality.
- It could help give a clean all-orders character-level definition of $r^8_n$ and $r^9_n$.

### C5. Recent covariant spectrum technology (lower priority)
- Markou–Skvortsov, *An excursion into the string spectrum*, [arXiv:2309.15988](https://arxiv.org/abs/2309.15988) · [INSPIRE 2703943](https://inspirehep.net/literature/2703943)
- Basile–Markou, *On the deep superstring spectrum*, [arXiv:2405.18467](https://arxiv.org/abs/2405.18467) · [INSPIRE 2791237](https://inspirehep.net/literature/2791237)

These use Howe duality: the SO(D−1) little group is paired with a symplectic or orthosymplectic algebra, whose lowest-weight states give whole Regge trajectories. That yields a *structural* decomposition of the spectrum into trajectories, a possible way to organise $r^9_n$ uniformly in $n$. The superstring paper is the relevant one.

Also: Pesando, *The bosonic string spectrum … up to level 10*, [arXiv:2405.09987](https://arxiv.org/abs/2405.09987) · [INSPIRE 2787469](https://inspirehep.net/literature/2787469). It builds covariant irreps from explicit light-cone states.

---

## D. Citers of the seed judged not relevant (for completeness)

These cite Bars 1995 only in passing, for general 11D/M-theory context:
- Witten-era duality papers: Sen [hep-th/9504027](https://inspirehep.net/literature/394043); Duff–Liu–Minasian [hep-th/9506126](https://inspirehep.net/literature/396364); Duff [hep-th/9509106](https://inspirehep.net/literature/399682).
- Cvetič–Youm [hep-th/9505045](https://inspirehep.net/literature/394856); Dasgupta–Mukhi [hep-th/9512196](https://inspirehep.net/literature/403938).
- Youm's review and thesis ([hep-th/9710046](https://inspirehep.net/literature/449382), [1484232](https://inspirehep.net/literature/1484232)).
- Bars 1996 black-hole entropy and 12th dimension ([hep-th/9610074](https://inspirehep.net/literature/424430)); Bars et al. hidden 12D super-Poincaré ([hep-th/0312147](https://inspirehep.net/literature/635576)); Bars twistors and 2T review ([hep-th/0502065](https://inspirehep.net/literature/676276)).
- Boya on F-theory ([833212](https://inspirehep.net/literature/833212)); Di Vecchia–Maskalaniec on D0-brane eikonal scattering ([arXiv:2609.07788](https://inspirehep.net/literature/3201111)).

Seed references (all pre-1995 membrane/duality background): Witten; Hull–Townsend; Townsend; and Bars, Bars–Pope–Sezgin on supermembrane spectra ([INSPIRE refs](https://inspirehep.net/search?p=citedby:recid:393910)). The supermembrane massive-level papers (Bars 1987, 1989) could be worth a look for the original motivation of III.C.

---

## Suggested next steps

1. **Compare with A1/A2 data.** Take their SO(10)×SO(2)$_\Delta$ tables (levels ≤ 5) and check:
   - that they agree with our $R^{10}_n$;
   - whether their $\Delta$ rule ($\Delta-J=\nu$) selects a unique uplift at levels where ours is not unique ($n\ge7$).

   If so, it is a candidate canonical construction with a physical definition.
2. **Try an SO(10)-covariant Hanany–Forcella–Troost formula (A4).** Write an analogue for $(1-q)^{-1}\sum_n r^9_n q^n$ in plethystic form over SO(10) fugacities with $\chi_{8_v}=\chi_{10}-2$ (or similar), and check whether it reproduces a valid uplift.
3. **Use the A3 rules.** Use the empirical support rule and shift relations to predict the fibre sizes $m_n(a,a,b,b)$ and the large-$l_1$ behaviour of $R^9_n$. Check the shift relations against our stored data to $n=24$.
4. **Extract the chiral factor from free SYM.** A1's Pólya-counting derivation of $Z_{SO(10,2)}$ is an all-orders formula. Investigate whether it can be factorised to give the chiral $\mathrm{vac}_\ell$, which is the most promising route to an all-orders *construction*. Positivity would remain to be shown.
5. **Get the Curtright–Thorn–Goldstone papers (C1).** Check whether their SO(D−1) generating-function derivation generalises to a cumulative SO(D) statement.

---

## Follow-up (2026-10-12): arc spaces, βγ systems and pure-spinor characters

Context: `uplift_structure/` S6. Regge monotonicity of the SO(9) spectrum is proved for the bosonic analogue, using the arc space of the SO(9) null cone. We wanted to know whether the literature supplies the fermionic (superstring) extension.

Seeds:
- Berkovits–Nekrasov, *The character of pure spinors*, [hep-th/0503075](https://arxiv.org/abs/hep-th/0503075), [INSPIRE 678015](https://inspirehep.net/literature/678015).
- Aisaka–Arroyo–Berkovits–Nekrasov, *Pure spinor partition function and the massive superstring spectrum*, [arXiv:0806.0584](https://arxiv.org/abs/0806.0584), [INSPIRE 787223](https://inspirehep.net/literature/787223).

### Directly relevant
- **AABN (0806.0584).** They compute the (λ, ω) pure-spinor partition function, with spin characters, up to the fifth mass level. Adding the (x, θ, p) contribution reproduces the massive spectrum.
  - **Lesson.** The naive Hilbert space, gauge-invariant polynomials in λ, ω and derivatives, is *incomplete*. From the second mass level there are extra **fermionic** states, from functions not globally defined on pure-spinor space (higher Čech/Dolbeault cohomology), tied to the $b$ ghost.
  - Their spectrum is SO(10)-covariant (10D Lorentz) plus a light-cone twist (§6). It is not an SO(9)-covariant non-negative model.
- **Aisaka–Arroyo, *Hilbert space of curved βγ systems on quadric cones*** ([arXiv:0806.0586](https://arxiv.org/abs/0806.0586), [INSPIRE 787224](https://inspirehep.net/literature/787224)). This is the βγ system on exactly our type of cone (quadric).
  - BRST cohomology lives at ghost numbers 0 and 1, which are in bijection.
  - The ghost-number-1 states are the non-globally-defined ones.
- **Chandrasekhar–Ray, *Beta-gamma system, pure spinors and Hilbert series of arc spaces*** ([arXiv:1407.3762](https://arxiv.org/abs/1407.3762), [INSPIRE 1306332](https://inspirehep.net/literature/1306332)). They obtain field generating functions as Hilbert series of **arc spaces** of the constraint variety. This is the same object as our $B=\mathbb C[J_\infty(Q)]$, so it supports S6's use of it.

### Background / possibly useful later
- Nekrasov, *Lectures on curved βγ systems, pure spinors, and anomalies*, [hep-th/0511008](https://arxiv.org/abs/hep-th/0511008), [INSPIRE 696666](https://inspirehep.net/literature/696666).
- Movshev, *A formula for the partition function of the βγ system on the cone of pure spinors*, [arXiv:1602.04673](https://arxiv.org/abs/1602.04673), [INSPIRE 1421675](https://inspirehep.net/literature/1421675). It gives a closed formula via theta functions and Eisenstein series.
- Cederwall–Palmkvist, *Superalgebras, constraints and partition functions*, [arXiv:1503.06215](https://arxiv.org/abs/1503.06215), [INSPIRE 1355238](https://inspirehep.net/literature/1355238). Partition functions of constrained bosonic variables (pure spinors, null vectors) are encoded by Borcherds superalgebras.
- Boffo–Grassi–Hulik–Sachs, *On superparticles and their partition functions*, [arXiv:2402.09868](https://arxiv.org/abs/2402.09868), [INSPIRE 2758691](https://inspirehep.net/literature/2758691). Twisted world-line partition functions for spinning/superparticles.
- Grassi–Policastro–Scheidegger, *Partition functions, localization, and the chiral de Rham complex*, [hep-th/0702044](https://arxiv.org/abs/hep-th/0702044), [INSPIRE 744007](https://inspirehep.net/literature/744007).
- Cederwall, *Pure spinor superfields — an overview*, [arXiv:1307.1762](https://arxiv.org/abs/1307.1762), [INSPIRE 1241859](https://inspirehep.net/literature/1241859).

**Upshot.**
- The bosonic step of S6 (the null-cone arc space as a domain with the right Hilbert series) is standard ground.
- No paper found gives an SO(9)-covariant, manifestly non-negative model of the massive superstring spectrum.
- The pure-spinor experience warns that a correct model needs more than functions on a (super)variety: fermionic states arise from non-global functions.
