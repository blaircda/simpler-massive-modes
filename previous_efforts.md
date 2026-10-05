Human-authored summary as of 05/10/26.

I will use "baseline uplift" to refer to that previously established in earlier efforts. In our current status reports, this is referred to as "canonical". This would be a tentative definition of "canonical" - ultimately there may be some different selection criteria (physical or mathematical) that prefers some other uplift.  

One weakness intrinsic to previous efforts was that we focused too closely on the baseline uplift, and wasted some time trying to establish uniqueness.

Another trap that I have noticed is that there is a lot of possible structure in the proliferation of SO(9) representations level by level. This means that promising investigations can get sucked into a cycle of aggregating refinements as patterns change or fail to hold at high levels. Let us keep this in mind. 

Previous efforts/observations:

1) Plethystic Logarithm
The plethystic log (both bosonic and fermionic) was computed for the baseline uplift up to level 10. This somewhat brute force approach did not seem to lead to an especially simple structure, however most of our attention went to the bosonic log. There was at least one nicer cancellation at low levels in the fermionic log, which we did not follow-up. It might be useful to see this revisited more carefully.

2) Graded uplifts
Inpsired by the paper https://arxiv.org/abs/hep-th/0310292 by Beisert, Bianchi, Morales, Samtleben (BBMS), we considered whether the additional SO(2) grading selected a preferred uplift. We found that the uplift selected is not the baseline uplift. Note the BBMS approach is linked directly to type IIB rather than type IIA. We also experimented with some more IIA-style gradings, but did not establish a clear link to the baseline uplift. This led to us refocusing usefully on the non-uniqueness of the uplifts, but perhaps we did not exploit the grading approach as fully as we could have.

More details supplied by my collaborator:
BBMS does provide a selection criterion, but only within its own \(SO(10)\times SO(2)_\Delta\)-graded, type-IIB/\(AdS_5\times S^5\) problem. We found no reason to regard that grading as physically preferred for Bars’s type-IIA/M-theory proposal.
The final status was:
- We independently reconstructed the BBMS oscillator grading through level 7, retaining oscillator number \(\nu\) and imposing
\[
\Delta-J=\nu.
\]
Its refined branching agrees exactly with the independently constructed oscillator spectrum.
- After forgetting \(\Delta\) and cumulatively summing, the BBMS construction gives a positive integral \(SO(10)\) character at every tested horizon through level 7, with exactly the required cumulative \(SO(9)\) restriction.
- It agrees with the baseline uplift through level 5.
- Beginning at level 6 it differs from the protected continuation by outer-chirality kernel terms:
  \[
  [0,0,0,1,0]+[1,0,0,1,0]
  -
  [0,0,0,0,1]-[1,0,0,0,1].
  \]The analogous phenomenon continues at level 7. These differences vanish after restriction to \(SO(9)\), so both representatives lift the same target.
Inside the finite BBMS-graded domain, the BBMS answer was actually unique. Exact rank calculations showed that the baseline uplift could not be given the same \(\Delta\)-grading at levels 6 and 7. The obstruction persisted under several enlargements of the allowed \(\Delta\)-range. Thus the two representatives reflect genuinely different refined structures, rather than two solutions within one fixed refined problem.
But that does not establish that the BBMS refinement is better for Bars. Its physical justification comes from type IIB on \(AdS_5\times S^5\). We did not find an argument that the same \(\Delta\)-grading should govern the type-IIA strong-coupling/M-theory interpretation.
We also tried to obtain a more intrinsic selector:
- Combining opposite left- and right-moving chiralities did not remove the difference. It converted the chiral outer-odd ambiguity into an outer-even closed-string kernel element.
- Positivity, exact restriction and closed outer invariance still allowed multiple level-6 lifts.
- The tentative identification of Bars’s additional label \(\ell'\) with interlacing depth worked at low levels but failed starting at level 5.
- That depth was blind to the level-6 chirality difference anyway: the BBMS and protected representatives had identical depth data.
So the final interpretation was:
BBMS supplies a coherent, independently motivated alternative positive uplift through the tested range, but it does not establish a preferred Bars representative.
It remains potentially useful for the existence theorem. If its cumulative construction could be proved positive and correctly restricting at every level, it would prove Bars’s existence claim—despite producing a different representative from the protected continuation. What it would not prove is that its chirality choice is canonical, physically superior, or uniquely appropriate to the IIA/M-theory setting.

3) Regge trajectories and positive fibres
More recent effort has focused on understanding all possible positive uplifts at a given level. Taking into account non-uniqueness we realised we just want to show non-emptiness of the "positive fibre" (i.e. the set of uplifts which have positive coefficients). To do this, we make use of Regge shifts, the empirical support formulae of https://arxiv.org/abs/2609.04107 by Georgoudis, Minahan, Ström, Zoumis (GMSZ) and observations in https://arxiv.org/abs/1007.2622 by Hanany, Fourcella, Troost (HFT) about the factorisation of the massive mode partition function via plethystic exponentials. This approach is I think morally similar to what you have already started investigating in the uplift_structure directory. 

4) Closed string structure
One curious point is that the SO(10) uplift seems to occur naturally in a single left- or right-moving sector of the type IIA string, after dividing out by the universal massive supermultiplet. A key point for the full construction is that the universal massive supermultiplet itself does not admit a positive SO(10) uplift. Rather its tensor product, which is what is actually present in the closed IIA string, does admit such an uplift. For the rest of the structure, it does seem curious that we can directly uplift while forgetting about the final closed string nature. One complication that can arise is attempting to do the divison by the universal massive supermultiplet in SO(10), treating it as a virtual character. This led to difficulties in the past. 
