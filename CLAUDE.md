# Project Goal

The core reference for this project is the following paper:

**Itzhak Bars** (1995)

*Stringy evidence for $D = 11$ structure in strongly coupled type IIA superstring*

[hep-th/9503228](https://arxiv.org/abs/hep-th/9503228) | [INSPIRE-HEP recid: 393910](https://inspirehep.net/literature/393910)

This proposed both the minimal cumulative scheme (Section III.B) and the higher multiplicity scheme (Section III.C).

The principal goal is an all-orders **existence/admissibility proof** for the minimal Bars III.B scheme.

# Using SageMath

SageMath can be run using bin/sage-run 

# Literature search skill
When doing literature searches or exploring arXiv papers, read and follow the procedure in:
.agents/skills/hepth-arxiv-inspire-literature-exploration/SKILL.md

# LaTeX formatting in .md files
Do not use \begin{aligned} ... \end{aligned}, instead just use individual equations.

# Scientific Integrity & Derivation Rules
1. **No Circular Seeding:**
   Scripts that claim to derive, uplift, or reconstruct representations from first principles or generating functions must NEVER import baseline verification files (once the latter have been established). Comparison against baseline must always happen in a separate verification step.
2. **Obstructions Over Artificial Success:**
   Do not patch failing algorithms with ground-truth lookups to make milestone checks pass. If an algorithm encounters an infinite loop, negative multiplicities, or algebraic non-divisibility, stop and document the mathematical obstruction. Discovering why an uplift fails algebraically is just as valuable as a successful proof.
3. **No Dead-Code Claims:**
   Never declare a milestone passed based on a function that is bypassed or commented out in the actual execution flow. Every claimed computation must be actively executed in the script.
