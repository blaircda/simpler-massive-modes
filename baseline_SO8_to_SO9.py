from sage.all import *
D4 = WeylCharacterRing("D4", style="coroots")
B4 = WeylCharacterRing("B4", style="coroots")

########################################################################
# SO8 -> SO9
########################################################################

def greedy_so8_to_so9(r):
    """
    Greedy SO(8) -> SO(9) lift.

    Input:
        r  -- a D4 Weyl character

    Output:
        R  -- a B4 Weyl character whose restriction to D4
              reproduces r, provided the greedy procedure succeeds.

    At each step:
      1. Find the highest remaining D4 irreducible.
      2. Convert its highest weight to ambient coordinates.
      3. Use those ambient coordinates as the B4 highest weight.
      4. Subtract the complete D4 branching of that B4 irrep.
      5. Repeat.
    """

    rem = r
    R = B4(0)

    while rem != D4(0):

        # Find the highest remaining D4 constituent.
        coeffs = rem.monomial_coefficients()

        candidates = []

        for w, coeff in coeffs.items():
            if coeff <= 0:
                continue

            irrep = D4(w)
            hw = irrep.highest_weight()

            candidates.append((hw, irrep, coeff))

        if not candidates:
            raise ValueError(
                "No positive D4 constituent found in remainder."
            )

        # Highest ambient D4 weight.
        hw, irrep, coeff = max(
            candidates,
            key=lambda x: tuple(x[0][i] for i in range(4))
        )

        # The SO(9) highest weight is the same ambient tuple,
        # with the D4 chirality sign removed.
        mu = tuple(hw[i] for i in range(4))
        mu = (mu[0], mu[1], mu[2], abs(mu[3]))

        # Convert ambient B4 weight -> Dynkin labels.
        l1, l2, l3, l4 = mu

        B4_hw = B4(
            l1-l2,
            l2-l3,
            l3-l4,
            2*l4
        )

        # Branch the candidate SO(9) irrep to SO(8).
        branch = B4_hw.branch(D4)

        # We need the whole branch to fit in the remainder.
        difference = rem - coeff * branch

        if any(c < 0 for c in difference.monomial_coefficients().values()):
            raise ValueError(
                "Greedy lift failed: candidate branch does not fit.\n"
                f"D4 highest weight = {hw}\n"
                f"candidate SO(9) irrep = {B4_hw}\n"
                f"coefficient = {coeff}\n"
                f"branch = {branch}\n"
                f"remainder = {rem}"
            )

        # Record the SO(9) representation.
        R += coeff * B4_hw

        # Remove its SO(8) content.
        rem = difference

    return R





