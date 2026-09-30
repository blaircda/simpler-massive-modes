from sage.all import *

D4 = WeylCharacterRing("D4", style="coroots")

one = D4(0,0,0,0)
v8  = D4(1,0,0,0)       # 8_v
s8  = D4(0,0,0,1)       # 8_s
c8 = D4(0,0,1,0)
susy_factor = v8 + s8
# ------------------------------------------------------------
# Basic truncated q-series arithmetic
#
# A series is represented by
#     [a_0, a_1, ..., a_N]
# where a_n is a D4 character.
# ------------------------------------------------------------

def zero_series(N):
    return [0*one for _ in range(N+1)]


def add_series(A, B):
    N = min(len(A), len(B)) - 1
    return [A[n] + B[n] for n in range(N+1)]


def multiply_series(A, B, N):
    C = zero_series(N)

    for i in range(min(len(A), N+1)):
        for j in range(min(len(B), N-i+1)):
            C[i+j] += A[i] * B[j]

    return C

# ------------------------------------------------------------
# Factor belonging to oscillator mode m
#
#   bosons:
#       sum_k Sym^k(8_v) q^(m k)
#
#   fermions:
#       sum_k Lambda^k(8_s) q^(m k)
#
# We only need terms up to q^N.
# ------------------------------------------------------------

def oscillator_factor(m, N):

    F = zero_series(N)

    # At oscillator mode m:
    #
    #   bosonic occupation a contributes
    #       Sym^a(8_v) q^(m a)
    #
    #   fermionic occupation b contributes
    #       Lambda^b(8_s) q^(m b)
    #
    # The total occupation is n = a+b.

    for n in range(N//m + 1):

        coeff = 0 * one

        for a in range(n + 1):

            b = n - a

            if b <= 8:
                bos = v8.symmetric_power(a)
                ferm = s8.exterior_power(b)

                coeff += bos * ferm

        F[m*n] += coeff

    return F

# ------------------------------------------------------------
# Full Fock generating function through level N
# ------------------------------------------------------------

def fock_series(N):

    F = zero_series(N)

    # Start with 1
    F[0] = one

    # Multiply in oscillator modes m = 1,...,N
    for m in range(1, N+1):

        factor = oscillator_factor(m, N)

        F = multiply_series(F, factor, N)

    return F

# ------------------------------------------------------------
# factor out SUSY multiplet
# ------------------------------------------------------------

def dominant_D4(w):
    x = w.to_vector()
    return (
        x[0] >= x[1] >= x[2]
        and x[2] >= abs(x[3])
    )


def factor_susy_d4(fock_char):

    if isinstance(fock_char, list):
        fock_char = sum(fock_char, D4(0))

    rem = fock_char
    r_decomp = D4(0)

    # Highest weights of the two pieces of the SUSY factor.
    omega_v = v8.highest_weight()
    omega_s = s8.highest_weight()

    while rem != D4(0):

        coeffs = rem.monomial_coefficients()

        # Find dominant weights in the remaining character.
        dominant_weights = [
            w for w in coeffs
            if dominant_D4(w)
        ]

        if not dominant_weights:
            raise ValueError(
                "No dominant weight in remainder:\n"
                + str(rem)
            )

        # Highest dominant weight.
        hw = max(
            dominant_weights,
            key=lambda w: tuple(w.to_vector())
        )

        coeff = coeffs[hw]

        candidates = []

        # The highest weight can arise from either 8v or 8s.
        for omega in [omega_v, omega_s]:

            q_hw = hw - omega

            try:
                q_irrep = D4(q_hw)

                test = rem - coeff * (
                    susy_factor * q_irrep
                )

                # Require a genuine nonnegative weight character.
                test_coeffs = test.monomial_coefficients()

                if all(c >= 0 for c in test_coeffs.values()):
                    candidates.append(q_irrep)

            except (ValueError, TypeError):
                pass

        if len(candidates) != 1:
            raise ValueError(
                "Could not uniquely determine quotient term.\n"
                "Highest weight = " + str(hw) + "\n"
                "Coefficient = " + str(coeff) + "\n"
                "Candidates = " + str(candidates)
            )

        q_irrep = candidates[0]

        r_decomp += coeff * q_irrep

        rem -= coeff * (
            susy_factor * q_irrep
        )

    return r_decomp
