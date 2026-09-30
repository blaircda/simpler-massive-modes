# Code to compute the massive states of the IIA superstring
# in either left or right moving sector
# firstly as SO8 oscillator states F_l in lightcone gauge
# F_l = r_l^8 otimes (8_v + 8_s)
# then factoring out the universal SUSY part and focusing on the massive states
# we find the uplift
# r_l^8 -> r_l^9
# where r_l^9 are SO9 states
# as required by the fact SO9 is the little group for massive states in 10d
# so that the physical state at level l is
#  r_l^9 otimes ( 44 + 84 + 128 )
# where the common massive supermultiplet ( 44 + 84 + 128 )
# then we want to verify that the representations
# R_l^9 = sum_{i = 1}^l r_i^9
# form SO(10) multiplets
# i.e. there exists
# R_n^10 s.t. R_n^10 -> R_n^9 under branching
# this code verifies this using a previously found set of SO(10) multiplets up to level 10

from sage.all import *
from baseline_SO8_construct import (
    fock_series, factor_susy_d4, susy_factor
)
from baseline_SO8_to_SO9 import greedy_so8_to_so9
from baseline_hardcoded_level10_SO10 import (
    SO10_SUPERMULTIPLETS_DYNKIN, get_so10_character
)
D5 = WeylCharacterRing("D5", style="coroots")
D4 = WeylCharacterRing("D4", style="coroots")
B4 = WeylCharacterRing("B4", style="coroots")

# ------------------------------------------------------------
# Calculate first N levels of the massive state SO8 Fock space
# ------------------------------------------------------------

N = 10

print("\n============================================================")
print(f"SO8 Fock space up to level {N}")
print("============================================================")

F = fock_series(N)

for n in range(N+1):
    print()
    print("F_%d =" % n)
    print(F[n])
    print("dimension =", F[n].degree())
    
# ------------------------------------------------------------
# factor out SUSY multiplet
# ------------------------------------------------------------
r_list = []

print("\n============================================================")
print("Result of factoring out SUSY multiplet}")
print("============================================================")

for n in range(1, N+1):
    r = factor_susy_d4(F[n])
    r_list.append(r)

    print(f"r_{n}^8 =", r)
    assert susy_factor * r == F[n], f"SUSY factor reconstruction does not match at level {i+1}"
    print()

# ------------------------------------------------------------
# SO8 -> SO9
# ------------------------------------------------------------

print("\n============================================================")
print("Uplifting to SO(9)")
print("============================================================")

R_list = []

for n, r in enumerate(r_list, start=1):

    print()
    print("LEVEL", n)
    print(f"r^8_{n} =", r)

    R = greedy_so8_to_so9(r)

    print(f"r^9_{n} =", R)
    print("degree =", R.degree())
    print("branch =", R.branch(D4))
    assert R.branch(D4) == r, f"SO(9) reconstruction does not match at level {i+1}"
    R_list.append(R)

# ------------------------------------------------------------
# SO9 -> SO10
# ------------------------------------------------------------

print("\n============================================================")
print("Verifying hardcoded previously computed SO(10) uplift}")
print("============================================================")

proposed_SO10 = [ get_so10_character(n) for n in range(1,N+1) ]
proposed_SO9 = [ sum(R_list[:n]) for n in range(1,N+1) ]

for i, S in enumerate(proposed_SO10):
    R = S.branch(B4, rule="default")
    target = proposed_SO9[i]
    print()
    print(f"Level {i+1}")
    print(f"R^{10}_{i+1}\n", S)
    print(f"Branched R^{10}_{i+1}\n", R)
    print(f"Target R^9_{i+1}\n", target)
    assert R == target, f"SO(10) restriction does not match at level {i+1}"
