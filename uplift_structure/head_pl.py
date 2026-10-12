"""
S8: super-plethystic log (SO(10)) of the head series E(q) = sum_n E_n q^{n-1}
of the Regge induction, i.e. generators minus relations of the would-be head
algebra H = A/(y).  Class-level heads are resolved into chiralities degree by
degree so as to match the products predicted from lower degrees (the algebra
fixes chiralities; S4 allows any choice).  Usage: head_pl.py [DMAX] [rule]
Output: results/head_pl_<rule>.txt
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, SO10
sys.path.insert(0, SO10)
from regge_heads import run
from regge_stress import DEFAULT_EXT, load_R9
from weights import dynkin_D5
from sage.all import WeylCharacterRing, QQ
D5 = WeylCharacterRing("D5", style="coroots", base_ring=QQ)
N, r9, R9 = load_R9(DEFAULT_EXT)
DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
rule = sys.argv[2] if len(sys.argv) > 2 else "max-height"
H = run(rule, r9, DMAX + 1)
Ecl = [H[n] for n in range(1, DMAX + 2)]          # class-level heads, degree d = n-1
def amb(w): return tuple(int(2*x) for x in w.to_vector())
def spin(w): return any((2*x) % 2 for x in w.to_vector())
def split(ch):
    t, s = D5.zero(), D5.zero()
    for w, m in ch.monomial_coefficients().items():
        if spin(w): s += m*D5(w)
        else: t += m*D5(w)
    return t, s
def psi(ch, k):
    if ch == 0: return D5.zero()
    wm = {}
    for w, m in ch.weight_multiplicities().items():
        wm[k*w] = wm.get(k*w, 0) + m
    return D5.char_from_weights(wm)
def ch_of(lam): return D5(*dynkin_D5(lam))
def pe_degree(P, d):
    """degree-d coefficient of PE_super[sum_{j<d} P_j q^j]"""
    Gp = [None]
    for j in range(1, d + 1):
        acc = D5.zero()
        for k in range(1, j + 1):
            if j % k or j // k >= d or P[j // k] == 0: continue
            t, s = split(P[j // k])
            acc += (j // k) * (psi(t, k) + (-1) ** (k + 1) * psi(s, k))
        Gp.append(acc)
    Ep = [D5.one()]
    for m in range(1, d + 1):
        acc = D5.zero()
        for j in range(1, m + 1):
            if Gp[j] != 0 and Ep[m - j] != 0: acc += Gp[j] * Ep[m - j]
        Ep.append(acc / m)
    return Ep[d]
P = [None] + [D5.zero()] * DMAX
lines = []
for d in range(1, DMAX + 1):
    pred = pe_degree(P, d)
    pc = {amb(w): c for w, c in pred.monomial_coefficients().items()}
    E = D5.zero()
    for lam, c in Ecl[d].items():
        if lam[4] == 0:
            E += c * ch_of(lam); continue
        bar = lam[:4] + (-lam[4],)
        pp, pm = max(0, pc.get(lam, 0)), max(0, pc.get(bar, 0))
        n_plus = min(c, pp); n_minus = min(c - n_plus, pm); n_plus += c - n_plus - n_minus
        E += n_plus * ch_of(lam) + n_minus * ch_of(bar)
    P[d] = E - pred
    terms = []
    for w, cc in sorted(P[d].monomial_coefficients().items(), key=lambda t: (D5(t[0]).degree(), amb(t[0]))):
        cc = int(cc) if cc == int(cc) else cc
        terms.append(f"{'+' if cc > 0 else ''}{cc}x{D5(w).degree()}{'s' if spin(w) else ''}{list(dynkin_D5(amb(w)))}")
    s = f"[{rule}] PL(heads)_{d} (level {d+1}): " + (" ".join(terms) or "0")
    print(s, flush=True)
    lines.append(s)
os.makedirs(RESULTS, exist_ok=True)
with open(os.path.join(RESULTS, f"head_pl_{rule}.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")
