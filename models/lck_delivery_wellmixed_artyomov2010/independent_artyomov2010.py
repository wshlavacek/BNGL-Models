"""Level-1 verification: the fate of one agonist pMHC, solved exactly.

In the well-mixed limit, and while the agonist pMHC are a small tracer population
against the TCR and coreceptor pools, each agonist molecule is independent of the
others, so the fate of one of them is a continuous-time Markov chain with a single
absorbing state and can be solved by matrix exponential rather than simulated.

The state is the bond pattern of the little complex around one agonist MHC: which of
MHC-TCR, MHC-coreceptor and TCR-coreceptor are present. Six transient states and one
absorbing one:

    m_t m_c t_c
    1   0   0     MHC on TCR                 1 1 0   both, ring open
    0   1   0     MHC on coreceptor          1 1 1   the closed ring
    0   1   1     coreceptor bridges         1 0 1   TCR bridges
    0   0   0     MHC holds nothing -> destroyed at once, so this is absorption

The TCR-coreceptor bond has to be in the chain even though it touches neither of the
MHC's own sites. It is what lets the ring close, and a ring closure is INTRAMOLECULAR:
it happens at the full local rate -- 150/s for MHC-TCR -- where the intermolecular
route is the mean-field 150/10^4 * 200 = 3/s, fifty times slower. That is the whole
cooperativity the paper is about. Leaving the bond out gives a chain that reports
almost no coreceptor enhancement, which is how the omission was caught.

The pseudo-first-order on-rates use the free pool sizes, which barely move because
the tracer is small; `verify_artyomov2010.ipynb` checks that assumption against the
model's own free-pool observables rather than asserting it.

Nothing here reads the .bngl, the .net or any simulator output. The rate constants
are Table 1 of the paper, and the mean-field division by the chamber count is the
paper's own SI Methods Eq. S1b.
"""
import numpy as np
from scipy.linalg import expm

# Table 1, per second. The three "local" constants are per chamber.
KON_LOCAL = 150.0        # TCR-MHC on
KOFF_AG = 0.02           # TCR-MHC off, agonist
KONCD_LOCAL = 1000.0     # coreceptor-MHC on
KONLCK_LOCAL = 1.0       # coreceptor-TCR on, via MHC
KONLCK1_LOCAL = 1.0      # coreceptor-TCR on, neither on MHC
KOFF_LCK = 1.0           # coreceptor-TCR off
KDIFF = 50.0             # hop rate, used only by the encounter correction
N_CHAMBERS = 1e4

# (m_t, m_c, t_c); "000" is absorbing
STATES = ("100", "010", "110", "111", "011", "101", "000")


def generator(koff_cd, n_tcr_free, n_cd4_free, encounter=False):
    """Rate matrix of the seven-state chain, in the mean-field limit."""
    # Intermolecular constants are per chamber, so the mean-field value divides by
    # the chamber count (SI Methods Eq. S1b). Ring closures are intramolecular and
    # keep the local constant.
    scale = 1.0 / N_CHAMBERS
    kon, kon_cd = KON_LOCAL * scale, KONCD_LOCAL * scale
    klck, klck1 = KONLCK_LOCAL * scale, KONLCK1_LOCAL * scale
    if encounter:
        # The paper's own encounter probability k/(k + m_A + m_B).
        kon *= KON_LOCAL / (KON_LOCAL + 2 * KDIFF)
        kon_cd *= KONCD_LOCAL / (KONCD_LOCAL + 2 * KDIFF)
        klck *= KONLCK_LOCAL / (KONLCK_LOCAL + 2 * KDIFF)
        klck1 *= KONLCK1_LOCAL / (KONLCK1_LOCAL + 2 * KDIFF)

    a = kon * n_tcr_free           # a free TCR finds this MHC
    b = kon_cd * n_cd4_free        # a free coreceptor finds this MHC
    c = klck1 * n_cd4_free         # a free coreceptor finds the bound TCR
    d = klck * n_tcr_free          # a free TCR finds the bound coreceptor

    i = {s: k for k, s in enumerate(STATES)}
    q = np.zeros((len(STATES), len(STATES)))

    q[i["100"], i["000"]] = KOFF_AG          # last bond goes: lost
    q[i["100"], i["110"]] = b
    q[i["100"], i["101"]] = c

    q[i["010"], i["000"]] = koff_cd
    q[i["010"], i["110"]] = a
    q[i["010"], i["011"]] = d

    q[i["110"], i["010"]] = KOFF_AG
    q[i["110"], i["100"]] = koff_cd
    q[i["110"], i["111"]] = KONLCK_LOCAL     # ring closes, intramolecular

    q[i["111"], i["011"]] = KOFF_AG
    q[i["111"], i["101"]] = koff_cd
    q[i["111"], i["110"]] = KOFF_LCK

    q[i["011"], i["000"]] = koff_cd          # coreceptor lets go, MHC held by nothing
    q[i["011"], i["010"]] = KOFF_LCK
    q[i["011"], i["111"]] = KON_LOCAL        # ring closes, intramolecular

    q[i["101"], i["000"]] = KOFF_AG
    q[i["101"], i["100"]] = KOFF_LCK
    q[i["101"], i["111"]] = KONCD_LOCAL      # ring closes, intramolecular

    np.fill_diagonal(q, -q.sum(axis=1))
    return q


def bound_to_tcr(t, koff_cd, n_tcr_free=200.0, n_cd4_free=100.0, encounter=False):
    """P(MHC still bound to TCR) at each time in `t`, starting from S10.

    This is the SSC listing's `record MHC(p="ag",t#1) TCR(m#1)`, as a fraction.
    """
    q = generator(koff_cd, n_tcr_free, n_cd4_free, encounter)
    p0 = np.zeros(len(STATES))
    p0[STATES.index("100")] = 1.0                  # seeded already on a TCR
    keep = np.array([1.0 if s[0] == "1" else 0.0 for s in STATES])   # m_t present
    return np.array([p0 @ expm(q * ti) @ keep for ti in np.atleast_1d(t)])


def half_time(koff_cd, **kw):
    """Time at which the TCR-bound fraction falls to one half."""
    t = np.linspace(0.0, 600.0, 6001)
    y = bound_to_tcr(t, koff_cd, **kw)
    i = int(np.argmax(y <= 0.5))
    if i == 0:
        return None
    t0, t1, v0, v1 = t[i - 1], t[i], y[i - 1], y[i]
    return float(t0 + (v0 - 0.5) * (t1 - t0) / (v0 - v1))


def half_time_no_coreceptor():
    """With no coreceptor the chain collapses to a single exponential."""
    return float(np.log(2.0) / KOFF_AG)


if __name__ == "__main__":
    print(f"no coreceptor (closed form ln2/koffAg): {half_time_no_coreceptor():.2f} s\n")
    print(f"{'koffCD':>7}{'mean field':>12}{'encounter':>11}{'published':>11}")
    published = {20: 45.0, 40: 43.0, 60: 40.5, 80: 38.5, 100: 37.0, 110: 36.2, 150: 35.0}
    for k, pub in published.items():
        print(f"{k:>7}{half_time(k):>12.2f}{half_time(k, encounter=True):>11.2f}{pub:>11.1f}")
