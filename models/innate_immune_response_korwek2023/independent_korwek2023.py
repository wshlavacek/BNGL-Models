"""Independent implementation of the Korwek et al. (2023) innate-immunity ODE system.

Written by hand from the *reaction rules* of the curated BNGL file — never from the
BioNetGen-generated `.net` file — so that integrating it and comparing against
BioNetGen independently checks three things BioNetGen does for us:

1. **Network generation.** The 53 species and 96 reactions the rules imply are
   enumerated here explicitly, by name, so a missing or spurious species shows up.
2. **Functional rate laws.** The model's transcription, transcript-degradation and
   translation rates are saturating functions of observables rather than constants,
   and BioNetGen multiplies a rate law by the reactant amounts automatically. The
   competitive NF-kB/IkBa promoter terms, the RNase L transcript-degradation factor,
   the phospho-eIF2alpha translation-arrest factor and the STAT saturation terms are
   all reconstructed here from table S1 of Korwek et al. (2023).
3. **Observable patterns.** Several observables match more than the obvious species:
   `IFNb_ext` counts the receptor-bound ligand as well as the free extracellular
   pool, `STAT1_p` and `STAT1_total` count the STAT1/2 heterodimer, `RIG_I_total`
   counts RIG-I inside both the binary poly(I:C) complex and the ternary
   poly(I:C):RIG-I:MAVS complex, and the poly(I:C) that switches on PKR and OAS3 is
   the *total* cytoplasmic pool, free and complexed. Each is written out below.

Amounts are dimensionless (the published model is nondimensionalized) and time is in
seconds. Used by `verify_korwek2023.ipynb`.
"""
import numpy as np
from scipy.integrate import solve_ivp

# --------------------------------------------------------------------------- #
# Parameters (table S1 of Korwek et al. 2023; values as in the curated BNGL)
# --------------------------------------------------------------------------- #
P = dict(
    k_v=5.0,
    h_Mavs=1.0, h_Pkr_gene=1.0, h_Rnasel_gene=1.0, h_A20_gene=1.0,
    EPSILON=1.0e-6,
    k_FAST=1.0,
    k_POLYIC=0.00024390060510912293,
    tg_TRANSCRIPT=7.259610076954435e-05,
    sg_PROTEIN=1.532955086663982e-05,
    m_Rnasel=0.031164323577208185,
    m_Eif2a=0.015830114342212334,
    p_Eif2a_basal=7.723454616604743e-06,
    ma_Rigi_gene_basal=0.05930147549419551,
    ma_Pkr_Oas3_gene_basal=0.7245345686782008,
    ma_Rnasel_gene_basal=1.020172036452892,
    sg_Pkr=4.2805996953166866e-05,
    a_Tak1_by_Tnfa=0.8903603794354649,
    a_Ikk=0.010904467109555071,
    d_Ikk_1=0.0008829700327658654,
    d_Ikk_2=0.05982062403935003,
    d_Ikk_3=0.00026826640310831975,
    p_Ikba_by_Ikk=0.008458552410073877,
    g_Ikba_u_with_Nfkb=3.944402357840543e-05,
    i_Ikba=0.0005635774134351319,
    e_Ikba=0.02051107491820277,
    tg_Ikba_mrna=0.0002565157989928973,
    a_Ikba_gene_by_Nfkb__=0.0684701851861643,
    tg_A20_mrna=0.0009406825155497426,
    a_A20_gene_by_Nfkb__=0.025686603943611396,
    s_Ikba=0.004255786898310512,
    sg_A20=5.014115674793064e-06,
    q_Tbk1_by_A20__=3.875729946075596e2,
    b_Ifnar_Ifnb_ext=0.000750779497463291,
    m_Rnasel_Ifnar_mrna=0.0039023754286057065,
    s_Ifnar=2.8778555528751826e-05,
    g_Ifnar=9.77848544160331e-05,
    g_Ifnar_w_Ifnb=0.00031054256866899256,
    m_Ifnb_mrna_NfkbIrf3=3.162948461272814e-05,
    q_Stat=0.0018294332440515334,
    m_Ifnar_a=0.020141727438512898,
    qu_Stat1_Stat2=0.052370747192221206,
    ma_Stat1_gene_basal=0.18314838022720048,
    ma_Stat2_gene_basal=0.07023843113385617,
    a_gene_by_Stat12dim=13126.669661929025,
)


def derived(p):
    """Parameters the BNGL file defines in terms of others."""
    d = dict(p)
    kp, kf, tg, sg = p["k_POLYIC"], p["k_FAST"], p["tg_TRANSCRIPT"], p["sg_PROTEIN"]
    d.update(
        i_Polyic=kp, b_Rigi_Polyic=kp, b_RigiPolyic_Mavs=kp,
        a_Pkr_by_Polyic=kp, d_Pkr=kp, a_Oas3_by_Polyic=kp, d_Oas3=kp,
        a_Rnasel_by_Oas3=kp, d_Rnasel=kp,
        p_Eif2a_by_Pkr=kp, q_Eif2a=kp, p_Irf3_by_Tbk1=kp, q_Irf3=kp,
        tg_Isg_mrna=tg, tg_Ifnar_mrna=tg, tg_Ifnb_mrna=tg, tg_Stat_mrna=tg,
        sg_Rigi=sg, sg_Oas3=sg, sg_Rnasel=sg, sg_Ifnb=sg, sg_Stat=sg,
        a_Tak1_by_RigiMavs=kf, d_Tak1=kf,
        b_Nfkb_Ikba_cyt=kf, b_Nfkb_Ikba_nuc=p["k_v"] * kf,
        g_Ikba_p_any=kf, g_Ikba_u_free=sg,
        i_Nfkb=kf, e_Nfkb_with_Ikba=kf,
        a_Ikba_gene_by_Nfkb=kf * p["a_Ikba_gene_by_Nfkb__"],
        d_Ikba_gene_by_Ikba=kf,
        a_A20_gene_by_Nfkb=kf * p["a_A20_gene_by_Nfkb__"],
        d_A20_gene_by_Ikba=kf,
        p_Tbk1_by_RigiMavs=kf, q_Tbk1=kf,
        q_Tbk1_by_A20=p["q_Tbk1_by_A20__"] * kf,
        b_Ifnar_Ifnb_cyt=kf, p_Stat=kf, b_Stat1_Stat2=kf,
    )
    return d


SPECIES = [
    # poly(I:C) module (18)
    "PolyIC_ext", "PolyIC_cyt", "RIGI", "RIGI_PolyIC", "RIGI_PolyIC_MAVS", "MAVS",
    "PKR_i", "PKR_a", "OAS3_i", "OAS3_a", "RNaseL_i", "RNaseL_a",
    "eIF2a_0", "eIF2a_p",
    "RIGI_mRNA", "PKR_mRNA", "OAS3_mRNA", "RNaseL_mRNA",
    # NF-kB module (18)
    "TNFa", "TAK1_i", "TAK1_a", "IKK_n", "IKK_a", "IKK_i", "IKK_ii",
    "C_cyt0", "C_cytpp", "NFkB_cyt", "NFkB_nuc", "C_nuc0",
    "IkBa_cyt0", "IkBa_cytpp", "IkBa_nuc0", "A20", "IkBa_mRNA", "A20_mRNA",
    # IRF3 module (4)
    "TBK1_0", "TBK1_p", "IRF3_0", "IRF3_p",
    # IFN-beta module (6)
    "IFNAR", "IFNAR_IFNb", "IFNb_cyt", "IFNb_ext", "IFNAR_mRNA", "IFNb_mRNA",
    # STAT1/2 module (7)
    "STAT1_0", "STAT1_p", "STAT2_0", "STAT2_p", "STAT12", "STAT1_mRNA", "STAT2_mRNA",
]
IX = {s: i for i, s in enumerate(SPECIES)}
assert len(SPECIES) == 53


def observables(y):
    """The BNGL observables, evaluated from the species vector."""
    g = lambda n: y[IX[n]]
    polyic_cyt_tot = g("PolyIC_cyt") + g("RIGI_PolyIC") + g("RIGI_PolyIC_MAVS")
    return dict(
        RIG_I_total=g("RIGI") + g("RIGI_PolyIC") + g("RIGI_PolyIC_MAVS"),
        PKR_total=g("PKR_i") + g("PKR_a"),
        OAS3_total=g("OAS3_i") + g("OAS3_a"),
        RNaseL_total=g("RNaseL_i") + g("RNaseL_a"),
        RNaseL_a=g("RNaseL_a"),
        eIF2a_total=g("eIF2a_0") + g("eIF2a_p"),
        eIF2a_p=g("eIF2a_p"),
        RIGI_mRNA=g("RIGI_mRNA"), PKR_mRNA=g("PKR_mRNA"),
        OAS3_mRNA=g("OAS3_mRNA"), RNaseL_mRNA=g("RNaseL_mRNA"),
        TAK1_a=g("TAK1_a"),
        NFkB_nuc_free=g("NFkB_nuc"),
        NFkB_nuc_total=g("NFkB_nuc") + g("C_nuc0"),
        NFkB_total=g("NFkB_cyt") + g("NFkB_nuc") + g("C_cyt0") + g("C_cytpp") + g("C_nuc0"),
        IkBa_total=(g("IkBa_cyt0") + g("IkBa_cytpp") + g("IkBa_nuc0")
                    + g("C_cyt0") + g("C_cytpp") + g("C_nuc0")),
        IkBa_nuc_total=g("IkBa_nuc0") + g("C_nuc0"),
        IkBa_cyt_total=g("IkBa_cyt0") + g("IkBa_cytpp") + g("C_cyt0") + g("C_cytpp"),
        IkBa_cyt_free=g("IkBa_cyt0") + g("IkBa_cytpp"),
        IkBa_p_cyt=g("IkBa_cytpp") + g("C_cytpp"),
        IkBa_nuc_free=g("IkBa_nuc0"),
        A20=g("A20"), IKK_a=g("IKK_a"),
        IkBa_mRNA=g("IkBa_mRNA"), A20_mRNA=g("A20_mRNA"),
        IRF3_total=g("IRF3_0") + g("IRF3_p"), IRF3_p=g("IRF3_p"),
        IFNAR_total=g("IFNAR") + g("IFNAR_IFNb"), IFNAR_a=g("IFNAR_IFNb"),
        IFNb_ext=g("IFNb_ext") + g("IFNAR_IFNb"), IFNb_cyt=g("IFNb_cyt"),
        IFNAR_mRNA=g("IFNAR_mRNA"), IFNb_mRNA=g("IFNb_mRNA"),
        STAT1_total=g("STAT1_0") + g("STAT1_p") + g("STAT12"),
        STAT2_total=g("STAT2_0") + g("STAT2_p") + g("STAT12"),
        STAT1_p=g("STAT1_p") + g("STAT12"),
        STAT2_p=g("STAT2_p") + g("STAT12"),
        STAT1_u=g("STAT1_0"), STAT2_u=g("STAT2_0"),
        STAT12_dimer=g("STAT12"),
        STAT1_mRNA=g("STAT1_mRNA"), STAT2_mRNA=g("STAT2_mRNA"),
        _polyic_cyt_tot=polyic_cyt_tot,
    )


def rhs(t, y, p):
    o = observables(y)
    g = lambda n: y[IX[n]]
    d = np.zeros_like(y)

    def flux(rate, *consumed, produce=()):
        for s in consumed:
            d[IX[s]] -= rate
        for s in produce:
            d[IX[s]] += rate

    # ---------------- poly(I:C) module ----------------
    flux(p["i_Polyic"] * g("PolyIC_ext"), "PolyIC_ext", produce=("PolyIC_cyt",))
    flux(p["b_Rigi_Polyic"] * g("RIGI") * g("PolyIC_cyt"),
         "RIGI", "PolyIC_cyt", produce=("RIGI_PolyIC",))
    flux(p["sg_Rigi"] * g("RIGI_PolyIC"), "RIGI_PolyIC", produce=("PolyIC_cyt",))
    flux(p["b_RigiPolyic_Mavs"] * g("RIGI_PolyIC") * g("MAVS"),
         "RIGI_PolyIC", "MAVS", produce=("RIGI_PolyIC_MAVS",))
    flux(p["sg_Rigi"] * g("RIGI_PolyIC_MAVS"),
         "RIGI_PolyIC_MAVS", produce=("PolyIC_cyt", "MAVS"))

    pic = o["_polyic_cyt_tot"]
    flux(p["a_Pkr_by_Polyic"] * pic * g("PKR_i"), "PKR_i", produce=("PKR_a",))
    flux(p["d_Pkr"] * g("PKR_a"), "PKR_a", produce=("PKR_i",))
    flux(p["a_Oas3_by_Polyic"] * pic * g("OAS3_i"), "OAS3_i", produce=("OAS3_a",))
    flux(p["d_Oas3"] * g("OAS3_a"), "OAS3_a", produce=("OAS3_i",))

    flux(p["p_Eif2a_by_Pkr"] * g("PKR_a") * g("eIF2a_0"), "eIF2a_0", produce=("eIF2a_p",))
    flux(p["p_Eif2a_basal"] * g("eIF2a_0"), "eIF2a_0", produce=("eIF2a_p",))
    flux(p["q_Eif2a"] * g("eIF2a_p"), "eIF2a_p", produce=("eIF2a_0",))

    flux(p["a_Rnasel_by_Oas3"] * g("OAS3_a") * g("RNaseL_i"), "RNaseL_i", produce=("RNaseL_a",))
    flux(p["d_Rnasel"] * g("RNaseL_a"), "RNaseL_a", produce=("RNaseL_i",))

    dim = p["a_gene_by_Stat12dim"] * o["STAT12_dimer"]
    isg = p["tg_Isg_mrna"]
    for mrna, basal, gate in [("RIGI_mRNA", p["ma_Rigi_gene_basal"], 1.0),
                              ("PKR_mRNA", p["ma_Pkr_Oas3_gene_basal"], p["h_Pkr_gene"]),
                              ("OAS3_mRNA", p["ma_Pkr_Oas3_gene_basal"], 1.0),
                              ("RNaseL_mRNA", p["ma_Rnasel_gene_basal"], p["h_Rnasel_gene"])]:
        d[IX[mrna]] += gate * isg * (basal + dim) / (basal + dim + 1.0)
    deg = isg * (p["m_Rnasel"] + o["RNaseL_a"]) / p["m_Rnasel"]
    for mrna in ("RIGI_mRNA", "PKR_mRNA", "OAS3_mRNA", "RNaseL_mRNA"):
        d[IX[mrna]] -= deg * g(mrna)

    d[IX["RIGI"]] += p["sg_Rigi"] * g("RIGI_mRNA")
    d[IX["PKR_i"]] += p["sg_Pkr"] * g("PKR_mRNA")
    d[IX["OAS3_i"]] += p["sg_Oas3"] * g("OAS3_mRNA")
    d[IX["RNaseL_i"]] += p["sg_Rnasel"] * g("RNaseL_mRNA")

    d[IX["RIGI"]] -= p["sg_Rigi"] * g("RIGI")
    for s in ("PKR_i", "PKR_a"):
        d[IX[s]] -= p["sg_Pkr"] * g(s)
    for s in ("OAS3_i", "OAS3_a"):
        d[IX[s]] -= p["sg_Oas3"] * g(s)
    for s in ("RNaseL_i", "RNaseL_a"):
        d[IX[s]] -= p["sg_Rnasel"] * g(s)

    # ---------------- NF-kB module ----------------
    flux(p["a_Tak1_by_RigiMavs"] * g("RIGI_PolyIC_MAVS") * g("TAK1_i"),
         "TAK1_i", produce=("TAK1_a",))
    flux(p["a_Tak1_by_Tnfa"] * g("TNFa") * g("TAK1_i"), "TAK1_i", produce=("TAK1_a",))
    flux(p["d_Tak1"] * g("TAK1_a"), "TAK1_a", produce=("TAK1_i",))

    flux(p["a_Ikk"] * o["TAK1_a"] ** 2 * g("IKK_n"), "IKK_n", produce=("IKK_a",))
    flux(p["d_Ikk_1"] / p["d_Ikk_2"] * (p["d_Ikk_2"] + o["A20"]) * g("IKK_a"),
         "IKK_a", produce=("IKK_i",))
    flux(p["d_Ikk_3"] * g("IKK_i"), "IKK_i", produce=("IKK_ii",))
    flux(p["d_Ikk_3"] * g("IKK_ii"), "IKK_ii", produce=("IKK_n",))

    flux(p["b_Nfkb_Ikba_cyt"] * g("IkBa_cyt0") * g("NFkB_cyt"),
         "IkBa_cyt0", "NFkB_cyt", produce=("C_cyt0",))
    flux(p["b_Nfkb_Ikba_nuc"] * g("IkBa_nuc0") * g("NFkB_nuc"),
         "IkBa_nuc0", "NFkB_nuc", produce=("C_nuc0",))

    flux(p["p_Ikba_by_Ikk"] * g("IKK_a") * g("IkBa_cyt0"), "IkBa_cyt0",
         produce=("IkBa_cytpp",))
    flux(p["p_Ikba_by_Ikk"] * g("IKK_a") * g("C_cyt0"), "C_cyt0", produce=("C_cytpp",))

    flux(p["g_Ikba_p_any"] * g("IkBa_cytpp"), "IkBa_cytpp")
    flux(p["g_Ikba_p_any"] * g("C_cytpp"), "C_cytpp", produce=("NFkB_cyt",))
    flux(p["g_Ikba_u_free"] * g("IkBa_cyt0"), "IkBa_cyt0")
    flux(p["g_Ikba_u_with_Nfkb"] * g("C_cyt0"), "C_cyt0", produce=("NFkB_cyt",))

    flux(p["i_Nfkb"] * g("NFkB_cyt"), "NFkB_cyt", produce=("NFkB_nuc",))
    flux(p["e_Nfkb_with_Ikba"] * g("C_nuc0"), "C_nuc0", produce=("C_cyt0",))
    flux(p["i_Ikba"] * g("IkBa_cyt0"), "IkBa_cyt0", produce=("IkBa_nuc0",))
    flux(p["e_Ikba"] * g("IkBa_nuc0"), "IkBa_nuc0", produce=("IkBa_cyt0",))

    nf, ik = o["NFkB_nuc_free"], o["IkBa_nuc_free"]
    d[IX["IkBa_mRNA"]] += (p["tg_Ikba_mrna"] * p["a_Ikba_gene_by_Nfkb"] * nf
                           / (p["a_Ikba_gene_by_Nfkb"] * nf
                              + p["d_Ikba_gene_by_Ikba"] * ik + p["EPSILON"]))
    d[IX["IkBa_mRNA"]] -= (p["tg_Ikba_mrna"] * (p["m_Rnasel"] + o["RNaseL_a"])
                           / p["m_Rnasel"] * g("IkBa_mRNA"))
    d[IX["A20_mRNA"]] += (p["h_A20_gene"] * p["tg_A20_mrna"]
                          * p["a_A20_gene_by_Nfkb"] * nf
                          / (p["a_A20_gene_by_Nfkb"] * nf
                             + p["d_A20_gene_by_Ikba"] * ik + p["EPSILON"]))
    d[IX["A20_mRNA"]] -= (p["tg_A20_mrna"] * (p["m_Rnasel"] + o["RNaseL_a"])
                          / p["m_Rnasel"] * g("A20_mRNA"))

    arrest = p["m_Eif2a"] / (p["m_Eif2a"] + o["eIF2a_p"])
    d[IX["IkBa_cyt0"]] += p["s_Ikba"] * arrest * g("IkBa_mRNA")
    d[IX["A20"]] += p["sg_A20"] * arrest * g("A20_mRNA")
    d[IX["A20"]] -= p["sg_A20"] * g("A20")

    # ---------------- IRF3 module ----------------
    flux(p["p_Tbk1_by_RigiMavs"] * g("RIGI_PolyIC_MAVS") * g("TBK1_0"),
         "TBK1_0", produce=("TBK1_p",))
    flux(p["q_Tbk1"] * g("TBK1_p"), "TBK1_p", produce=("TBK1_0",))
    flux(p["q_Tbk1_by_A20"] * g("A20") * g("TBK1_p"), "TBK1_p", produce=("TBK1_0",))
    flux(p["p_Irf3_by_Tbk1"] * g("TBK1_p") * g("IRF3_0"), "IRF3_0", produce=("IRF3_p",))
    flux(p["q_Irf3"] * g("IRF3_p"), "IRF3_p", produce=("IRF3_0",))

    # ---------------- IFN-beta module ----------------
    flux(p["b_Ifnar_Ifnb_cyt"] * g("IFNAR") * g("IFNb_cyt"),
         "IFNAR", "IFNb_cyt", produce=("IFNAR_IFNb",))
    flux(p["b_Ifnar_Ifnb_ext"] * g("IFNAR") * g("IFNb_ext"),
         "IFNAR", "IFNb_ext", produce=("IFNAR_IFNb",))
    d[IX["IFNAR_mRNA"]] += p["tg_Ifnar_mrna"]
    d[IX["IFNAR_mRNA"]] -= (p["tg_Ifnar_mrna"]
                            * (p["m_Rnasel_Ifnar_mrna"] + o["RNaseL_a"])
                            / p["m_Rnasel_Ifnar_mrna"] * g("IFNAR_mRNA"))
    d[IX["IFNAR"]] += p["s_Ifnar"] * arrest * g("IFNAR_mRNA")
    d[IX["IFNAR"]] -= p["g_Ifnar"] * g("IFNAR")
    d[IX["IFNAR_IFNb"]] -= p["g_Ifnar_w_Ifnb"] * g("IFNAR_IFNb")
    d[IX["IFNb_mRNA"]] += (p["tg_Ifnb_mrna"] * nf * o["IRF3_p"]
                           / (p["m_Ifnb_mrna_NfkbIrf3"] + nf * o["IRF3_p"]))
    d[IX["IFNb_mRNA"]] -= p["tg_Ifnb_mrna"] * g("IFNb_mRNA")
    d[IX["IFNb_cyt"]] += p["sg_Ifnb"] * g("IFNb_mRNA")
    d[IX["IFNb_cyt"]] -= p["sg_Ifnb"] * g("IFNb_cyt")

    # ---------------- STAT1/2 module ----------------
    kph = p["p_Stat"] * o["IFNAR_a"] * p["m_Ifnar_a"]
    flux(kph / (p["m_Ifnar_a"] + o["STAT1_u"]) * g("STAT1_0"), "STAT1_0",
         produce=("STAT1_p",))
    flux(p["q_Stat"] * g("STAT1_p"), "STAT1_p", produce=("STAT1_0",))
    flux(kph / (p["m_Ifnar_a"] + o["STAT2_u"]) * g("STAT2_0"), "STAT2_0",
         produce=("STAT2_p",))
    flux(p["q_Stat"] * g("STAT2_p"), "STAT2_p", produce=("STAT2_0",))
    flux(p["b_Stat1_Stat2"] * g("STAT1_p") * g("STAT2_p"),
         "STAT1_p", "STAT2_p", produce=("STAT12",))
    flux(p["qu_Stat1_Stat2"] * g("STAT12"), "STAT12", produce=("STAT1_0", "STAT2_0"))

    tgs = p["tg_Stat_mrna"]
    d[IX["STAT1_mRNA"]] += tgs * (p["ma_Stat1_gene_basal"] + dim) / (p["ma_Stat1_gene_basal"] + dim + 1.0)
    d[IX["STAT2_mRNA"]] += tgs * (p["ma_Stat2_gene_basal"] + dim) / (p["ma_Stat2_gene_basal"] + dim + 1.0)
    sdeg = tgs * (p["m_Rnasel"] + o["RNaseL_a"]) / p["m_Rnasel"]
    d[IX["STAT1_mRNA"]] -= sdeg * g("STAT1_mRNA")
    d[IX["STAT2_mRNA"]] -= sdeg * g("STAT2_mRNA")
    d[IX["STAT1_0"]] += p["sg_Stat"] * g("STAT1_mRNA")
    d[IX["STAT2_0"]] += p["sg_Stat"] * g("STAT2_mRNA")
    for s in ("STAT1_0", "STAT1_p", "STAT2_0", "STAT2_p", "STAT12"):
        d[IX[s]] -= p["sg_Stat"] * g(s)
    return d


def seed(p):
    y = np.zeros(len(SPECIES))
    y[IX["MAVS"]] = p["h_Mavs"]
    y[IX["eIF2a_0"]] = 1.0
    y[IX["C_cyt0"]] = 1.0
    y[IX["TAK1_i"]] = 1.0
    y[IX["IKK_n"]] = 1.0
    y[IX["TBK1_0"]] = 1.0
    y[IX["IRF3_0"]] = 1.0
    return y


def integrate(phases, times, overrides=None, equilibrate=2592000.0,
              rtol=1e-11, atol=1e-16):
    """Integrate one protocol and sample it on BioNetGen's own output grid.

    Parameters
    ----------
    phases : list of (t_start, t_end, sets)
        The recorded phases of the protocol, in seconds on the recorded clock.
        `sets` maps species names to amounts applied at `t_start`, mirroring the
        `setConcentration` calls in the BNGL actions block.
    times : array
        The BioNetGen output times to sample, covering all the phases.
    overrides : dict, optional
        Parameter overrides, e.g. ``{"h_A20_gene": 0.0}`` for the A20 KO cell line.
    equilibrate : float
        Duration of the unrecorded relaxation to the resting state that every
        protocol begins with (30 days, as in the BNGL files).

    Returns
    -------
    (t, Y) : the sampled times and the corresponding 53-species state matrix.
    """
    p = derived({**P, **(overrides or {})})
    y = seed(p)
    y = solve_ivp(rhs, (0.0, equilibrate), y, args=(p,), method="LSODA",
                  rtol=rtol, atol=atol).y[:, -1]
    times = np.asarray(times, float)
    ts, ys = [], []
    for i, (t0, t1, sets) in enumerate(phases):
        for name, value in sets.items():
            y[IX[name]] = value
        lo = times > t0 + 1e-9 if i else times >= t0 - 1e-9
        grid = times[lo & (times <= t1 + 1e-9)]
        sol = solve_ivp(rhs, (0.0, t1 - t0), y, args=(p,), method="LSODA",
                        t_eval=grid - t0, rtol=rtol, atol=atol)
        y = sol.y[:, -1]
        ts.append(grid)
        ys.append(sol.y.T)
    return np.concatenate(ts), np.concatenate(ys, axis=0)


def observable_frame(Y):
    """Evaluate every BNGL observable over a state matrix; returns a DataFrame."""
    import pandas as pd
    recs = [observables(y) for y in Y]
    names = [k for k in recs[0] if not k.startswith("_")]
    return pd.DataFrame({n: [r[n] for r in recs] for n in names})
