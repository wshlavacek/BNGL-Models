#!/usr/bin/env python3
"""Independent isotopomer balance equations for the Calvin cycle of Mu et al. (2010).

BioNetGen turns the BNGL carbon fate maps of Sec. 15.4 into isotopomer balance
equations by rule-based graph rewriting: it matches molecule patterns, applies the
component-label mapping, and enumerates the reachable species. This module reaches
the same equations by the other route the chapter contrasts with it (Sec. 15.2.1),
enumerating every isotopomer of every metabolite up front and mapping labeling
patterns with integer bit arithmetic. Nothing is shared with BioNetGen except the
published network and its fate maps, so agreement between the two is a real check on
the curated BNGL file.

Representation
--------------
A metabolite with n carbons has 2**n isotopomers. An isotopomer is an integer whose
bit n-1-i is the labeling state of carbon i, so carbon C1 is the most significant bit
and the ordering matches the (C1 C2) convention of Eq. 15.14. A fate map is given as,
for each product molecule, the list of (reactant slot, carbon index) that each of its
carbons comes from; the map of a reversible reaction's backward direction is derived
by inverting the forward map rather than transcribed a second time.

Expanding one network reaction over all combinations of reactant isotopomers gives the
elementary isotopomer reactions. Each carries rate constant flux/(product of the pool
sizes of its reactants), so the total flux through the reaction is the specified flux
whenever the pools are at their specified sizes, and the balance equations are mass
action on the isotopomer amounts.

The expansion produces 4042 elementary reactions over 562 isotopomers, which are the
counts BioNetGen reports for the same model and the 562 isotopomers quoted in Sec. 15.4.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

# Carbon count of each metabolite, in the InChI atom ordering used by Mu et al. (2007).
CARBONS = {
    "CO2": 1,
    "G3P": 3, "P13G": 3, "T3P": 3, "DHAP": 3,
    "E4P": 4,
    "R5P": 5, "R5DP": 5, "Ri5P": 5, "X5P": 5,
    "F16P": 6, "F6P": 6,
    "S17P": 7, "S7P": 7,
}

# Metabolite pool sizes; Mu et al. (2010) use the default pool size of one.
POOLS = {name: 1.0 for name in CARBONS}

# Published fluxes of Sec. 15.4. `v_R01529_b` and `vs_R5P` are printed as 0.82 and
# 0.05; both are lowered by one unit in their last digit here, as in the primary
# BNGL file, so that every steady-state balance closes exactly. Set
# `PUBLISHED_FLUXES` instead to reproduce the washout of the variant BNGL file.
FLUXES = {
    "v_R00024": 1.00,
    "v_R01512": 1.95,
    "v_R01063_f": 0.05, "v_R01063_b": 2.00,
    "v_R01015_f": 0.86, "v_R01015_b": 0.10,
    "v_R01068_f": 0.15, "v_R01068_b": 0.58,
    "v_R00762": 0.43,
    "v_R01067_f": 0.50, "v_R01067_b": 0.12,
    "v_R01829_f": 0.20, "v_R01829_b": 0.53,
    "v_R01845": 0.33,
    "v_R01641_f": 0.43, "v_R01641_b": 0.10,
    "v_R01529_f": 0.10, "v_R01529_b": 0.81,
    "v_R01056_f": 0.50, "v_R01056_b": 0.17,
    "v_R01523": 1.00,
    "vs_T3P": 0.05, "vs_F6P": 0.05, "vs_R5P": 0.04,
    "vs_E4P": 0.05, "vs_G3P": 0.05,
}
PUBLISHED_FLUXES = {**FLUXES, "v_R01529_b": 0.82, "vs_R5P": 0.05}

# Carbon atom fate maps, forward direction only. Each entry lists the reactants, then
# for every product molecule the source (reactant slot, carbon index) of each of its
# carbons. `reversible` marks the reactions for which a `_b` flux is also specified;
# their backward fate map is the inverse of the map given here.
FATE_MAPS = [
    # Ribulose-bisphosphate carboxylase. Three carbons of R5DP plus the CO2 carbon
    # make one G3P; the remaining two R5DP carbons make the other.
    ("R00024", ["R5DP", "CO2"],
     [("G3P", [(0, 1), (0, 3), (1, 0)]), ("G3P", [(0, 0), (0, 2), (0, 4)])], False),
    # Phosphoglycerate kinase.
    ("R01512", ["G3P"], [("P13G", [(0, 0), (0, 1), (0, 2)])], False),
    # Glyceraldehyde-3-phosphate dehydrogenase.
    ("R01063", ["T3P"], [("P13G", [(0, 1), (0, 2), (0, 0)])], True),
    # Triose-phosphate isomerase.
    ("R01015", ["T3P"], [("DHAP", [(0, 0), (0, 1), (0, 2)])], True),
    # Fructose-bisphosphate aldolase.
    ("R01068", ["F16P"],
     [("DHAP", [(0, 4), (0, 1), (0, 5)]), ("T3P", [(0, 3), (0, 0), (0, 2)])], True),
    # Fructose-1,6-bisphosphatase.
    ("R00762", ["F16P"],
     [("F6P", [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5)])], False),
    # Transketolase, F6P + T3P.
    ("R01067", ["F6P", "T3P"],
     [("E4P", [(0, 4), (0, 0), (0, 3), (0, 2)]),
      ("X5P", [(0, 1), (1, 1), (0, 5), (1, 2), (1, 0)])], True),
    # Sedoheptulose-bisphosphate aldolase.
    ("R01829", ["S17P"],
     [("DHAP", [(0, 5), (0, 1), (0, 3)]),
      ("E4P", [(0, 6), (0, 0), (0, 4), (0, 2)])], True),
    # Sedoheptulose-bisphosphatase.
    ("R01845", ["S17P"],
     [("S7P", [(0, 1), (0, 0), (0, 3), (0, 2), (0, 5), (0, 4), (0, 6)])], False),
    # Transketolase, S7P + T3P.
    ("R01641", ["S7P", "T3P"],
     [("Ri5P", [(0, 1), (0, 3), (0, 5), (0, 6), (0, 4)]),
      ("X5P", [(0, 0), (1, 1), (0, 2), (1, 2), (1, 0)])], True),
    # Ribulose-phosphate 3-epimerase.
    ("R01529", ["R5P"], [("X5P", [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4)])], True),
    # Ribose-5-phosphate isomerase.
    ("R01056", ["Ri5P"], [("R5P", [(0, 4), (0, 0), (0, 3), (0, 1), (0, 2)])], True),
    # Phosphoribulokinase.
    ("R01523", ["R5P"], [("R5DP", [(0, 1), (0, 0), (0, 3), (0, 2), (0, 4)])], False),
]

# Pools drained to central metabolism; see Fig. 15.5 in Mu et al. (2010).
SINKS = ["T3P", "F6P", "R5P", "E4P", "G3P"]


def invert_map(reactants, products):
    """Return the fate map of the reverse reaction.

    The forward map says where each product carbon comes from; the reverse map says
    where each reactant carbon comes from, which is the same bijection read backwards.
    """
    inverse: dict[tuple[int, int], tuple[int, int]] = {}
    for p_slot, (_, sources) in enumerate(products):
        for p_carbon, src in enumerate(sources):
            inverse[src] = (p_slot, p_carbon)
    rev_reactants = [name for name, _ in products]
    rev_products = []
    for r_slot, name in enumerate(reactants):
        rev_products.append((name, [inverse[(r_slot, c)] for c in range(CARBONS[name])]))
    n_fwd = sum(CARBONS[name] for name in reactants)
    assert len(inverse) == n_fwd, "fate map is not a bijection on carbons"
    return rev_reactants, rev_products


def stoichiometry() -> tuple[list[str], list[str], np.ndarray]:
    """Return (metabolites, flux names, N) with N @ v the vector of pool balances.

    Fluxes are counted in the direction they are specified, so every entry of v is
    nonnegative and a reversible reaction contributes two columns. `vi_CO2` is the
    total CO2 input, the sum of its unlabeled and labeled parts.
    """
    metabolites = list(CARBONS)
    row = {name: i for i, name in enumerate(metabolites)}
    columns: list[str] = []
    entries: list[np.ndarray] = []

    def column(label, reactants, products):
        col = np.zeros(len(metabolites))
        for name in reactants:
            col[row[name]] -= 1.0
        for name, _ in products:
            col[row[name]] += 1.0
        columns.append(label)
        entries.append(col)

    for name, reactants, products, reversible in FATE_MAPS:
        column(f"v_{name}_f" if reversible else f"v_{name}", reactants, products)
        if reversible:
            rr, rp = invert_map(reactants, products)
            column(f"v_{name}_b", rr, rp)
    for name in SINKS:
        column(f"vs_{name}", [name], [])
    column("vi_CO2", [], [("CO2", None)])
    return metabolites, columns, np.column_stack(entries)


def flux_vector(fluxes: dict, columns: list[str]) -> np.ndarray:
    """Order a flux dictionary to match the columns of `stoichiometry`."""
    return np.array([1.0 if c == "vi_CO2" else fluxes[c] for c in columns])


class Network:
    """Isotopomer network expanded from the fate maps of Fig. 15.5."""

    def __init__(self, fluxes=None, pools=None):
        self.fluxes = dict(FLUXES if fluxes is None else fluxes)
        self.pools = dict(POOLS if pools is None else pools)

        # Global index of every (metabolite, isotopomer) pair.
        self.offset: dict[str, int] = {}
        n = 0
        for name, nc in CARBONS.items():
            self.offset[name] = n
            n += 1 << nc
        self.n_species = n

        reac: list[list[int]] = []
        prod: list[list[int]] = []
        rate: list[float] = []

        def expand(label, reactants, products, flux):
            k = flux / np.prod([self.pools[r] for r in reactants])
            sizes = [1 << CARBONS[r] for r in reactants]
            for combo in np.ndindex(*sizes):
                bits = [
                    [(combo[s] >> (CARBONS[r] - 1 - c)) & 1 for c in range(CARBONS[r])]
                    for s, r in enumerate(reactants)
                ]
                r_idx = [self.offset[r] + combo[s] for s, r in enumerate(reactants)]
                p_idx = []
                for pname, sources in products:
                    v = 0
                    for c, (slot, carbon) in enumerate(sources):
                        v |= bits[slot][carbon] << (CARBONS[pname] - 1 - c)
                    p_idx.append(self.offset[pname] + v)
                reac.append(r_idx)
                prod.append(p_idx)
                rate.append(k)
            self.per_reaction[label] = int(np.prod(sizes))

        self.per_reaction: dict[str, int] = {}
        for name, reactants, products, reversible in FATE_MAPS:
            fwd = self.fluxes.get(f"v_{name}", self.fluxes.get(f"v_{name}_f"))
            expand(name if not reversible else f"{name}_f", reactants, products, fwd)
            if reversible:
                rr, rp = invert_map(reactants, products)
                expand(f"{name}_b", rr, rp, self.fluxes[f"v_{name}_b"])

        # Sinks: one elementary reaction per isotopomer of the drained pool.
        for name in SINKS:
            k = self.fluxes[f"vs_{name}"] / self.pools[name]
            for i in range(1 << CARBONS[name]):
                reac.append([self.offset[name] + i])
                prod.append([])
                rate.append(k)
            self.per_reaction[f"sink_{name}"] = 1 << CARBONS[name]

        # CO2 feed, split into its unlabeled and labeled parts.
        self.co2_unlabeled = len(rate)
        reac.append([])
        prod.append([self.offset["CO2"] + 0])
        rate.append(0.0)
        self.co2_labeled = len(rate)
        reac.append([])
        prod.append([self.offset["CO2"] + 1])
        rate.append(0.0)
        self.per_reaction["CO2_feed"] = 2

        self.reac = reac
        self.prod = prod
        self.rate = np.array(rate)
        self.n_reactions = len(rate)

        # Flat arrays for a vectorized right-hand side.
        self._r1 = np.array([r[0] if len(r) > 0 else -1 for r in reac])
        self._r2 = np.array([r[1] if len(r) > 1 else -1 for r in reac])
        self._p_rxn = np.array([i for i, p in enumerate(prod) for _ in p])
        self._p_sp = np.array([s for p in prod for s in p])
        self._c_rxn = np.array([i for i, r in enumerate(reac) for _ in r])
        self._c_sp = np.array([s for r in reac for s in r])

    def set_feed(self, unlabeled: float, labeled: float) -> None:
        self.rate[self.co2_unlabeled] = unlabeled
        self.rate[self.co2_labeled] = labeled

    def initial_state(self) -> np.ndarray:
        """Every pool at its specified size and fully unlabeled."""
        x = np.zeros(self.n_species)
        for name in CARBONS:
            x[self.offset[name] + 0] = self.pools[name]
        return x

    def rhs(self, _t: float, x: np.ndarray) -> np.ndarray:
        v = self.rate.copy()
        m1 = self._r1 >= 0
        v[m1] *= x[self._r1[m1]]
        m2 = self._r2 >= 0
        v[m2] *= x[self._r2[m2]]
        dx = np.zeros_like(x)
        np.add.at(dx, self._c_sp, -v[self._c_rxn])
        np.add.at(dx, self._p_sp, v[self._p_rxn])
        return dx

    def totals(self, x: np.ndarray) -> dict[str, float]:
        """Pool size of each metabolite, summed over its isotopomers."""
        return {
            name: float(np.sum(x[self.offset[name]:self.offset[name] + (1 << nc)]))
            for name, nc in CARBONS.items()
        }

    def mass_isotopomers(self, name: str, x: np.ndarray) -> np.ndarray:
        """Amounts grouped by number of 13C atoms, i.e. the mass isotopomer spectrum."""
        nc = CARBONS[name]
        block = x[self.offset[name]:self.offset[name] + (1 << nc)]
        weight = np.array([bin(i).count("1") for i in range(1 << nc)])
        return np.array([block[weight == m].sum() for m in range(nc + 1)])

    def labeled_carbons(self, name: str, x: np.ndarray) -> float:
        """Total number of 13C atoms carried by a metabolite pool."""
        nc = CARBONS[name]
        block = x[self.offset[name]:self.offset[name] + (1 << nc)]
        weight = np.array([bin(i).count("1") for i in range(1 << nc)])
        return float(block @ weight)

    def run(self, t_eval: np.ndarray, x0: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        sol = solve_ivp(
            self.rhs, (float(t_eval[0]), float(t_eval[-1])), x0, t_eval=t_eval,
            method="LSODA", rtol=1e-10, atol=1e-12,
        )
        assert sol.success, sol.message
        return sol.t, sol.y


def binomial_spectrum(n: int, p: float) -> np.ndarray:
    """Mass isotopomer spectrum of n independently labeled carbons at enrichment p.

    At isotopic steady state the labeling of every carbon in the network is an
    independent Bernoulli draw at the enrichment of the CO2 feed, because CO2 is the
    only carbon input and every fate map is a bijection on carbons; the product
    distribution is therefore a fixed point of the balance equations, and it is the
    unique one.
    """
    from math import comb
    return np.array([comb(n, m) * p**m * (1 - p) ** (n - m) for m in range(n + 1)])


if __name__ == "__main__":
    net = Network()
    print(f"{net.n_species} isotopomers, {net.n_reactions} elementary reactions")

    net.set_feed(1.0, 0.0)
    x0 = net.initial_state()
    _, y = net.run(np.array([0.0, 10000.0]), x0)
    drift = max(abs(v - 1.0) for v in net.totals(y[:, -1]).values())
    print(f"equilibration on unlabeled feed: max pool drift = {drift:.3e}")

    net.set_feed(0.9, 0.1)
    t, y = net.run(np.linspace(0.0, 200.0, 201), y[:, -1])
    spectrum = net.mass_isotopomers("T3P", y[:, -1]) / net.pools["T3P"]
    print("T3P mass isotopomer fractions at t = 200: "
          + ", ".join(f"{v:.6f}" for v in spectrum))
    print("binomial(3, 0.1) limit:               "
          + ", ".join(f"{v:.6f}" for v in binomial_spectrum(3, 0.1)))
