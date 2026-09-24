"""Gene-drive population genetics: deterministic recursion + stochastic Wright-Fisher.

Model class: homing drive with end-joining resistance (Burt 2003 Proc R Soc B;
Unckless et al. 2017 PNAS functional-form recursions). Three alleles at the drive
locus:
  W = wild-type, D = drive, R = resistant (non-functional) allele.
Each generation, in D-heterozygote germlines, the W allele is converted to D with
probability h (homing rate); of the failed conversions, fraction e becomes R.
Fitness: homozygous cost c_hom for D/D (sterility class), s_het for D/W, R neutral.
"""
from __future__ import annotations

import numpy as np


def homing_recurse(w: float, d: float, r: float, h: float, e: float,
                   c_hom: float = 1.0, s_het: float = 0.0) -> tuple[float, float, float]:
    """One discrete generation of drive spread. Allele frequencies sum to 1.

    Genotype fitnesses: WW=1, WD=1-s_het, DD=1-c_hom, WR=1, DR=1-s_het, RR=1.
    Gametes from WD heterozygotes: W->D at rate h (homing); of non-homed W,
    fraction e becomes R (end-joining), rest stay W.
    """
    assert h >= 0 and e >= 0
    p = {"WW": w * w, "WD": 2 * w * d, "DD": d * d, "WR": 2 * w * r, "DR": 2 * d * r, "RR": r * r}
    fit = {"WW": 1.0, "WD": 1.0 - s_het, "DD": 1.0 - c_hom, "WR": 1.0, "DR": 1.0 - s_het, "RR": 1.0}
    mean_fit = sum(p[g] * fit[g] for g in p)
    # allele contributions after selection + homing in WD/DR germlines
    def gametes(g: str) -> tuple[float, float, float]:
        if g == "WW": return 1.0, 0.0, 0.0
        if g == "DD": return 0.0, 1.0, 0.0
        if g == "RR": return 0.0, 0.0, 1.0
        if g == "WR": return 0.5, 0.0, 0.5
        if g == "WD":
            # drive-bearing homolog converts W with prob h; failures: e -> R, else W
            return 0.5 * (1 - h) * (1 - e), 0.5 * (1 + h), 0.5 * (1 - h) * e
        if g == "DR":
            return 0.0, 0.5, 0.5
        raise ValueError(g)
    gw = sum(p[g] * fit[g] * gametes(g)[0] for g in p) / mean_fit
    gd = sum(p[g] * fit[g] * gametes(g)[1] for g in p) / mean_fit
    gr = sum(p[g] * fit[g] * gametes(g)[2] for g in p) / mean_fit
    return gw, gd, gr


def simulate(w0: float, d0: float, r0: float, h: float, e: float,
             c_hom: float = 1.0, s_het: float = 0.0, generations: int = 40) -> np.ndarray:
    """Deterministic trajectory; returns (generations+1, 3) array of (W, D, R)."""
    traj = [(w0, d0, r0)]
    for _ in range(generations):
        traj.append(homing_recurse(*traj[-1], h, e, c_hom, s_het))
    return np.array(traj)


def wright_fisher(n: int, w0: float, d0: float, r0: float, h: float, e: float,
                  c_hom: float = 1.0, s_het: float = 0.0, generations: int = 40,
                  seed: int = 0) -> np.ndarray:
    """Stochastic Wright-Fisher with the same expected recursions (multinomial draws)."""
    rng = np.random.default_rng(seed)
    freqs = np.array([w0, d0, r0])
    traj = [freqs.copy()]
    for _ in range(generations):
        expected = np.array(homing_recurse(*freqs, h, e, c_hom, s_het))
        expected = np.clip(expected, 0.0, 1.0)
        expected /= expected.sum()
        freqs = rng.multinomial(n, expected) / n
        traj.append(freqs.copy())
    return np.array(traj)
