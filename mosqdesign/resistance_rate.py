"""Per-candidate resistance-rate (e) prediction for dsx homing-drive guides.

Honesty contract
----------------
This module does NOT claim an absolute laboratory resistance rate. It
computes, per guide, three fully derived quantities and combines them into
an interval-valued e estimate with every assumption stated:

  1. r0  - standing compromised-allele fraction in the wild population
           (from results/ag1000g/allele_freq.json when present; the
           Ag1000G phase-3 all-sites VCF scan). This is DATA.
  2. L   - mutational target size: the number of single-nucleotide
           substitutions inside the 23-nt target window that abolish or
           strongly reduce cleavage (PAM G2/G3 positions x 3 alternative
           bases; all 8 PAM-proximal seed positions x 3, counting every
           seed mismatch as potentially resistance-causing - an upper
           bound, stated as such). Pure combinatorics.
  3. MH  - microhomology inventory around the predicted cut site: every
           deletion of length 1..MH_MAX_DEL with >=MH_MIN homologous bases
           at its breakpoints, enumerated from the AgamP4 reference
           sequence, classified resistance-generating when the deleted
           interval removes >=1 PAM-GG base or >=1 seed base. Counts and
           length histogram are reported; no claim to reproduce Bae et al.
           2014's exact scoring formula is made.

The de-novo per-generation resistance rate is parameterized as

      e(g) = p_EJ * w(g),   w(g) = n_MH_res(g) / n_MH_total(g)

where p_EJ (probability that a failed homing event is repaired as a
resistant allele) is swept over a documented grid because published
estimates vary by locus and construct:

  * Kyrou et al. 2018 (PMC6871539, verified): the dsx Ag(QFS)1 drive
    showed no resistant-allele selection in caged populations - the
    observed-zero anchor at the dsx splice-junction target.
  * Hammond et al. 2016 (PMC4913862, verified abstract): transmission
    rates 91.4-99.6% across three An. gambiae fertility-gene drives,
    i.e. failed-homing fractions 0.4-8.6%, within which end-joining
    alleles were documented. The fraction of failed events becoming
    RESISTANT alleles is construct-dependent and is NOT read from the
    abstract; it is swept, not asserted.

Outputs the interval e in [p_EJ_min*w, p_EJ_max*w] plus drive-dynamic
consequences (generations to 95% drive frequency, resistance-fixation
generation) under each grid point via mosqdesign.drive_sim.

References: Burt 2003 Proc R Soc B; Unckless et al. 2017 PNAS;
Kyrou et al. 2018 Nat Biotechnol (PMC6871539); Hammond et al. 2016
Nat Biotechnol (PMC4913862); Bae et al. 2014 NAR (microhomology concept).
"""
import json
import math
import os

from mosqdesign.drive_sim import simulate

COMP = str.maketrans("ACGT", "TGCA")
MH_MIN = 2          # minimum breakpoint homology (bp)
MH_MAX_DEL = 40     # maximum deletion length scanned (bp)
FLANK = 60          # sequence context each side of the cut (bp)
SEED_N = 8          # PAM-proximal seed length (nt)
P_EJ_GRID = [1e-4, 1e-3, 1e-2, 5e-2, 0.2]   # documented sweep, not a fit

# 1-based inclusive AgamP4 windows (verified twice; see
# scripts/ag1000g_allele_freq.py GUIDES).
GUIDE_WINDOWS = {
    "dsx-v3-1": {"proto": (48711501, 48711520), "pam": (48711498, 48711500)},
    "dsx-v3-2": {"proto": (48712765, 48712784), "pam": (48712762, 48712764)},
    "kyrou":    {"proto": (48714640, 48714659), "pam": (48714637, 48714639)},
}

REGION_OFFSET = 48700000  # AgamP4 2R coordinate of FASTA base 1


def load_region(path):
    seq = "".join(l.strip() for l in open(path) if not l.startswith(">"))
    return seq.upper()


def cut_site_plus(gname):
    """Plus-strand cut position for a minus-strand NGG target: Cas9 cuts
    3 nt into the protospacer from the PAM-proximal end. PAM is upstream
    (lower coords) for these minus-strand guides, so the cut falls between
    proto_start+2 and proto_start+3; we report the left base."""
    pa, _ = GUIDE_WINDOWS[gname]["proto"]
    return pa + 2


def mutational_target_size(gname):
    """L: upper-bound count of cleavage-abolishing single substitutions."""
    # PAM G2/G3 (2 positions x 3 alt bases) + seed (SEED_N x 3)
    return 2 * 3 + SEED_N * 3


def microhomology_deletions(seq, cut, pam_pos, seed_pos,
                            mh_min=MH_MIN, max_del=MH_MAX_DEL, flank=FLANK):
    """Enumerate breakpoint-homology deletions around the cut site.

    A deletion [s, s+l) (0-based, seq coordinates) is microhomology-
    mediated if the m bases immediately 5' of the left breakpoint equal
    the m bases immediately 5' of the right breakpoint (the classic MMEJ
    signature: one copy is deleted with the interval). Resistance-
    generating if the deleted interval contains a PAM-GG or seed base.

    pam_pos, seed_pos: sets of 0-based seq indices to test overlap.
    Returns list of dicts.
    """
    n = len(seq)
    lo = max(mh_min, cut - flank)
    hi = min(n - max_del - mh_min, cut + flank)
    out = []
    for s in range(lo, hi + 1):
        for l in range(1, max_del + 1):
            e = s + l
            if e + mh_min > n:
                break
            m = 0
            while (m < l and s - 1 - m >= 0 and e - 1 - m >= 0
                   and seq[s - 1 - m] == seq[e - 1 - m]):
                m += 1
            if m < mh_min:
                continue
            deleted = set(range(s, e))
            res = bool(deleted & pam_pos) or bool(deleted & seed_pos)
            out.append({"start": s, "del_len": l, "mh_len": m,
                        "resistance_generating": res})
    return out


def window_pos_sets(gname):
    """0-based seq-index sets for PAM-GG and seed positions of a guide."""
    w = GUIDE_WINDOWS[gname]
    pa, pb = w["proto"]
    sa, sb = w["pam"]
    pam_gg = {sa - REGION_OFFSET, (sa + 1) - REGION_OFFSET}  # G3, G2
    seed = set(range(pa - REGION_OFFSET,
                     pa + SEED_N - REGION_OFFSET))
    return pam_gg, seed


def predict(gname, seq, af_sites=None):
    """Full per-guide prediction record."""
    pam_pos, seed_pos = window_pos_sets(gname)
    cut = cut_site_plus(gname) - REGION_OFFSET
    dels = microhomology_deletions(seq, cut, pam_pos, seed_pos)
    n_tot = len(dels)
    n_res = sum(1 for d in dels if d["resistance_generating"])
    w = n_res / n_tot if n_tot else 0.0
    L = mutational_target_size(gname)
    r0 = None
    if af_sites:
        # standing compromised AF = sum AF of alleles flagged compromised
        r0 = af_sites.get(gname)
    e_interval = [P_EJ_GRID[0] * w, P_EJ_GRID[-1] * w]
    e_grid = {f"{p:g}": p * w for p in P_EJ_GRID}
    return {
        "guide": gname,
        "cut_site_plus_1based": cut + REGION_OFFSET,
        "mutational_target_size_L": L,
        "mh_deletions_total": n_tot,
        "mh_deletions_resistance_generating": n_res,
        "mh_resistance_fraction_w": w,
        "mh_del_len_hist_res": _hist([d["del_len"] for d in dels
                                      if d["resistance_generating"]]),
        "standing_compromised_AF_r0": r0,
        "e_interval": e_interval,
        "e_grid_pEJ": e_grid,
    }


def _hist(lengths):
    h = {}
    for l in lengths:
        h[l] = h.get(l, 0) + 1
    return dict(sorted(h.items()))


def drive_consequences(pred, h=0.99, c_hom=1.0, generations=200):
    """Drive outcomes at the predicted e interval endpoints.

    Returns generations to 95% drive frequency and to 50% resistance
    (None if not reached) for each endpoint of the e interval, with
    r0 taken from standing variation when available.
    """
    r0 = pred["standing_compromised_AF_r0"] or 0.0
    out = {}
    for tag, e in zip(("e_low", "e_high"), pred["e_interval"]):
        w0 = 0.99 - r0
        d0 = 0.01
        traj = simulate(w0, d0, r0, h, e, c_hom=c_hom,
                        generations=generations)
        g95 = next((i for i, (w, d, r) in enumerate(traj)
                    if d >= 0.95), None)
        r50 = next((i for i, (w, d, r) in enumerate(traj)
                    if r >= 0.5), None)
        out[tag] = {"e": e, "gen_to_95pct_drive": g95,
                    "gen_to_50pct_resistance": r50}
    return out


def run_all(fasta_path, af_json=None, out_json="results/resistance_rates.json",
            out_md="results/resistance_rates.md"):
    seq = load_region(fasta_path)
    af_sites = {}
    if af_json and os.path.exists(af_json):
        with open(af_json) as fh:
            af = json.load(fh)
        for g, s in af.get("guide_summary", {}).items():
            af_sites[g] = s["compromised_allele_fraction"]
    preds = {g: predict(g, seq, af_sites) for g in GUIDE_WINDOWS}
    for g in preds:
        preds[g]["drive_consequences"] = drive_consequences(preds[g])
    with open(out_json, "w") as fh:
        json.dump({"model": "e = p_EJ * w; w = resistance-generating "
                            "microhomology fraction; L = upper-bound "
                            "mutational target size",
                   "p_EJ_grid": P_EJ_GRID, "guides": preds}, fh, indent=2)
    lines = [
        "## Per-candidate resistance-rate (e) prediction", "",
        "Model: e = p_EJ x w, w = fraction of cut-site microhomology "
        "deletions that remove PAM-GG or seed sequence. p_EJ swept "
        "(not fitted) over " +
        ", ".join(f"{p:g}" for p in P_EJ_GRID) +
        ". L = upper-bound single-substitution resistance target size "
        "(PAM GG 2x3 + seed 8x3). r0 = standing compromised-allele "
        "fraction (Ag1000G).", "",
        "| guide | L | MH dels (res/total) | w | r0 | e interval | "
        "gen to 95% drive (e_low/e_high) | gen to 50% R (e_low/e_high) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for g, p in preds.items():
        dc = p["drive_consequences"]
        def fmt(x):
            return str(x) if x is not None else ">200"
        r0 = p["standing_compromised_AF_r0"]
        lines.append(
            f"| {g} | {p['mutational_target_size_L']} | "
            f"{p['mh_deletions_resistance_generating']}/"
            f"{p['mh_deletions_total']} | {p['mh_resistance_fraction_w']:.3f} "
            f"| {r0 if r0 is not None else 'n/a'} | "
            f"[{p['e_interval'][0]:.2e}, {p['e_interval'][1]:.2e}] | "
            f"{fmt(dc['e_low']['gen_to_95pct_drive'])}/"
            f"{fmt(dc['e_high']['gen_to_95pct_drive'])} | "
            f"{fmt(dc['e_low']['gen_to_50pct_resistance'])}/"
            f"{fmt(dc['e_high']['gen_to_50pct_resistance'])} |")
    with open(out_md, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return preds
