"""Independent bcftools recheck of the Ag1000G dsx-window scan (tool 42).

The scan of record (scripts/ag1000g_scan_all.py) used pysam.TabixFile: raw
tabix text slicing, Python-side column splitting, manual GT/GQ/DP indexing.
This recheck uses bcftools 1.21 (htslib full VCF parse, header validation,
record re-emission) on a stratified subsample of the same 4,693 per-sample
VCFs and compares row-for-row inside the six target windows
(3 guides x protospacer+PAM).

Honesty note: both engines link htslib at the transport layer (bgzf/tabix/
libcurl). The independence is at the parse layer: pysam path did text slicing
with Python splits; bcftools does a full VCF parse + re-emit. A column-index
or region-bounds bug in either path shows up as discordance here.

Arms (seed 42, recorded in results):
  carrier: N_C random samples with >=1 non-ref call inside the six windows
  zero:    N_Z random samples with zero target-window calls and no error

Usage:
  python3 scripts/bcftools_afcheck.py --scan     # fetch via bcftools (checkpointed)
  python3 scripts/bcftools_afcheck.py --compare  # offline comparison -> results/bcftools_afcheck.{json,md}
"""
import json, os, random, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

BASE = "https://vo_agam_output.cog.sanger.ac.uk/"
REGION = "2R:48711450-48714700"
FETCH = (48711450, 48714700)
TARGETS = {
    "dsx-v3-1_proto": (48711501, 48711520),
    "dsx-v3-1_PAM":   (48711498, 48711500),
    "dsx-v3-2_proto": (48712765, 48712784),
    "dsx-v3-2_PAM":   (48712762, 48712764),
    "kyrou_proto":    (48714640, 48714659),
    "kyrou_PAM":      (48714637, 48714639),
}
BCFTOOLS = os.environ.get("BCFTOOLS_BIN", "/tmp/bcftools-1.21/bcftools")
SCAN_IN = "results/ag1000g/on_target.jsonl"
RECHECK_OUT = "results/ag1000g/bcftools_recheck.jsonl"
RESULT_JSON = "results/bcftools_afcheck.json"
RESULT_MD = "results/bcftools_afcheck.md"
N_CARRIER, N_ZERO, SEED, WORKERS = 400, 300, 42, 8

def in_targets(pos):
    return [n for n, (a, b) in TARGETS.items() if a <= pos <= b]

def select_subsample(path=SCAN_IN, n_carrier=N_CARRIER, n_zero=N_ZERO, seed=SEED):
    carriers, zeros = [], []
    for line in open(path):
        r = json.loads(line)
        if "error" in r:
            continue
        th = [h for h in r.get("target_hits", r.get("hits", [])) if h.get("targets")]
        (carriers if th else zeros).append(r["sample"])
    rng = random.Random(seed)
    rng.shuffle(carriers); rng.shuffle(zeros)
    return carriers[:n_carrier], zeros[:n_zero]

def fetch_bcftools(sid):
    """Extract rows in REGION via bcftools view; parse like the original scan."""
    t0 = time.time()
    env = dict(os.environ, LD_LIBRARY_PATH="/tmp/htslib-install/lib")
    url = BASE + sid + ".vcf.gz"
    p = None
    for _attempt in range(3):
        try:
            p = subprocess.run([BCFTOOLS, "view", "-r", REGION, url],
                               capture_output=True, text=True, timeout=180, env=env)
            break
        except subprocess.TimeoutExpired:
            p = None
    if p is None:
        raise RuntimeError("bcftools: timed out 3x (180s each)")
    if p.returncode != 0:
        raise RuntimeError("bcftools: " + p.stderr.strip()[:150])
    hits, n_rows = [], 0
    for line in p.stdout.splitlines():
        if line.startswith("#"):
            continue
        n_rows += 1
        f = line.split("\t")
        gt_field = f[9].split(":"); fmt = f[8].split(":")
        gt = gt_field[0]
        if gt in ("0/0", "./.", "0|0", ".|."):
            continue
        try:
            gq = int(gt_field[fmt.index("GQ")]); dp = int(gt_field[fmt.index("DP")])
        except (ValueError, IndexError):
            gq, dp = None, None
        pos = int(f[1])
        hits.append({"pos": pos, "ref": f[3], "alt": f[4], "gt": gt,
                     "gq": gq, "dp": dp, "targets": in_targets(pos)})
    return {"sample": sid, "n_rows": n_rows, "n_nonref": len(hits),
            "target_hits": hits, "seconds": round(time.time() - t0, 2)}

def run_scan():
    carriers, zeros = select_subsample()
    todo = set(carriers) | set(zeros)
    done = set()
    if os.path.exists(RECHECK_OUT):
        for line in open(RECHECK_OUT):
            done.add(json.loads(line)["sample"])
    todo = [s for s in (carriers + zeros) if s not in done]
    print(f"subsample: {len(carriers)} carriers + {len(zeros)} zeros; {len(todo)} to fetch", flush=True)
    def _safe(sid):
        try:
            return fetch_bcftools(sid)
        except Exception as e:
            return {"sample": sid, "error": str(e)[:200]}
    with open(RECHECK_OUT, "a") as fh, ThreadPoolExecutor(max_workers=WORKERS) as ex:
        count = 0
        for rec in ex.map(_safe, todo, chunksize=1):
            fh.write(json.dumps(rec) + "\n"); fh.flush()
            count += 1
            if count % 25 == 0:
                print(f"fetched {count}/{len(todo)}", flush=True)

def _key(h):
    return (h["pos"], h["ref"], h["alt"], h["gt"].replace("|", "/"))

def _hits(rec):
    # scan of record has two schemas: later rows use target_hits+n_rows,
    # earlier rows use hits (no n_rows). Same hit dict shape otherwise.
    return rec.get("target_hits", rec.get("hits", []))

def compare_records(pysam_rec, bcf_rec):
    """Row-level comparison inside the six target windows."""
    pa = [h for h in _hits(pysam_rec) if h.get("targets")]
    pb = [h for h in bcf_rec.get("target_hits", []) if h.get("targets")]
    sa, sb = set(map(_key, pa)), set(map(_key, pb))
    return {
        "n_pysam": len(pa), "n_bcftools": len(pb),
        "only_pysam": sorted(sa - sb), "only_bcftools": sorted(sb - sa),
        "concordant": sa == sb,
        "n_rows_pysam": pysam_rec.get("n_rows"), "n_rows_bcftools": bcf_rec.get("n_rows"),
    }

def run_compare():
    pysam_recs = {}
    for line in open(SCAN_IN):
        r = json.loads(line)
        pysam_recs[r["sample"]] = r
    rows = []
    for line in open(RECHECK_OUT):
        r = json.loads(line)
        if "error" in r:
            rows.append({"sample": r["sample"], "error": r["error"]}); continue
        p = pysam_recs.get(r["sample"])
        if p is None or "error" in p:
            rows.append({"sample": r["sample"], "error": "no pysam record"}); continue
        c = compare_records(p, r); c["sample"] = r["sample"]
        rows.append(c)
    n_ok = sum(1 for r in rows if r.get("concordant"))
    rowcount_cmp = [r for r in rows if not r.get("error") and r["n_rows_pysam"] is not None]
    n_rows_ok = sum(1 for r in rowcount_cmp if r["n_rows_pysam"] == r["n_rows_bcftools"])
    n_rows_offby1 = sum(1 for r in rowcount_cmp
                        if r["n_rows_bcftools"] == r["n_rows_pysam"] + 1)
    n_err = sum(1 for r in rows if r.get("error"))
    discordant = [r for r in rows if not r.get("error") and not r["concordant"]]
    # Vintage split: legacy rows (ag1000g_on_target.py, GQ>=20/DP>=5 filter,
    # "hits" key, no n_rows) vs scan-of-record rows (ag1000g_scan_all.py,
    # unfiltered, "target_hits" + n_rows). Legacy discordance is expected
    # wherever bcftools recovers calls the legacy filter dropped.
    pysam_recs_all = {}
    for line in open(SCAN_IN):
        r = json.loads(line)
        pysam_recs_all[r["sample"]] = r
    new_rows = [r for r in rows if not r.get("error")
                and "n_rows" in pysam_recs_all.get(r["sample"], {})]
    legacy_rows = [r for r in rows if not r.get("error")
                   and "n_rows" not in pysam_recs_all.get(r["sample"], {})]
    n_new_ok = sum(1 for r in new_rows if r["concordant"])
    bcf_by_sample = {}
    for line in open(RECHECK_OUT):
        r = json.loads(line)
        if "error" not in r:
            bcf_by_sample[r["sample"]] = r
    def legacy_explained(r):
        hits = {(h["pos"], h["ref"], h["alt"], h["gt"].replace("|", "/")): h
                for h in bcf_by_sample[r["sample"]]["target_hits"]}
        for k in map(tuple, r["only_bcftools"]):
            h = hits.get(k)
            if h is None:
                return False
            gq, dp = h.get("gq"), h.get("dp")
            if gq is not None and dp is not None and gq >= 20 and dp >= 5:
                return False
        return not r["only_pysam"]
    legacy_disc = [r for r in legacy_rows if not r["concordant"]]
    n_legacy_explained = sum(1 for r in legacy_disc if legacy_explained(r))
    out = {
        "engine": "bcftools 1.21 (htslib full VCF parse) vs pysam.TabixFile (raw text slice)",
        "region": REGION, "targets": TARGETS, "seed": SEED,
        "n_compared": len(rows) - n_err, "n_errors": n_err,
        "n_target_rows_concordant": n_ok,
        "n_rowcount_comparable": len(rowcount_cmp),
        "n_scan_vintage_new": len(new_rows),
        "n_scan_vintage_legacy": len(legacy_rows),
        "n_new_schema_target_concordant": n_new_ok,
        "new_schema_target_concordance": round(n_new_ok / max(1, len(new_rows)), 6),
        "n_legacy_discordant": len(legacy_disc),
        "n_legacy_discordant_explained_by_legacy_filter": n_legacy_explained,
        "n_window_rowcount_concordant": n_rows_ok,
        "n_window_rowcount_offby1": n_rows_offby1,
        "target_row_concordance": round(n_ok / max(1, len(rows) - n_err), 6),
        "window_rowcount_concordance": round(n_rows_ok / max(1, len(rowcount_cmp)), 6),
        "window_rowcount_offby1_fraction": round(n_rows_offby1 / max(1, len(rowcount_cmp)), 6),
        "discordant_samples": [
            {"sample": r["sample"], "only_pysam": r["only_pysam"], "only_bcftools": r["only_bcftools"]}
            for r in discordant[:50]],
    }
    with open(RESULT_JSON, "w") as fh:
        json.dump(out, fh, indent=2)
    md = [
        "## bcftools independent recheck of the Ag1000G dsx-window scan",
        "",
        f"Engine: {out['engine']}",
        f"Subsample: seed {SEED}, {N_CARRIER} carriers + {N_ZERO} zero-hit samples of 4,693.",
        f"Compared: {out['n_compared']} samples ({n_err} fetch errors).",
        f"Target-window row concordance (pos/ref/alt/GT): **{n_ok}/{out['n_compared']} = {out['target_row_concordance']}**",
        f"Whole-window row-count concordance (n_rows-carrying schemas only, "
        f"{len(rowcount_cmp)} samples): {n_rows_ok}/{len(rowcount_cmp)} "
        f"= {out['window_rowcount_concordance']}",
        f"Off-by-one row-count delta (bcftools = pysam + 1): {n_rows_offby1}/{len(rowcount_cmp)} "
        f"- coordinate convention: pysam fetch() is 0-based half-open, bcftools -r is 1-based "
        f"inclusive, so bcftools also returns the row at 2R:48,711,450 (outside all six targets).",
        "",
        "Shared transport (both link htslib for bgzf/tabix/HTTPS); independent parse layer",
        "(raw text slicing vs full VCF parse + re-emit).",
        "",
        f"Scan-vintage split: {len(new_rows)} scan-of-record rows (unfiltered) vs "
        f"{len(legacy_rows)} legacy rows (ag1000g_on_target.py, GQ>=20/DP>=5 filter).",
        f"Target concordance on scan-of-record rows: **{n_new_ok}/{len(new_rows)} = "
        f"{out['new_schema_target_concordance']}**",
        f"Legacy-row discordance: {len(legacy_disc)} samples, of which "
        f"{n_legacy_explained} are fully explained by the legacy quality filter "
        f"(bcftools recovers exactly the calls that filter drops).",
        "",
        f"Discordant samples (first 50): {len(discordant)} total",
    ]
    for d in out["discordant_samples"]:
        md.append(f"- {d['sample']}: only_pysam={d['only_pysam']} only_bcftools={d['only_bcftools']}")
    with open(RESULT_MD, "w") as fh:
        fh.write("\n".join(md) + "\n")
    print(json.dumps({k: out[k] for k in ("n_compared", "n_errors", "target_row_concordance", "window_rowcount_concordance")}, indent=2))

if __name__ == "__main__":
    if "--scan" in sys.argv:
        run_scan()
    elif "--compare" in sys.argv:
        run_compare()
    else:
        print(__doc__)
