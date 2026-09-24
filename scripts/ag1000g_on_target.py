"""Ag1000G on-target polymorphism check for dsx guide targets (AgamP4 coords).

For each per-sample all-sites VCF on the Sanger S3 mirror, fetch the tabix
window covering the three guide protospacers + PAMs and record every site
with a non-reference, reasonably-confident genotype inside each target window.

Targets (AgamP4 2R, minus-strand guides; PAM = 3 bp upstream genomic):
  dsx-v3-1: proto 48711501-48711520, PAM 48711498-48711500
  dsx-v3-2: proto 48712765-48712784, PAM 48712762-48712764
  kyrou   : proto 48714640-48714659, PAM 48714637-48714639
Checkpointed: appends one JSONL row per sample; skips samples already done.
"""
import json, os, sys, time
import pysam

BASE = "https://vo_agam_output.cog.sanger.ac.uk/"
OUT = "results/ag1000g/on_target.jsonl"
FETCH = (48711450, 48714700)
TARGETS = {
    "dsx-v3-1_proto": (48711501, 48711520),
    "dsx-v3-1_PAM":   (48711498, 48711500),
    "dsx-v3-2_proto": (48712765, 48712784),
    "dsx-v3-2_PAM":   (48712762, 48712764),
    "kyrou_proto":    (48714640, 48714659),
    "kyrou_PAM":      (48714637, 48714639),
}
GQ_MIN, DP_MIN = 20, 5

def done_set():
    done = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            done.add(json.loads(line)["sample"])
    return done

def fetch_sample(sid):
    t0 = time.time()
    tf = pysam.TabixFile(BASE + sid + ".vcf.gz")
    hits, n_rows = [], 0
    for row in tf.fetch("2R", FETCH[0], FETCH[1]):
        n_rows += 1
        f = row.rstrip("\n").split("\t")
        pos = int(f[1])
        gt_field = f[9].split(":")
        fmt = f[8].split(":")
        gt = gt_field[0]
        if gt in ("0/0", "./.") or gt.startswith("0/0"):
            continue
        try:
            gq = int(gt_field[fmt.index("GQ")])
            dp = int(gt_field[fmt.index("DP")])
        except (ValueError, IndexError):
            continue
        if gq < GQ_MIN or dp < DP_MIN:
            continue
        in_targets = [name for name, (a, b) in TARGETS.items() if a <= pos <= b]
        hits.append({"pos": pos, "ref": f[3], "alt": f[4], "gt": gt,
                     "gq": gq, "dp": dp, "targets": in_targets})
    tf.close()
    return {"sample": sid, "n_rows": n_rows, "n_nonref": len(hits),
            "target_hits": [h for h in hits if h["targets"]],
            "seconds": round(time.time() - t0, 2)}

def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    ids = [l.strip() for l in open("/tmp/ag3_samples.txt") if l.strip()]
    done = done_set()
    todo = [s for s in ids if s not in done][:n]
    with open(OUT, "a") as fh:
        for sid in todo:
            try:
                rec = fetch_sample(sid)
            except Exception as e:
                rec = {"sample": sid, "error": str(e)[:200]}
            fh.write(json.dumps(rec) + "\n"); fh.flush()
            tgt = len(rec.get("target_hits", []))
            print(f"{sid}: rows={rec.get('n_rows')} nonref={rec.get('n_nonref')} "
                  f"target_hits={tgt} t={rec.get('seconds')}s", flush=True)

if __name__ == "__main__":
    main()
