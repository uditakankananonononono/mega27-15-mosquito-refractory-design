"""Full Ag1000G dsx-window scan (AgamP4 2R:48,711,450-48,714,700).

Threaded (WORKERS below - 8 triggered server-side throttling at sample ~590),
checkpointed. Appends one JSONL row per sample to
results/ag1000g/on_target.jsonl; skips samples already recorded.
Keeps ALL non-ref calls in the window with their GQ/DP recorded (no quality
filtering here; downstream analysis re-filters), flagging sites inside guide
protospacer/PAM targets. Schema matches scripts/ag1000g_on_target.py.
Failed samples are recorded as {"sample": ..., "error": ...} rows so resume
skips them; drop error rows to retry.
"""
import json, os, sys, time
from concurrent.futures import ThreadPoolExecutor
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
WORKERS = 4  # 4->2 at 22:58, reverted 2->4 at 23:28 IST 2026-09-24 per Main 23:25 directive: 4w kill-and-relaunch ~17 rows/min effective beats 2w lumpy ~4-13 (Appendix Q)

def fetch_sample(sid):
    t0 = time.time()
    tf = pysam.TabixFile(BASE + sid + ".vcf.gz")
    hits, n_rows = [], 0
    for row in tf.fetch("2R", FETCH[0], FETCH[1]):
        n_rows += 1
        f = row.rstrip("\n").split("\t")
        gt_field = f[9].split(":"); fmt = f[8].split(":")
        gt = gt_field[0]
        if gt == "0/0" or gt == "./.":
            continue
        try:
            gq = int(gt_field[fmt.index("GQ")]); dp = int(gt_field[fmt.index("DP")])
        except (ValueError, IndexError):
            gq, dp = None, None
        pos = int(f[1])
        tgt = [n for n, (a, b) in TARGETS.items() if a <= pos <= b]
        hits.append({"pos": pos, "ref": f[3], "alt": f[4], "gt": gt,
                     "gq": gq, "dp": dp, "targets": tgt})
    tf.close()
    return {"sample": sid, "n_rows": n_rows, "n_nonref": len(hits),
            "target_hits": hits, "seconds": round(time.time() - t0, 2)}

def main():
    ids = [l.strip() for l in open("/tmp/ag3_samples.txt") if l.strip()]
    done = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            done.add(json.loads(line)["sample"])
    todo = [s for s in ids if s not in done]
    print(f"{len(todo)} samples to go", flush=True)
    lock_fh = open(OUT, "a")
    count = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for rec in ex.map(_safe, todo, chunksize=1):
            lock_fh.write(json.dumps(rec) + "\n"); lock_fh.flush()
            count += 1
            if count % 50 == 0:
                print(f"done {count}/{len(todo)}", flush=True)

def _safe(sid):
    try:
        return fetch_sample(sid)
    except Exception as e:
        return {"sample": sid, "error": str(e)[:200]}

if __name__ == "__main__":
    main()
