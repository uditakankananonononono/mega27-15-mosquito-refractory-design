"""mosqdesign command-line interface.

Subcommands
-----------
run-v1      Scan the dsx locus, score with the CNN, rank (v1), simulate drive grid.
run-v2      Constraint-aware re-ranking (v2).
figures     Regenerate the figure set.
rank-v3     Regenerate the v3 specificity-fused ranking artifacts.
offtargets  Enumerate <=N-mismatch near-matches of protospacer(s) in a FASTA.
"""
import argparse
import json
import sys


def _cmd_run_v1(_args):
    from mosqdesign import run_analysis
    summary = run_analysis.run()
    print(json.dumps({"top": summary["top_candidates"][:3]}, indent=1))


def _cmd_run_v2(_args):
    from mosqdesign import run_analysis_v2
    summary = run_analysis_v2.run()
    print(json.dumps({k: v for k, v in summary.items() if k != "candidates"} , indent=1)[:2000])


def _cmd_figures(_args):
    from mosqdesign import make_figures
    for p in (make_figures.make(), make_figures.make_fig6(),
              make_figures.make_fig7(), make_figures.make_fig8()):
        print(p)


def _cmd_rank_v3(_args):
    from mosqdesign import v3_rank
    summary = v3_rank.run()
    import json as _json
    print(_json.dumps({"n_pass": summary["n_pass"], "n_fail": summary["n_fail"],
                       "control": summary["control_kyrou"]}, indent=1))


def _cmd_offtargets(args):
    from mosqdesign.genome_scan import scan_genome_offtargets
    exclude = {}
    if args.exclude_self:
        for spec in args.exclude_self:
            proto, chrom, pos = spec.split(":")
            exclude[proto] = (chrom, int(pos))
    hits = scan_genome_offtargets(args.proto, args.fasta,
                                  max_mismatches=args.max_mismatches,
                                  exclude=exclude)
    counts = {p: {m: 0 for m in range(args.max_mismatches + 1)} for p in args.proto}
    for h in hits:
        counts[h.guide][h.mismatches] += 1
    out = {"fasta": args.fasta, "max_mismatches": args.max_mismatches,
           "counts": counts,
           "hits": [h.__dict__ for h in hits] if args.show_hits else "use --show-hits"}
    print(json.dumps(out, indent=1))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="mosqdesign", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("run-v1", help="locus scan + CNN efficacy + v1 ranking + drive grid")
    sub.add_parser("run-v2", help="constraint-aware v2 ranking")
    sub.add_parser("figures", help="regenerate figures 1, 6, 7, 8")
    sub.add_parser("rank-v3", help="v3 specificity-fused ranking (regenerates ranked_designs_v3.csv + v3_summary.json)")
    ot = sub.add_parser("offtargets", help="enumerate near-matches in a FASTA")
    ot.add_argument("--proto", action="append", required=True,
                    help="20-nt protospacer; repeat for several")
    ot.add_argument("--fasta", required=True, help="FASTA to scan")
    ot.add_argument("--max-mismatches", type=int, default=3)
    ot.add_argument("--exclude-self", action="append", default=[],
                    metavar="PROTO:CHROM:POS",
                    help="drop the on-target locus (repeatable)")
    ot.add_argument("--show-hits", action="store_true")
    args = ap.parse_args(argv)
    {"run-v1": _cmd_run_v1, "run-v2": _cmd_run_v2, "rank-v3": _cmd_rank_v3,
     "figures": _cmd_figures, "offtargets": _cmd_offtargets}[args.cmd](args)


if __name__ == "__main__":
    main()
