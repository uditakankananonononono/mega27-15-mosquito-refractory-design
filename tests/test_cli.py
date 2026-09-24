"""Pin the CLI: dispatch works and the offtargets subcommand reproduces the
paper's locus-scale control result on synthetic and archived substrates."""
import json
import os

from mosqdesign.cli import main


def test_cli_offtargets_synthetic(capsys, tmp_path):
    fasta = tmp_path / "synth.fasta"
    # one exact occurrence of the guide plus one 2-mismatch copy, both plus strand
    guide = "ACGTACGTACGTACGTACGT"
    near = "ACGTACGTACGTACGTACGG".replace("GG", "GA")  # 1 mismatch at pos 19
    seq = "TTTT" + guide + "A" * 50 + near + "TTTT"
    fasta.write_text(">synth\n" + seq + "\n")
    main(["offtargets", "--proto", guide, "--fasta", str(fasta),
          "--max-mismatches", "3"])
    out = json.loads(capsys.readouterr().out)
    assert out["counts"][guide]["0"] == 1
    assert out["counts"][guide]["1"] == 1
    assert out["counts"][guide]["2"] == 0


def test_cli_offtargets_kyrou_locus_control(capsys):
    # archived AgamP5 dsx locus; Kyrou guide must show exactly its own site and
    # nothing else at <=2 mismatches (paper Section 3.1 control result)
    fasta = os.path.join(os.path.dirname(__file__), "..", "data",
                         "dsx_locus_region.fasta")
    guide = "GTTTAACACAGGTCAAGCGG"
    main(["offtargets", "--proto", guide, "--fasta", fasta,
          "--max-mismatches", "2"])
    out = json.loads(capsys.readouterr().out)
    assert out["counts"][guide]["0"] == 1  # its own locus
    assert out["counts"][guide]["1"] == 0
    assert out["counts"][guide]["2"] == 0


def test_cli_rank_v3_regenerates_committed_artifacts(tmp_path):
    """rank-v3 must reproduce the committed v3 artifacts byte-for-byte."""
    import hashlib
    import shutil
    before = {}
    for f in ("results/ranked_designs_v3.csv", "results/v3_summary.json"):
        before[f] = hashlib.sha256(open(f, "rb").read()).hexdigest()
    main(["rank-v3"])
    for f, h in before.items():
        assert hashlib.sha256(open(f, "rb").read()).hexdigest() == h, f"{f} changed"
