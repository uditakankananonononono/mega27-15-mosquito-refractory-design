"""Offline fixture tests for the SnpEff window annotation (tool 44).

Exercises distinct_sites GT allele-index resolution (including multiallelic
rows and both scan vintages), VCF writing, ANN parsing, and summarize's
worst-impact ranking - no SnpEff run, no network.
"""
import json
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "snpeff_window_ann", Path(__file__).parent.parent / "scripts" / "snpeff_window_ann.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_distinct_sites_multiallelic_and_vintages(tmp_path, monkeypatch):
    scan = tmp_path / "on_target.jsonl"
    rows = [
        # new vintage, multiallelic: GT 0/2 -> second ALT
        {"sample": "S1", "n_rows": 10, "target_hits": [
            {"pos": 48711500, "ref": "C", "alt": "A,T,G", "gt": "0/2",
             "gq": 90, "dp": 20, "targets": ["dsx-v3-1_PAM"]}]},
        # legacy vintage, GT 2/2 -> two copies of second ALT
        {"sample": "S2", "hits": [
            {"pos": 48711500, "ref": "C", "alt": "A,T,G", "gt": "2/2",
             "gq": 90, "dp": 20, "targets": ["dsx-v3-1_PAM"]},
            {"pos": 48711600, "ref": "G", "alt": "T", "gt": "1/1",
             "gq": 90, "dp": 20, "targets": []}]},
        {"sample": "S3", "error": "fetch failed"},
    ]
    scan.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    monkeypatch.setattr(m, "SCAN_IN", str(scan))
    sites = m.distinct_sites(str(scan))
    assert sites[(48711500, "C", "T")] == 3   # 0/2 + 2/2 = 1 + 2 copies
    assert sites[(48711600, "G", "T")] == 2   # 1/1
    assert len(sites) == 2                    # error row skipped


def test_write_vcf_and_parse_ann_roundtrip(tmp_path):
    vcf = tmp_path / "sites.vcf"
    m.write_vcf({(48711500, "C", "T"): 3, (48711507, "C", "A"): 1}, str(vcf))
    lines = vcf.read_text().splitlines()
    assert lines[0].startswith("##fileformat") and lines[2].startswith("#CHROM")
    assert lines[3].split("\t")[1] == "48711500"
    # fabricate the corresponding ANN output and parse it
    ann = tmp_path / "sites.ann.vcf"
    ann.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "2R\t48711500\t.\tC\tT\t.\t.\tAC=3;ANN=T|intron_variant|MODIFIER|AGAP004050|"
        "transcript|AGAP004050-RA|.\t\n"
        "2R\t48711507\t.\tC\tA\t.\t.\tAC=1;ANN=A|missense_variant|MODERATE|AGAP004050|"
        "transcript|AGAP004050-RA|.\t\n")
    recs = m.parse_ann(str(ann))
    assert len(recs) == 2
    assert recs[0]["ac"] == 3 and recs[0]["effects"][0]["impact"] == "MODIFIER"
    assert recs[1]["targets"] == ["dsx-v3-1_proto"]


def test_summarize_worst_impact_ranking():
    recs = [
        {"pos": 48712765, "ref": "A", "alt": "G", "ac": 2,
         "targets": ["dsx-v3-2_proto"],
         "effects": [{"effect": "synonymous_variant", "impact": "LOW", "gene": "AGAP004050", "transcript": "t1"},
                     {"effect": "missense_variant", "impact": "MODERATE", "gene": "AGAP004050", "transcript": "t2"}]},
        {"pos": 48711500, "ref": "C", "alt": "T", "ac": 9,
         "targets": ["dsx-v3-1_PAM"],
         "effects": [{"effect": "intron_variant", "impact": "MODIFIER", "gene": "AGAP004050", "transcript": "t1"}]},
        {"pos": 48712000, "ref": "G", "alt": "A", "ac": 1, "targets": [], "effects": []},
    ]
    out = m.summarize(recs)
    assert out["n_sites"] == 3
    assert out["n_protein_altering_target_sites"] == 1
    assert out["protein_altering_target_sites"][0]["pos"] == 48712765
    assert out["effect_classes"]["intron_variant"] == 1
    assert out["effect_classes"]["intergenic_or_unannotated"] == 1
    assert out["per_guide_impact_of_variant_sites"]["dsx-v3-1|MODIFIER"] == 1
