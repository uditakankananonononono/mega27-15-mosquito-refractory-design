"""Hermetic tests for the Biopython window re-derivation."""
import importlib.util
import os

spec = importlib.util.spec_from_file_location(
    "verify_windows_biopython",
    os.path.join(os.path.dirname(__file__), "..", "scripts", "verify_windows_biopython.py"))
vw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vw)


def test_all_windows_pass():
    report = vw.derive()
    assert set(report) == {"dsx-v3-1", "dsx-v3-2", "kyrou"}
    for name, r in report.items():
        assert r["pass"], f"{name} failed independent re-derivation: {r}"


def test_kyrou_exact_match_to_agamp4():
    r = vw.derive()["kyrou"]
    assert r["ref_minus_derived"] == "GTTTAACACAGGTCAAGCGG"
    assert r["pam_minus_derived"] == "TGG"


def test_pam_strand_convention():
    # minus-strand PAMs must read 5'-NGG-3' on the minus strand
    report = vw.derive()
    for r in report.values():
        assert r["pam_minus_derived"].endswith("GG")
