"""Hermetic tests for ViennaRNA folding and Primer3 amplicon design."""
import importlib.util
import os

spec = importlib.util.spec_from_file_location(
    "guide_fold_primers",
    os.path.join(os.path.dirname(__file__), "..", "scripts", "guide_fold_primers.py"))
gp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gp)

KYROU = "GTTTAACACAGGTCAAGCGG"


def test_fold_is_deterministic():
    s1, m1 = gp.fold_spacer(KYROU)
    s2, m2 = gp.fold_spacer(KYROU)
    assert (s1, m1) == (s2, m2)
    assert len(s1) == 20 and m1 <= 0.0


def test_amplicon_spans_window():
    for name in ("dsx-v3-1", "dsx-v3-2", "kyrou"):
        a = gp.design_amplicon(name)
        assert a["designed"], f"no primer pair for {name}: {a}"
        lo, hi = a["amplicon_spans"]
        assert 250 <= a["product_size"] <= 450
        # primers must sit outside the targeted window
        assert a["left"] and a["right"]


def test_full_report_consistent():
    report = gp.run()
    for name, r in report.items():
        assert len(r["spacer"]) == 20
        assert r["mfe_kcal_mol"] <= 0.0
        assert r["amplicon"]["designed"]
