from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
STUDY_ROOT = REPO_ROOT / "evidence-book" / "studies" / "primate-pgls-and-signal"
ALLOWED_STUDY_ENTRIES = {
    "README.md",
    "datasets",
    "provenance",
    "reference",
    *{f"evidence-{index:03d}" for index in range(1, 11)},
}
ALLOWED_BUNDLE_ENTRIES = {
    "README.md",
    "analysis.py",
    "checks.json",
    "claims.json",
    "inputs.manifest.json",
    "manifest.json",
    "provenance.json",
    "reference.R",
    "report.md",
    "results",
}


def test_primate_pgls_and_signal_study_root_is_human_first() -> None:
    assert {path.name for path in STUDY_ROOT.iterdir()} == ALLOWED_STUDY_ENTRIES


def test_primate_pgls_and_signal_bundles_use_one_repeatable_layout() -> None:
    for evidence_id in sorted(
        name for name in ALLOWED_STUDY_ENTRIES if name.startswith("evidence-")
    ):
        bundle_root = STUDY_ROOT / evidence_id
        assert {path.name for path in bundle_root.iterdir()} == ALLOWED_BUNDLE_ENTRIES
        assert (bundle_root / "results" / "manifest.json").is_file()
        assert (bundle_root / "results" / "README.md").is_file()


def test_primate_pgls_and_signal_scalar_outputs_live_under_results() -> None:
    bundle_root = STUDY_ROOT / "evidence-001" / "results"
    assert (bundle_root / "scalar-parity-table.json").is_file()
    assert (bundle_root / "scalar-parity-table.md").is_file()
