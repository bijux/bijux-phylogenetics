from __future__ import annotations

import json
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


def test_primate_pgls_and_signal_reference_uses_governed_dataset_paths() -> None:
    reference_script = (
        STUDY_ROOT / "reference" / "primate_pgls_and_signal_reference_r.R"
    ).read_text(encoding="utf-8")
    assert '"datasets"' in reference_script
    assert '"reference_primate.csv"' in reference_script
    assert '"reference_trimmed_primatetree.nwk"' in reference_script


def test_primate_pgls_and_signal_pagel_lambda_report_records_fixed_and_estimated_modes() -> (
    None
):
    payload = json.loads(
        (
            STUDY_ROOT
            / "evidence-003"
            / "results"
            / "pagel-lambda-regression-parity.json"
        ).read_text(encoding="utf-8")
    )
    assert "r_fixed_reference_lambda" in payload
    assert "bijux_fixed_reference_lambda" in payload
    assert "r_estimated_lambda" in payload
    assert "bijux_estimated_lambda" in payload
    assert (
        payload["r_fixed_reference_lambda"]["aic"]
        < payload["r_estimated_lambda"]["aic"]
    )
