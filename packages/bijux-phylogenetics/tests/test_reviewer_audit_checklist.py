from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.comparative.reporting.analysis_package import (
    build_comparative_report_package,
)
from bijux_phylogenetics.reports import (
    build_reviewer_audit_checklist,
    write_reviewer_audit_checklist_from_manifest,
)
from bijux_phylogenetics.reports.publication.tree import build_tree_report_package

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_build_reviewer_audit_checklist_supports_tree_package_manifests(
    tmp_path: Path,
) -> None:
    result = build_tree_report_package(
        tree_fixture("example_tree_support_invalid.nwk"),
        out_dir=tmp_path / "tree-report-package",
    )
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))

    checklist = build_reviewer_audit_checklist(manifest)

    statuses = {item.section: item.status for item in checklist.items}
    assert checklist.report_kind == "tree_package"
    assert statuses["inputs"] == "pass"
    assert statuses["methods"] == "pass"
    assert statuses["validity"] == "blocked"
    assert statuses["support_surface"] == "blocked"
    assert statuses["interpretation_limits"] == "risk"


def test_write_reviewer_audit_checklist_from_manifest_writes_tsv(
    tmp_path: Path,
) -> None:
    result = build_comparative_report_package(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        out_dir=tmp_path / "comparative-report-package",
        response="response",
        predictors=["predictor_one"],
        lambda_value=0.0,
    )

    checklist_result = write_reviewer_audit_checklist_from_manifest(
        tmp_path / "reviewer-audit-checklist.tsv",
        result.manifest_path,
    )

    assert checklist_result.output_path.exists()
    lines = checklist_result.output_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "section\tstatus\tsummary\tevidence\tartifact_paths"
    assert any(line.startswith("model_selection\t") for line in lines[1:])
    assert checklist_result.checklist.report_kind == "comparative_package"
