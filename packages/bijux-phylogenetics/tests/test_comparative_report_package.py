from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.comparative.reporting import build_comparative_method_report
from bijux_phylogenetics.comparative.reporting.analysis_package import (
    build_comparative_report_package,
    summarize_comparative_interpretation,
)

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


def test_build_comparative_report_package_writes_html_and_ledgers(
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

    assert result.report_path.exists()
    assert result.methods_summary_path.exists()
    assert result.reviewer_audit_checklist_path.exists()
    assert result.summary_table_path.exists()
    assert result.coefficient_table_path.exists()
    assert result.residual_table_path.exists()
    assert result.signal_table_path.exists()
    assert result.model_comparison_table_path.exists()
    assert result.interpretation_table_path.exists()
    assert result.audit_table_path.exists()
    assert result.contrast_table_path.exists()
    assert result.manifest_path.exists()

    html = result.report_path.read_text(encoding="utf-8")
    assert "Bijux Comparative Analysis Report" in html
    assert "Method Tier" in html
    assert "supported" in html
    assert "Methods Summary" in html
    assert "Reviewer Audit Checklist" in html
    assert "Coefficient Table" in html
    assert "Residual Summary" in html
    assert "Phylogenetic Signal" in html
    assert "Model Comparison" in html
    assert "Biological Interpretation" in html

    summary_rows = result.summary_table_path.read_text(encoding="utf-8").splitlines()
    methods_summary = result.methods_summary_path.read_text(encoding="utf-8")
    coefficient_rows = result.coefficient_table_path.read_text(
        encoding="utf-8"
    ).splitlines()
    residual_rows = result.residual_table_path.read_text(encoding="utf-8").splitlines()
    signal_rows = result.signal_table_path.read_text(encoding="utf-8").splitlines()
    comparison_rows = result.model_comparison_table_path.read_text(
        encoding="utf-8"
    ).splitlines()
    interpretation_rows = result.interpretation_table_path.read_text(
        encoding="utf-8"
    ).splitlines()
    checklist_rows = result.reviewer_audit_checklist_path.read_text(
        encoding="utf-8"
    ).splitlines()
    contrast_rows = result.contrast_table_path.read_text(encoding="utf-8").splitlines()

    assert summary_rows[0].startswith("response\tformula\tpredictor_count")
    assert "Comparative Analysis Methods Summary" in methods_summary
    assert "- predictor terms: `predictor_one`" in methods_summary
    assert coefficient_rows[0].startswith("term\testimate\tstandard_error")
    assert residual_rows[0].startswith("analysis\tresidual_variance")
    assert signal_rows[0].startswith("trait\ttaxon_count\tblombergs_k")
    assert comparison_rows[0].startswith("model\tparameter_count\tlog_likelihood")
    assert interpretation_rows[0].startswith("topic\tclaim\tevidence\tcaution")
    assert checklist_rows[0] == "section\tstatus\tsummary\tevidence\tartifact_paths"
    assert any(line.startswith("model_selection\t") for line in checklist_rows[1:])
    assert contrast_rows[0].startswith("node\tleft_taxa\tright_taxa")

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest["report_kind"] == "comparative_package"
    assert manifest["metrics"]["analysis_taxa"] == 4
    assert manifest["metrics"]["selected_model"] in {"brownian", "ou"}
    assert manifest["metrics"]["methods_summary_warning_count"] >= 0
    assert manifest["outputs"]["methods_summary_path"].endswith(
        "comparative-methods-summary.md"
    )
    assert manifest["outputs"]["reviewer_audit_checklist_path"].endswith(
        "reviewer-audit-checklist.tsv"
    )
    assert len(manifest["reviewer_audit_checklist"]["items"]) == 5
    assert result.method_tier.tier == "supported"


def test_summarize_comparative_interpretation_includes_signal_model_and_coefficients() -> (
    None
):
    report = build_comparative_method_report(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value=0.0,
    )

    rows = summarize_comparative_interpretation(report)

    topics = {row.topic for row in rows}
    assert "formula" in topics
    assert "phylogenetic-signal" in topics
    assert "model-comparison" in topics
    assert "coefficient" in topics
