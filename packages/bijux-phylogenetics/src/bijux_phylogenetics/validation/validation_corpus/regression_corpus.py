from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.core.dataset import DatasetAuditReport, audit_dataset_inputs

from .contracts import RegressionDatasetCaseResult, RegressionDatasetCorpusReport
from .dataset_corpora import default_fixtures_root, fixture


def regression_summary(report: DatasetAuditReport) -> dict[str, object]:
    return {
        "readiness_decision": report.readiness_decision,
        "analysis_taxa": report.analysis_taxa,
        "allowed_analyses": report.allowed_analyses,
        "blocked_analyses": report.blocked_analyses,
        "blocker_count": len(report.blockers),
        "warning_count": len(report.warnings),
        "risk_level": report.risk_score.risk_level,
        "risk_score": round(report.risk_score.total_score, 15),
    }


def build_regression_dataset_corpus(
    *,
    fixtures_root: Path | None = None,
    expected_path: Path | None = None,
) -> RegressionDatasetCorpusReport:
    """Compare current benchmark summaries against checked-in regression expectations."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    expectation_path = (
        fixture(root, "expected", "benchmark_regression_dataset_corpus.json")
        if expected_path is None
        else expected_path
    )
    expected_payload = json.loads(expectation_path.read_text(encoding="utf-8"))
    current_reports = {
        "core_inference_ready_dataset": audit_dataset_inputs(
            fixture(root, "trees", "example_tree.nwk"),
            fixture(root, "metadata", "example_metadata.tsv"),
            fixture(root, "metadata", "example_traits_validate.tsv"),
            alignment_path=fixture(root, "alignments", "example_alignment.fasta"),
        ),
        "warning_rich_dataset": audit_dataset_inputs(
            fixture(root, "trees", "example_tree.nwk"),
            fixture(root, "metadata", "example_metadata_reordered.tsv"),
            fixture(root, "metadata", "example_traits_validate.tsv"),
            alignment_path=fixture(
                root, "alignments", "example_alignment_extra_taxon.fasta"
            ),
            tip_dates_path=fixture(root, "metadata", "example_tip_dates_invalid.tsv"),
            calibration_path=fixture(
                root, "metadata", "example_calibrations_invalid.tsv"
            ),
        ),
    }

    results: list[RegressionDatasetCaseResult] = []
    for name, expected in expected_payload.items():
        observed = regression_summary(current_reports[name])
        passed = expected == observed
        notes = (
            []
            if passed
            else ["observed regression summary drifted from the checked-in expectation"]
        )
        results.append(
            RegressionDatasetCaseResult(
                name=name,
                passed=passed,
                expected=expected,
                observed=observed,
                notes=notes,
            )
        )

    passed_case_count = sum(1 for case in results if case.passed)
    return RegressionDatasetCorpusReport(
        goal_id=245,
        corpus="regression-dataset-corpus",
        passed=passed_case_count == len(results),
        case_count=len(results),
        passed_case_count=passed_case_count,
        failed_case_count=len(results) - passed_case_count,
        cases=results,
        limitations=[
            "the regression corpus currently tracks stable summary fields rather than every nested report detail",
        ],
    )
