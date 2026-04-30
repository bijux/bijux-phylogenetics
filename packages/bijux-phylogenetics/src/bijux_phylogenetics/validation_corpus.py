from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.benchmark import (
    BenchmarkObservation,
    benchmark_alignment_site_scaling,
    benchmark_tree_comparison,
    benchmark_tree_set_consensus,
    benchmark_tree_validation,
)
from bijux_phylogenetics.core.dataset import DatasetAuditReport, audit_dataset_inputs
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.errors import PhylogeneticsError
from bijux_phylogenetics.io.fasta import summarise_fasta
from bijux_phylogenetics.reference_validation import build_core_workflow_validation_report


@dataclass(frozen=True, slots=True)
class CorpusDatasetCase:
    """One checked-in dataset case used in a benchmark corpus."""

    name: str
    tree_path: Path
    metadata_path: Path
    traits_path: Path
    alignment_path: Path | None = None
    tip_dates_path: Path | None = None
    calibration_path: Path | None = None
    required_allowed_analyses: tuple[str, ...] = ()
    forbidden_blockers: tuple[str, ...] = ()
    allowed_warning_prefixes: tuple[str, ...] = ()


@dataclass(slots=True)
class CorpusDatasetCaseResult:
    """Observed result for one evaluated dataset case."""

    name: str
    passed: bool
    readiness_decision: str
    analysis_taxa: list[str]
    allowed_analyses: list[str]
    blocked_analyses: list[str]
    blockers: list[str]
    warnings: list[str]
    notes: list[str]
    observed_code: str | None = None


@dataclass(slots=True)
class BenchmarkCorpusReport:
    """Reviewer-facing summary for one benchmark corpus."""

    goal_id: int
    corpus: str
    passed: bool
    case_count: int
    passed_case_count: int
    failed_case_count: int
    cases: list[CorpusDatasetCaseResult]
    limitations: list[str]


@dataclass(slots=True)
class RegressionDatasetCaseResult:
    """Observed-versus-expected summary for one regression dataset case."""

    name: str
    passed: bool
    expected: dict[str, object]
    observed: dict[str, object]
    notes: list[str]


@dataclass(slots=True)
class RegressionDatasetCorpusReport:
    """Stable biological summary snapshots tracked across releases."""

    goal_id: int
    corpus: str
    passed: bool
    case_count: int
    passed_case_count: int
    failed_case_count: int
    cases: list[RegressionDatasetCaseResult]
    limitations: list[str]


@dataclass(slots=True)
class MethodAccuracyRow:
    """One validation surface summarized for accuracy, error, and coverage."""

    surface: str
    accuracy: float
    passed_count: int
    failed_count: int
    coverage_count: int
    bias_notes: list[str]
    error_notes: list[str]


@dataclass(slots=True)
class MethodAccuracyDashboard:
    """Goal 246 dashboard across the main checked-in validation surfaces."""

    goal_id: int
    rows: list[MethodAccuracyRow]
    limitations: list[str]


@dataclass(slots=True)
class BenchmarkDashboardRow:
    """One workflow scaling curve in runtime or memory dashboards."""

    workflow: str
    scaling_axis: str
    observations: list[BenchmarkObservation]


@dataclass(slots=True)
class RuntimeBenchmarkDashboard:
    """Goal 247 runtime scaling summary across major benchmark axes."""

    goal_id: int
    rows: list[BenchmarkDashboardRow]
    limitations: list[str]


@dataclass(slots=True)
class MemoryBenchmarkDashboard:
    """Goal 248 memory scaling summary across major benchmark axes."""

    goal_id: int
    rows: list[BenchmarkDashboardRow]
    limitations: list[str]


def _default_fixtures_root() -> Path:
    return Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _fixture(root: Path, *parts: str) -> Path:
    return root.joinpath(*parts)


def _evaluate_dataset_case(case: CorpusDatasetCase) -> tuple[DatasetAuditReport, list[str]]:
    report = audit_dataset_inputs(
        case.tree_path,
        case.metadata_path,
        case.traits_path,
        alignment_path=case.alignment_path,
        tip_dates_path=case.tip_dates_path,
        calibration_path=case.calibration_path,
    )
    notes: list[str] = []
    missing_allowed = sorted(set(case.required_allowed_analyses) - set(report.allowed_analyses))
    if missing_allowed:
        notes.append(f"missing required allowed analyses: {', '.join(missing_allowed)}")
    forbidden = [blocker for blocker in report.blockers if blocker in case.forbidden_blockers]
    if forbidden:
        notes.append(f"observed forbidden blockers: {', '.join(forbidden)}")
    disallowed_warnings = [
        warning
        for warning in report.warnings
        if case.allowed_warning_prefixes
        and not any(warning.startswith(prefix) for prefix in case.allowed_warning_prefixes)
    ]
    if disallowed_warnings:
        notes.append(f"unexpected warnings: {', '.join(disallowed_warnings)}")
    return report, notes


def build_clean_benchmark_corpus(*, fixtures_root: Path | None = None) -> BenchmarkCorpusReport:
    """Validate that the checked-in clean corpus stays usable for core workflows."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    cases = [
        CorpusDatasetCase(
            name="core_inference_ready_dataset",
            tree_path=_fixture(root, "trees", "example_tree.nwk"),
            metadata_path=_fixture(root, "metadata", "example_metadata.tsv"),
            traits_path=_fixture(root, "metadata", "example_traits_validate.tsv"),
            alignment_path=_fixture(root, "alignments", "example_alignment.fasta"),
            required_allowed_analyses=(
                "inspection",
                "distance",
                "maximum_likelihood",
                "bayesian",
                "comparative",
                "publication",
            ),
            forbidden_blockers=(
                "metadata table is missing one or more tree taxa",
                "trait table is missing one or more tree taxa",
                "alignment is missing one or more tree taxa",
                "fewer than two taxa remain after intersecting all requested dataset surfaces",
            ),
            allowed_warning_prefixes=(
                "equal-length ungapped FASTA may be aligned",
            ),
        )
    ]

    results: list[CorpusDatasetCaseResult] = []
    for case in cases:
        report, notes = _evaluate_dataset_case(case)
        passed = not notes
        results.append(
            CorpusDatasetCaseResult(
                name=case.name,
                passed=passed,
                readiness_decision=report.readiness_decision,
                analysis_taxa=report.analysis_taxa,
                allowed_analyses=report.allowed_analyses,
                blocked_analyses=report.blocked_analyses,
                blockers=report.blockers,
                warnings=report.warnings,
                notes=notes,
            )
        )
    passed_case_count = sum(1 for case in results if case.passed)
    return BenchmarkCorpusReport(
        goal_id=242,
        corpus="clean-benchmark-corpus",
        passed=passed_case_count == len(results),
        case_count=len(results),
        passed_case_count=passed_case_count,
        failed_case_count=len(results) - passed_case_count,
        cases=results,
        limitations=[
            "the clean corpus currently focuses on core inference and publication readiness rather than time-tree readiness",
        ],
    )


def build_broken_benchmark_corpus(*, fixtures_root: Path | None = None) -> BenchmarkCorpusReport:
    """Validate that intentionally broken fixtures still fail for the expected reason."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    results: list[CorpusDatasetCaseResult] = []

    def record_error_case(name: str, callback, expected_code: str) -> None:
        observed_code: str | None = None
        notes: list[str] = []
        try:
            callback()
            notes.append(f"expected error code {expected_code} but the fixture completed successfully")
        except PhylogeneticsError as error:
            observed_code = error.code
            if observed_code != expected_code:
                notes.append(f"expected error code {expected_code} but observed {observed_code}")
        results.append(
            CorpusDatasetCaseResult(
                name=name,
                passed=not notes,
                readiness_decision="error" if observed_code is not None else "unexpected_success",
                analysis_taxa=[],
                allowed_analyses=[],
                blocked_analyses=[],
                blockers=[],
                warnings=[],
                notes=notes,
                observed_code=observed_code,
            )
        )

    record_error_case(
        "duplicate_tip_tree",
        lambda: validate_tree_path(_fixture(root, "trees", "example_tree_duplicate.nwk")),
        "duplicate_taxon_error",
    )
    record_error_case(
        "invalid_alignment_lengths",
        lambda: summarise_fasta(_fixture(root, "alignments", "example_alignment_invalid_lengths.fasta")),
        "invalid_alignment_error",
    )

    broken_dataset = CorpusDatasetCase(
        name="dataset_missing_metadata_taxon",
        tree_path=_fixture(root, "trees", "example_tree.nwk"),
        metadata_path=_fixture(root, "metadata", "example_metadata_missing_taxon.csv"),
        traits_path=_fixture(root, "metadata", "example_traits_validate.tsv"),
        alignment_path=_fixture(root, "alignments", "example_alignment.fasta"),
    )
    report, notes = _evaluate_dataset_case(broken_dataset)
    if report.readiness_decision != "blocked":
        notes.append(f"expected blocked readiness but observed {report.readiness_decision}")
    if "metadata table is missing one or more tree taxa" not in report.blockers:
        notes.append("expected metadata-missing blocker was not observed")
    results.append(
        CorpusDatasetCaseResult(
            name=broken_dataset.name,
            passed=not notes,
            readiness_decision=report.readiness_decision,
            analysis_taxa=report.analysis_taxa,
            allowed_analyses=report.allowed_analyses,
            blocked_analyses=report.blocked_analyses,
            blockers=report.blockers,
            warnings=report.warnings,
            notes=notes,
            observed_code="dataset_blocked" if report.readiness_decision == "blocked" else None,
        )
    )

    passed_case_count = sum(1 for case in results if case.passed)
    return BenchmarkCorpusReport(
        goal_id=243,
        corpus="broken-benchmark-corpus",
        passed=passed_case_count == len(results),
        case_count=len(results),
        passed_case_count=passed_case_count,
        failed_case_count=len(results) - passed_case_count,
        cases=results,
        limitations=[
            "the broken corpus emphasizes stable failure signatures rather than every possible malformed input surface",
        ],
    )


def build_messy_benchmark_corpus(*, fixtures_root: Path | None = None) -> BenchmarkCorpusReport:
    """Validate that realistic multi-problem datasets surface the expected warning mix."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    cases: list[tuple[CorpusDatasetCase, tuple[str, ...], tuple[str, ...]]] = [
        (
            CorpusDatasetCase(
                name="reordered_alignment_extra_taxa_invalid_dates_and_calibrations",
                tree_path=_fixture(root, "trees", "example_tree.nwk"),
                metadata_path=_fixture(root, "metadata", "example_metadata_reordered.tsv"),
                traits_path=_fixture(root, "metadata", "example_traits_validate.tsv"),
                alignment_path=_fixture(root, "alignments", "example_alignment_extra_taxon.fasta"),
                tip_dates_path=_fixture(root, "metadata", "example_tip_dates_invalid.tsv"),
                calibration_path=_fixture(root, "metadata", "example_calibrations_invalid.tsv"),
            ),
            (
                "calibration table contains invalid fossil calibration targets or ages",
                "tip-date metadata contains invalid or missing tree-tip dates",
            ),
            (
                "alignment contains taxa absent from the tree",
                "one or more dataset surfaces silently reorder shared taxa relative to the tree",
                "tip-date metadata contains taxa absent from the tree or alignment",
            ),
        ),
        (
            CorpusDatasetCase(
                name="low_information_alignment_with_trait_mismatch",
                tree_path=_fixture(root, "trees", "example_tree.nwk"),
                metadata_path=_fixture(root, "metadata", "example_metadata.tsv"),
                traits_path=_fixture(root, "metadata", "example_traits.tsv"),
                alignment_path=_fixture(root, "alignments", "example_alignment_ambiguity.fasta"),
            ),
            (
                "alignment is missing one or more tree taxa",
                "alignment is not currently safe for core inference workflows",
                "trait table is missing one or more tree taxa",
            ),
            (
                "alignment contains near-duplicate sequences",
                "alignment has fewer parsimony-informative sites than the minimum threshold for defensible inference",
                "trait table contains taxa absent from the tree",
            ),
        ),
    ]

    results: list[CorpusDatasetCaseResult] = []
    for case, expected_blockers, expected_warnings in cases:
        report, notes = _evaluate_dataset_case(case)
        missing_blockers = [message for message in expected_blockers if message not in report.blockers]
        if missing_blockers:
            notes.append(f"missing expected blockers: {', '.join(missing_blockers)}")
        missing_warnings = [message for message in expected_warnings if message not in report.warnings]
        if missing_warnings:
            notes.append(f"missing expected warnings: {', '.join(missing_warnings)}")
        if report.readiness_decision != "blocked":
            notes.append(f"expected blocked readiness but observed {report.readiness_decision}")
        results.append(
            CorpusDatasetCaseResult(
                name=case.name,
                passed=not notes,
                readiness_decision=report.readiness_decision,
                analysis_taxa=report.analysis_taxa,
                allowed_analyses=report.allowed_analyses,
                blocked_analyses=report.blocked_analyses,
                blockers=report.blockers,
                warnings=report.warnings,
                notes=notes,
                observed_code="dataset_blocked",
            )
        )

    passed_case_count = sum(1 for case in results if case.passed)
    return BenchmarkCorpusReport(
        goal_id=244,
        corpus="messy-benchmark-corpus",
        passed=passed_case_count == len(results),
        case_count=len(results),
        passed_case_count=passed_case_count,
        failed_case_count=len(results) - passed_case_count,
        cases=results,
        limitations=[
            "the messy corpus currently emphasizes warning-rich integration failures rather than every possible domain-specific artifact",
        ],
    )


def _regression_summary(report: DatasetAuditReport) -> dict[str, object]:
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
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    expectation_path = (
        _fixture(root, "expected", "benchmark_regression_dataset_corpus.json")
        if expected_path is None
        else expected_path
    )
    expected_payload = json.loads(expectation_path.read_text(encoding="utf-8"))
    current_reports = {
        "core_inference_ready_dataset": audit_dataset_inputs(
            _fixture(root, "trees", "example_tree.nwk"),
            _fixture(root, "metadata", "example_metadata.tsv"),
            _fixture(root, "metadata", "example_traits_validate.tsv"),
            alignment_path=_fixture(root, "alignments", "example_alignment.fasta"),
        ),
        "warning_rich_dataset": audit_dataset_inputs(
            _fixture(root, "trees", "example_tree.nwk"),
            _fixture(root, "metadata", "example_metadata_reordered.tsv"),
            _fixture(root, "metadata", "example_traits_validate.tsv"),
            alignment_path=_fixture(root, "alignments", "example_alignment_extra_taxon.fasta"),
            tip_dates_path=_fixture(root, "metadata", "example_tip_dates_invalid.tsv"),
            calibration_path=_fixture(root, "metadata", "example_calibrations_invalid.tsv"),
        ),
    }

    results: list[RegressionDatasetCaseResult] = []
    for name, expected in expected_payload.items():
        observed = _regression_summary(current_reports[name])
        passed = expected == observed
        notes = [] if passed else ["observed regression summary drifted from the checked-in expectation"]
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


def build_method_accuracy_dashboard(*, fixtures_root: Path | None = None) -> MethodAccuracyDashboard:
    """Summarize validation accuracy, error counts, and coverage across benchmark surfaces."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    core = build_core_workflow_validation_report(fixtures_root=root)
    clean = build_clean_benchmark_corpus(fixtures_root=root)
    broken = build_broken_benchmark_corpus(fixtures_root=root)
    messy = build_messy_benchmark_corpus(fixtures_root=root)
    regression = build_regression_dataset_corpus(fixtures_root=root)

    def row(
        surface: str,
        passed_count: int,
        failed_count: int,
        coverage_count: int,
        bias_notes: list[str],
        error_notes: list[str],
    ) -> MethodAccuracyRow:
        total = max(coverage_count, 1)
        return MethodAccuracyRow(
            surface=surface,
            accuracy=round(passed_count / total, 15),
            passed_count=passed_count,
            failed_count=failed_count,
            coverage_count=coverage_count,
            bias_notes=bias_notes,
            error_notes=error_notes,
        )

    rows = [
        row(
            "level1-reference-validation",
            core.passed_fixture_count,
            core.failed_fixture_count,
            core.total_fixture_count,
            core.limitations,
            [case.fixture_name for case in core.failure_gallery if not case.passed],
        ),
        row(
            "clean-benchmark-corpus",
            clean.passed_case_count,
            clean.failed_case_count,
            clean.case_count,
            clean.limitations,
            [case.name for case in clean.cases if not case.passed],
        ),
        row(
            "broken-benchmark-corpus",
            broken.passed_case_count,
            broken.failed_case_count,
            broken.case_count,
            broken.limitations,
            [case.name for case in broken.cases if not case.passed],
        ),
        row(
            "messy-benchmark-corpus",
            messy.passed_case_count,
            messy.failed_case_count,
            messy.case_count,
            messy.limitations,
            [case.name for case in messy.cases if not case.passed],
        ),
        row(
            "regression-dataset-corpus",
            regression.passed_case_count,
            regression.failed_case_count,
            regression.case_count,
            regression.limitations,
            [case.name for case in regression.cases if not case.passed],
        ),
    ]
    return MethodAccuracyDashboard(
        goal_id=246,
        rows=rows,
        limitations=[
            "accuracy currently summarizes checked-in fixture and corpus pass rates; it does not yet replace external software comparison studies",
        ],
    )


def build_runtime_benchmark_dashboard(*, replicates: int = 1) -> RuntimeBenchmarkDashboard:
    """Summarize runtime scaling across taxa, sites, tree counts, and posterior-like samples."""
    rows = [
        BenchmarkDashboardRow(
            workflow="tree-validation",
            scaling_axis="taxa",
            observations=benchmark_tree_validation(replicates=replicates).observations,
        ),
        BenchmarkDashboardRow(
            workflow="tree-comparison",
            scaling_axis="taxa",
            observations=benchmark_tree_comparison(replicates=replicates).observations,
        ),
        BenchmarkDashboardRow(
            workflow="alignment-diagnostics",
            scaling_axis="sites",
            observations=benchmark_alignment_site_scaling(replicates=replicates).observations,
        ),
        BenchmarkDashboardRow(
            workflow="tree-set-consensus",
            scaling_axis="posterior_samples",
            observations=benchmark_tree_set_consensus(replicates=replicates).observations,
        ),
    ]
    return RuntimeBenchmarkDashboard(
        goal_id=247,
        rows=rows,
        limitations=[
            "runtime summaries measure local benchmark fixtures and should be re-run on target hardware before operational promises are made",
        ],
    )


def build_memory_benchmark_dashboard(*, replicates: int = 1) -> MemoryBenchmarkDashboard:
    """Summarize peak memory scaling across the main benchmark axes."""
    rows = [
        BenchmarkDashboardRow(
            workflow="tree-validation",
            scaling_axis="taxa",
            observations=benchmark_tree_validation(replicates=replicates).observations,
        ),
        BenchmarkDashboardRow(
            workflow="tree-comparison",
            scaling_axis="taxa",
            observations=benchmark_tree_comparison(replicates=replicates).observations,
        ),
        BenchmarkDashboardRow(
            workflow="alignment-diagnostics",
            scaling_axis="sites",
            observations=benchmark_alignment_site_scaling(replicates=replicates).observations,
        ),
        BenchmarkDashboardRow(
            workflow="tree-set-consensus",
            scaling_axis="posterior_samples",
            observations=benchmark_tree_set_consensus(replicates=replicates).observations,
        ),
    ]
    return MemoryBenchmarkDashboard(
        goal_id=248,
        rows=rows,
        limitations=[
            "memory summaries capture Python-side peak allocations during benchmark runs and do not model every external-engine workflow",
        ],
    )
