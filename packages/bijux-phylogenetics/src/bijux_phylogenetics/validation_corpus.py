from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.dataset import DatasetAuditReport, audit_dataset_inputs


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
