from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.core.dataset import DatasetAuditReport, audit_dataset_inputs
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .contracts import BenchmarkCorpusReport, CorpusDatasetCase, CorpusDatasetCaseResult


def default_fixtures_root() -> Path:
    return Path(__file__).resolve().parents[3] / "tests" / "fixtures"


def fixture(root: Path, *parts: str) -> Path:
    return root.joinpath(*parts)


def evaluate_dataset_case(
    case: CorpusDatasetCase,
) -> tuple[DatasetAuditReport, list[str]]:
    report = audit_dataset_inputs(
        case.tree_path,
        case.metadata_path,
        case.traits_path,
        alignment_path=case.alignment_path,
        tip_dates_path=case.tip_dates_path,
        calibration_path=case.calibration_path,
    )
    notes: list[str] = []
    missing_allowed = sorted(
        set(case.required_allowed_analyses) - set(report.allowed_analyses)
    )
    if missing_allowed:
        notes.append(f"missing required allowed analyses: {', '.join(missing_allowed)}")
    forbidden = [
        blocker for blocker in report.blockers if blocker in case.forbidden_blockers
    ]
    if forbidden:
        notes.append(f"observed forbidden blockers: {', '.join(forbidden)}")
    disallowed_warnings = [
        warning
        for warning in report.warnings
        if case.allowed_warning_prefixes
        and not any(
            warning.startswith(prefix) for prefix in case.allowed_warning_prefixes
        )
    ]
    if disallowed_warnings:
        notes.append(f"unexpected warnings: {', '.join(disallowed_warnings)}")
    return report, notes


def build_clean_benchmark_corpus(
    *, fixtures_root: Path | None = None
) -> BenchmarkCorpusReport:
    """Validate that the checked-in clean corpus stays usable for core workflows."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    cases = [
        CorpusDatasetCase(
            name="core_inference_ready_dataset",
            tree_path=fixture(root, "trees", "example_tree.nwk"),
            metadata_path=fixture(root, "metadata", "example_metadata.tsv"),
            traits_path=fixture(root, "metadata", "example_traits_validate.tsv"),
            alignment_path=fixture(root, "alignments", "example_alignment.fasta"),
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
            allowed_warning_prefixes=("equal-length ungapped FASTA may be aligned",),
        )
    ]

    results: list[CorpusDatasetCaseResult] = []
    for case in cases:
        report, notes = evaluate_dataset_case(case)
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


def build_broken_benchmark_corpus(
    *, fixtures_root: Path | None = None
) -> BenchmarkCorpusReport:
    """Validate that intentionally broken fixtures still fail for the expected reason."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    results: list[CorpusDatasetCaseResult] = []

    def record_error_case(name: str, callback, expected_code: str) -> None:
        observed_code: str | None = None
        notes: list[str] = []
        try:
            callback()
            notes.append(
                f"expected error code {expected_code} but the fixture completed successfully"
            )
        except PhylogeneticsError as error:
            observed_code = error.code
            if observed_code != expected_code:
                notes.append(
                    f"expected error code {expected_code} but observed {observed_code}"
                )
        results.append(
            CorpusDatasetCaseResult(
                name=name,
                passed=not notes,
                readiness_decision="error"
                if observed_code is not None
                else "unexpected_success",
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
        lambda: validate_tree_path(
            fixture(root, "trees", "example_tree_duplicate.nwk")
        ),
        "duplicate_taxon_error",
    )
    record_error_case(
        "invalid_alignment_lengths",
        lambda: summarise_fasta(
            fixture(root, "alignments", "example_alignment_invalid_lengths.fasta")
        ),
        "invalid_alignment_error",
    )

    broken_dataset = CorpusDatasetCase(
        name="dataset_missing_metadata_taxon",
        tree_path=fixture(root, "trees", "example_tree.nwk"),
        metadata_path=fixture(root, "metadata", "example_metadata_missing_taxon.csv"),
        traits_path=fixture(root, "metadata", "example_traits_validate.tsv"),
        alignment_path=fixture(root, "alignments", "example_alignment.fasta"),
    )
    report, notes = evaluate_dataset_case(broken_dataset)
    if report.readiness_decision != "blocked":
        notes.append(
            f"expected blocked readiness but observed {report.readiness_decision}"
        )
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
            observed_code="dataset_blocked"
            if report.readiness_decision == "blocked"
            else None,
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


def build_messy_benchmark_corpus(
    *, fixtures_root: Path | None = None
) -> BenchmarkCorpusReport:
    """Validate that realistic multi-problem datasets surface the expected warning mix."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    cases: list[tuple[CorpusDatasetCase, tuple[str, ...], tuple[str, ...]]] = [
        (
            CorpusDatasetCase(
                name="reordered_alignment_extra_taxa_invalid_dates_and_calibrations",
                tree_path=fixture(root, "trees", "example_tree.nwk"),
                metadata_path=fixture(
                    root, "metadata", "example_metadata_reordered.tsv"
                ),
                traits_path=fixture(root, "metadata", "example_traits_validate.tsv"),
                alignment_path=fixture(
                    root, "alignments", "example_alignment_extra_taxon.fasta"
                ),
                tip_dates_path=fixture(
                    root, "metadata", "example_tip_dates_invalid.tsv"
                ),
                calibration_path=fixture(
                    root, "metadata", "example_calibrations_invalid.tsv"
                ),
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
                tree_path=fixture(root, "trees", "example_tree.nwk"),
                metadata_path=fixture(root, "metadata", "example_metadata.tsv"),
                traits_path=fixture(root, "metadata", "example_traits.tsv"),
                alignment_path=fixture(
                    root, "alignments", "example_alignment_ambiguity.fasta"
                ),
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
        report, notes = evaluate_dataset_case(case)
        missing_blockers = [
            message for message in expected_blockers if message not in report.blockers
        ]
        if missing_blockers:
            notes.append(f"missing expected blockers: {', '.join(missing_blockers)}")
        missing_warnings = [
            message for message in expected_warnings if message not in report.warnings
        ]
        if missing_warnings:
            notes.append(f"missing expected warnings: {', '.join(missing_warnings)}")
        if report.readiness_decision != "blocked":
            notes.append(
                f"expected blocked readiness but observed {report.readiness_decision}"
            )
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
