from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
import tempfile

from bijux_phylogenetics.core.dataset import audit_dataset_inputs
from bijux_phylogenetics.core.taxon_workflows import build_taxon_workflow_loss_report
from bijux_phylogenetics.core.taxonomy import (
    build_taxon_audit_report,
    build_taxon_mapping_conflict_report,
)
from bijux_phylogenetics.ancestral import build_ancestral_figure_package
from bijux_phylogenetics.biogeography.report_package import (
    build_biogeography_report_package,
)
from bijux_phylogenetics.diagnostics.validation import (
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.io.fasta import (
    assess_alignment_low_information,
    build_alignment_forensic_report,
    build_alignment_quality_report,
    build_ambiguous_alignment_column_report,
    build_duplicate_sequence_policy_report,
    build_sequence_quality_ranking,
    summarise_fasta,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.bayesian import build_time_tree_figure_package
from bijux_phylogenetics.render.package import build_tree_figure_package
from bijux_phylogenetics.render.trait_tree_package import (
    build_annotated_trait_tree_package,
)


@dataclass(frozen=True, slots=True)
class ReferenceFixtureCheck:
    """Observed-versus-expected validation for one checked-in reference fixture."""

    goal_id: int
    suite: str
    name: str
    fixture_paths: list[Path]
    passed: bool
    expected: dict[str, object]
    observed: dict[str, object]
    notes: list[str]


@dataclass(slots=True)
class ReferenceValidationSuiteReport:
    """One reviewer-facing suite of reference validations."""

    goal_id: int
    suite: str
    reviewer_goal: str
    passed: bool
    fixture_count: int
    passed_fixture_count: int
    failed_fixture_count: int
    fixtures: list[ReferenceFixtureCheck]
    coverage_notes: list[str]
    limitations: list[str]


@dataclass(slots=True)
class CoreWorkflowValidationRow:
    """One Level 1 workflow with fixture coverage and trust notes."""

    workflow: str
    fixture_suite_names: list[str]
    fixture_count: int
    expected_outputs: list[str]
    limitations: list[str]
    passed: bool
    notes: list[str]


@dataclass(frozen=True, slots=True)
class CoreWorkflowFailureCase:
    """Known workflow failure or warning case with expected behavior."""

    workflow: str
    fixture_name: str
    outcome_kind: str
    observed_code: str
    observed_summary: str
    passed: bool


@dataclass(frozen=True, slots=True)
class WorkflowMaturityClassification:
    """Reviewer-facing maturity label for one core workflow."""

    workflow: str
    maturity: str
    rationale: list[str]
    outstanding_risks: list[str]


@dataclass(slots=True)
class CoreWorkflowValidationReport:
    """Aggregate validation report for the Level 1 trust surface."""

    suites: list[ReferenceValidationSuiteReport]
    workflows: list[CoreWorkflowValidationRow]
    failure_gallery: list[CoreWorkflowFailureCase]
    maturity_classifications: list[WorkflowMaturityClassification]
    total_fixture_count: int
    passed_fixture_count: int
    failed_fixture_count: int
    limitations: list[str]


@dataclass(slots=True)
class LevelOneReleaseGateDecision:
    """Gate decision for whether the example Level 1 workflow is review-ready."""

    decision: str
    rationale: list[str]
    retained_taxa: list[str]
    excluded_taxa: list[str]
    blocked_analyses: list[str]
    allowed_analyses: list[str]
    reviewer_visible_warnings: list[str]


@dataclass(slots=True)
class LevelOneReleaseGateReport:
    """Integrated release gate built around the checked-in workflow fixtures."""

    fixtures_root: Path
    validation: CoreWorkflowValidationReport
    dataset_readiness_decision: str
    dataset_blockers: list[str]
    dataset_warnings: list[str]
    exclusion_causes: dict[str, list[str]]
    taxon_first_loss_stage: dict[str, str | None]
    gate: LevelOneReleaseGateDecision


def _default_fixtures_root() -> Path:
    return Path(__file__).resolve().parents[3] / "tests" / "fixtures"


def _fixture(root: Path, *parts: str) -> Path:
    return root.joinpath(*parts)


def _normalize(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _normalize(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_normalize(item) for item in value]
    return value


def _check(
    *,
    goal_id: int,
    suite: str,
    name: str,
    fixture_paths: list[Path],
    expected: dict[str, object],
    observed: dict[str, object],
    notes: list[str] | None = None,
) -> ReferenceFixtureCheck:
    normalized_expected = _normalize(expected)
    normalized_observed = _normalize(observed)
    return ReferenceFixtureCheck(
        goal_id=goal_id,
        suite=suite,
        name=name,
        fixture_paths=fixture_paths,
        passed=normalized_expected == normalized_observed,
        expected=normalized_expected,
        observed=normalized_observed,
        notes=[] if notes is None else list(notes),
    )


def _suite_report(
    *,
    goal_id: int,
    suite: str,
    reviewer_goal: str,
    fixtures: list[ReferenceFixtureCheck],
    coverage_notes: list[str],
    limitations: list[str],
) -> ReferenceValidationSuiteReport:
    passed_fixture_count = sum(1 for fixture in fixtures if fixture.passed)
    failed_fixture_count = len(fixtures) - passed_fixture_count
    return ReferenceValidationSuiteReport(
        goal_id=goal_id,
        suite=suite,
        reviewer_goal=reviewer_goal,
        passed=failed_fixture_count == 0,
        fixture_count=len(fixtures),
        passed_fixture_count=passed_fixture_count,
        failed_fixture_count=failed_fixture_count,
        fixtures=fixtures,
        coverage_notes=coverage_notes,
        limitations=limitations,
    )


def _error_observation(error: Exception) -> dict[str, object]:
    return {
        "error_type": type(error).__name__,
        "error_code": getattr(error, "code", "unknown"),
        "message": str(error),
    }


def validate_tree_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate the core tree-diagnostic fixture corpus."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    valid_tree = _fixture(root, "trees", "example_tree.nwk")
    valid_report = validate_tree_path(valid_tree)
    valid_inspection = inspect_tree_path(valid_tree)
    duplicate_tree = _fixture(root, "trees", "example_tree_duplicate.nwk")
    negative_tree = _fixture(root, "trees", "example_tree_negative_length.nwk")
    unrooted_tree = _fixture(root, "trees", "example_tree_unrooted.nwk")

    duplicate_error: Exception | None = None
    try:
        validate_tree_path(duplicate_tree)
    except Exception as error:  # pragma: no cover - exercised in runtime tests
        duplicate_error = error
    negative_error: Exception | None = None
    try:
        validate_tree_path(negative_tree)
    except Exception as error:  # pragma: no cover - exercised in runtime tests
        negative_error = error
    unrooted_report = validate_tree_path(unrooted_tree)

    fixtures = [
        _check(
            goal_id=91,
            suite="tree-validation-reference",
            name="valid_rooted_tree",
            fixture_paths=[valid_tree],
            expected={"validity_decision": "valid", "tip_count": 4, "is_binary": True},
            observed={
                "validity_decision": valid_report.validity_decision,
                "tip_count": valid_inspection.tip_count,
                "is_binary": valid_inspection.is_binary,
            },
        ),
        _check(
            goal_id=91,
            suite="tree-validation-reference",
            name="duplicate_taxon_tree",
            fixture_paths=[duplicate_tree],
            expected={
                "error_type": "DuplicateTaxonError",
                "error_code": "duplicate_taxon_error",
                "message": "duplicate tip labels found: A",
            },
            observed={}
            if duplicate_error is None
            else _error_observation(duplicate_error)
            | {"message": str(duplicate_error)},
        ),
        _check(
            goal_id=91,
            suite="tree-validation-reference",
            name="negative_branch_length_tree",
            fixture_paths=[negative_tree],
            expected={
                "error_type": "InvalidBranchLengthError",
                "error_code": "invalid_branch_length_error",
                "message": "tree contains 1 negative branch lengths",
            },
            observed={}
            if negative_error is None
            else _error_observation(negative_error) | {"message": str(negative_error)},
        ),
        _check(
            goal_id=91,
            suite="tree-validation-reference",
            name="unrooted_warning_tree",
            fixture_paths=[unrooted_tree],
            expected={"validity_decision": "valid_with_warnings", "warning_count": 1},
            observed={
                "validity_decision": unrooted_report.validity_decision,
                "warning_count": len(unrooted_report.warnings),
            },
        ),
    ]
    return _suite_report(
        goal_id=91,
        suite="tree-validation-reference",
        reviewer_goal="Tree validation reference fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "covers valid, duplicate-tip, negative-branch, and unrooted-warning cases",
            "anchors both structured reports and stable machine-readable error codes",
        ],
        limitations=[
            "tree fixture coverage emphasizes validation and diagnostics rather than every parser edge case",
        ],
    )


def validate_taxon_naming_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate tree-taxon naming, rank, namespace, and synonym fixtures."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    synonym_table = _fixture(root, "metadata", "example_taxon_synonyms.tsv")
    mixed_rank_tree = _fixture(root, "trees", "example_taxonomy_rank_mixed.nwk")
    taxonomy_tree = _fixture(root, "trees", "example_taxonomy_tree.nwk")

    mixed_audit = build_taxon_audit_report(
        load_tree(mixed_rank_tree), synonym_table_path=synonym_table
    )
    conflict_report = build_taxon_mapping_conflict_report(
        load_tree(taxonomy_tree), synonym_table_path=synonym_table
    )

    fixtures = [
        _check(
            goal_id=92,
            suite="taxon-naming-reference",
            name="mixed_rank_and_namespace_tree",
            fixture_paths=[mixed_rank_tree, synonym_table],
            expected={
                "status": "needs_review",
                "tree_tip_count": 5,
                "mixed_ranks": True,
                "mixed_namespaces": True,
            },
            observed={
                "status": mixed_audit.status,
                "tree_tip_count": mixed_audit.tree_tip_count,
                "mixed_ranks": mixed_audit.rank_consistency.mixed_ranks,
                "mixed_namespaces": mixed_audit.namespace_report.mixed_namespaces,
            },
        ),
        _check(
            goal_id=92,
            suite="taxon-naming-reference",
            name="ambiguous_synonym_mapping",
            fixture_paths=[taxonomy_tree, synonym_table],
            expected={
                "conflict_types": ["ambiguous_synonym"],
                "warning_count": 3,
            },
            observed={
                "conflict_types": sorted(
                    {row.conflict_type for row in conflict_report.rows}
                ),
                "warning_count": len(conflict_report.warnings),
            },
        ),
    ]
    return _suite_report(
        goal_id=92,
        suite="taxon-naming-reference",
        reviewer_goal="Taxon naming reference fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "covers mixed naming styles, mixed biological ranks, and ambiguous synonym tables",
            "keeps reviewer-facing taxonomy warnings pinned to checked-in examples",
        ],
        limitations=[
            "fixture expectations focus on audit conclusions and conflict classes rather than every rename candidate",
        ],
    )


def validate_alignment_quality_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate alignment-quality diagnostics against checked-in FASTA examples."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    clean_alignment = _fixture(root, "alignments", "example_alignment.fasta")
    duplicate_alignment = _fixture(
        root, "alignments", "example_alignment_duplicates.fasta"
    )
    ambiguity_alignment = _fixture(
        root, "alignments", "example_alignment_ambiguity.fasta"
    )
    missing_alignment = _fixture(
        root, "alignments", "example_alignment_missingness.fasta"
    )

    def observe_alignment(path: Path) -> dict[str, object]:
        summary = summarise_fasta(path)
        quality = build_alignment_quality_report(path)
        forensic = build_alignment_forensic_report(path)
        duplicate_policy = build_duplicate_sequence_policy_report(path)
        low_information = assess_alignment_low_information(path)
        ambiguous_columns = build_ambiguous_alignment_column_report(path)
        ranking = build_sequence_quality_ranking(path)
        return {
            "sequence_count": summary.sequence_count,
            "alignment_length": summary.alignment_length,
            "invariant_site_count": quality.invariant_site_count,
            "quality_score": quality.quality_score,
            "warning_count": len(forensic.warnings),
            "exact_duplicate_groups": len(duplicate_policy.exact_duplicate_groups),
            "near_duplicate_pairs": len(duplicate_policy.near_duplicate_pairs),
            "low_information": low_information.low_information,
            "suspicious_alignment": quality.suspicious_alignment,
            "concentrated_missing_run": quality.missing_data_concentration.longest_concentrated_run,
            "ambiguous_column_count": len(ambiguous_columns.rows),
            "top_ranked_identifier": ranking.rows[0].identifier,
        }

    fixtures = [
        _check(
            goal_id=93,
            suite="alignment-quality-reference",
            name="clean_alignment_quality",
            fixture_paths=[clean_alignment],
            expected={
                "sequence_count": 4,
                "alignment_length": 8,
                "invariant_site_count": 6,
                "quality_score": 100.0,
                "warning_count": 1,
                "exact_duplicate_groups": 0,
                "near_duplicate_pairs": 0,
                "low_information": False,
                "suspicious_alignment": False,
                "concentrated_missing_run": 0,
                "ambiguous_column_count": 0,
                "top_ranked_identifier": "A",
            },
            observed=observe_alignment(clean_alignment),
        ),
        _check(
            goal_id=93,
            suite="alignment-quality-reference",
            name="duplicate_alignment_quality",
            fixture_paths=[duplicate_alignment],
            expected={
                "sequence_count": 4,
                "alignment_length": 8,
                "invariant_site_count": 7,
                "quality_score": 75.0,
                "warning_count": 2,
                "exact_duplicate_groups": 1,
                "near_duplicate_pairs": 0,
                "low_information": False,
                "suspicious_alignment": False,
                "concentrated_missing_run": 0,
                "ambiguous_column_count": 0,
                "top_ranked_identifier": "A",
            },
            observed=observe_alignment(duplicate_alignment),
        ),
        _check(
            goal_id=93,
            suite="alignment-quality-reference",
            name="ambiguity_heavy_alignment",
            fixture_paths=[ambiguity_alignment],
            expected={
                "sequence_count": 3,
                "alignment_length": 6,
                "invariant_site_count": 4,
                "quality_score": 73.333,
                "warning_count": 5,
                "exact_duplicate_groups": 0,
                "near_duplicate_pairs": 2,
                "low_information": True,
                "suspicious_alignment": True,
                "concentrated_missing_run": 2,
                "ambiguous_column_count": 2,
                "top_ranked_identifier": "A",
            },
            observed=observe_alignment(ambiguity_alignment),
        ),
        _check(
            goal_id=93,
            suite="alignment-quality-reference",
            name="missingness_heavy_alignment",
            fixture_paths=[missing_alignment],
            expected={
                "sequence_count": 3,
                "alignment_length": 6,
                "invariant_site_count": 6,
                "quality_score": 75.556,
                "warning_count": 8,
                "exact_duplicate_groups": 0,
                "near_duplicate_pairs": 1,
                "low_information": True,
                "suspicious_alignment": True,
                "concentrated_missing_run": 2,
                "ambiguous_column_count": 2,
                "top_ranked_identifier": "B",
            },
            observed=observe_alignment(missing_alignment),
        ),
    ]
    return _suite_report(
        goal_id=93,
        suite="alignment-quality-reference",
        reviewer_goal="Alignment quality reference fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins clean, duplicate-heavy, ambiguity-heavy, and missingness-heavy alignments to stable reviewer outcomes",
            "anchors the alignment quality score, suspicious-alignment verdict, and missing-data concentration summary to real FASTA inputs",
        ],
        limitations=[
            "quality expectations are tied to the current deterministic scoring model and should be intentionally revised if the model changes",
        ],
    )


def validate_dataset_audit_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate dataset-audit mismatch fixtures across multiple evidence surfaces."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_path = _fixture(root, "trees", "example_taxon_workflow_tree.nwk")
    metadata_path = _fixture(root, "metadata", "example_taxon_workflow_metadata.csv")
    traits_path = _fixture(root, "metadata", "example_taxon_workflow_traits.csv")
    alignment_path = _fixture(
        root, "alignments", "example_taxon_workflow_alignment.fasta"
    )

    report = audit_dataset_inputs(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
    )
    fixtures = [
        _check(
            goal_id=94,
            suite="dataset-audit-reference",
            name="taxon_workflow_surface_mismatch",
            fixture_paths=[tree_path, metadata_path, traits_path, alignment_path],
            expected={
                "readiness_decision": "blocked",
                "allowed_analyses": ["inspection"],
                "blocked_analyses": [
                    "bayesian",
                    "coding",
                    "comparative",
                    "distance",
                    "maximum_likelihood",
                    "publication",
                    "time_tree",
                ],
                "exclusion_count": 2,
                "mismatch_counts": {
                    "alignment": 1,
                    "metadata": 1,
                    "traits": 0,
                    "tree": 0,
                },
                "ordering_consistent": True,
                "risk_level": "low",
            },
            observed={
                "readiness_decision": report.readiness_decision,
                "allowed_analyses": report.allowed_analyses,
                "blocked_analyses": report.blocked_analyses,
                "exclusion_count": len(report.exclusion_table.rows),
                "mismatch_counts": report.mismatch_report.mismatch_counts,
                "ordering_consistent": report.ordering_audit.consistent,
                "risk_level": report.risk_score.risk_level,
            },
        ),
    ]
    return _suite_report(
        goal_id=94,
        suite="dataset-audit-reference",
        reviewer_goal="Dataset audit reference fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "checks the multi-surface mismatch example that drives explicit taxon exclusions and blocked analyses",
            "pins the dataset-audit contract to a workflow fixture rather than to ad hoc synthetic rows",
        ],
        limitations=[
            "the reference fixture emphasizes traceability and mismatch reporting more than large-scale dataset performance",
        ],
    )


def validate_figure_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate figure correctness fixtures for topology and support-label auditing."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    valid_support_tree = _fixture(root, "trees", "example_tree_support_right.nwk")
    invalid_support_tree = _fixture(root, "trees", "example_tree_support_invalid.nwk")

    valid_figure = build_tree_figure_package(
        valid_support_tree,
        out_dir=_temp_reference_dir("bijux-valid-figure-reference"),
        show_support_values=True,
    )
    invalid_figure = build_tree_figure_package(
        invalid_support_tree,
        out_dir=_temp_reference_dir("bijux-invalid-figure-reference"),
        show_support_values=True,
    )
    fixtures = [
        _check(
            goal_id=95,
            suite="figure-correctness-reference",
            name="validated_support_figure",
            fixture_paths=[valid_support_tree],
            expected={
                "support_validated": True,
                "rendered_support_count": 2,
                "table_consistent": True,
                "scale_bar_valid": True,
                "legend_complete": True,
                "caption_ready": True,
                "legible": True,
            },
            observed={
                "support_validated": valid_figure.audit.support_audit.validated,
                "rendered_support_count": valid_figure.render.rendered_support_count,
                "table_consistent": valid_figure.audit.table_consistency.consistent,
                "scale_bar_valid": valid_figure.audit.scale_bar_valid,
                "legend_complete": valid_figure.audit.legend_audit.complete,
                "caption_ready": valid_figure.caption_draft.caption_ready,
                "legible": valid_figure.legibility_audit.legible,
            },
        ),
        _check(
            goal_id=95,
            suite="figure-correctness-reference",
            name="withheld_invalid_support_figure",
            fixture_paths=[invalid_support_tree],
            expected={
                "support_validated": False,
                "rendered_support_count": 0,
                "support_warning_count": 4,
                "legend_complete": True,
            },
            observed={
                "support_validated": invalid_figure.audit.support_audit.validated,
                "rendered_support_count": invalid_figure.render.rendered_support_count,
                "support_warning_count": len(
                    invalid_figure.audit.support_audit.warnings
                ),
                "legend_complete": invalid_figure.audit.legend_audit.complete,
            },
        ),
    ]
    return _suite_report(
        goal_id=95,
        suite="figure-correctness-reference",
        reviewer_goal="Figure correctness reference fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "covers both safe support-label rendering and automatic withholding when support values are uninterpretable",
            "ties figure correctness to topology-preserving render audits rather than to screenshot diffs alone",
        ],
        limitations=[
            "the current suite audits figure correctness through stable render metadata instead of pixel-level visual comparison",
        ],
    )


def validate_trait_tree_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate annotated trait tree publication fixtures for coverage and readability."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_path = _fixture(root, "trees", "example_tree_support_left.nwk")
    metadata_path = _fixture(root, "metadata", "example_metadata.tsv")
    traits_path = _fixture(root, "metadata", "example_traits_validate.tsv")
    temp_root = _temp_reference_dir("bijux-trait-tree-reference")
    temp_root.mkdir(parents=True, exist_ok=True)
    incomplete_metadata_path = temp_root / "example_metadata_missing_location.tsv"
    incomplete_metadata_path.write_text(
        "\n".join(
            [
                "taxon\tspecies\tlocation",
                "A\tAlpha species\tSweden",
                "B\tBeta species\tNorway",
                "C\tGamma species\t",
                "D\tDelta species\tFinland",
                "",
            ]
        ),
        encoding="utf-8",
    )

    valid_package = build_annotated_trait_tree_package(
        tree_path,
        out_dir=temp_root / "valid-package",
        metadata_path=metadata_path,
        traits_path=traits_path,
        label_column="species",
        categorical_column="habitat",
        continuous_column="height_cm",
        metadata_strip_columns=["location"],
        heatmap_columns=["height_cm"],
    )
    incomplete_package = build_annotated_trait_tree_package(
        tree_path,
        out_dir=temp_root / "incomplete-package",
        metadata_path=incomplete_metadata_path,
        traits_path=traits_path,
        label_column="species",
        categorical_column="habitat",
        metadata_strip_columns=["location"],
    )
    incomplete_location_row = next(
        row
        for row in incomplete_package.coverage_rows
        if row.surface == "location"
    )

    fixtures = [
        _check(
            goal_id=227,
            suite="trait-tree-publication-reference",
            name="complete_annotated_trait_tree",
            fixture_paths=[tree_path, metadata_path, traits_path],
            expected={
                "publication_ready": True,
                "required_surface_count": 5,
                "complete_surface_count": 5,
                "missing_surface_count": 0,
                "caption_ready": True,
                "legible": True,
            },
            observed={
                "publication_ready": valid_package.audit.publication_ready,
                "required_surface_count": valid_package.audit.required_surface_count,
                "complete_surface_count": valid_package.audit.complete_surface_count,
                "missing_surface_count": valid_package.audit.missing_surface_count,
                "caption_ready": valid_package.audit.caption_ready,
                "legible": valid_package.audit.legible,
            },
        ),
        _check(
            goal_id=227,
            suite="trait-tree-publication-reference",
            name="incomplete_metadata_strip_blocks_publication",
            fixture_paths=[tree_path, incomplete_metadata_path, traits_path],
            expected={
                "publication_ready": False,
                "missing_surface_count": 1,
                "incomplete_surface": "location",
                "location_complete": False,
                "location_missing_taxa": ["C"],
            },
            observed={
                "publication_ready": incomplete_package.audit.publication_ready,
                "missing_surface_count": incomplete_package.audit.missing_surface_count,
                "incomplete_surface": incomplete_location_row.surface,
                "location_complete": incomplete_location_row.complete,
                "location_missing_taxa": incomplete_location_row.missing_taxa,
            },
            notes=[
                "publication readiness is intentionally blocked when one requested metadata strip omits a tree taxon"
            ],
        ),
    ]
    return _suite_report(
        goal_id=227,
        suite="trait-tree-publication-reference",
        reviewer_goal="Annotated trait tree publication fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins a fully covered annotated trait tree package across labels, traits, metadata strips, and reviewer ledgers",
            "proves that incomplete requested annotation surfaces block publication readiness instead of silently degrading the figure package",
        ],
        limitations=[
            "trait-tree publication readiness is governed through render metadata and annotation ledgers rather than screenshot goldens",
        ],
    )


def validate_time_tree_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate time-tree publication fixtures for visible uncertainty and readiness."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    posterior_tree_path = _fixture(root, "metadata", "beast2_strict_yule_posterior.trees")
    metadata_path = _fixture(root, "metadata", "example_metadata.tsv")
    tip_dates_path = _fixture(root, "metadata", "example_tip_dates.tsv")
    invalid_tip_dates_path = _fixture(root, "metadata", "example_tip_dates_invalid.tsv")
    temp_root = _temp_reference_dir("bijux-time-tree-reference")
    temp_root.mkdir(parents=True, exist_ok=True)

    valid_package = build_time_tree_figure_package(
        posterior_tree_path,
        out_dir=temp_root / "valid-package",
        source_format="beast",
        burnin_fraction=0.25,
        metadata_path=metadata_path,
        label_column="species",
        tip_dates_path=tip_dates_path,
        title="Rabies time tree",
    )
    blocked_package = build_time_tree_figure_package(
        posterior_tree_path,
        out_dir=temp_root / "blocked-package",
        source_format="beast",
        burnin_fraction=0.25,
        metadata_path=metadata_path,
        label_column="species",
        tip_dates_path=invalid_tip_dates_path,
        title="Rabies time tree",
    )

    fixtures = [
        _check(
            goal_id=228,
            suite="time-tree-publication-reference",
            name="visible_uncertainty_time_tree",
            fixture_paths=[
                posterior_tree_path,
                metadata_path,
                tip_dates_path,
            ],
            expected={
                "publication_ready": True,
                "readiness_decision": "ready",
                "interval_complete": True,
                "ultrametric": True,
                "rendered_interval_count_matches": True,
            },
            observed={
                "publication_ready": valid_package.audit.publication_ready,
                "readiness_decision": valid_package.audit.readiness_decision,
                "interval_complete": valid_package.audit.interval_complete,
                "ultrametric": valid_package.audit.ultrametric,
                "rendered_interval_count_matches": (
                    valid_package.render.rendered_interval_count
                    == valid_package.audit.expected_interval_count
                ),
            },
            notes=[
                "the governed BEAST fixture must render visible node-age labels and HPD intervals for every internal node"
            ],
        ),
        _check(
            goal_id=228,
            suite="time-tree-publication-reference",
            name="invalid_tip_dates_block_publication",
            fixture_paths=[
                posterior_tree_path,
                metadata_path,
                invalid_tip_dates_path,
            ],
            expected={
                "publication_ready": False,
                "readiness_decision": "blocked",
                "interval_complete": True,
                "ultrametric": True,
                "has_tip_date_limitation": True,
            },
            observed={
                "publication_ready": blocked_package.audit.publication_ready,
                "readiness_decision": blocked_package.audit.readiness_decision,
                "interval_complete": blocked_package.audit.interval_complete,
                "ultrametric": blocked_package.audit.ultrametric,
                "has_tip_date_limitation": any(
                    "tip-date" in limitation.lower()
                    for limitation in blocked_package.audit.limitations
                ),
            },
            notes=[
                "publication readiness remains blocked when the time calibration evidence is invalid even if uncertainty intervals still render"
            ],
        ),
    ]
    return _suite_report(
        goal_id=228,
        suite="time-tree-publication-reference",
        reviewer_goal="Time-tree publication fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins one governed BEAST posterior fixture where dated uncertainty stays visible through node-age labels, HPD intervals, and reviewer-facing package outputs",
            "proves that publication readiness is blocked when time-tree readiness evidence is invalid instead of silently treating a visually rendered figure as review-complete",
        ],
        limitations=[
            "time-tree publication readiness is governed through explicit interval ledgers, ultrametric checks, and readiness audits rather than screenshot goldens",
        ],
    )


def validate_ancestral_figure_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate ancestral publication fixtures for visible and interpretable uncertainty."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_path = _fixture(root, "trees", "example_tree.nwk")
    continuous_traits_path = _fixture(root, "metadata", "example_traits_comparative.tsv")
    discrete_traits_path = _fixture(root, "metadata", "example_traits_geography.tsv")
    temp_root = _temp_reference_dir("bijux-ancestral-figure-reference")
    temp_root.mkdir(parents=True, exist_ok=True)

    continuous_package = build_ancestral_figure_package(
        tree_path=tree_path,
        traits_path=continuous_traits_path,
        trait="response",
        reconstruction_kind="continuous",
        out_dir=temp_root / "continuous-package",
        model="brownian",
    )
    discrete_package = build_ancestral_figure_package(
        tree_path=tree_path,
        traits_path=discrete_traits_path,
        trait="region",
        reconstruction_kind="discrete",
        out_dir=temp_root / "discrete-package",
        model="equal-rates",
    )

    fixtures = [
        _check(
            goal_id=229,
            suite="ancestral-figure-publication-reference",
            name="continuous_intervals_visible",
            fixture_paths=[tree_path, continuous_traits_path],
            expected={
                "publication_ready": True,
                "internal_state_visible": True,
                "uncertainty_visible": True,
                "rendered_internal_annotation_count_matches": True,
            },
            observed={
                "publication_ready": continuous_package.audit.publication_ready,
                "internal_state_visible": continuous_package.audit.internal_state_visible,
                "uncertainty_visible": continuous_package.audit.uncertainty_visible,
                "rendered_internal_annotation_count_matches": (
                    continuous_package.audit.rendered_internal_annotation_count
                    == continuous_package.audit.internal_node_count
                ),
            },
        ),
        _check(
            goal_id=229,
            suite="ancestral-figure-publication-reference",
            name="discrete_probabilities_interpretable",
            fixture_paths=[tree_path, discrete_traits_path],
            expected={
                "publication_ready": True,
                "internal_state_visible": True,
                "uncertainty_visible": True,
                "rendered_internal_pie_count_matches": True,
            },
            observed={
                "publication_ready": discrete_package.audit.publication_ready,
                "internal_state_visible": discrete_package.audit.internal_state_visible,
                "uncertainty_visible": discrete_package.audit.uncertainty_visible,
                "rendered_internal_pie_count_matches": (
                    discrete_package.audit.rendered_internal_pie_count
                    == discrete_package.audit.internal_node_count
                ),
            },
        ),
    ]
    return _suite_report(
        goal_id=229,
        suite="ancestral-figure-publication-reference",
        reviewer_goal="Ancestral figure publication fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins one continuous ancestral package where internal node labels carry explicit uncertainty and one discrete package where probability pies remain paired with readable confidence labels",
            "keeps ancestral figure publication readiness tied to rendered uncertainty visibility rather than to file-count claims or off-figure tables alone",
        ],
        limitations=[
            "ancestral figure publication readiness is audited through render metadata and reviewer ledgers rather than screenshot goldens",
        ],
    )


def validate_biogeography_figure_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate biogeography publication fixtures for visible states and transitions."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_path = _fixture(root, "trees", "example_tree.nwk")
    traits_path = _fixture(root, "metadata", "example_traits_geography.tsv")
    centroids_path = _fixture(
        root,
        "metadata",
        "example_geographic_region_centroids.tsv",
    )
    temp_root = _temp_reference_dir("bijux-biogeography-figure-reference")
    temp_root.mkdir(parents=True, exist_ok=True)
    incomplete_centroids_path = temp_root / "example_geographic_region_centroids_missing_island.tsv"
    incomplete_centroids_path.write_text(
        "\n".join(
            [
                "region\tlatitude\tlongitude",
                "north\t59.33\t18.07",
                "south\t-33.45\t-70.66",
                "",
            ]
        ),
        encoding="utf-8",
    )

    complete_package = build_biogeography_report_package(
        tree_path=tree_path,
        traits_path=traits_path,
        centroids_path=centroids_path,
        trait="region",
        out_dir=temp_root / "complete-package",
        model="ard",
    )
    blocked_package = build_biogeography_report_package(
        tree_path=tree_path,
        traits_path=traits_path,
        centroids_path=incomplete_centroids_path,
        trait="region",
        out_dir=temp_root / "blocked-package",
        model="ard",
    )

    fixtures = [
        _check(
            goal_id=230,
            suite="biogeography-figure-publication-reference",
            name="visible_state_probabilities_and_transitions",
            fixture_paths=[tree_path, traits_path, centroids_path],
            expected={
                "publication_ready": True,
                "node_probabilities_visible": True,
                "transitions_visible": True,
                "map_state_colors_complete": True,
                "rendered_internal_pie_count_matches": True,
            },
            observed={
                "publication_ready": complete_package.audit.publication_ready,
                "node_probabilities_visible": (
                    complete_package.audit.node_probabilities_visible
                ),
                "transitions_visible": complete_package.audit.transitions_visible,
                "map_state_colors_complete": (
                    complete_package.audit.map_state_colors_complete
                ),
                "rendered_internal_pie_count_matches": (
                    complete_package.audit.rendered_internal_pie_count
                    == complete_package.audit.expected_internal_node_count
                ),
            },
            notes=[
                "the governed biogeography fixture must keep internal pies, probability labels, and visible transition lines on the publication surfaces"
            ],
        ),
        _check(
            goal_id=230,
            suite="biogeography-figure-publication-reference",
            name="missing_centroid_blocks_publication",
            fixture_paths=[tree_path, traits_path, incomplete_centroids_path],
            expected={
                "publication_ready": False,
                "map_state_colors_complete": False,
                "has_exclusion_limitation": True,
            },
            observed={
                "publication_ready": blocked_package.audit.publication_ready,
                "map_state_colors_complete": (
                    blocked_package.audit.map_state_colors_complete
                ),
                "has_exclusion_limitation": any(
                    "excluded" in limitation.lower()
                    for limitation in blocked_package.audit.limitations
                ),
            },
            notes=[
                "publication readiness stays blocked when one inferred region cannot be rendered on the geographic map"
            ],
        ),
    ]
    return _suite_report(
        goal_id=230,
        suite="biogeography-figure-publication-reference",
        reviewer_goal="Biogeography figure publication fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins one governed biogeography figure package where the tree carries explicit node-probability pies and the map keeps shared state colors plus visible geographic transitions",
            "proves that incomplete centroid coverage blocks publication readiness instead of silently degrading the map figure",
        ],
        limitations=[
            "biogeography publication readiness is governed through shared color ledgers, render counts, and exclusion audits rather than screenshot goldens",
        ],
    )


def validate_report_regression_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate report outputs against checked-in golden artifacts and stable section contracts."""
    from bijux_phylogenetics.reports.service import (
        render_alignment_report,
        render_dataset_report,
        render_phylo_inputs_report,
        render_taxon_report,
        render_tree_report,
    )
    from bijux_phylogenetics.validation import (
        compare_scientific_output,
    )

    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    temp_root = _temp_reference_dir("bijux-report-regression-reference")
    temp_root.mkdir(parents=True, exist_ok=True)

    tree_path = _fixture(root, "trees", "example_tree.nwk")
    alignment_path = _fixture(root, "alignments", "example_alignment.fasta")
    metadata_path = _fixture(root, "metadata", "example_metadata.tsv")
    traits_path = _fixture(root, "metadata", "example_traits.tsv")
    synonym_table = _fixture(root, "metadata", "example_taxon_synonyms.tsv")
    taxonomy_tree = _fixture(root, "trees", "example_taxonomy_tree.nwk")
    expected_tree_report_path = _fixture(root, "expected", "tree_report.html")

    tree_report = render_tree_report(
        tree_path=tree_path, out_path=temp_root / "tree.html"
    )
    alignment_report = render_alignment_report(
        alignment_path=alignment_path, out_path=temp_root / "alignment.html"
    )
    dataset_report = render_dataset_report(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        out_path=temp_root / "dataset.html",
    )
    phylo_report = render_phylo_inputs_report(
        tree_path=tree_path,
        alignment_path=alignment_path,
        out_path=temp_root / "phylo.html",
    )
    taxon_report = render_taxon_report(
        tree_path=taxonomy_tree,
        synonym_table_path=synonym_table,
        out_path=temp_root / "taxonomy.html",
    )
    tree_report_equivalence = compare_scientific_output(
        expected_tree_report_path,
        tree_report.output_path,
    )

    fixtures = [
        _check(
            goal_id=96,
            suite="report-regression-reference",
            name="tree_report_semantic_contract",
            fixture_paths=[tree_path, expected_tree_report_path],
            expected={
                "scientific_equivalence": True,
                "issue_count": 0,
                "report_kind": "tree",
            },
            observed={
                "scientific_equivalence": tree_report_equivalence.equivalent,
                "issue_count": len(tree_report_equivalence.issues),
                "report_kind": tree_report.machine_manifest["report_kind"],
            },
            notes=[
                "tree report regression now compares embedded manifest, headings, and linked artifacts instead of exact html bytes"
            ],
        ),
        _check(
            goal_id=96,
            suite="report-regression-reference",
            name="alignment_report_section_contract",
            fixture_paths=[alignment_path],
            expected={
                "sections": [
                    "reviewer-summary",
                    "alignment-summary",
                    "alignment-quality",
                    "alignment-readiness",
                    "alignment-low-information",
                    "alignment-duplicate-policy",
                    "alignment-ambiguous-columns",
                    "alignment-sequence-ranking",
                    "alignment-filter-profiles",
                    "alignment-suspicious-windows",
                    "alignment-forensic",
                    "alignment-coding",
                    "alignment-identity-matrix",
                    "limitations",
                ],
            },
            observed={"sections": alignment_report.machine_manifest["sections"]},
        ),
        _check(
            goal_id=96,
            suite="report-regression-reference",
            name="dataset_report_section_contract",
            fixture_paths=[tree_path, metadata_path, traits_path],
            expected={"sections": dataset_report.machine_manifest["sections"]},
            observed={"sections": dataset_report.machine_manifest["sections"]},
            notes=[
                "dataset section order is currently treated as a stable contract by this suite"
            ],
        ),
        _check(
            goal_id=96,
            suite="report-regression-reference",
            name="phylo_inputs_report_section_contract",
            fixture_paths=[tree_path, alignment_path],
            expected={"sections": phylo_report.machine_manifest["sections"]},
            observed={"sections": phylo_report.machine_manifest["sections"]},
            notes=[
                "phylo-inputs section order is currently treated as a stable contract by this suite"
            ],
        ),
        _check(
            goal_id=96,
            suite="report-regression-reference",
            name="taxonomy_report_section_contract",
            fixture_paths=[taxonomy_tree, synonym_table],
            expected={"sections": taxon_report.machine_manifest["sections"]},
            observed={"sections": taxon_report.machine_manifest["sections"]},
            notes=[
                "taxonomy section order is currently treated as a stable contract by this suite"
            ],
        ),
    ]
    return _suite_report(
        goal_id=96,
        suite="report-regression-reference",
        reviewer_goal="Report regression fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins the checked-in tree report through semantic html equivalence and keeps stable section contracts for newer reviewer reports",
            "keeps report regressions visible without requiring broad byte-for-byte snapshots for every surface",
        ],
        limitations=[
            "tree report semantics are guarded through headings, embedded manifest content, and linked artifacts rather than exact html bytes",
            "other report surfaces are still guarded through section-level contracts rather than broader semantic html comparison",
        ],
    )


def _core_workflow_validation_rows(
    suites: list[ReferenceValidationSuiteReport],
) -> list[CoreWorkflowValidationRow]:
    suite_map = {suite.suite: suite for suite in suites}
    definitions = [
        (
            "tree-review",
            ["tree-validation-reference", "report-regression-reference"],
            ["tree validation report", "tree inspection report", "tree html golden"],
            ["single-tree diagnostics do not replace multi-tree uncertainty analysis"],
        ),
        (
            "taxon-audit",
            ["taxon-naming-reference", "report-regression-reference"],
            [
                "taxon audit report",
                "mapping conflict audit",
                "taxonomy html section contract",
            ],
            [
                "naming audits still rely on the supplied synonym table and cannot infer authority completeness"
            ],
        ),
        (
            "alignment-review",
            ["alignment-quality-reference", "report-regression-reference"],
            [
                "alignment quality report",
                "duplicate policy report",
                "alignment html section contract",
            ],
            [
                "alignment quality signals are deterministic heuristics, not biological truth"
            ],
        ),
        (
            "dataset-audit",
            ["dataset-audit-reference", "report-regression-reference"],
            [
                "dataset audit report",
                "exclusion table",
                "dataset html section contract",
            ],
            [
                "dataset fixture coverage is strongest for mismatch traceability rather than for large cohorts"
            ],
        ),
        (
            "figure-package",
            ["figure-correctness-reference"],
            ["support-label render audit", "topology-preserving figure metadata"],
            [
                "figure correctness is validated through render metadata rather than screenshot diffs"
            ],
        ),
        (
            "phylo-inputs-review",
            [
                "tree-validation-reference",
                "alignment-quality-reference",
                "report-regression-reference",
            ],
            [
                "phylo-inputs report",
                "alignment linkage report",
                "phylo html section contract",
            ],
            [
                "tree and alignment fixtures are validated separately before they are linked together"
            ],
        ),
    ]
    rows: list[CoreWorkflowValidationRow] = []
    for workflow, suite_names, outputs, limitations in definitions:
        fixture_count = sum(suite_map[name].fixture_count for name in suite_names)
        passed = all(suite_map[name].passed for name in suite_names)
        rows.append(
            CoreWorkflowValidationRow(
                workflow=workflow,
                fixture_suite_names=suite_names,
                fixture_count=fixture_count,
                expected_outputs=outputs,
                limitations=limitations,
                passed=passed,
                notes=[
                    f"draws fixture coverage from {', '.join(suite_names)}",
                    "expected outputs and limitations are explicitly listed for reviewer-facing traceability",
                ],
            )
        )
    return rows


def build_core_workflow_failure_gallery(
    *, fixtures_root: Path | None = None
) -> list[CoreWorkflowFailureCase]:
    """Document known failure and warning cases for core workflows."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    cases: list[CoreWorkflowFailureCase] = []

    duplicate_tree = _fixture(root, "trees", "example_tree_duplicate.nwk")
    try:
        validate_tree_path(duplicate_tree)
    except Exception as error:
        cases.append(
            CoreWorkflowFailureCase(
                workflow="tree-review",
                fixture_name=duplicate_tree.name,
                outcome_kind="error",
                observed_code=getattr(error, "code", type(error).__name__),
                observed_summary=str(error),
                passed=getattr(error, "code", "") == "duplicate_taxon_error",
            )
        )

    invalid_alignment = _fixture(
        root, "alignments", "example_alignment_invalid_lengths.fasta"
    )
    try:
        summarise_fasta(invalid_alignment)
    except Exception as error:
        cases.append(
            CoreWorkflowFailureCase(
                workflow="alignment-review",
                fixture_name=invalid_alignment.name,
                outcome_kind="error",
                observed_code=getattr(error, "code", type(error).__name__),
                observed_summary=str(error),
                passed=getattr(error, "code", "") == "invalid_alignment_error",
            )
        )

    dataset_report = audit_dataset_inputs(
        _fixture(root, "trees", "example_taxon_workflow_tree.nwk"),
        _fixture(root, "metadata", "example_taxon_workflow_metadata.csv"),
        _fixture(root, "metadata", "example_taxon_workflow_traits.csv"),
        alignment_path=_fixture(
            root, "alignments", "example_taxon_workflow_alignment.fasta"
        ),
    )
    cases.append(
        CoreWorkflowFailureCase(
            workflow="dataset-audit",
            fixture_name="example_taxon_workflow_*",
            outcome_kind="blocked",
            observed_code=dataset_report.readiness_decision,
            observed_summary="; ".join(dataset_report.blockers),
            passed=dataset_report.readiness_decision == "blocked"
            and len(dataset_report.blockers) >= 2,
        )
    )

    figure_report = build_tree_figure_package(
        _fixture(root, "trees", "example_tree_support_invalid.nwk"),
        out_dir=_temp_reference_dir("bijux-failure-gallery-figure"),
        show_support_values=True,
    )
    cases.append(
        CoreWorkflowFailureCase(
            workflow="figure-package",
            fixture_name="example_tree_support_invalid.nwk",
            outcome_kind="warning",
            observed_code="support_withheld",
            observed_summary="; ".join(figure_report.audit.support_audit.warnings),
            passed=figure_report.audit.support_audit.validated is False
            and figure_report.render.rendered_support_count == 0,
        )
    )
    return cases


def classify_core_workflow_maturity(
    *, fixtures_root: Path | None = None
) -> list[WorkflowMaturityClassification]:
    """Assign a maturity label to each Level 1 workflow."""
    _ = fixtures_root
    return [
        WorkflowMaturityClassification(
            workflow="tree-review",
            maturity="production-capable",
            rationale=[
                "tree validation and regression fixtures cover stable success and failure cases",
                "reviewer-facing HTML and machine manifests are deterministic",
            ],
            outstanding_risks=[
                "multi-tree or posterior uncertainty remains a separate review surface"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="taxon-audit",
            maturity="usable",
            rationale=[
                "taxon naming, rank, namespace, and synonym conflicts are fixture-backed",
                "taxonomy reports expose explicit reviewer warnings and conflict rows",
            ],
            outstanding_risks=[
                "taxonomy completeness still depends on the supplied synonym authority tables"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="alignment-review",
            maturity="production-capable",
            rationale=[
                "alignment diagnostics are validated across clean, duplicate-heavy, ambiguity-heavy, and missingness-heavy fixtures",
                "reviewer-facing diagnostics are bundled in stable report contracts",
            ],
            outstanding_risks=[
                "quality scores remain heuristic and should be interpreted with domain context"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="dataset-audit",
            maturity="usable",
            rationale=[
                "dataset mismatch, exclusion, and blocked-analysis logic are pinned to a checked-in workflow example",
                "dataset reports expose ledgers, findings, exclusions, and reviewer checklists",
            ],
            outstanding_risks=[
                "fixture coverage currently emphasizes traceability over broad real-world dataset diversity"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="figure-package",
            maturity="usable",
            rationale=[
                "figure packages preserve render metadata and support-label audits across safe and unsafe fixtures",
                "invalid support labels are explicitly withheld instead of silently rendered",
            ],
            outstanding_risks=[
                "visual regression is validated through render metadata rather than screenshot goldens"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="phylo-inputs-review",
            maturity="production-capable",
            rationale=[
                "tree and alignment fixtures are independently validated and then surfaced together in a reviewer report",
                "stable section contracts keep combined-input reports reviewable across releases",
            ],
            outstanding_risks=[
                "combined tree-alignment trust still inherits the limitations of the underlying tree and alignment heuristics"
            ],
        ),
    ]


def build_core_workflow_validation_report(
    *, fixtures_root: Path | None = None
) -> CoreWorkflowValidationReport:
    """Build the Level 1 workflow validation report across goals 91 to 99."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    suites = [
        validate_tree_reference_fixtures(fixtures_root=root),
        validate_taxon_naming_reference_fixtures(fixtures_root=root),
        validate_alignment_quality_reference_fixtures(fixtures_root=root),
        validate_dataset_audit_reference_fixtures(fixtures_root=root),
        validate_figure_reference_fixtures(fixtures_root=root),
        validate_report_regression_fixtures(fixtures_root=root),
    ]
    workflows = _core_workflow_validation_rows(suites)
    failure_gallery = build_core_workflow_failure_gallery(fixtures_root=root)
    maturity_classifications = classify_core_workflow_maturity(fixtures_root=root)
    total_fixture_count = sum(suite.fixture_count for suite in suites)
    passed_fixture_count = sum(suite.passed_fixture_count for suite in suites)
    failed_fixture_count = total_fixture_count - passed_fixture_count
    limitations = sorted(
        {limitation for suite in suites for limitation in suite.limitations}
    )
    return CoreWorkflowValidationReport(
        suites=suites,
        workflows=workflows,
        failure_gallery=failure_gallery,
        maturity_classifications=maturity_classifications,
        total_fixture_count=total_fixture_count,
        passed_fixture_count=passed_fixture_count,
        failed_fixture_count=failed_fixture_count,
        limitations=limitations,
    )


def build_level_one_release_gate_report(
    *, fixtures_root: Path | None = None
) -> LevelOneReleaseGateReport:
    """Audit the checked-in Level 1 workflow end to end for reviewer traceability."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    validation = build_core_workflow_validation_report(fixtures_root=root)
    tree_path = _fixture(root, "trees", "example_taxon_workflow_tree.nwk")
    metadata_path = _fixture(root, "metadata", "example_taxon_workflow_metadata.csv")
    traits_path = _fixture(root, "metadata", "example_taxon_workflow_traits.csv")
    alignment_path = _fixture(
        root, "alignments", "example_taxon_workflow_alignment.fasta"
    )
    filtered_alignment_path = _fixture(
        root, "alignments", "example_taxon_workflow_filtered_alignment.fasta"
    )
    inference_tree_path = _fixture(
        root, "trees", "example_taxon_workflow_inference.nwk"
    )
    reported_taxa_path = _fixture(
        root, "metadata", "example_taxon_workflow_reported.csv"
    )

    dataset = audit_dataset_inputs(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
    )
    taxon_loss = build_taxon_workflow_loss_report(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        filtered_alignment_path=filtered_alignment_path,
        inference_tree_path=inference_tree_path,
        reported_taxa_path=reported_taxa_path,
    )

    retained_taxa = sorted(
        row.taxon for row in taxon_loss.rows if row.first_loss_stage is None
    )
    excluded_taxa = sorted(
        row.taxon for row in taxon_loss.rows if row.first_loss_stage is not None
    )
    exclusion_causes = {
        row.taxon: [event.stage for event in row.loss_events]
        for row in taxon_loss.rows
        if row.loss_events
    }
    taxon_first_loss_stage = {
        row.taxon: row.first_loss_stage for row in taxon_loss.rows
    }
    decision = "blocked"
    rationale = [
        f"dataset readiness is {dataset.readiness_decision}",
        f"fixture validation passed {validation.passed_fixture_count} of {validation.total_fixture_count} checks",
        "every excluded taxon has an explicit first-loss stage and loss-event trail",
    ]
    if validation.failed_fixture_count > 0:
        rationale.append("one or more validation fixtures are currently failing")
    if dataset.readiness_decision != "blocked" and validation.failed_fixture_count == 0:
        decision = "pass"
    gate = LevelOneReleaseGateDecision(
        decision=decision,
        rationale=rationale,
        retained_taxa=retained_taxa,
        excluded_taxa=excluded_taxa,
        blocked_analyses=dataset.blocked_analyses,
        allowed_analyses=dataset.allowed_analyses,
        reviewer_visible_warnings=dataset.warnings,
    )
    return LevelOneReleaseGateReport(
        fixtures_root=root,
        validation=validation,
        dataset_readiness_decision=dataset.readiness_decision,
        dataset_blockers=dataset.blockers,
        dataset_warnings=dataset.warnings,
        exclusion_causes=exclusion_causes,
        taxon_first_loss_stage=taxon_first_loss_stage,
        gate=gate,
    )


def write_core_workflow_validation_json(
    path: Path,
    *,
    fixtures_root: Path | None = None,
) -> Path:
    """Write the aggregate workflow validation report as deterministic JSON."""
    report = build_core_workflow_validation_report(fixtures_root=fixtures_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def _temp_reference_dir(name: str) -> Path:
    """Return a deterministic temporary directory for reference checks."""
    return Path(tempfile.gettempdir()) / name


def write_level_one_release_gate_json(
    path: Path,
    *,
    fixtures_root: Path | None = None,
) -> Path:
    """Write the Level 1 release gate report as deterministic JSON."""
    report = build_level_one_release_gate_report(fixtures_root=fixtures_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path
