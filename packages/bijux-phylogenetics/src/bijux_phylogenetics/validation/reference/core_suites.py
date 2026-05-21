from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.core.dataset import audit_dataset_inputs
from bijux_phylogenetics.diagnostics.validation import (
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.io.fasta.quality import (
    assess_alignment_low_information,
    build_alignment_forensic_report,
    build_alignment_quality_report,
    build_ambiguous_alignment_column_report,
    build_duplicate_sequence_policy_report,
    build_sequence_quality_ranking,
)
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.taxa import (
    build_taxon_audit_report,
    build_taxon_mapping_conflict_report,
)

from .models import ReferenceValidationSuiteReport
from .shared import (
    check,
    default_fixtures_root,
    error_observation,
    fixture,
    suite_report,
)


def validate_tree_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate the core tree-diagnostic fixture corpus."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    valid_tree = fixture(root, "trees", "example_tree.nwk")
    valid_report = validate_tree_path(valid_tree)
    valid_inspection = inspect_tree_path(valid_tree)
    duplicate_tree = fixture(root, "trees", "example_tree_duplicate.nwk")
    negative_tree = fixture(root, "trees", "example_tree_negative_length.nwk")
    unrooted_tree = fixture(root, "trees", "example_tree_unrooted.nwk")

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
        check(
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
        check(
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
            else error_observation(duplicate_error) | {"message": str(duplicate_error)},
        ),
        check(
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
            else error_observation(negative_error) | {"message": str(negative_error)},
        ),
        check(
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
    return suite_report(
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
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    synonym_table = fixture(root, "metadata", "example_taxon_synonyms.tsv")
    mixed_rank_tree = fixture(root, "trees", "example_taxonomy_rank_mixed.nwk")
    taxonomy_tree = fixture(root, "trees", "example_taxonomy_tree.nwk")

    mixed_audit = build_taxon_audit_report(
        load_tree(mixed_rank_tree), synonym_table_path=synonym_table
    )
    conflict_report = build_taxon_mapping_conflict_report(
        load_tree(taxonomy_tree), synonym_table_path=synonym_table
    )

    fixtures = [
        check(
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
        check(
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
    return suite_report(
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
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    clean_alignment = fixture(root, "alignments", "example_alignment.fasta")
    duplicate_alignment = fixture(
        root, "alignments", "example_alignment_duplicates.fasta"
    )
    ambiguity_alignment = fixture(
        root, "alignments", "example_alignment_ambiguity.fasta"
    )
    missing_alignment = fixture(
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
        check(
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
        check(
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
        check(
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
        check(
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
    return suite_report(
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
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_path = fixture(root, "trees", "example_taxon_workflow_tree.nwk")
    metadata_path = fixture(root, "metadata", "example_taxon_workflow_metadata.csv")
    traits_path = fixture(root, "metadata", "example_taxon_workflow_traits.csv")
    alignment_path = fixture(
        root, "alignments", "example_taxon_workflow_alignment.fasta"
    )

    report = audit_dataset_inputs(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
    )
    fixtures = [
        check(
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
    return suite_report(
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
