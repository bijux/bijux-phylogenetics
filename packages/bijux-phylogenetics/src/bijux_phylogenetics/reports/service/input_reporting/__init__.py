from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.core.dataset import (
    audit_dataset_inputs,
    summarize_dataset_readiness,
)
from bijux_phylogenetics.datasets.study_inputs import detect_missing_trait_values
from bijux_phylogenetics.diagnostics.validation import (
    forensic_tree_path,
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.io.fasta.cleaning import (
    compute_pairwise_sequence_identity_matrix,
    list_alignment_filter_profiles,
)
from bijux_phylogenetics.io.fasta.coding import inspect_coding_alignment
from bijux_phylogenetics.io.fasta.quality import (
    assess_alignment_low_information,
    build_alignment_forensic_report,
    build_alignment_quality_report,
    build_ambiguous_alignment_column_report,
    build_duplicate_sequence_policy_report,
    build_sequence_quality_ranking,
)
from bijux_phylogenetics.io.fasta.records import link_alignment_to_tree, summarise_fasta
from bijux_phylogenetics.render.html import write_html_report

from ..artifacts import report_sidecar_path, section, write_machine_manifest
from ..ledger import build_input_ledger, serialize_input_ledger, sha256
from ..linkage import annotate_tree_against_table
from ..models import (
    AlignmentReportBuildResult,
    ReportBuildResult,
    ReportInputLedgerEntry,
)
from ..summary import build_machine_manifest, report_summary_and_limitations
from .tree_report import render_tree_report


def render_dataset_report(
    *,
    tree_path: Path,
    metadata_path: Path,
    out_path: Path,
    traits_path: Path | None = None,
    alignment_path: Path | None = None,
    tip_dates_path: Path | None = None,
    calibration_path: Path | None = None,
) -> ReportBuildResult:
    """Build the explicit tree plus table dataset report contract."""
    validation = validate_tree_path(tree_path)
    inspection = inspect_tree_path(tree_path)
    forensic = forensic_tree_path(tree_path)
    metadata_linkage = annotate_tree_against_table(tree_path, metadata_path)
    traits_linkage = (
        annotate_tree_against_table(tree_path, traits_path)
        if traits_path is not None
        else None
    )
    trait_missing_values = (
        detect_missing_trait_values(traits_path) if traits_path is not None else None
    )
    dataset_readiness = (
        summarize_dataset_readiness(tree_path, metadata_path, traits_path)
        if traits_path is not None
        else None
    )
    dataset_audit = (
        audit_dataset_inputs(
            tree_path,
            metadata_path,
            traits_path,
            alignment_path=alignment_path,
            tip_dates_path=tip_dates_path,
            calibration_path=calibration_path,
        )
        if traits_path is not None
        else None
    )
    title = "Bijux Dataset Report"
    reviewer_summary, limitations = report_summary_and_limitations(
        report_kind="dataset",
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        dataset_audit=dataset_audit,
    )
    sections = [
        section("reviewer-summary", reviewer_summary),
        section("tree-validation", asdict(validation)),
        section("tree-inspection", asdict(inspection)),
        section("tree-forensic", asdict(forensic)),
        section("metadata-linkage", asdict(metadata_linkage)),
    ]
    if traits_linkage is not None:
        sections.append(section("traits-linkage", asdict(traits_linkage)))
    if trait_missing_values is not None:
        sections.append(section("trait-missing-values", asdict(trait_missing_values)))
    if dataset_readiness is not None:
        sections.append(section("dataset-readiness", asdict(dataset_readiness)))
    if dataset_audit is not None:
        sections.append(section("dataset-audit", asdict(dataset_audit)))
        sections.append(
            section("dataset-findings", [asdict(row) for row in dataset_audit.findings])
        )
        sections.append(
            section(
                "dataset-analysis-decisions",
                [asdict(row) for row in dataset_audit.analysis_decisions],
            )
        )
        sections.append(
            section(
                "dataset-readiness-levels",
                [asdict(row) for row in dataset_audit.readiness_levels],
            )
        )
        sections.append(section("dataset-crosswalk", asdict(dataset_audit.crosswalk)))
        sections.append(
            section("dataset-completeness", asdict(dataset_audit.completeness_matrix))
        )
        sections.append(
            section("dataset-exclusions", asdict(dataset_audit.exclusion_table))
        )
        sections.append(
            section("dataset-mismatches", asdict(dataset_audit.mismatch_report))
        )
        sections.append(section("dataset-risk-score", asdict(dataset_audit.risk_score)))
        sections.append(
            section("dataset-minimal-fix-plan", asdict(dataset_audit.minimal_fix_plan))
        )
        sections.append(
            section(
                "dataset-reviewer-checklist", asdict(dataset_audit.reviewer_checklist)
            )
        )
        sections.append(
            section("dataset-ordering", asdict(dataset_audit.ordering_audit))
        )
        sections.append(
            section(
                "dataset-pruning", [asdict(row) for row in dataset_audit.pruning_steps]
            )
        )
        sections.append(
            section(
                "dataset-group-imbalance",
                [asdict(row) for row in dataset_audit.group_imbalance_warnings],
            )
        )
    input_ledger_entries: list[tuple[Path, str, list[str]]] = [
        (
            tree_path,
            "tree",
            ["tree_validation", "tree_inspection", "tree_forensic", "dataset_audit"],
        ),
        (metadata_path, "metadata", ["metadata_linkage", "dataset_audit"]),
    ]
    if traits_path is not None:
        input_ledger_entries.append(
            (
                traits_path,
                "traits",
                ["traits_linkage", "trait_missing_values", "dataset_audit"],
            )
        )
    if alignment_path is not None:
        input_ledger_entries.append(
            (alignment_path, "alignment", ["alignment_forensic", "dataset_audit"])
        )
    if tip_dates_path is not None:
        input_ledger_entries.append((tip_dates_path, "tip_dates", ["dataset_audit"]))
    if calibration_path is not None:
        input_ledger_entries.append(
            (calibration_path, "calibrations", ["dataset_audit"])
        )
    input_ledger = build_input_ledger(input_ledger_entries)
    sections.append(
        section("dataset-input-ledger", serialize_input_ledger(input_ledger))
    )
    sections.append(section("limitations", limitations))
    input_paths = [tree_path, metadata_path]
    if traits_path is not None:
        input_paths.append(traits_path)
    if alignment_path is not None:
        input_paths.append(alignment_path)
    if tip_dates_path is not None:
        input_paths.append(tip_dates_path)
    if calibration_path is not None:
        input_paths.append(calibration_path)
    machine_manifest = build_machine_manifest(
        report_kind="dataset",
        title=title,
        input_paths=input_paths,
        sections=sections,
        inspection=inspection,
    )
    machine_manifest["reviewer_summary"] = reviewer_summary
    machine_manifest["limitations"] = limitations
    machine_manifest["input_ledger"] = serialize_input_ledger(input_ledger)
    machine_manifest_path = write_machine_manifest(
        report_sidecar_path(out_path), machine_manifest
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return ReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="dataset",
        title=title,
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        metadata_linkage=metadata_linkage,
        traits_linkage=traits_linkage,
        trait_missing_values=trait_missing_values,
        alignment=None,
        alignment_quality=None,
        alignment_forensic=None,
        alignment_low_information=None,
        alignment_duplicate_policy=None,
        alignment_ambiguous_columns=None,
        alignment_sequence_ranking=None,
        alignment_coding=None,
        alignment_identity_matrix=None,
        alignment_linkage=None,
        dataset_readiness=dataset_readiness,
        dataset_audit=dataset_audit,
        input_ledger=input_ledger,
        machine_manifest=machine_manifest,
    )


def render_phylo_inputs_report(
    *,
    tree_path: Path,
    alignment_path: Path,
    out_path: Path,
) -> ReportBuildResult:
    """Build the explicit tree plus alignment input report contract."""
    validation = validate_tree_path(tree_path)
    inspection = inspect_tree_path(tree_path)
    forensic = forensic_tree_path(tree_path)
    alignment = summarise_fasta(alignment_path)
    alignment_quality = build_alignment_quality_report(alignment_path)
    alignment_forensic = build_alignment_forensic_report(alignment_path)
    alignment_low_information = assess_alignment_low_information(alignment_path)
    alignment_duplicate_policy = build_duplicate_sequence_policy_report(alignment_path)
    alignment_ambiguous_columns = build_ambiguous_alignment_column_report(
        alignment_path
    )
    alignment_sequence_ranking = build_sequence_quality_ranking(alignment_path)
    alignment_coding = (
        inspect_coding_alignment(alignment_path)
        if alignment.inferred_alphabet in {"dna", "rna"}
        else None
    )
    alignment_identity_matrix = compute_pairwise_sequence_identity_matrix(
        alignment_path
    )
    alignment_linkage = link_alignment_to_tree(tree_path, alignment_path)
    title = "Bijux Phylo Inputs Report"
    reviewer_summary, limitations = report_summary_and_limitations(
        report_kind="phylo-inputs",
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        alignment_forensic=alignment_forensic,
    )
    sections = [
        section("reviewer-summary", reviewer_summary),
        section("tree-validation", asdict(validation)),
        section("tree-inspection", asdict(inspection)),
        section("tree-forensic", asdict(forensic)),
        section("alignment-summary", asdict(alignment)),
        section("alignment-quality", asdict(alignment_quality)),
        section("alignment-low-information", asdict(alignment_low_information)),
        section("alignment-duplicate-policy", asdict(alignment_duplicate_policy)),
        section("alignment-ambiguous-columns", asdict(alignment_ambiguous_columns)),
        section("alignment-sequence-ranking", asdict(alignment_sequence_ranking)),
        section("alignment-forensic", asdict(alignment_forensic)),
        *(
            [section("alignment-coding", asdict(alignment_coding))]
            if alignment_coding is not None
            else []
        ),
        section("alignment-identity-matrix", asdict(alignment_identity_matrix)),
        section("alignment-linkage", asdict(alignment_linkage)),
        section("limitations", limitations),
    ]
    input_ledger = build_input_ledger(
        [
            (
                tree_path,
                "tree",
                [
                    "tree_validation",
                    "tree_inspection",
                    "tree_forensic",
                    "alignment_linkage",
                ],
            ),
            (
                alignment_path,
                "alignment",
                [
                    "alignment_summary",
                    "alignment_quality",
                    "alignment_low_information",
                    "alignment_duplicate_policy",
                    "alignment_ambiguous_columns",
                    "alignment_sequence_ranking",
                    "alignment_forensic",
                    "alignment_identity_matrix",
                    "alignment_linkage",
                ],
            ),
        ]
    )
    machine_manifest = build_machine_manifest(
        report_kind="phylo-inputs",
        title=title,
        input_paths=[tree_path, alignment_path],
        sections=sections,
        inspection=inspection,
    )
    machine_manifest["reviewer_summary"] = reviewer_summary
    machine_manifest["limitations"] = limitations
    machine_manifest["input_ledger"] = serialize_input_ledger(input_ledger)
    machine_manifest_path = write_machine_manifest(
        report_sidecar_path(out_path), machine_manifest
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return ReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="phylo-inputs",
        title=title,
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        metadata_linkage=None,
        traits_linkage=None,
        trait_missing_values=None,
        alignment=alignment,
        alignment_quality=alignment_quality,
        alignment_forensic=alignment_forensic,
        alignment_low_information=alignment_low_information,
        alignment_duplicate_policy=alignment_duplicate_policy,
        alignment_ambiguous_columns=alignment_ambiguous_columns,
        alignment_sequence_ranking=alignment_sequence_ranking,
        alignment_coding=alignment_coding,
        alignment_identity_matrix=alignment_identity_matrix,
        alignment_linkage=alignment_linkage,
        dataset_readiness=None,
        dataset_audit=None,
        input_ledger=input_ledger,
        machine_manifest=machine_manifest,
    )


def render_alignment_report(
    *, alignment_path: Path, out_path: Path
) -> AlignmentReportBuildResult:
    """Build a reviewer-facing alignment-only report."""
    alignment = summarise_fasta(alignment_path)
    alignment_quality = build_alignment_quality_report(alignment_path)
    alignment_low_information = assess_alignment_low_information(alignment_path)
    alignment_duplicate_policy = build_duplicate_sequence_policy_report(alignment_path)
    alignment_ambiguous_columns = build_ambiguous_alignment_column_report(
        alignment_path
    )
    alignment_sequence_ranking = build_sequence_quality_ranking(alignment_path)
    alignment_forensic = build_alignment_forensic_report(alignment_path)
    alignment_coding = (
        inspect_coding_alignment(alignment_path)
        if alignment.inferred_alphabet in {"dna", "rna"}
        else None
    )
    alignment_identity_matrix = compute_pairwise_sequence_identity_matrix(
        alignment_path
    )
    title = "Bijux Alignment Report"
    reviewer_summary = [
        f"alignment quality score: {alignment_quality.quality_score}",
        (
            "alignment suspicious diagnostics: flagged"
            if alignment_quality.suspicious_alignment
            else "alignment suspicious diagnostics: clear"
        ),
        (
            "alignment remains suitable for at least one inference family"
            if any(
                (
                    alignment_forensic.safe_for_distance_analysis,
                    alignment_forensic.safe_for_maximum_likelihood,
                    alignment_forensic.safe_for_bayesian_inference,
                    alignment_forensic.safe_for_coding_analysis,
                )
            )
            else "alignment is currently blocked for the main inference families reviewed here"
        ),
        f"reviewer-facing warnings: {len(alignment_forensic.warnings)}",
    ]
    if alignment_quality.suspicious_alignment:
        reviewer_summary.append(
            "longest concentrated missing-data run: "
            f"{alignment_quality.missing_data_concentration.longest_concentrated_run}"
        )
    limitations = sorted(
        dict.fromkeys([*alignment_forensic.limitations, *alignment_forensic.warnings])
    )
    sections = [
        section("reviewer-summary", reviewer_summary),
        section("alignment-summary", asdict(alignment)),
        section("alignment-quality", asdict(alignment_quality)),
        section("alignment-readiness", asdict(alignment_forensic.readiness)),
        section("alignment-low-information", asdict(alignment_low_information)),
        section("alignment-duplicate-policy", asdict(alignment_duplicate_policy)),
        section("alignment-ambiguous-columns", asdict(alignment_ambiguous_columns)),
        section("alignment-sequence-ranking", asdict(alignment_sequence_ranking)),
        section(
            "alignment-filter-profiles",
            [asdict(profile) for profile in list_alignment_filter_profiles()],
        ),
        section(
            "alignment-suspicious-windows",
            {
                "over_aligned_regions": [
                    asdict(row) for row in alignment_forensic.over_aligned_regions
                ],
                "under_aligned_regions": [
                    asdict(row) for row in alignment_forensic.under_aligned_regions
                ],
            },
        ),
        section("alignment-forensic", asdict(alignment_forensic)),
        *(
            [section("alignment-coding", asdict(alignment_coding))]
            if alignment_coding is not None
            else []
        ),
        section("alignment-identity-matrix", asdict(alignment_identity_matrix)),
        section("limitations", limitations),
    ]
    machine_manifest = {
        "report_kind": "alignment",
        "title": title,
        "input_paths": [str(alignment_path)],
        "input_checksums": {str(alignment_path): sha256(alignment_path)},
        "sections": [name for name, _ in sections],
        "metrics": {
            "sequence_count": alignment.sequence_count,
            "alignment_length": alignment.alignment_length,
            "quality_score": alignment_quality.quality_score,
        },
        "reviewer_summary": reviewer_summary,
        "limitations": limitations,
    }
    machine_manifest_path = write_machine_manifest(
        report_sidecar_path(out_path), machine_manifest
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return AlignmentReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="alignment",
        title=title,
        alignment=alignment,
        alignment_quality=alignment_quality,
        alignment_forensic=alignment_forensic,
        alignment_low_information=alignment_low_information,
        alignment_duplicate_policy=alignment_duplicate_policy,
        alignment_ambiguous_columns=alignment_ambiguous_columns,
        alignment_sequence_ranking=alignment_sequence_ranking,
        alignment_coding=alignment_coding,
        alignment_identity_matrix=alignment_identity_matrix,
        machine_manifest=machine_manifest,
    )


def render_phylogenetics_report(
    *,
    tree_path: Path,
    out_path: Path,
    alignment_path: Path | None = None,
    traits_path: Path | None = None,
    metadata_path: Path | None = None,
) -> ReportBuildResult:
    """Build an HTML report around a tree and optional evidence tables."""
    validation = validate_tree_path(tree_path)
    inspection = inspect_tree_path(tree_path)
    forensic = forensic_tree_path(tree_path)
    alignment = summarise_fasta(alignment_path) if alignment_path else None
    alignment_quality = (
        build_alignment_quality_report(alignment_path) if alignment_path else None
    )
    alignment_forensic = (
        build_alignment_forensic_report(alignment_path) if alignment_path else None
    )
    alignment_coding = (
        inspect_coding_alignment(alignment_path)
        if alignment_path is not None
        and alignment is not None
        and alignment.inferred_alphabet in {"dna", "rna"}
        else None
    )
    alignment_identity_matrix = (
        compute_pairwise_sequence_identity_matrix(alignment_path)
        if alignment_path is not None
        else None
    )
    traits_linkage = (
        annotate_tree_against_table(tree_path, traits_path) if traits_path else None
    )
    trait_missing_values = (
        detect_missing_trait_values(traits_path) if traits_path else None
    )
    metadata_linkage = (
        annotate_tree_against_table(tree_path, metadata_path) if metadata_path else None
    )
    dataset_readiness = (
        summarize_dataset_readiness(tree_path, metadata_path, traits_path)
        if traits_path and metadata_path
        else None
    )
    dataset_audit = (
        audit_dataset_inputs(
            tree_path,
            metadata_path,
            traits_path,
            alignment_path=alignment_path,
        )
        if traits_path and metadata_path
        else None
    )
    reviewer_summary, limitations = report_summary_and_limitations(
        report_kind="dataset" if dataset_audit is not None else "tree",
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        dataset_audit=dataset_audit,
        alignment_forensic=alignment_forensic,
    )

    sections = [
        section("reviewer-summary", reviewer_summary),
        section("tree-validation", asdict(validation)),
        section("tree-inspection", asdict(inspection)),
        section("tree-forensic", asdict(forensic)),
    ]
    if alignment is not None:
        sections.append(section("alignment-summary", asdict(alignment)))
    if alignment_quality is not None:
        sections.append(section("alignment-quality", asdict(alignment_quality)))
    if alignment_forensic is not None:
        sections.append(section("alignment-forensic", asdict(alignment_forensic)))
    if alignment_coding is not None:
        sections.append(section("alignment-coding", asdict(alignment_coding)))
    if alignment_identity_matrix is not None:
        sections.append(
            section("alignment-identity-matrix", asdict(alignment_identity_matrix))
        )
    if traits_linkage is not None:
        sections.append(section("traits-linkage", asdict(traits_linkage)))
    if trait_missing_values is not None:
        sections.append(section("trait-missing-values", asdict(trait_missing_values)))
    if metadata_linkage is not None:
        sections.append(section("metadata-linkage", asdict(metadata_linkage)))
    if dataset_readiness is not None:
        sections.append(section("dataset-readiness", asdict(dataset_readiness)))
    if dataset_audit is not None:
        sections.append(section("dataset-audit", asdict(dataset_audit)))
        sections.append(section("dataset-crosswalk", asdict(dataset_audit.crosswalk)))
        sections.append(
            section("dataset-completeness", asdict(dataset_audit.completeness_matrix))
        )
    sections.append(section("limitations", limitations))

    title = "Bijux Phylogenetics Report"
    input_paths = [tree_path]
    if alignment_path is not None:
        input_paths.append(alignment_path)
    if traits_path is not None:
        input_paths.append(traits_path)
    if metadata_path is not None:
        input_paths.append(metadata_path)
    machine_manifest = build_machine_manifest(
        report_kind="phylogenetics",
        title=title,
        input_paths=input_paths,
        sections=sections,
        inspection=inspection,
    )
    machine_manifest["reviewer_summary"] = reviewer_summary
    machine_manifest["limitations"] = limitations
    input_ledger: list[ReportInputLedgerEntry] = []
    machine_manifest["input_ledger"] = serialize_input_ledger(input_ledger)
    machine_manifest_path = write_machine_manifest(
        report_sidecar_path(out_path), machine_manifest
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return ReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="phylogenetics",
        title=title,
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        metadata_linkage=metadata_linkage,
        traits_linkage=traits_linkage,
        trait_missing_values=trait_missing_values,
        alignment=alignment,
        alignment_quality=alignment_quality,
        alignment_forensic=alignment_forensic,
        alignment_low_information=None,
        alignment_duplicate_policy=None,
        alignment_ambiguous_columns=None,
        alignment_sequence_ranking=None,
        alignment_coding=alignment_coding,
        alignment_identity_matrix=alignment_identity_matrix,
        alignment_linkage=None,
        dataset_readiness=dataset_readiness,
        dataset_audit=dataset_audit,
        input_ledger=input_ledger,
        machine_manifest=machine_manifest,
    )
