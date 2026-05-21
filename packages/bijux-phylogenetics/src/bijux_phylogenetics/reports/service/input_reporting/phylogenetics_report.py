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
)
from bijux_phylogenetics.io.fasta.coding import inspect_coding_alignment
from bijux_phylogenetics.io.fasta.quality import (
    build_alignment_forensic_report,
    build_alignment_quality_report,
)
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.render.html import write_html_report

from ..artifacts import report_sidecar_path, section, write_machine_manifest
from ..ledger import serialize_input_ledger
from ..linkage import annotate_tree_against_table
from ..models import ReportBuildResult, ReportInputLedgerEntry
from ..summary import build_machine_manifest, report_summary_and_limitations


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
