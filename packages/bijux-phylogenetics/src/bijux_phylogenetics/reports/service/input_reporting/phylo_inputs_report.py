from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

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
from ..ledger import build_input_ledger, serialize_input_ledger
from ..models import ReportBuildResult
from ..summary import build_machine_manifest, report_summary_and_limitations


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
