from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import _load_tree
from bijux_phylogenetics.phylo.taxa import build_taxon_audit_report
from bijux_phylogenetics.render.html import write_html_report

from ..artifacts import report_sidecar_path, write_machine_manifest
from ..ledger import build_input_ledger, serialize_input_ledger, sha256
from ..models import TaxonReportBuildResult
from .linked_evidence import build_taxon_linked_evidence
from .presentation import (
    TAXON_REPORT_TITLE,
    build_taxon_reviewer_summary,
    build_taxon_sections,
)


def render_taxon_report(
    *,
    tree_path: Path,
    out_path: Path,
    synonym_table_path: Path | None = None,
    metadata_path: Path | None = None,
    traits_path: Path | None = None,
    alignment_path: Path | None = None,
    filtered_alignment_path: Path | None = None,
    inference_tree_path: Path | None = None,
    reported_taxa_path: Path | None = None,
) -> TaxonReportBuildResult:
    """Build a reviewer-facing taxon audit report."""
    tree = _load_tree(tree_path)
    audit = build_taxon_audit_report(tree, synonym_table_path=synonym_table_path)
    linked_evidence = build_taxon_linked_evidence(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        alignment_path=alignment_path,
        filtered_alignment_path=filtered_alignment_path,
        inference_tree_path=inference_tree_path,
        reported_taxa_path=reported_taxa_path,
    )
    taxon_crosswalk = linked_evidence.taxon_crosswalk
    taxon_exclusions = linked_evidence.taxon_exclusions
    taxon_workflow_loss = linked_evidence.taxon_workflow_loss
    taxon_stability = linked_evidence.taxon_stability
    title = TAXON_REPORT_TITLE
    reviewer_summary = build_taxon_reviewer_summary(
        audit=audit,
        linked_evidence=linked_evidence,
    )
    limitations = sorted(dict.fromkeys(audit.warnings))
    sections = build_taxon_sections(
        audit=audit,
        linked_evidence=linked_evidence,
        reviewer_summary=reviewer_summary,
        limitations=limitations,
    )
    input_paths = [
        tree_path,
        *([synonym_table_path] if synonym_table_path is not None else []),
        *([metadata_path] if metadata_path is not None else []),
        *([traits_path] if traits_path is not None else []),
        *([alignment_path] if alignment_path is not None else []),
        *([filtered_alignment_path] if filtered_alignment_path is not None else []),
        *([inference_tree_path] if inference_tree_path is not None else []),
        *([reported_taxa_path] if reported_taxa_path is not None else []),
    ]
    machine_manifest = {
        "report_kind": "taxonomy",
        "title": title,
        "input_paths": [str(path) for path in input_paths],
        "input_checksums": {str(path): sha256(path) for path in input_paths},
        "sections": [name for name, _ in sections],
        "metrics": {
            "tree_tip_count": audit.tree_tip_count,
            "warning_count": len(audit.warnings),
            "conflict_count": len(audit.mapping_conflicts.rows),
            "crosswalk_rows": 0
            if taxon_crosswalk is None
            else len(taxon_crosswalk.rows),
            "excluded_taxa": 0
            if taxon_exclusions is None
            else len(taxon_exclusions.rows),
            "loss_stage_count": 0
            if taxon_workflow_loss is None
            else len(taxon_workflow_loss.loss_stage_counts),
            "unstable_taxa": 0
            if taxon_stability is None
            else len(taxon_stability.unstable_taxa),
        },
        "reviewer_summary": reviewer_summary,
        "limitations": limitations,
    }
    input_ledger_entries: list[tuple[Path, str, list[str]]] = [
        (tree_path, "tree", ["taxon_audit", "taxon_stability"]),
        *(
            [
                (
                    synonym_table_path,
                    "synonym_table",
                    ["taxon_synonyms", "accepted_name_export"],
                )
            ]
            if synonym_table_path is not None
            else []
        ),
        *(
            [(metadata_path, "metadata", ["taxon_crosswalk", "taxon_exclusions"])]
            if metadata_path is not None
            else []
        ),
        *(
            [(traits_path, "traits", ["taxon_crosswalk", "taxon_exclusions"])]
            if traits_path is not None
            else []
        ),
        *(
            [(alignment_path, "alignment", ["taxon_crosswalk", "taxon_loss"])]
            if alignment_path is not None
            else []
        ),
        *(
            [
                (
                    filtered_alignment_path,
                    "filtered_alignment",
                    ["taxon_loss", "taxon_stability"],
                )
            ]
            if filtered_alignment_path is not None
            else []
        ),
        *(
            [
                (
                    inference_tree_path,
                    "inference_tree",
                    ["taxon_loss", "taxon_stability"],
                )
            ]
            if inference_tree_path is not None
            else []
        ),
        *(
            [
                (
                    reported_taxa_path,
                    "reported_taxa",
                    ["taxon_loss", "taxon_stability"],
                )
            ]
            if reported_taxa_path is not None
            else []
        ),
    ]
    machine_manifest["input_ledger"] = serialize_input_ledger(
        build_input_ledger(input_ledger_entries)
    )
    machine_manifest_path = write_machine_manifest(
        report_sidecar_path(out_path), machine_manifest
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return TaxonReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="taxonomy",
        title=title,
        tree_path=tree_path,
        taxon_audit=audit,
        taxon_crosswalk=taxon_crosswalk,
        taxon_exclusions=taxon_exclusions,
        taxon_workflow_loss=taxon_workflow_loss,
        taxon_stability=taxon_stability,
        machine_manifest=machine_manifest,
    )
