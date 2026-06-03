from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.diagnostics.validation.structure import _load_tree
from bijux_phylogenetics.phylo.taxa import build_taxon_audit_report
from bijux_phylogenetics.render.html import write_html_report

from ..artifacts import report_sidecar_path, write_machine_manifest
from ..models import TaxonReportBuildResult
from .linked_evidence import build_taxon_linked_evidence
from .machine_manifest import build_taxon_machine_manifest
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
    machine_manifest = build_taxon_machine_manifest(
        audit=audit,
        linked_evidence=linked_evidence,
        tree_path=tree_path,
        synonym_table_path=synonym_table_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        alignment_path=alignment_path,
        filtered_alignment_path=filtered_alignment_path,
        inference_tree_path=inference_tree_path,
        reported_taxa_path=reported_taxa_path,
        title=title,
        sections=sections,
        reviewer_summary=reviewer_summary,
        limitations=limitations,
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
        taxon_crosswalk=linked_evidence.taxon_crosswalk,
        taxon_exclusions=linked_evidence.taxon_exclusions,
        taxon_workflow_loss=linked_evidence.taxon_workflow_loss,
        taxon_stability=linked_evidence.taxon_stability,
        machine_manifest=machine_manifest,
    )
