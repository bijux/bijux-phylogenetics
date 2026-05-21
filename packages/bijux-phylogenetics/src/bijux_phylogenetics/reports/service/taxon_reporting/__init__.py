from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import _load_tree
from bijux_phylogenetics.phylo.taxa import build_taxon_audit_report
from bijux_phylogenetics.render.html import write_html_report

from ..artifacts import report_sidecar_path, section, write_machine_manifest
from ..ledger import build_input_ledger, serialize_input_ledger, sha256
from ..models import TaxonReportBuildResult
from .linked_evidence import build_taxon_linked_evidence


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
    title = "Bijux Taxon Audit Report"
    reviewer_summary = [
        f"taxon audit status: {audit.status}",
        f"tree tip count: {audit.tree_tip_count}",
        *audit.summary,
    ]
    if taxon_crosswalk is not None:
        reviewer_summary.append(
            f"crosswalk rows: {len(taxon_crosswalk.rows)} across linked dataset surfaces"
        )
    if taxon_exclusions is not None:
        reviewer_summary.append(
            f"excluded taxa with explicit causes: {len(taxon_exclusions.rows)}"
        )
    if taxon_workflow_loss is not None:
        reviewer_summary.append(
            f"workflow loss stages observed: {len(taxon_workflow_loss.loss_stage_counts)}"
        )
    if taxon_stability is not None:
        reviewer_summary.append(
            f"unstable taxa across linked sources: {len(taxon_stability.unstable_taxa)}"
        )
    limitations = sorted(dict.fromkeys(audit.warnings))
    sections = [
        section("reviewer-summary", reviewer_summary),
        section("taxon-audit", asdict(audit)),
        section("taxon-identity", asdict(audit.identity_audit)),
        section("taxon-safety", asdict(audit.safety_report)),
        section("taxon-namespaces", asdict(audit.namespace_report)),
        section("taxon-rank-consistency", asdict(audit.rank_consistency)),
        *(
            [section("taxon-synonyms", asdict(audit.synonym_audit))]
            if audit.synonym_audit is not None
            else []
        ),
        section("taxon-duplicate-identities", asdict(audit.duplicate_identities)),
        section("taxon-mapping-conflicts", asdict(audit.mapping_conflicts)),
        *(
            [section("taxon-accepted-names", asdict(audit.accepted_name_export))]
            if audit.accepted_name_export is not None
            else []
        ),
        *(
            [section("taxon-crosswalk", asdict(taxon_crosswalk))]
            if taxon_crosswalk is not None
            else []
        ),
        *(
            [section("taxon-exclusions", asdict(taxon_exclusions))]
            if taxon_exclusions is not None
            else []
        ),
        *(
            [section("taxon-loss", asdict(taxon_workflow_loss))]
            if taxon_workflow_loss is not None
            else []
        ),
        *(
            [
                section(
                    "taxon-loss-events",
                    [
                        {
                            "taxon": row.taxon,
                            "first_loss_stage": row.first_loss_stage,
                            "loss_events": [asdict(event) for event in row.loss_events],
                        }
                        for row in taxon_workflow_loss.rows
                        if row.loss_events
                    ],
                )
            ]
            if taxon_workflow_loss is not None
            else []
        ),
        *(
            [section("taxon-stability", asdict(taxon_stability))]
            if taxon_stability is not None
            else []
        ),
        section("limitations", limitations),
    ]
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
