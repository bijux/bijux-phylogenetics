from __future__ import annotations

from pathlib import Path
from typing import Any

from ..ledger import build_input_ledger, serialize_input_ledger, sha256
from .linked_evidence import TaxonLinkedEvidence


def build_taxon_machine_manifest(
    *,
    audit: Any,
    linked_evidence: TaxonLinkedEvidence,
    tree_path: Path,
    synonym_table_path: Path | None,
    metadata_path: Path | None,
    traits_path: Path | None,
    alignment_path: Path | None,
    filtered_alignment_path: Path | None,
    inference_tree_path: Path | None,
    reported_taxa_path: Path | None,
    title: str,
    sections: list[tuple[str, Any]],
    reviewer_summary: list[str],
    limitations: list[str],
) -> dict[str, Any]:
    """Build the taxon report machine manifest and input ledger."""
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
            if linked_evidence.taxon_crosswalk is None
            else len(linked_evidence.taxon_crosswalk.rows),
            "excluded_taxa": 0
            if linked_evidence.taxon_exclusions is None
            else len(linked_evidence.taxon_exclusions.rows),
            "loss_stage_count": 0
            if linked_evidence.taxon_workflow_loss is None
            else len(linked_evidence.taxon_workflow_loss.loss_stage_counts),
            "unstable_taxa": 0
            if linked_evidence.taxon_stability is None
            else len(linked_evidence.taxon_stability.unstable_taxa),
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
    return machine_manifest
