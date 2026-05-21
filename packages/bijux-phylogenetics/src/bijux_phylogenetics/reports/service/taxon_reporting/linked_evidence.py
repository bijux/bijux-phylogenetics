from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from bijux_phylogenetics.core.dataset import audit_dataset_inputs
from bijux_phylogenetics.phylo.taxa import (
    build_taxon_stability_report,
    build_taxon_workflow_loss_report,
    load_taxon_run_source,
)


@dataclass(frozen=True)
class TaxonLinkedEvidence:
    taxon_crosswalk: Any | None
    taxon_exclusions: Any | None
    taxon_workflow_loss: Any | None
    taxon_stability: Any | None


def build_taxon_linked_evidence(
    *,
    tree_path: Path,
    metadata_path: Path | None,
    traits_path: Path | None,
    alignment_path: Path | None,
    filtered_alignment_path: Path | None,
    inference_tree_path: Path | None,
    reported_taxa_path: Path | None,
) -> TaxonLinkedEvidence:
    """Load linked taxon evidence across dataset and workflow surfaces."""
    dataset_audit = (
        None
        if metadata_path is None or traits_path is None
        else audit_dataset_inputs(
            tree_path,
            metadata_path,
            traits_path,
            alignment_path=alignment_path,
        )
    )
    taxon_workflow_loss = (
        None
        if metadata_path is None or traits_path is None
        else build_taxon_workflow_loss_report(
            tree_path,
            metadata_path,
            traits_path,
            alignment_path=alignment_path,
            filtered_alignment_path=filtered_alignment_path,
            inference_tree_path=inference_tree_path,
            reported_taxa_path=reported_taxa_path,
        )
    )
    stability_sources = [
        load_taxon_run_source(label="tree", path=tree_path),
        *(
            [load_taxon_run_source(label="metadata", path=metadata_path)]
            if metadata_path is not None
            else []
        ),
        *(
            [load_taxon_run_source(label="traits", path=traits_path)]
            if traits_path is not None
            else []
        ),
        *(
            [load_taxon_run_source(label="alignment", path=alignment_path)]
            if alignment_path is not None
            else []
        ),
        *(
            [
                load_taxon_run_source(
                    label="filtered_alignment", path=filtered_alignment_path
                )
            ]
            if filtered_alignment_path is not None
            else []
        ),
        *(
            [load_taxon_run_source(label="inference_tree", path=inference_tree_path)]
            if inference_tree_path is not None
            else []
        ),
        *(
            [load_taxon_run_source(label="reported_taxa", path=reported_taxa_path)]
            if reported_taxa_path is not None
            else []
        ),
    ]
    taxon_stability = (
        build_taxon_stability_report(stability_sources)
        if len(stability_sources) >= 2
        else None
    )
    return TaxonLinkedEvidence(
        taxon_crosswalk=None if dataset_audit is None else dataset_audit.crosswalk,
        taxon_exclusions=None
        if dataset_audit is None
        else dataset_audit.exclusion_table,
        taxon_workflow_loss=taxon_workflow_loss,
        taxon_stability=taxon_stability,
    )
