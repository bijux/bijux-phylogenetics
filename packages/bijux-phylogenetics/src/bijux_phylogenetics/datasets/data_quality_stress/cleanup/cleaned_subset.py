from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import (
    TraitLinkageReport,
    TraitValidationReport,
    link_tree_to_traits,
    validate_traits_table,
    write_taxon_rows,
)
from bijux_phylogenetics.diagnostics.validation import TreeValidationReport
from bijux_phylogenetics.io.fasta import (
    AlignmentRecord,
    load_fasta_alignment,
    write_fasta_alignment,
)
from bijux_phylogenetics.io.fasta.records import FastaInputValidationReport
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.phylo.pruning import drop_tree_taxa

from ..models import CatarrhineDataQualityStressPanelDataset
from ..panel import TREE_BRANCH_FLOOR
from .shared import apply_branch_length_floor


@dataclass(slots=True)
class CleanedSubsetSurfaces:
    cleaned_traits_path: Path
    cleaned_tree_path: Path
    cleaned_alignment_path: Path
    cleaned_trait_validation: TraitValidationReport
    cleaned_tree_validation: TreeValidationReport
    cleaned_alignment_validation: FastaInputValidationReport
    cleaned_linkage: TraitLinkageReport
    cleaned_taxa: list[str]
    cleaned_alignment_records: list[AlignmentRecord]
    repaired_branch_nodes: list[str]


def assemble_cleaned_subset(
    *,
    dataset: CatarrhineDataQualityStressPanelDataset,
    assembled_root: Path,
    excluded_taxa: list[str],
    cleaned_trait_candidates: list[dict[str, str]],
    validate_cleaned_tree: Callable[[Path], TreeValidationReport],
    validate_cleaned_alignment: Callable[[Path], FastaInputValidationReport],
) -> CleanedSubsetSurfaces:
    cleaned_trait_rows = [
        dict(row)
        for row in cleaned_trait_candidates
        if row["taxon"] not in excluded_taxa
        and all(row[trait] for trait in dataset.required_traits)
    ]
    cleaned_traits_path = write_taxon_rows(
        assembled_root / "cleaned-traits.csv",
        columns=list(cleaned_trait_rows[0].keys()),
        rows=cleaned_trait_rows,
    )

    cleaned_tree, _ = drop_tree_taxa(dataset.raw_tree_path, excluded_taxa)
    repaired_branch_nodes = apply_branch_length_floor(
        cleaned_tree.root,
        floor=TREE_BRANCH_FLOOR,
    )
    cleaned_tree_path = write_newick(assembled_root / "cleaned-tree.nwk", cleaned_tree)

    cleaned_alignment_records = [
        record
        for record in load_fasta_alignment(dataset.raw_alignment_path)
        if record.identifier not in excluded_taxa
    ]
    cleaned_alignment_path = write_fasta_alignment(
        assembled_root / "cleaned-alignment.fasta",
        cleaned_alignment_records,
    )

    cleaned_trait_validation = validate_traits_table(cleaned_traits_path)
    cleaned_tree_validation = validate_cleaned_tree(cleaned_tree_path)
    cleaned_alignment_validation = validate_cleaned_alignment(cleaned_alignment_path)
    cleaned_linkage = link_tree_to_traits(
        cleaned_tree_path,
        cleaned_traits_path,
        strict=True,
    )
    return CleanedSubsetSurfaces(
        cleaned_traits_path=cleaned_traits_path,
        cleaned_tree_path=cleaned_tree_path,
        cleaned_alignment_path=cleaned_alignment_path,
        cleaned_trait_validation=cleaned_trait_validation,
        cleaned_tree_validation=cleaned_tree_validation,
        cleaned_alignment_validation=cleaned_alignment_validation,
        cleaned_linkage=cleaned_linkage,
        cleaned_taxa=sorted(cleaned_linkage.usable_taxa),
        cleaned_alignment_records=cleaned_alignment_records,
        repaired_branch_nodes=repaired_branch_nodes,
    )
