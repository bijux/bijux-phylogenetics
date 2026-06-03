from __future__ import annotations

from pathlib import Path
import tempfile

from Bio import Phylo

from bijux_phylogenetics.ancestral.common import (
    load_continuous_dataset,
    load_discrete_dataset,
)
from bijux_phylogenetics.io.biopython import tree_from_biophylo
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import detect_tree_format
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    AncestralReconstructionError,
    InvalidAlignmentError,
)
from bijux_phylogenetics.trees import (
    TreeSetReport,
    load_tree_set,
)

from .models import AncestralTreeSetExclusion


def validate_burnin_fraction(burnin_fraction: float) -> None:
    """Validate the configured tree-set burn-in fraction."""
    if not 0.0 <= burnin_fraction < 1.0:
        raise AncestralReconstructionError(
            "ancestral tree-set burnin fraction must be between 0 inclusive and 1 exclusive"
        )


def load_tree_set_trees(path: Path) -> tuple[str, list[PhyloTree]]:
    """Load tree-set records as native phylogenetic trees."""
    if not path.exists():
        raise FileNotFoundError(f"tree-set file not found: {path}")
    source_format = detect_tree_format(path)
    bio_trees = list(Phylo.parse(path, source_format))
    if not bio_trees:
        raise InvalidAlignmentError(f"tree set contains no trees: {path}")
    return source_format, [
        tree_from_biophylo(tree, source_format=source_format) for tree in bio_trees
    ]


def shared_taxa(trees: list[PhyloTree]) -> set[str]:
    """Return taxa shared across every retained tree."""
    shared = set(trees[0].tip_names)
    for tree in trees[1:]:
        shared &= set(tree.tip_names)
    return shared


def prune_tree_to_taxa(
    tree: PhyloTree,
    requested_taxa: list[str],
    *,
    scratch_path: Path,
) -> PhyloTree:
    """Prune one tree down to the requested taxa using the canonical workflow."""
    scratch_path.write_text(dumps_newick(tree) + "\n", encoding="utf-8")
    pruned_tree, _report = prune_tree_to_requested_taxa(scratch_path, requested_taxa)
    return pruned_tree


def prepare_analysis_tree_set(
    *,
    traits_path: Path,
    taxon_column: str | None,
    trait: str,
    kept_tree_entries: list[tuple[int, PhyloTree]],
    shared_tree_taxa: list[str],
    dataset_kind: str,
) -> tuple[
    list[tuple[int, PhyloTree]],
    TreeSetReport,
    list[str],
    list[AncestralTreeSetExclusion],
    list[str],
    str,
]:
    """Build the pruned analysis tree set and aligned dataset view used by reviews."""
    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-ancestral-tree-set-prepare-"
    ) as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        reference_tree_path = tmp_dir_path / "ancestral-tree-set-reference.nwk"
        reference_tree = prune_tree_to_taxa(
            kept_tree_entries[0][1],
            shared_tree_taxa,
            scratch_path=reference_tree_path,
        )
        reference_tree_path.write_text(
            dumps_newick(reference_tree) + "\n",
            encoding="utf-8",
        )
        if dataset_kind == "continuous":
            dataset = load_continuous_dataset(
                reference_tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
            )
            exclusions = [
                *[
                    AncestralTreeSetExclusion(
                        taxon=taxon,
                        reason="missing_trait_value",
                    )
                    for taxon in dataset.dropped_missing_taxa
                ],
                *[
                    AncestralTreeSetExclusion(
                        taxon=taxon,
                        reason="non_numeric_trait_value",
                    )
                    for taxon in dataset.dropped_non_numeric_taxa
                ],
            ]
        else:
            dataset = load_discrete_dataset(
                reference_tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
            )
            exclusions = [
                AncestralTreeSetExclusion(
                    taxon=taxon,
                    reason="missing_discrete_trait_state",
                )
                for taxon in dataset.dropped_missing_taxa
            ]
        analysis_taxa = list(dataset.taxa)
        analysis_tree_set_path = tmp_dir_path / "ancestral-tree-set-analysis.nwk"
        analysis_trees = [
            (
                source_tree_index,
                prune_tree_to_taxa(
                    tree,
                    analysis_taxa,
                    scratch_path=tmp_dir_path
                    / f"ancestral-tree-set-pruned-{source_tree_index}.nwk",
                ),
            )
            for source_tree_index, tree in kept_tree_entries
        ]
        analysis_tree_set_path.write_text(
            "".join(dumps_newick(tree) + "\n" for _, tree in analysis_trees),
            encoding="utf-8",
        )
        topology_summary = load_tree_set(analysis_tree_set_path)
        return (
            analysis_trees,
            topology_summary,
            analysis_taxa,
            sorted(exclusions, key=lambda row: (row.taxon, row.reason)),
            list(dataset.warnings),
            dataset.taxon_column,
        )
