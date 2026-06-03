from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.io.trees import load_tree

from .contracts import MetadataClusteringReport, MetadataClusterObservation


def compare_inferred_tree_to_taxon_metadata(
    tree_path: Path,
    metadata_path: Path,
    *,
    group_column: str,
    taxon_column: str | None = None,
) -> MetadataClusteringReport:
    """Report whether metadata-defined biological groups cluster monophyletically in one inferred tree."""
    from bijux_phylogenetics.ancestral.common import node_descendant_taxa

    tree = load_tree(tree_path)
    table = load_taxon_table(metadata_path, taxon_column=taxon_column)
    if group_column not in table.columns:
        raise ValueError(f"metadata table does not contain column '{group_column}'")
    tree_taxa = set(tree.tip_names)
    group_taxa: dict[str, set[str]] = {}
    for row in table.rows:
        taxon = row[table.taxon_column]
        group = row[group_column].strip()
        if group and taxon in tree_taxa:
            group_taxa.setdefault(group, set()).add(taxon)
    node_taxa_sets = [set(node_descendant_taxa(node)) for node in tree.iter_nodes()]
    observations: list[MetadataClusterObservation] = []
    for group, taxa in sorted(group_taxa.items()):
        ordered_taxa = sorted(taxa)
        if len(ordered_taxa) < 2:
            observations.append(
                MetadataClusterObservation(
                    group=group,
                    tree_taxa=ordered_taxa,
                    monophyletic=None,
                    status="not_evaluable",
                    note="group has fewer than two taxa in the inferred tree",
                )
            )
            continue
        monophyletic = any(node_taxa == taxa for node_taxa in node_taxa_sets)
        observations.append(
            MetadataClusterObservation(
                group=group,
                tree_taxa=ordered_taxa,
                monophyletic=monophyletic,
                status="clusters_as_expected" if monophyletic else "split_unexpectedly",
                note=(
                    "all observed group members collapse to one internal clade"
                    if monophyletic
                    else "group members are distributed across multiple clades in the inferred tree"
                ),
            )
        )
    return MetadataClusteringReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        taxon_column=table.taxon_column,
        group_column=group_column,
        group_count=len(observations),
        monophyletic_group_count=sum(
            1 for row in observations if row.monophyletic is True
        ),
        split_group_count=sum(1 for row in observations if row.monophyletic is False),
        observations=observations,
    )
