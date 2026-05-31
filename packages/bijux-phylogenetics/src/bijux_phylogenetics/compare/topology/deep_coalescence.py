from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .models import (
    DeepCoalescenceBranchRow,
    DeepCoalescenceReport,
    DeepCoalescenceTaxonMapRow,
)


def count_deep_coalescences(
    species_tree_path: Path,
    gene_tree_path: Path,
    *,
    taxon_map_path: Path | None = None,
) -> DeepCoalescenceReport:
    """Count extra lineages by reconciling a rooted gene tree onto a rooted species tree."""
    species_tree = load_tree(species_tree_path)
    gene_tree = load_tree(gene_tree_path)
    _validate_species_tree(species_tree)
    _validate_gene_tree(gene_tree)

    species_descendants_by_key = _species_descendant_sets(species_tree)
    tip_map = _resolve_gene_to_species_tip_map(
        gene_tree=gene_tree,
        species_tree=species_tree,
        taxon_map_path=taxon_map_path,
    )
    event_count_by_branch: dict[str, int] = defaultdict(int)

    def map_gene_node(node: TreeNode) -> set[str]:
        if node.is_leaf():
            species_taxon = tip_map[node.name or ""]
            return {species_taxon}
        descendant_species: set[str] = set()
        for child in node.children:
            descendant_species.update(map_gene_node(child))
        mapped_branch = _find_species_lca_key(
            descendant_species,
            species_tree=species_tree,
            species_descendants_by_key=species_descendants_by_key,
        )
        event_count_by_branch[mapped_branch] += 1
        return descendant_species

    map_gene_node(gene_tree.root)
    branch_rows: list[DeepCoalescenceBranchRow] = []

    def summarize_species_node(node: TreeNode) -> int:
        branch_key = _species_branch_key(node)
        if node.is_leaf():
            lineage_count_entering = sum(
                1
                for species_taxon in tip_map.values()
                if species_taxon == (node.name or "")
            )
        else:
            lineage_count_entering = sum(
                summarize_species_node(child) for child in node.children
            )
        coalescent_event_count = event_count_by_branch.get(branch_key, 0)
        lineage_count_exiting = lineage_count_entering - coalescent_event_count
        if lineage_count_exiting < 0:
            raise ValueError(
                "gene-tree reconciliation produced more in-branch coalescences than entering lineages"
            )
        included_in_total = node.parent is not None
        branch_rows.append(
            DeepCoalescenceBranchRow(
                species_branch=branch_key,
                branch_role=_species_branch_role(node),
                descendant_species=sorted(species_descendants_by_key[branch_key]),
                lineage_count_entering=lineage_count_entering,
                coalescent_event_count=coalescent_event_count,
                lineage_count_exiting=lineage_count_exiting,
                extra_lineage_count=(
                    max(lineage_count_exiting - 1, 0) if included_in_total else 0
                ),
                included_in_deep_coalescence_total=included_in_total,
            )
        )
        return lineage_count_exiting

    root_lineage_count = summarize_species_node(species_tree.root)
    if root_lineage_count != 1:
        raise ValueError(
            "deep-coalescence counting requires the reconciled gene tree to reduce to one root lineage"
        )
    branch_rows.reverse()
    observed_species_taxa = sorted(set(tip_map.values()))
    species_only_taxa = sorted(set(species_tree.tip_names) - set(observed_species_taxa))
    return DeepCoalescenceReport(
        species_tree_path=species_tree_path,
        gene_tree_path=gene_tree_path,
        observed_species_taxa=observed_species_taxa,
        species_only_taxa=species_only_taxa,
        gene_tip_count=gene_tree.tip_count,
        deep_coalescence_total=sum(
            row.extra_lineage_count
            for row in branch_rows
            if row.included_in_deep_coalescence_total
        ),
        mapping_rows=[
            DeepCoalescenceTaxonMapRow(
                gene_taxon=gene_taxon,
                species_taxon=species_taxon,
            )
            for gene_taxon, species_taxon in sorted(tip_map.items())
        ],
        branch_rows=branch_rows,
    )


def _resolve_gene_to_species_tip_map(
    *,
    gene_tree: PhyloTree,
    species_tree: PhyloTree,
    taxon_map_path: Path | None,
) -> dict[str, str]:
    species_taxa = set(species_tree.tip_names)
    if taxon_map_path is None:
        if any(
            (tip_name or "") not in species_taxa for tip_name in gene_tree.tip_names
        ):
            raise ValueError(
                "deep-coalescence counting requires --taxon-map when gene tips do not exactly match species-tree taxa"
            )
        return {tip_name: tip_name for tip_name in gene_tree.tip_names}

    table = load_taxon_table(taxon_map_path)
    if "species_taxon" not in table.columns:
        raise ValueError("taxon map must include a 'species_taxon' column")
    tip_map = {row[table.taxon_column]: row["species_taxon"] for row in table.rows}
    for gene_taxon in gene_tree.tip_names:
        if gene_taxon not in tip_map:
            raise ValueError(f"taxon map is missing gene tip '{gene_taxon}'")
        if tip_map[gene_taxon] not in species_taxa:
            raise ValueError(
                f"taxon map assigns gene tip '{gene_taxon}' to unknown species '{tip_map[gene_taxon]}'"
            )
    extra_gene_taxa = sorted(set(tip_map) - set(gene_tree.tip_names))
    if extra_gene_taxa:
        raise ValueError(
            "taxon map contains gene taxa absent from the gene tree: "
            + ", ".join(extra_gene_taxa)
        )
    return tip_map


def _validate_species_tree(tree: PhyloTree) -> None:
    if len(tree.root.children) != 2:
        raise ValueError("deep-coalescence counting requires a rooted species tree")
    if any(node.name is None for node in tree.iter_leaves()):
        raise ValueError("deep-coalescence counting requires named species-tree tips")
    if len(set(tree.tip_names)) != tree.tip_count:
        raise ValueError(
            "deep-coalescence counting requires unique species-tree tip labels"
        )
    if any(len(node.children) != 2 for node in tree.iter_internal_nodes()):
        raise ValueError(
            "deep-coalescence counting requires a strictly binary species tree"
        )


def _validate_gene_tree(tree: PhyloTree) -> None:
    if len(tree.root.children) != 2:
        raise ValueError("deep-coalescence counting requires a rooted gene tree")
    if any(node.name is None for node in tree.iter_leaves()):
        raise ValueError("deep-coalescence counting requires named gene-tree tips")
    if len(set(tree.tip_names)) != tree.tip_count:
        raise ValueError(
            "deep-coalescence counting requires unique gene-tree tip labels"
        )
    if any(len(node.children) != 2 for node in tree.iter_internal_nodes()):
        raise ValueError(
            "deep-coalescence counting requires a strictly binary gene tree"
        )


def _species_descendant_sets(tree: PhyloTree) -> dict[str, set[str]]:
    descendant_sets: dict[str, set[str]] = {}
    for node in tree.iter_nodes(order="postorder"):
        descendant_sets[_species_branch_key(node)] = set(node.descendant_taxa)
    return descendant_sets


def _species_branch_key(node: TreeNode) -> str:
    if node.is_leaf():
        return node.name or "<unnamed-species>"
    return "|".join(node.descendant_taxa)


def _species_branch_role(node: TreeNode) -> str:
    if node.parent is None:
        return "root-population"
    if node.is_leaf():
        return "tip-branch"
    return "internal-branch"


def _find_species_lca_key(
    descendant_species: set[str],
    *,
    species_tree: PhyloTree,
    species_descendants_by_key: dict[str, set[str]],
) -> str:
    node = species_tree.root
    while True:
        containing_child = None
        for child in node.children:
            child_key = _species_branch_key(child)
            if descendant_species <= species_descendants_by_key[child_key]:
                containing_child = child
                break
        if containing_child is None:
            return _species_branch_key(node)
        node = containing_child
