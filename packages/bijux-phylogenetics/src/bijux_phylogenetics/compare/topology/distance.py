from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation.structure import _load_tree
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_clade_id,
    informative_rooted_clades,
    informative_unrooted_splits,
    robinson_foulds_metrics,
    split_sort_key,
    tree_has_polytomy,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .comparison import (
    _resolve_shared_taxa,
    _validate_rf_mode,
    _validate_taxon_overlap_policy,
)
from .models import (
    RobinsonFouldsMode,
    TaxonOverlapPolicy,
)


@dataclass(frozen=True, slots=True)
class TopologyDistanceSplitRow:
    """One clade or split contributing to a topology-distance comparison."""

    split_id: str
    split_kind: str
    comparison_status: str
    taxon_count: int
    descendant_taxa: tuple[str, ...]
    left_present: bool
    right_present: bool


@dataclass(slots=True)
class TopologyDistanceReport:
    """Owned RF-style topology distance summary aligned to `ape::dist.topo` cases."""

    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    taxon_overlap_policy: str
    rf_mode: str
    rooted_left: bool | None
    rooted_right: bool | None
    polytomy_present_left: bool
    polytomy_present_right: bool
    left_split_count: int
    right_split_count: int
    shared_split_count: int
    left_only_split_count: int
    right_only_split_count: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    topology_equal: bool
    split_rows: list[TopologyDistanceSplitRow]


def _topology_signatures(
    tree: PhyloTree,
    shared_taxa: set[str],
    *,
    rf_mode: RobinsonFouldsMode,
) -> set[frozenset[str]]:
    if rf_mode == "rooted":
        return informative_rooted_clades(tree, shared_taxa)
    return informative_unrooted_splits(tree, shared_taxa)


def compare_topology_distance_trees(
    left: PhyloTree,
    right: PhyloTree,
    *,
    left_path: Path,
    right_path: Path,
    rf_mode: RobinsonFouldsMode = "rooted",
    taxon_overlap_policy: TaxonOverlapPolicy = "require-identical",
) -> TopologyDistanceReport:
    _validate_rf_mode(rf_mode)
    _validate_taxon_overlap_policy(taxon_overlap_policy)

    left_taxa = set(left.tip_names)
    right_taxa = set(right.tip_names)
    shared_taxa, left_only_taxa, right_only_taxa = _resolve_shared_taxa(
        left_taxa,
        right_taxa,
        taxon_overlap_policy=taxon_overlap_policy,
    )
    left_signatures = _topology_signatures(left, shared_taxa, rf_mode=rf_mode)
    right_signatures = _topology_signatures(right, shared_taxa, rf_mode=rf_mode)
    metrics = robinson_foulds_metrics(
        left,
        right,
        shared_taxa,
        rf_mode=rf_mode,
    )

    shared_signatures = left_signatures & right_signatures
    left_only_signatures = left_signatures - right_signatures
    right_only_signatures = right_signatures - left_signatures
    split_kind = "clade" if rf_mode == "rooted" else "split"
    all_signatures = sorted(
        left_signatures | right_signatures,
        key=split_sort_key,
    )
    split_rows = [
        TopologyDistanceSplitRow(
            split_id=canonical_clade_id(signature),
            split_kind=split_kind,
            comparison_status=(
                "shared"
                if signature in shared_signatures
                else "left_only"
                if signature in left_only_signatures
                else "right_only"
            ),
            taxon_count=len(signature),
            descendant_taxa=tuple(sorted(signature)),
            left_present=signature in left_signatures,
            right_present=signature in right_signatures,
        )
        for signature in all_signatures
    ]
    return TopologyDistanceReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=sorted(shared_taxa),
        left_only_taxa=left_only_taxa,
        right_only_taxa=right_only_taxa,
        taxon_overlap_policy=taxon_overlap_policy,
        rf_mode=rf_mode,
        rooted_left=left.rooted,
        rooted_right=right.rooted,
        polytomy_present_left=tree_has_polytomy(left),
        polytomy_present_right=tree_has_polytomy(right),
        left_split_count=metrics.left_count,
        right_split_count=metrics.right_count,
        shared_split_count=len(shared_signatures),
        left_only_split_count=len(left_only_signatures),
        right_only_split_count=len(right_only_signatures),
        robinson_foulds_distance=metrics.distance,
        normalized_robinson_foulds=metrics.normalized_distance,
        topology_equal=metrics.distance == 0,
        split_rows=split_rows,
    )


def compare_topology_distance(
    left_path: Path,
    right_path: Path,
    *,
    rf_mode: RobinsonFouldsMode = "rooted",
    taxon_overlap_policy: TaxonOverlapPolicy = "require-identical",
) -> TopologyDistanceReport:
    """Compare tree topology by explicit RF-style split ledgers."""

    left = _load_tree(left_path)
    right = _load_tree(right_path)
    return compare_topology_distance_trees(
        left,
        right,
        left_path=left_path,
        right_path=right_path,
        rf_mode=rf_mode,
        taxon_overlap_policy=taxon_overlap_policy,
    )


def write_topology_distance_split_table(
    path: Path,
    left_path: Path,
    right_path: Path,
    *,
    rf_mode: RobinsonFouldsMode = "rooted",
    taxon_overlap_policy: TaxonOverlapPolicy = "require-identical",
) -> Path:
    """Write the split ledger behind one topology-distance comparison."""

    report = compare_topology_distance(
        left_path,
        right_path,
        rf_mode=rf_mode,
        taxon_overlap_policy=taxon_overlap_policy,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "split_id",
                "split_kind",
                "comparison_status",
                "taxon_count",
                "descendant_taxa",
                "left_present",
                "right_present",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.split_rows:
            writer.writerow(
                {
                    "split_id": row.split_id,
                    "split_kind": row.split_kind,
                    "comparison_status": row.comparison_status,
                    "taxon_count": row.taxon_count,
                    "descendant_taxa": "|".join(row.descendant_taxa),
                    "left_present": str(row.left_present).lower(),
                    "right_present": str(row.right_present).lower(),
                }
            )
    return path
