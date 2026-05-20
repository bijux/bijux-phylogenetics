from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


@dataclass(frozen=True, slots=True)
class StructuralTreeParityReport:
    """Structural equality result for one tree."""

    equivalent: bool
    mismatch_reason: str | None


@dataclass(frozen=True, slots=True)
class StructuralTreeSetParityReport:
    """Structural equality result for one ordered tree set."""

    equivalent: bool
    mismatch_reason: str | None


def _leaf_names(node: TreeNode) -> tuple[str, ...]:
    if node.is_leaf():
        if node.name is None:
            raise ValueError(
                "tree contains unnamed tips and cannot be compared structurally"
            )
        return (node.name,)
    labels: list[str] = []
    for child in node.children:
        labels.extend(_leaf_names(child))
    return tuple(sorted(labels))


def _rooted_flag(tree: PhyloTree) -> bool | None:
    return tree.rooted


def _edge_key(
    descendant_taxa: tuple[str, ...],
    *,
    all_taxa: tuple[str, ...],
    rooted: bool | None,
    is_tip: bool,
) -> tuple[str, tuple[str, ...] | str]:
    if is_tip:
        return ("tip", descendant_taxa[0])
    if rooted is False:
        descendant_set = set(descendant_taxa)
        complement_taxa = tuple(
            taxon for taxon in all_taxa if taxon not in descendant_set
        )
        canonical = min(
            (descendant_taxa, complement_taxa),
            key=lambda value: (len(value), value),
        )
        return ("split", canonical)
    return ("clade", descendant_taxa)


def _branch_lengths_match(
    expected: float | None,
    observed: float | None,
    *,
    tolerance: float,
) -> bool:
    if expected is None or observed is None:
        return expected is observed
    return abs(expected - observed) <= tolerance


def _internal_label(node: TreeNode) -> str:
    return "" if node.name is None else node.name


def _edge_records(
    tree: PhyloTree,
) -> tuple[
    Counter[str],
    dict[tuple[str, tuple[str, ...] | str], float | None],
    dict[tuple[str, tuple[str, ...] | str], str],
    str,
]:
    tip_counts = Counter(tree.tip_names)
    all_taxa = tuple(sorted(tree.tip_names))
    branch_lengths: dict[tuple[str, tuple[str, ...] | str], float | None] = {}
    internal_labels: dict[tuple[str, tuple[str, ...] | str], str] = {}

    def visit(node: TreeNode) -> tuple[str, ...]:
        descendant_taxa = _leaf_names(node)
        for child in node.children:
            visit(child)
        if node is not tree.root:
            edge_key = _edge_key(
                descendant_taxa,
                all_taxa=all_taxa,
                rooted=_rooted_flag(tree),
                is_tip=node.is_leaf(),
            )
            branch_lengths[edge_key] = node.branch_length
            if not node.is_leaf():
                internal_labels[edge_key] = _internal_label(node)
        return descendant_taxa

    visit(tree.root)
    root_label = _internal_label(tree.root)
    return tip_counts, branch_lengths, internal_labels, root_label


def compare_tree_structurally(
    expected: PhyloTree,
    observed: PhyloTree,
    *,
    tolerance: float = 0.0,
    compare_internal_labels: bool = True,
) -> StructuralTreeParityReport:
    """Compare two trees by rootedness, edge structure, branch lengths, and labels."""

    if _rooted_flag(expected) != _rooted_flag(observed):
        return StructuralTreeParityReport(
            equivalent=False,
            mismatch_reason=(
                f"tree rootedness differs: expected {_rooted_flag(expected)}, observed {_rooted_flag(observed)}"
            ),
        )

    expected_tips, expected_branch_lengths, expected_labels, expected_root_label = (
        _edge_records(expected)
    )
    observed_tips, observed_branch_lengths, observed_labels, observed_root_label = (
        _edge_records(observed)
    )

    if expected_tips != observed_tips:
        return StructuralTreeParityReport(
            equivalent=False,
            mismatch_reason=(
                f"tip labels differ: expected {sorted(expected_tips.elements())}, observed {sorted(observed_tips.elements())}"
            ),
        )

    if set(expected_branch_lengths) != set(observed_branch_lengths):
        return StructuralTreeParityReport(
            equivalent=False,
            mismatch_reason="clades or splits differ between trees",
        )

    for edge_key in sorted(expected_branch_lengths, key=str):
        expected_length = expected_branch_lengths[edge_key]
        observed_length = observed_branch_lengths[edge_key]
        if not _branch_lengths_match(
            expected_length, observed_length, tolerance=tolerance
        ):
            return StructuralTreeParityReport(
                equivalent=False,
                mismatch_reason=(
                    f"branch lengths differ for {edge_key[0]} {edge_key[1]}: expected {expected_length}, observed {observed_length}"
                ),
            )

    if compare_internal_labels:
        if (
            _rooted_flag(expected) is not False
            and expected_root_label != observed_root_label
        ):
            return StructuralTreeParityReport(
                equivalent=False,
                mismatch_reason=(
                    f"root label differs: expected {expected_root_label!r}, observed {observed_root_label!r}"
                ),
            )
        if expected_labels != observed_labels:
            differing_keys = sorted(
                set(expected_labels) | set(observed_labels),
                key=str,
            )
            for edge_key in differing_keys:
                if expected_labels.get(edge_key, "") != observed_labels.get(
                    edge_key, ""
                ):
                    return StructuralTreeParityReport(
                        equivalent=False,
                        mismatch_reason=(
                            f"internal labels differ for {edge_key[0]} {edge_key[1]}: expected {expected_labels.get(edge_key, '')!r}, observed {observed_labels.get(edge_key, '')!r}"
                        ),
                    )

    return StructuralTreeParityReport(equivalent=True, mismatch_reason=None)


def compare_tree_sets_structurally(
    expected: list[PhyloTree],
    observed: list[PhyloTree],
    *,
    tolerance: float = 0.0,
    compare_internal_labels: bool = True,
) -> StructuralTreeSetParityReport:
    """Compare two ordered tree sets by count and per-tree structure."""

    if len(expected) != len(observed):
        return StructuralTreeSetParityReport(
            equivalent=False,
            mismatch_reason=(
                f"tree counts differ: expected {len(expected)}, observed {len(observed)}"
            ),
        )
    for index, (expected_tree, observed_tree) in enumerate(
        zip(expected, observed, strict=True),
        start=1,
    ):
        report = compare_tree_structurally(
            expected_tree,
            observed_tree,
            tolerance=tolerance,
            compare_internal_labels=compare_internal_labels,
        )
        if not report.equivalent:
            return StructuralTreeSetParityReport(
                equivalent=False,
                mismatch_reason=f"tree {index}: {report.mismatch_reason}",
            )
    return StructuralTreeSetParityReport(equivalent=True, mismatch_reason=None)
