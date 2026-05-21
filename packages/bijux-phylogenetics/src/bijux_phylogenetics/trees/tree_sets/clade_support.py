from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_bipartition,
    informative_unrooted_splits,
)
from bijux_phylogenetics.phylo.topology.node_identity import build_ape_internal_node_map
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

from .contracts import (
    CladeFrequency,
    CladeFrequencyReport,
    TreeSetCladeSupportReport,
    TreeSetCladeSupportRow,
)
from .inventory import _analyze_tree_set, _require_exact_taxa, _TreeSetAnalysis
from .topology import _format_clade


def _support_classification(frequency: float, conflict_count: int) -> str:
    if frequency >= 0.9 and conflict_count == 0:
        return "robust"
    if 0.3 <= frequency <= 0.7:
        return "intermediate-support"
    if conflict_count > 0:
        return "credibility-conflicted"
    return "weak-support"


def _build_clade_frequency_report(analysis: _TreeSetAnalysis) -> CladeFrequencyReport:
    exact_taxa = _require_exact_taxa(analysis)
    counts = analysis.clade_counts or {}
    total = len(analysis.trees)
    return CladeFrequencyReport(
        path=analysis.path,
        tree_count=total,
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        clade_frequencies=[
            CladeFrequency(
                clade=_format_clade(clade),
                tree_count=count,
                frequency=round(count / total, 15),
            )
            for clade, count in sorted(
                counts.items(), key=lambda item: _format_clade(item[0])
            )
        ],
    )


def _descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_descendant_taxa(child))
    return sorted(taxa)


def _reference_tree_root_depths(tree: PhyloTree) -> dict[str, float | None]:
    depths: dict[str, float | None] = {tree.root.node_id or "": 0.0}

    def visit(node: TreeNode) -> None:
        base_depth = depths[node.node_id or ""]
        for child in node.children:
            if base_depth is None or child.branch_length is None:
                depths[child.node_id or ""] = None
            else:
                depths[child.node_id or ""] = round(
                    base_depth + float(child.branch_length), 15
                )
            visit(child)

    visit(tree.root)
    return depths


def _split_counts(
    trees: list[PhyloTree],
    shared_taxa: set[str],
) -> dict[frozenset[str], int]:
    counts: dict[frozenset[str], int] = {}
    for tree in trees:
        for split in informative_unrooted_splits(tree, shared_taxa):
            counts[split] = counts.get(split, 0) + 1
    return counts


def _clade_support_status(
    *,
    supporting_tree_count: int | None,
    tree_count: int,
    node_kind: str,
    unscored_reason: str | None = None,
) -> tuple[str, str]:
    if node_kind == "root":
        return (
            "fixed",
            "the root spans the full compatible taxon set and is present in every comparison tree",
        )
    if supporting_tree_count is None:
        if unscored_reason == "absent-root-split":
            return (
                "not-counted",
                "ape::prop.clades leaves this root-adjacent split unscored when the comparison tree set never realizes the matching bipartition",
            )
        return (
            "not-counted",
            "ape::prop.clades leaves this root-adjacent clade unscored because its complement is a singleton tip",
        )
    if supporting_tree_count == 0:
        return (
            "absent",
            "the reference clade is absent from the comparison tree set",
        )
    if supporting_tree_count == tree_count:
        return (
            "fixed",
            "the reference clade is present in every comparison tree",
        )
    return (
        "partial-support",
        "the reference clade is present in only a subset of comparison trees",
    )


def _build_reference_tree_clade_support_report(
    *,
    reference_tree_path: Path,
    reference_tree: PhyloTree,
    analysis: _TreeSetAnalysis,
) -> TreeSetCladeSupportReport:
    exact_taxa = _require_exact_taxa(analysis)
    reference_taxa = set(exact_taxa)
    clade_counts = analysis.clade_counts or {}
    split_counts = _split_counts(analysis.trees, reference_taxa)
    tree_count = len(analysis.trees)
    root_depths = _reference_tree_root_depths(reference_tree)
    rows: list[TreeSetCladeSupportRow] = []
    supported_clade_count = 0
    absent_clade_count = 0
    unscored_clade_count = 0

    for node_id, node in build_ape_internal_node_map(reference_tree).items():
        descendant_taxa = _descendant_taxa(node)
        clade = frozenset(descendant_taxa)
        node_kind = "root" if node is reference_tree.root else "internal"
        unscored_reason: str | None = None
        if node_kind == "root":
            supporting_tree_count = tree_count
            clade_frequency = 1.0
            support_percent = 100.0
        elif len(clade) == len(reference_taxa) - 1:
            supporting_tree_count = None
            clade_frequency = None
            support_percent = None
            unscored_reason = "singleton-complement"
            unscored_clade_count += 1
        elif node in reference_tree.root.children:
            split = canonical_bipartition(set(descendant_taxa), reference_taxa)
            split_support = split_counts.get(split, 0)
            if split_support == 0:
                supporting_tree_count = None
                clade_frequency = None
                support_percent = None
                unscored_reason = "absent-root-split"
                unscored_clade_count += 1
            else:
                supporting_tree_count = split_support
                clade_frequency = round(supporting_tree_count / tree_count, 15)
                support_percent = round(clade_frequency * 100.0, 15)
                supported_clade_count += 1
        else:
            supporting_tree_count = clade_counts.get(clade, 0)
            clade_frequency = round(supporting_tree_count / tree_count, 15)
            support_percent = round(clade_frequency * 100.0, 15)
            if supporting_tree_count == 0:
                absent_clade_count += 1
            else:
                supported_clade_count += 1
        support_status, explanation = _clade_support_status(
            supporting_tree_count=supporting_tree_count,
            tree_count=tree_count,
            node_kind=node_kind,
            unscored_reason=unscored_reason,
        )
        rows.append(
            TreeSetCladeSupportRow(
                node_id=node_id,
                node_kind=node_kind,
                node_label=node.name,
                descendant_taxa=descendant_taxa,
                supporting_tree_count=supporting_tree_count,
                clade_frequency=clade_frequency,
                support_percent=support_percent,
                support_status=support_status,
                explanation=explanation,
                reference_branch_length=node.branch_length,
                reference_root_depth=root_depths[node.node_id or ""],
            )
        )
    return TreeSetCladeSupportReport(
        reference_tree_path=reference_tree_path,
        comparison_tree_set_path=analysis.path,
        tree_count=tree_count,
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        supported_clade_count=supported_clade_count,
        absent_clade_count=absent_clade_count,
        unscored_clade_count=unscored_clade_count,
        rows=rows,
    )


def compute_clade_frequency_table(path: Path) -> CladeFrequencyReport:
    """Compute informative clade frequencies across a tree set with a shared taxon set."""
    return _build_clade_frequency_report(_analyze_tree_set(path))


def compute_reference_tree_clade_support(
    reference_tree_path: Path,
    comparison_tree_set_path: Path,
) -> TreeSetCladeSupportReport:
    """Map tree-set clade support onto one reference tree by descendant tip set."""
    reference_tree = load_tree(reference_tree_path)
    analysis = _analyze_tree_set(comparison_tree_set_path)
    if analysis.exact_taxa is None:
        raise InvalidAlignmentError(
            "reference tree support mapping requires all comparison trees to share the exact same taxon set"
        )
    exact_taxa = analysis.exact_taxa
    if sorted(reference_tree.tip_names) != exact_taxa:
        raise InvalidAlignmentError(
            "reference tree and comparison tree set must share the exact same taxon set"
        )
    return _build_reference_tree_clade_support_report(
        reference_tree_path=reference_tree_path,
        reference_tree=reference_tree,
        analysis=analysis,
    )


def write_clade_frequency_table(path: Path, report: CladeFrequencyReport) -> Path:
    """Write a clade-frequency table as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["clade", "tree_count", "frequency"], delimiter="\t"
        )
        writer.writeheader()
        for row in report.clade_frequencies:
            writer.writerow(
                {
                    "clade": row.clade,
                    "tree_count": row.tree_count,
                    "frequency": format(row.frequency, ".15g"),
                }
            )
    return path


def write_reference_tree_clade_support_table(
    path: Path,
    report: TreeSetCladeSupportReport,
) -> Path:
    """Write one reference-tree clade support table as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "node_id",
                "node_kind",
                "node_label",
                "descendant_taxa",
                "supporting_tree_count",
                "clade_frequency",
                "support_percent",
                "support_status",
                "explanation",
                "reference_branch_length",
                "reference_root_depth",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.rows:
            writer.writerow(
                {
                    "node_id": row.node_id,
                    "node_kind": row.node_kind,
                    "node_label": "" if row.node_label is None else row.node_label,
                    "descendant_taxa": "|".join(row.descendant_taxa),
                    "supporting_tree_count": (
                        ""
                        if row.supporting_tree_count is None
                        else row.supporting_tree_count
                    ),
                    "clade_frequency": (
                        ""
                        if row.clade_frequency is None
                        else format(row.clade_frequency, ".15g")
                    ),
                    "support_percent": (
                        ""
                        if row.support_percent is None
                        else format(row.support_percent, ".15g")
                    ),
                    "support_status": row.support_status,
                    "explanation": row.explanation,
                    "reference_branch_length": (
                        ""
                        if row.reference_branch_length is None
                        else format(row.reference_branch_length, ".15g")
                    ),
                    "reference_root_depth": (
                        ""
                        if row.reference_root_depth is None
                        else format(row.reference_root_depth, ".15g")
                    ),
                }
            )
    return path
