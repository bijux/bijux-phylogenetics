from __future__ import annotations

import csv
import math
from pathlib import Path

from bijux_phylogenetics.phylo.topology.node_identity import build_ape_internal_node_map
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_bipartition,
    informative_unrooted_splits,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidAlignmentError,
    WorkflowBudgetError,
)
from bijux_phylogenetics.io.newick import (
    dumps_newick,
    write_newick,
)
from bijux_phylogenetics.io.trees import load_tree
from .contracts import (
    CladeFrequency,
    CladeFrequencyReport,
    ConsensusTreeReport,
    TreeDistanceMatrixReport,
    TreeDistancePair,
    TreeSetCladeSupportReport,
    TreeSetCladeSupportRow,
    TreeSetProcessingSummary,
    TreeSetRecord,
    TreeSetReport,
    TreeSetWorkflowBudget,
    TreeSetWorkflowBudgetReport,
)
from .inventory import (
    _TreeSetAnalysis,
    _analyze_tree_set,
    _require_exact_taxa,
    _require_tree_set,
    _validate_same_taxa,
    load_tree_set,
)
from .topology import (
    _clade_counts,
    _clade_signature,
    _clades_conflict,
    _format_clade,
    _rooted_topology_id,
    _tree_distance,
)


def _maximal_nested_clades(
    parent: frozenset[str], clades: set[frozenset[str]]
) -> list[frozenset[str]]:
    nested = [clade for clade in clades if clade < parent]
    return sorted(
        [
            clade
            for clade in nested
            if not any(clade < other < parent for other in nested)
        ],
        key=lambda clade: (len(clade), sorted(clade)),
    )


def _mean(values: list[float]) -> float:
    return round(sum(values) / len(values), 15)


def _validate_budget_limit(
    value: int | None,
    *,
    name: str,
) -> int | None:
    if value is None:
        return None
    if value < 1:
        raise ValueError(f"{name} must be at least 1, got {value}")
    return value


def build_tree_set_workflow_budget(
    *,
    max_tree_count: int | None = None,
    max_report_table_rows: int | None = None,
    memory_warning_threshold_bytes: int | None = None,
) -> TreeSetWorkflowBudget:
    """Normalize one reviewer-facing resource budget for tree-set workflows."""
    validated_threshold = (
        None
        if memory_warning_threshold_bytes is None
        else _validate_budget_limit(
            memory_warning_threshold_bytes,
            name="memory_warning_threshold_bytes",
        )
    )
    return TreeSetWorkflowBudget(
        max_tree_count=_validate_budget_limit(
            max_tree_count,
            name="max_tree_count",
        ),
        max_report_table_rows=_validate_budget_limit(
            max_report_table_rows,
            name="max_report_table_rows",
        ),
        memory_warning_threshold_bytes=validated_threshold,
    )


def enforce_tree_set_tree_budget(
    *,
    tree_count: int,
    budget: TreeSetWorkflowBudget,
    workflow_name: str,
    source_path: Path,
) -> None:
    """Reject tree-set workflows that exceed an explicit input-size budget."""
    if budget.max_tree_count is None or tree_count <= budget.max_tree_count:
        return
    raise WorkflowBudgetError(
        (
            f"{workflow_name} budget allows at most {budget.max_tree_count} trees, "
            f"but {source_path} contains {tree_count}"
        ),
        code="tree_set_tree_budget_exceeded",
        details={
            "workflow_name": workflow_name,
            "source_path": str(source_path),
            "tree_count": tree_count,
            "max_tree_count": budget.max_tree_count,
        },
    )


def build_tree_set_budget_report(
    *,
    budget: TreeSetWorkflowBudget,
    peak_memory_bytes: int,
    truncated_section_names: list[str] | None = None,
) -> TreeSetWorkflowBudgetReport:
    """Summarize how one tree-set workflow budget was applied."""
    warning_messages: list[str] = []
    if (
        budget.memory_warning_threshold_bytes is not None
        and peak_memory_bytes > budget.memory_warning_threshold_bytes
    ):
        warning_messages.append(
            "peak memory exceeded the configured workflow warning threshold"
        )
    truncated_names = sorted(set(truncated_section_names or []))
    if truncated_names:
        warning_messages.append(
            "reviewer-facing sections were truncated to the configured row limit"
        )
    return TreeSetWorkflowBudgetReport(
        max_tree_count=budget.max_tree_count,
        max_report_table_rows=budget.max_report_table_rows,
        memory_warning_threshold_bytes=budget.memory_warning_threshold_bytes,
        truncated_section_names=truncated_names,
        warning_messages=warning_messages,
    )


def _support_classification(frequency: float, conflict_count: int) -> str:
    if frequency >= 0.9 and conflict_count == 0:
        return "robust"
    if 0.3 <= frequency <= 0.7:
        return "intermediate-support"
    if conflict_count > 0:
        return "credibility-conflicted"
    return "weak-support"


def _build_consensus_node(
    taxa: frozenset[str],
    *,
    majority_clades: set[frozenset[str]],
    clade_support: dict[frozenset[str], float],
    clade_lengths: dict[frozenset[str], float],
    terminal_lengths: dict[str, float],
    is_root: bool = False,
) -> TreeNode:
    child_clades = _maximal_nested_clades(taxa, majority_clades)
    covered: set[str] = set()
    children: list[TreeNode] = []
    for child_clade in child_clades:
        covered.update(child_clade)
        children.append(
            _build_consensus_node(
                child_clade,
                majority_clades=majority_clades,
                clade_support=clade_support,
                clade_lengths=clade_lengths,
                terminal_lengths=terminal_lengths,
            )
        )
    for taxon in sorted(taxa - covered):
        children.append(TreeNode(name=taxon, branch_length=terminal_lengths.get(taxon)))
    if len(children) == 1:
        return children[0]
    label = None if is_root else format(clade_support[taxa], ".15g")
    return TreeNode(
        name=label,
        branch_length=None if is_root else clade_lengths.get(taxa),
        children=children,
    )


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


def _build_consensus_tree_with_threshold(
    analysis: _TreeSetAnalysis,
    *,
    threshold: float,
) -> tuple[PhyloTree, ConsensusTreeReport]:
    shared_taxa = _require_exact_taxa(analysis)
    universe = frozenset(shared_taxa)
    counts = analysis.clade_counts or {}
    majority_clades = {
        clade
        for clade, count in counts.items()
        if count / len(analysis.trees) >= threshold
    }
    clade_support = {
        clade: round((counts[clade] / len(analysis.trees)) * 100.0, 15)
        for clade in majority_clades
    }
    clade_lengths = {
        clade: _mean(lengths)
        for clade, lengths in analysis.clade_branch_lengths.items()
        if clade in majority_clades and lengths
    }
    terminal_length_means = {
        taxon: _mean(lengths)
        for taxon, lengths in analysis.terminal_lengths.items()
        if lengths
    }
    tree = PhyloTree(
        root=_build_consensus_node(
            universe,
            majority_clades=majority_clades,
            clade_support=clade_support,
            clade_lengths=clade_lengths,
            terminal_lengths=terminal_length_means,
            is_root=True,
        ),
        source_format=analysis.source_format,
        rooted=True,
    )
    if math.isclose(threshold, 1.0):
        consensus_method = "strict"
    elif math.isclose(threshold, 0.5):
        consensus_method = "majority-rule"
    else:
        consensus_method = "thresholded"
    return tree, ConsensusTreeReport(
        path=analysis.path,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=shared_taxa,
        consensus_method=consensus_method,
        consensus_threshold=threshold,
        included_clade_count=len(majority_clades),
        consensus_newick=dumps_newick(tree),
    )


def _build_tree_distance_matrix_report(
    analysis: _TreeSetAnalysis,
) -> TreeDistanceMatrixReport:
    shared_taxa = set(_require_exact_taxa(analysis))
    pairs: list[TreeDistancePair] = []
    for left_index, left in enumerate(analysis.trees, start=1):
        for right_index, right in enumerate(
            analysis.trees[left_index - 1 :], start=left_index
        ):
            distance, normalized = _tree_distance(left, right, shared_taxa)
            pairs.append(
                TreeDistancePair(
                    left_index=left_index,
                    right_index=right_index,
                    robinson_foulds_distance=distance,
                    normalized_robinson_foulds=normalized,
                )
            )
    return TreeDistanceMatrixReport(
        path=analysis.path,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=sorted(shared_taxa),
        pairs=pairs,
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


def compute_consensus_tree(path: Path) -> tuple[PhyloTree, ConsensusTreeReport]:
    """Compute a majority-rule consensus tree from a tree set."""
    return compute_consensus_tree_with_threshold(path, threshold=0.5)


def compute_strict_consensus_tree(path: Path) -> tuple[PhyloTree, ConsensusTreeReport]:
    """Compute a strict consensus tree from a tree set."""
    return compute_consensus_tree_with_threshold(path, threshold=1.0)


def compute_consensus_tree_with_threshold(
    path: Path,
    *,
    threshold: float,
) -> tuple[PhyloTree, ConsensusTreeReport]:
    """Compute a deterministic consensus tree at a caller-supplied clade frequency threshold."""
    if not 0.0 < threshold <= 1.0:
        raise ValueError(
            f"consensus threshold must be greater than 0 and at most 1, got {threshold}"
        )
    return _build_consensus_tree_with_threshold(
        _analyze_tree_set(path), threshold=threshold
    )


def write_consensus_tree(path: Path, tree: PhyloTree) -> Path:
    """Write a consensus tree as canonical Newick."""
    return write_newick(path, tree)


def compute_tree_distance_matrix(path: Path) -> TreeDistanceMatrixReport:
    """Compute a pairwise RF-distance matrix across a tree set."""
    return _build_tree_distance_matrix_report(_analyze_tree_set(path))


def write_tree_distance_matrix(path: Path, report: TreeDistanceMatrixReport) -> Path:
    """Write a pairwise tree-distance matrix as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "left_index",
                "right_index",
                "robinson_foulds_distance",
                "normalized_robinson_foulds",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.pairs:
            writer.writerow(
                {
                    "left_index": row.left_index,
                    "right_index": row.right_index,
                    "robinson_foulds_distance": row.robinson_foulds_distance,
                    "normalized_robinson_foulds": format(
                        row.normalized_robinson_foulds, ".15g"
                    ),
                }
            )
    return path
