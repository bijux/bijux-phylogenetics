from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass
import json
import math
from pathlib import Path

from bijux_phylogenetics.io.newick import loads_newick, write_newick
from bijux_phylogenetics.phylo.topology.clades import canonical_clade_id
from bijux_phylogenetics.phylo.topology.models import (
    StepwiseAdditionCandidateScore,
    StepwiseAdditionTraceRow,
    StepwiseAdditionTreeReport,
)

from .tree import PhyloTree, TreeNode, descendant_taxa

STEPWISE_ADDITION_ROOT_BRANCH_ID = "root"
_SUPPORTED_STEPWISE_OBJECTIVE_DIRECTIONS = frozenset({"minimize", "maximize"})


@dataclass(frozen=True, slots=True)
class StepwiseAdditionEdgeCandidate:
    """One branch on which a new taxon may be inserted."""

    target_node_id: str | None
    branch_id: str
    descendant_taxa: tuple[str, ...]


def validate_stepwise_addition_taxa(taxa: list[str]) -> list[str]:
    """Require at least two distinct non-empty taxa while preserving insertion order."""
    if len(taxa) < 2:
        raise ValueError("stepwise addition requires at least two taxa")
    seen: set[str] = set()
    duplicates: list[str] = []
    for taxon in taxa:
        if not taxon.strip():
            raise ValueError("stepwise addition does not allow blank taxon labels")
        if taxon in seen and taxon not in duplicates:
            duplicates.append(taxon)
        seen.add(taxon)
    if duplicates:
        raise ValueError(
            "stepwise addition requires distinct taxa; duplicates: "
            + ", ".join(duplicates)
        )
    return list(taxa)


def validate_stepwise_objective_direction(objective_direction: str) -> str:
    """Validate whether the score objective is minimized or maximized."""
    normalized_direction = objective_direction.strip().lower()
    if normalized_direction not in _SUPPORTED_STEPWISE_OBJECTIVE_DIRECTIONS:
        raise ValueError(
            "objective_direction must be one of "
            + ", ".join(sorted(_SUPPORTED_STEPWISE_OBJECTIVE_DIRECTIONS))
        )
    return normalized_direction


def stepwise_addition_node_sort_key(node: TreeNode) -> tuple[int, tuple[str, ...]]:
    """Sort insertion edges deterministically by descendant taxa."""
    descendants = tuple(descendant_taxa(node))
    return (len(descendants), descendants)


def iter_stepwise_addition_edge_candidates(
    tree: PhyloTree,
) -> Iterator[StepwiseAdditionEdgeCandidate]:
    """Yield all legal rooted insertion edges exactly once in deterministic order."""
    current_taxa = tuple(sorted(tree.tip_names))
    yield StepwiseAdditionEdgeCandidate(
        target_node_id=None,
        branch_id=STEPWISE_ADDITION_ROOT_BRANCH_ID,
        descendant_taxa=current_taxa,
    )
    candidate_nodes = sorted(
        (node for node in tree.iter_nodes(order="preorder") if node is not tree.root),
        key=stepwise_addition_node_sort_key,
    )
    for node in candidate_nodes:
        if node.node_id is None:
            raise AssertionError(
                "stepwise addition candidates require refreshed node IDs"
            )
        yield StepwiseAdditionEdgeCandidate(
            target_node_id=node.node_id,
            branch_id=canonical_clade_id(frozenset(node.descendant_taxa)),
            descendant_taxa=tuple(sorted(node.descendant_taxa)),
        )


def apply_stepwise_addition_candidate(
    tree: PhyloTree,
    candidate: StepwiseAdditionEdgeCandidate,
    taxon: str,
) -> PhyloTree:
    """Return one copied rooted tree with the selected taxon inserted on one edge."""
    if not taxon.strip():
        raise ValueError("stepwise addition does not allow blank taxon labels")
    if taxon in tree.tip_names:
        raise ValueError(f"stepwise addition tree already contains taxon '{taxon}'")
    working_tree = tree.copy().refresh()
    new_leaf = TreeNode(name=taxon)
    tree_has_explicit_branch_lengths = any(
        child.branch_length is not None for _parent, child in working_tree.iter_edges()
    )
    if candidate.target_node_id is None:
        if tree_has_explicit_branch_lengths:
            if working_tree.root.branch_length is None:
                working_tree.root.branch_length = 0.0
            new_leaf.branch_length = 0.0
        return PhyloTree(
            root=TreeNode(children=[working_tree.root, new_leaf]),
            source_format=working_tree.source_format,
            rooted=True,
        ).refresh()

    target = working_tree.node_by_id(candidate.target_node_id)
    target_parent = target.parent
    if target_parent is None:
        raise AssertionError("non-root stepwise addition candidates require one parent")
    inserted = TreeNode(children=[target, new_leaf])
    if tree_has_explicit_branch_lengths:
        original_target_branch_length = target.branch_length
        if original_target_branch_length is not None:
            split_length = original_target_branch_length / 2.0
            inserted.branch_length = split_length
            target.branch_length = split_length
            new_leaf.branch_length = split_length
        else:
            inserted.branch_length = 0.0
            target.branch_length = 0.0
            new_leaf.branch_length = 0.0
    target_parent.replace_children(
        [inserted if child is target else child for child in target_parent.children]
    )
    return working_tree.refresh()


def build_greedy_stepwise_addition_tree(
    taxa: list[str],
    *,
    score_tree: Callable[[PhyloTree], float],
    objective_name: str,
    objective_direction: str = "minimize",
) -> tuple[PhyloTree, StepwiseAdditionTreeReport]:
    """Build one rooted tree by inserting each taxon at the best-scoring edge."""
    ordered_taxa = validate_stepwise_addition_taxa(taxa)
    resolved_direction = validate_stepwise_objective_direction(objective_direction)
    current_tree = PhyloTree(
        root=TreeNode(
            children=[
                TreeNode(name=ordered_taxa[0]),
                TreeNode(name=ordered_taxa[1]),
            ]
        ),
        rooted=True,
    ).refresh()
    trace_rows: list[StepwiseAdditionTraceRow] = []
    current_score = float(score_tree(current_tree))
    for step_index, taxon in enumerate(ordered_taxa[2:], start=1):
        tested_edge_rows: list[StepwiseAdditionCandidateScore] = []
        best_candidate: StepwiseAdditionEdgeCandidate | None = None
        best_tree: PhyloTree | None = None
        best_score: float | None = None
        best_tree_newick: str | None = None
        for candidate in iter_stepwise_addition_edge_candidates(current_tree):
            candidate_tree = apply_stepwise_addition_candidate(
                current_tree, candidate, taxon
            )
            candidate_score = float(score_tree(candidate_tree))
            candidate_tree_newick = candidate_tree.to_newick()
            tested_edge_rows.append(
                StepwiseAdditionCandidateScore(
                    branch_id=candidate.branch_id,
                    descendant_taxa=list(candidate.descendant_taxa),
                    score=candidate_score,
                    candidate_tree_newick=candidate_tree_newick,
                )
            )
            if best_score is None or _prefer_stepwise_score(
                candidate_score,
                candidate_tree_newick,
                best_score,
                best_tree_newick,
                objective_direction=resolved_direction,
            ):
                best_candidate = candidate
                best_tree = candidate_tree
                best_score = candidate_score
                best_tree_newick = candidate_tree_newick
        if (
            best_candidate is None
            or best_tree is None
            or best_score is None
            or best_tree_newick is None
        ):
            raise AssertionError("stepwise addition must evaluate at least one edge")
        current_tree = best_tree
        current_score = best_score
        trace_rows.append(
            StepwiseAdditionTraceRow(
                step_index=step_index,
                taxon=taxon,
                inserted_taxa=ordered_taxa[: step_index + 2],
                tested_edge_rows=tested_edge_rows,
                best_edge_id=best_candidate.branch_id,
                best_edge_descendant_taxa=list(best_candidate.descendant_taxa),
                best_score=best_score,
                selected_tree_newick=best_tree_newick,
            )
        )
    report = summarize_stepwise_addition_tree(
        current_tree,
        insertion_order=ordered_taxa,
        objective_name=objective_name,
        objective_direction=resolved_direction,
        final_score=current_score,
        trace_rows=trace_rows,
    )
    return current_tree, report


def summarize_stepwise_addition_tree(
    tree: PhyloTree,
    *,
    insertion_order: list[str],
    objective_name: str,
    objective_direction: str,
    final_score: float,
    trace_rows: list[StepwiseAdditionTraceRow],
) -> StepwiseAdditionTreeReport:
    """Summarize one greedy stepwise-addition build against the requested taxon order."""
    tip_order = tree.tip_names
    tip_name_counts: dict[str, int] = {}
    for taxon in tip_order:
        tip_name_counts[taxon] = tip_name_counts.get(taxon, 0) + 1
    duplicate_generated_taxa = sorted(
        taxon for taxon, count in tip_name_counts.items() if count > 1
    )
    generated_taxa = set(tip_order)
    requested_taxon_set = set(insertion_order)
    missing_requested_taxa = sorted(requested_taxon_set - generated_taxa)
    unexpected_generated_taxa = sorted(generated_taxa - requested_taxon_set)
    validation_errors = tree.validation_errors()
    return StepwiseAdditionTreeReport(
        algorithm="greedy-stepwise-addition-tree",
        objective_name=objective_name,
        objective_direction=objective_direction,
        insertion_order=list(insertion_order),
        starting_taxa=list(insertion_order[:2]),
        tip_order=tip_order,
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        rooted=tree.rooted,
        strictly_bifurcating=_tree_is_strictly_bifurcating(tree),
        all_requested_taxa_present_once=not (
            missing_requested_taxa
            or duplicate_generated_taxa
            or unexpected_generated_taxa
        ),
        missing_requested_taxa=missing_requested_taxa,
        duplicate_generated_taxa=duplicate_generated_taxa,
        unexpected_generated_taxa=unexpected_generated_taxa,
        validation_errors=validation_errors,
        final_score=final_score,
        trace_rows=trace_rows,
        tree_newick=tree.to_newick(),
    )


def write_stepwise_addition_trace_table(
    path: Path,
    report: StepwiseAdditionTreeReport,
) -> Path:
    """Write one row per tested insertion edge across the stepwise trace."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "step_index",
                "taxon",
                "inserted_taxa",
                "tested_edge_id",
                "tested_edge_descendant_taxa",
                "tested_edge_score",
                "best_edge_id",
                "best_edge_descendant_taxa",
                "best_score",
                "selected",
                "candidate_tree_newick",
            ]
        )
    ]
    for row in report.trace_rows:
        for tested_edge in row.tested_edge_rows:
            lines.append(
                "\t".join(
                    [
                        str(row.step_index),
                        row.taxon,
                        ",".join(row.inserted_taxa),
                        tested_edge.branch_id,
                        ",".join(tested_edge.descendant_taxa),
                        str(tested_edge.score),
                        row.best_edge_id,
                        ",".join(row.best_edge_descendant_taxa),
                        str(row.best_score),
                        (
                            "true"
                            if tested_edge.branch_id == row.best_edge_id
                            and math.isclose(tested_edge.score, row.best_score)
                            and tested_edge.candidate_tree_newick
                            == row.selected_tree_newick
                            else "false"
                        ),
                        tested_edge.candidate_tree_newick,
                    ]
                )
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_stepwise_addition_run_json(
    path: Path,
    report: StepwiseAdditionTreeReport,
) -> Path:
    """Write one machine-readable greedy stepwise-addition payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "algorithm": report.algorithm,
        "objective_name": report.objective_name,
        "objective_direction": report.objective_direction,
        "insertion_order": report.insertion_order,
        "starting_taxa": report.starting_taxa,
        "tip_order": report.tip_order,
        "tip_count": report.tip_count,
        "internal_node_count": report.internal_node_count,
        "rooted": report.rooted,
        "strictly_bifurcating": report.strictly_bifurcating,
        "all_requested_taxa_present_once": report.all_requested_taxa_present_once,
        "missing_requested_taxa": report.missing_requested_taxa,
        "duplicate_generated_taxa": report.duplicate_generated_taxa,
        "unexpected_generated_taxa": report.unexpected_generated_taxa,
        "validation_errors": report.validation_errors,
        "final_score": report.final_score,
        "trace_rows": [
            {
                "step_index": row.step_index,
                "taxon": row.taxon,
                "inserted_taxa": row.inserted_taxa,
                "tested_edge_rows": [
                    {
                        "branch_id": tested_edge.branch_id,
                        "descendant_taxa": tested_edge.descendant_taxa,
                        "score": tested_edge.score,
                        "candidate_tree_newick": tested_edge.candidate_tree_newick,
                    }
                    for tested_edge in row.tested_edge_rows
                ],
                "best_edge_id": row.best_edge_id,
                "best_edge_descendant_taxa": row.best_edge_descendant_taxa,
                "best_score": row.best_score,
                "selected_tree_newick": row.selected_tree_newick,
            }
            for row in report.trace_rows
        ],
        "tree_newick": report.tree_newick,
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def write_stepwise_addition_artifacts(
    out_dir: Path,
    report: StepwiseAdditionTreeReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one stepwise-addition build."""
    out_dir.mkdir(parents=True, exist_ok=True)
    tree_path = write_newick(out_dir / "tree.nwk", loads_newick(report.tree_newick))
    trace_path = write_stepwise_addition_trace_table(out_dir / "trace.tsv", report)
    run_json_path = write_stepwise_addition_run_json(out_dir / "run.json", report)
    return {
        "tree_path": tree_path,
        "trace_path": trace_path,
        "run_json_path": run_json_path,
    }


def _prefer_stepwise_score(
    left_score: float,
    left_newick: str,
    right_score: float,
    right_newick: str | None,
    *,
    objective_direction: str,
) -> bool:
    if objective_direction == "minimize":
        if left_score < right_score and not math.isclose(left_score, right_score):
            return True
        if right_score < left_score and not math.isclose(left_score, right_score):
            return False
        return left_newick < (right_newick or "")
    if left_score > right_score and not math.isclose(left_score, right_score):
        return True
    if right_score > left_score and not math.isclose(left_score, right_score):
        return False
    return left_newick < (right_newick or "")


def _tree_is_strictly_bifurcating(tree: PhyloTree) -> bool:
    return all(node.is_leaf() or len(node.children) == 2 for node in tree.iter_nodes())
