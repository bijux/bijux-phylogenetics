from __future__ import annotations

from dataclasses import asdict
import json
import math
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick, loads_newick, write_newick
from bijux_phylogenetics.phylo.topology.bionj import build_bionj_tree
from bijux_phylogenetics.phylo.topology.neighbor_joining import (
    build_neighbor_joining_tree,
)
from bijux_phylogenetics.phylo.topology.rooted_nni import (
    apply_rooted_nni_move,
    iter_rooted_nni_move_candidates,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .balanced_minimum_evolution import score_balanced_minimum_evolution
from .imported import (
    _distance_lookup_from_imported,
    load_imported_distance_matrix,
    validate_imported_distance_matrix,
)
from .models import (
    BalancedMinimumEvolutionNniSearchReport,
    BalancedMinimumEvolutionNniTraceRow,
)


def _strip_branch_lengths(tree: PhyloTree) -> PhyloTree:
    stripped = tree.copy().refresh()
    for node in stripped.iter_nodes(order="preorder"):
        if node is stripped.root:
            continue
        node.branch_length = None
    return stripped.refresh()


def _prefer_score(
    left_score: float,
    left_newick: str,
    right_score: float,
    right_newick: str,
) -> bool:
    if left_score < right_score and not math.isclose(left_score, right_score):
        return True
    if right_score < left_score and not math.isclose(left_score, right_score):
        return False
    return left_newick < right_newick


def search_balanced_minimum_evolution_nni(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
    *,
    start_method: str,
) -> BalancedMinimumEvolutionNniSearchReport:
    """Start from NJ or BIONJ and accept score-improving rooted NNI moves under BME."""
    resolved_start_method = start_method.strip().lower()
    if resolved_start_method not in {"neighbor-joining", "bionj"}:
        raise ValueError("start_method must be 'neighbor-joining' or 'bionj'")
    if resolved_start_method == "neighbor-joining":
        start_tree = build_neighbor_joining_tree(identifiers, distance_lookup)
    else:
        start_tree = build_bionj_tree(identifiers, distance_lookup)
    current_tree = _strip_branch_lengths(start_tree)
    current_score = score_balanced_minimum_evolution(
        current_tree,
        identifiers,
        distance_lookup,
    )
    start_tree_newick = dumps_newick(current_tree)
    trace_rows = [
        BalancedMinimumEvolutionNniTraceRow(
            event_index=1,
            event_kind="start",
            iteration=0,
            score_before=None,
            score_after=current_score,
            score_delta=None,
            tree_before_newick=None,
            tree_after_newick=start_tree_newick,
            pivot_branch_id=None,
            sibling_clade_id=None,
            exchanged_clade_id=None,
            stopping_reason=None,
        )
    ]
    accepted_move_count = 0
    evaluated_neighbor_count = 0
    while True:
        improving_tree: PhyloTree | None = None
        improving_score: float | None = None
        improving_newick: str | None = None
        improving_candidate = None
        for candidate in iter_rooted_nni_move_candidates(current_tree):
            neighbor_tree = apply_rooted_nni_move(current_tree, candidate)
            neighbor_score = score_balanced_minimum_evolution(
                neighbor_tree,
                identifiers,
                distance_lookup,
            )
            evaluated_neighbor_count += 1
            neighbor_newick = dumps_newick(neighbor_tree)
            if neighbor_score > current_score or math.isclose(
                neighbor_score,
                current_score,
            ):
                continue
            if improving_tree is None or _prefer_score(
                neighbor_score,
                neighbor_newick,
                improving_score if improving_score is not None else neighbor_score,
                improving_newick if improving_newick is not None else neighbor_newick,
            ):
                improving_tree = neighbor_tree
                improving_score = neighbor_score
                improving_newick = neighbor_newick
                improving_candidate = candidate
        if (
            improving_tree is None
            or improving_score is None
            or improving_newick is None
        ):
            stopping_reason = "no-improving-neighbor"
            break
        accepted_move_count += 1
        tree_before_newick = dumps_newick(current_tree)
        score_before = current_score
        current_tree = improving_tree.copy().refresh()
        current_score = improving_score
        trace_rows.append(
            BalancedMinimumEvolutionNniTraceRow(
                event_index=len(trace_rows) + 1,
                event_kind="accepted-move",
                iteration=accepted_move_count,
                score_before=score_before,
                score_after=current_score,
                score_delta=current_score - score_before,
                tree_before_newick=tree_before_newick,
                tree_after_newick=improving_newick,
                pivot_branch_id=improving_candidate.pivot_branch_id,
                sibling_clade_id=improving_candidate.sibling_clade_id,
                exchanged_clade_id=improving_candidate.exchanged_clade_id,
                stopping_reason=None,
            )
        )
    trace_rows.append(
        BalancedMinimumEvolutionNniTraceRow(
            event_index=len(trace_rows) + 1,
            event_kind="final",
            iteration=accepted_move_count,
            score_before=None,
            score_after=current_score,
            score_delta=None,
            tree_before_newick=None,
            tree_after_newick=dumps_newick(current_tree),
            pivot_branch_id=None,
            sibling_clade_id=None,
            exchanged_clade_id=None,
            stopping_reason=stopping_reason,
        )
    )
    pair_count = (len(identifiers) * (len(identifiers) - 1)) // 2
    return BalancedMinimumEvolutionNniSearchReport(
        algorithm="balanced-minimum-evolution-nni-search",
        matrix_path=None,
        start_method=resolved_start_method,
        taxon_count=len(identifiers),
        pair_count=pair_count,
        start_tree_newick=start_tree_newick,
        start_score=trace_rows[0].score_after,
        final_tree_newick=trace_rows[-1].tree_after_newick,
        final_score=trace_rows[-1].score_after,
        accepted_move_count=accepted_move_count,
        evaluated_neighbor_count=evaluated_neighbor_count,
        stopping_reason=stopping_reason,
        trace_rows=trace_rows,
    )


def search_balanced_minimum_evolution_nni_from_imported_distance_matrix(
    matrix_path: Path,
    *,
    start_method: str,
) -> BalancedMinimumEvolutionNniSearchReport:
    """Search one imported distance matrix by rooted NNI under the BME objective."""
    entries = load_imported_distance_matrix(matrix_path)
    validation = validate_imported_distance_matrix(matrix_path)
    distance_lookup, _missing_distance_policy_report = _distance_lookup_from_imported(
        validation,
        entries,
    )
    report = search_balanced_minimum_evolution_nni(
        validation.identifiers,
        distance_lookup,
        start_method=start_method,
    )
    report.matrix_path = matrix_path
    return report


def write_balanced_minimum_evolution_nni_trace_table(
    path: Path,
    report: BalancedMinimumEvolutionNniSearchReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "event_index",
        "event_kind",
        "iteration",
        "score_before",
        "score_after",
        "score_delta",
        "tree_before_newick",
        "tree_after_newick",
        "pivot_branch_id",
        "sibling_clade_id",
        "exchanged_clade_id",
        "stopping_reason",
    ]
    rows = ["\t".join(columns)]
    for row in report.trace_rows:
        payload = [
            row.event_index,
            row.event_kind,
            row.iteration,
            row.score_before,
            row.score_after,
            row.score_delta,
            row.tree_before_newick,
            row.tree_after_newick,
            row.pivot_branch_id,
            row.sibling_clade_id,
            row.exchanged_clade_id,
            row.stopping_reason,
        ]
        rows.append(
            "\t".join(
                ""
                if value is None
                else format(value, ".15g")
                if isinstance(value, float)
                else str(value)
                for value in payload
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_balanced_minimum_evolution_nni_run_json(
    path: Path,
    report: BalancedMinimumEvolutionNniSearchReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "algorithm": report.algorithm,
        "matrix_path": None if report.matrix_path is None else str(report.matrix_path),
        "start_method": report.start_method,
        "taxon_count": report.taxon_count,
        "pair_count": report.pair_count,
        "start_tree_newick": report.start_tree_newick,
        "start_score": report.start_score,
        "final_tree_newick": report.final_tree_newick,
        "final_score": report.final_score,
        "accepted_move_count": report.accepted_move_count,
        "evaluated_neighbor_count": report.evaluated_neighbor_count,
        "stopping_reason": report.stopping_reason,
        "trace_rows": [asdict(row) for row in report.trace_rows],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def write_balanced_minimum_evolution_nni_artifacts(
    out_dir: Path,
    report: BalancedMinimumEvolutionNniSearchReport,
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    start_tree_path = write_newick(
        out_dir / "start_tree.nwk",
        loads_newick(report.start_tree_newick),
    )
    final_tree_path = write_newick(
        out_dir / "final_tree.nwk",
        loads_newick(report.final_tree_newick),
    )
    trace_path = write_balanced_minimum_evolution_nni_trace_table(
        out_dir / "search_trace.tsv",
        report,
    )
    run_json_path = write_balanced_minimum_evolution_nni_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "start_tree_path": start_tree_path,
        "final_tree_path": final_tree_path,
        "trace_path": trace_path,
        "run_json_path": run_json_path,
    }
