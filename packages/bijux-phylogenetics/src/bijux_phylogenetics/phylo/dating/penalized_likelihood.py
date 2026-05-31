from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path

import numpy

from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.newick import dumps_newick, loads_newick, write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood.parameter_search import (
    run_bounded_coordinate_likelihood_search,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    PhylogeneticsError,
    UnrootedTreeError,
)

from .inputs import load_tip_dates_for_tree, validate_tip_dates_against_tree
from .least_squares import fit_least_squares_dating
from .models import (
    PenalizedLikelihoodDatingBranchRow,
    PenalizedLikelihoodDatingNodeRow,
    PenalizedLikelihoodDatingReport,
)

_DATE_TOLERANCE = 1e-9
_ALPHA_EPSILON = 1e-6
_ROOT_HEIGHT_EPSILON = 1e-6


@dataclass(frozen=True, slots=True)
class _PenalizedLikelihoodTreeLayout:
    all_nodes: list[TreeNode]
    internal_nodes: list[TreeNode]
    internal_non_root_nodes: list[TreeNode]
    node_index: dict[str, int]
    upper_date_by_node_id: dict[str, float]
    alpha_parameter_by_node_id: dict[str, str]
    tip_node_id_by_name: dict[str, str]
    fixed_node_dates: dict[str, float]


@dataclass(frozen=True, slots=True)
class _PenalizedLikelihoodCandidate:
    node_dates: dict[str, float]
    node_log_rates: dict[str, float]
    node_rates: dict[str, float]
    data_score: float
    penalty_score: float
    total_score: float
    condition_number: float
    node_rows: list[PenalizedLikelihoodDatingNodeRow]
    branch_rows: list[PenalizedLikelihoodDatingBranchRow]
    dated_tree_newick: str


def fit_penalized_likelihood_dating(
    tree: PhyloTree,
    tip_dates: Mapping[str, float],
    *,
    fixed_node_dates: Mapping[str, float] | None = None,
    smoothing_parameter: float = 1.0,
    max_coordinate_passes: int = 8,
    tree_path: Path | None = None,
    metadata_path: Path | None = None,
    taxon_column: str = "taxon",
    date_column: str = "date",
) -> PenalizedLikelihoodDatingReport:
    """Fit one rooted tree to tip dates and optional fixed node dates with penalized rate smoothing."""
    normalized_fixed_node_dates = _normalize_fixed_node_dates(
        tree,
        fixed_node_dates=fixed_node_dates,
    )
    _validate_penalized_likelihood_inputs(
        tree,
        tip_dates,
        fixed_node_dates=normalized_fixed_node_dates,
        smoothing_parameter=smoothing_parameter,
    )
    least_squares_report = fit_least_squares_dating(
        tree,
        tip_dates,
        tree_path=tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        date_column=date_column,
    )
    layout = _build_tree_layout(
        tree,
        tip_dates,
        fixed_node_dates=normalized_fixed_node_dates,
    )
    initial_values, bounds_by_name = _build_initial_parameters(
        tree,
        layout=layout,
        least_squares_report=least_squares_report,
        tip_dates=tip_dates,
    )

    def evaluate(
        parameter_values: dict[str, float],
    ) -> tuple[_PenalizedLikelihoodCandidate | None, float]:
        try:
            candidate = _evaluate_candidate(
                tree,
                layout=layout,
                tip_dates=tip_dates,
                parameter_values=parameter_values,
                smoothing_parameter=smoothing_parameter,
            )
        except (ArithmeticError, ValueError):
            return None, -math.inf
        return candidate, -candidate.total_score

    search_result = run_bounded_coordinate_likelihood_search(
        initial_values=initial_values,
        bounds_by_name=bounds_by_name,
        evaluate=evaluate,
        improvement_tolerance=1e-6,
        max_coordinate_passes=max_coordinate_passes,
    )
    best_candidate = search_result.payload
    if best_candidate is None:
        raise PhylogeneticsError(
            "penalized likelihood dating could not find one valid chronological fit",
            code="penalized_likelihood_dating_error",
        )
    return PenalizedLikelihoodDatingReport(
        tree_newick=dumps_newick(tree),
        dated_tree_newick=best_candidate.dated_tree_newick,
        taxa=sorted(tree.tip_names),
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        branch_count=len(best_candidate.branch_rows),
        parameter_count=len(search_result.parameter_values) + len(layout.all_nodes),
        tree_path=None if tree_path is None else str(tree_path),
        metadata_path=None if metadata_path is None else str(metadata_path),
        taxon_column=taxon_column,
        date_column=date_column,
        minimum_tip_date=min(tip_dates.values()),
        maximum_tip_date=max(tip_dates.values()),
        root_date=best_candidate.node_dates[tree.root.node_id or ""],
        smoothing_parameter=smoothing_parameter,
        data_score=best_candidate.data_score,
        penalty_score=best_candidate.penalty_score,
        total_score=best_candidate.total_score,
        condition_number=best_candidate.condition_number,
        optimizer_name="bounded-coordinate-search with closed-form penalized log-rate solve",
        optimization_pass_count=search_result.optimization_pass_count,
        function_evaluation_count=search_result.function_evaluation_count,
        converged=search_result.converged,
        node_rows=best_candidate.node_rows,
        branch_rows=best_candidate.branch_rows,
    )


def fit_penalized_likelihood_dating_from_metadata(
    tree_path: Path,
    metadata_path: Path,
    *,
    smoothing_parameter: float = 1.0,
    max_coordinate_passes: int = 8,
    taxon_column: str | None = None,
    date_column: str = "date",
) -> PenalizedLikelihoodDatingReport:
    """Fit one penalized likelihood dated-tree report from one rooted tree path and one tip-date table."""
    validate_tree_path(tree_path, require_rooted=True)
    tree = load_tree(tree_path)
    tree.rooted = True
    tip_dates, resolved_taxon_column = load_tip_dates_for_tree(
        metadata_path,
        tree_taxa=tree.tip_names,
        taxon_column=taxon_column,
        date_column=date_column,
    )
    return fit_penalized_likelihood_dating(
        tree,
        tip_dates,
        smoothing_parameter=smoothing_parameter,
        max_coordinate_passes=max_coordinate_passes,
        tree_path=tree_path,
        metadata_path=metadata_path,
        taxon_column=resolved_taxon_column,
        date_column=date_column,
    )


def _validate_penalized_likelihood_inputs(
    tree: PhyloTree,
    tip_dates: Mapping[str, float],
    *,
    fixed_node_dates: Mapping[str, float],
    smoothing_parameter: float,
) -> None:
    if tree.rooted is not True:
        raise UnrootedTreeError(
            "penalized likelihood dating requires one rooted tree",
            code="penalized_likelihood_dating_requires_rooted_tree",
        )
    if smoothing_parameter <= 0.0:
        raise PhylogeneticsError(
            "penalized likelihood dating requires a strictly positive smoothing parameter",
            code="penalized_likelihood_dating_error",
        )
    validate_tip_dates_against_tree(tree, tip_dates)
    if len(set(tip_dates.values())) < 2:
        raise PhylogeneticsError(
            "penalized likelihood dating requires variation in tip dates",
            code="penalized_likelihood_dating_error",
        )
    for node_id, fixed_date in fixed_node_dates.items():
        node = tree.node_by_id(node_id)
        if node is tree.root:
            continue
        parent = node.parent
        if parent is None:
            raise PhylogeneticsError(
                "penalized likelihood dating requires stable rooted parentage",
                code="penalized_likelihood_dating_error",
            )
        parent_node_id = parent.node_id or ""
        if parent_node_id in fixed_node_dates and fixed_date <= (
            fixed_node_dates[parent_node_id] + _DATE_TOLERANCE
        ):
            raise PhylogeneticsError(
                "fixed node dates must remain chronological along the rooted tree",
                code="penalized_likelihood_dating_error",
            )
    for branch_length in tree.branch_lengths():
        if branch_length is None:
            raise InvalidBranchLengthError(
                "penalized likelihood dating requires complete branch lengths"
            )
        if branch_length <= 0.0:
            raise InvalidBranchLengthError(
                "penalized likelihood dating requires strictly positive branch lengths"
            )


def _normalize_fixed_node_dates(
    tree: PhyloTree,
    *,
    fixed_node_dates: Mapping[str, float] | None,
) -> dict[str, float]:
    if fixed_node_dates is None:
        return {}
    normalized: dict[str, float] = {}
    for node_id, raw_date in fixed_node_dates.items():
        if tree.node_by_id(node_id) is None:
            raise PhylogeneticsError(
                f"fixed node date references unknown node_id '{node_id}'",
                code="penalized_likelihood_dating_error",
            )
        try:
            normalized[node_id] = float(raw_date)
        except (TypeError, ValueError) as error:
            raise PhylogeneticsError(
                f"fixed node date for '{node_id}' must be numeric",
                code="penalized_likelihood_dating_error",
            ) from error
    return normalized


def _build_tree_layout(
    tree: PhyloTree,
    tip_dates: Mapping[str, float],
    *,
    fixed_node_dates: Mapping[str, float],
) -> _PenalizedLikelihoodTreeLayout:
    all_nodes = list(tree.iter_nodes(order="preorder"))
    internal_nodes = list(tree.iter_internal_nodes(order="preorder"))
    internal_non_root_nodes = [node for node in internal_nodes if node is not tree.root]
    node_index = {node.node_id or "": index for index, node in enumerate(all_nodes)}
    upper_date_by_node_id = _build_upper_date_bounds(
        tree,
        tip_dates=tip_dates,
        fixed_node_dates=fixed_node_dates,
    )
    alpha_parameter_by_node_id = {
        node.node_id or "": f"alpha::{node.node_id}"
        for node in internal_non_root_nodes
        if (node.node_id or "") not in fixed_node_dates
    }
    tip_node_id_by_name = {}
    for node in tree.iter_leaves():
        if node.name is None or node.node_id is None:
            raise PhylogeneticsError(
                "penalized likelihood dating requires stable named tips",
                code="penalized_likelihood_dating_error",
            )
        tip_node_id_by_name[node.name] = node.node_id
    return _PenalizedLikelihoodTreeLayout(
        all_nodes=all_nodes,
        internal_nodes=internal_nodes,
        internal_non_root_nodes=internal_non_root_nodes,
        node_index=node_index,
        upper_date_by_node_id=upper_date_by_node_id,
        alpha_parameter_by_node_id=alpha_parameter_by_node_id,
        tip_node_id_by_name=tip_node_id_by_name,
        fixed_node_dates=dict(fixed_node_dates),
    )


def _build_upper_date_bounds(
    tree: PhyloTree,
    *,
    tip_dates: Mapping[str, float],
    fixed_node_dates: Mapping[str, float],
) -> dict[str, float]:
    upper_date_by_node_id: dict[str, float] = {}

    def visit(node: TreeNode) -> float:
        node_id = node.node_id or ""
        if node.is_leaf():
            if node.name is None:
                raise PhylogeneticsError(
                    "penalized likelihood dating requires stable named tips",
                    code="penalized_likelihood_dating_error",
                )
            return float(tip_dates[node.name])
        proper_descendant_upper = min(visit(child) for child in node.children)
        upper_date_by_node_id[node_id] = proper_descendant_upper
        if node_id in fixed_node_dates:
            return min(fixed_node_dates[node_id], proper_descendant_upper)
        return proper_descendant_upper

    visit(tree.root)
    return upper_date_by_node_id


def _build_initial_parameters(
    tree: PhyloTree,
    *,
    layout: _PenalizedLikelihoodTreeLayout,
    least_squares_report,
    tip_dates: Mapping[str, float],
) -> tuple[dict[str, float], dict[str, tuple[float, float]]]:
    minimum_tip_date = min(tip_dates.values())
    maximum_tip_date = max(tip_dates.values())
    tip_span = max(maximum_tip_date - minimum_tip_date, 1.0)
    max_root_to_tip_branch_length = max(_root_to_tip_branch_lengths(tree), default=1.0)
    rate_floor = max(least_squares_report.estimated_clock_rate * 0.1, 1e-6)
    initial_node_dates = _build_initial_node_dates(
        tree,
        layout=layout,
        least_squares_report=least_squares_report,
        tip_dates=tip_dates,
    )
    initial_values: dict[str, float] = {}
    bounds_by_name: dict[str, tuple[float, float]] = {}
    root_node_id = tree.root.node_id or ""
    if root_node_id not in layout.fixed_node_dates:
        initial_root_height = max(
            minimum_tip_date - initial_node_dates[root_node_id],
            _ROOT_HEIGHT_EPSILON * 10.0,
        )
        root_height_upper = max(
            initial_root_height * 20.0,
            tip_span + (max_root_to_tip_branch_length / rate_floor),
            initial_root_height + tip_span + 1.0,
        )
        initial_values["root_height"] = initial_root_height
        bounds_by_name["root_height"] = (_ROOT_HEIGHT_EPSILON, root_height_upper)
    for node in layout.internal_non_root_nodes:
        node_id = node.node_id or ""
        if node_id in layout.fixed_node_dates:
            continue
        parent = node.parent
        if parent is None or parent.node_id is None:
            raise PhylogeneticsError(
                "penalized likelihood dating requires stable rooted parentage",
                code="penalized_likelihood_dating_error",
            )
        parent_date = initial_node_dates[parent.node_id]
        upper_date = layout.upper_date_by_node_id[node_id]
        numerator = initial_node_dates[node_id] - parent_date
        denominator = max(upper_date - parent_date, _DATE_TOLERANCE)
        initial_alpha = min(
            max(numerator / denominator, _ALPHA_EPSILON),
            1.0 - _ALPHA_EPSILON,
        )
        parameter_name = layout.alpha_parameter_by_node_id[node_id]
        initial_values[parameter_name] = initial_alpha
        bounds_by_name[parameter_name] = (_ALPHA_EPSILON, 1.0 - _ALPHA_EPSILON)
    return initial_values, bounds_by_name


def _build_initial_node_dates(
    tree: PhyloTree,
    *,
    layout: _PenalizedLikelihoodTreeLayout,
    least_squares_report,
    tip_dates: Mapping[str, float],
) -> dict[str, float]:
    least_squares_node_dates = {
        row.node_id: row.estimated_date for row in least_squares_report.node_rows
    }
    minimum_tip_date = min(tip_dates.values())
    maximum_tip_date = max(tip_dates.values())
    tip_span = max(maximum_tip_date - minimum_tip_date, 1.0)
    node_dates: dict[str, float] = {}
    root_node_id = tree.root.node_id or ""
    root_upper_date = layout.upper_date_by_node_id[root_node_id]
    if root_node_id in layout.fixed_node_dates:
        root_date = layout.fixed_node_dates[root_node_id]
    else:
        root_date = min(
            least_squares_node_dates[root_node_id],
            root_upper_date - max((tip_span * 0.1), (_ROOT_HEIGHT_EPSILON * 10.0)),
        )
    _require_chronological_node_date(
        node_label="root",
        fixed_date=root_date,
        lower_date=None,
        upper_date=root_upper_date,
    )
    node_dates[root_node_id] = root_date
    for node in layout.internal_non_root_nodes:
        node_id = node.node_id or ""
        parent = node.parent
        if parent is None or parent.node_id is None:
            raise PhylogeneticsError(
                "penalized likelihood dating requires stable rooted parentage",
                code="penalized_likelihood_dating_error",
            )
        lower_date = node_dates[parent.node_id]
        upper_date = layout.upper_date_by_node_id[node_id]
        if node_id in layout.fixed_node_dates:
            node_date = layout.fixed_node_dates[node_id]
            _require_chronological_node_date(
                node_label=node_id,
                fixed_date=node_date,
                lower_date=lower_date,
                upper_date=upper_date,
            )
        else:
            node_date = _clamp_initial_node_date(
                guess=least_squares_node_dates[node_id],
                lower_date=lower_date,
                upper_date=upper_date,
            )
        node_dates[node_id] = node_date
    return node_dates


def _clamp_initial_node_date(
    *,
    guess: float,
    lower_date: float,
    upper_date: float,
) -> float:
    margin = _chronology_margin(lower_date, upper_date)
    clamped = min(max(guess, lower_date + margin), upper_date - margin)
    if not (lower_date + _DATE_TOLERANCE < clamped < upper_date - _DATE_TOLERANCE):
        raise PhylogeneticsError(
            "penalized likelihood dating could not place an initial chronological node date",
            code="penalized_likelihood_dating_error",
        )
    return clamped


def _require_chronological_node_date(
    *,
    node_label: str,
    fixed_date: float,
    lower_date: float | None,
    upper_date: float,
) -> None:
    if lower_date is not None and fixed_date <= (lower_date + _DATE_TOLERANCE):
        raise PhylogeneticsError(
            f"fixed node date for '{node_label}' is not younger than its parent",
            code="penalized_likelihood_dating_error",
        )
    if fixed_date >= (upper_date - _DATE_TOLERANCE):
        raise PhylogeneticsError(
            f"fixed node date for '{node_label}' is not older than its descendants",
            code="penalized_likelihood_dating_error",
        )


def _chronology_margin(lower_date: float, upper_date: float) -> float:
    return max((upper_date - lower_date) * 1e-6, _DATE_TOLERANCE * 10.0)


def _root_to_tip_branch_lengths(tree: PhyloTree) -> list[float]:
    path_lengths: list[float] = []

    def visit(node: TreeNode, running_length: float) -> None:
        if node.is_leaf():
            path_lengths.append(running_length)
            return
        for child in node.children:
            visit(child, running_length + float(child.branch_length or 0.0))

    visit(tree.root, 0.0)
    return path_lengths


def _evaluate_candidate(
    tree: PhyloTree,
    *,
    layout: _PenalizedLikelihoodTreeLayout,
    tip_dates: Mapping[str, float],
    parameter_values: Mapping[str, float],
    smoothing_parameter: float,
) -> _PenalizedLikelihoodCandidate:
    node_dates = _build_node_dates(
        tree,
        layout=layout,
        tip_dates=tip_dates,
        parameter_values=parameter_values,
    )
    node_log_rates, condition_number = _solve_penalized_node_log_rates(
        tree,
        layout=layout,
        node_dates=node_dates,
        smoothing_parameter=smoothing_parameter,
    )
    node_rates = {
        node_id: math.exp(log_rate) for node_id, log_rate in node_log_rates.items()
    }
    maximum_tip_date = max(tip_dates.values())
    node_rows = _build_penalized_node_rows(
        tree,
        node_dates=node_dates,
        node_log_rates=node_log_rates,
        node_rates=node_rates,
        tip_dates=tip_dates,
        maximum_tip_date=maximum_tip_date,
    )
    branch_rows = _build_penalized_branch_rows(
        tree,
        node_dates=node_dates,
        node_log_rates=node_log_rates,
        node_rates=node_rates,
        smoothing_parameter=smoothing_parameter,
    )
    data_score = sum(row.data_score_contribution for row in branch_rows)
    penalty_score = sum(row.smoothing_penalty_contribution for row in branch_rows)
    dated_tree = _build_dated_tree(tree, node_dates)
    return _PenalizedLikelihoodCandidate(
        node_dates=node_dates,
        node_log_rates=node_log_rates,
        node_rates=node_rates,
        data_score=data_score,
        penalty_score=penalty_score,
        total_score=data_score + penalty_score,
        condition_number=condition_number,
        node_rows=node_rows,
        branch_rows=branch_rows,
        dated_tree_newick=dumps_newick(dated_tree),
    )


def _build_node_dates(
    tree: PhyloTree,
    *,
    layout: _PenalizedLikelihoodTreeLayout,
    tip_dates: Mapping[str, float],
    parameter_values: Mapping[str, float],
) -> dict[str, float]:
    minimum_tip_date = min(tip_dates.values())
    root_node_id = tree.root.node_id or ""
    if root_node_id in layout.fixed_node_dates:
        root_date = layout.fixed_node_dates[root_node_id]
    else:
        root_date = minimum_tip_date - float(parameter_values["root_height"])
    node_dates = {root_node_id: root_date}
    for node in layout.internal_non_root_nodes:
        node_id = node.node_id or ""
        parent = node.parent
        if parent is None or parent.node_id is None:
            raise ValueError("internal node is missing one rooted parent")
        parent_date = node_dates[parent.node_id]
        upper_date = layout.upper_date_by_node_id[node_id]
        if node_id in layout.fixed_node_dates:
            node_date = layout.fixed_node_dates[node_id]
            if node_date <= (parent_date + _DATE_TOLERANCE) or node_date >= (
                upper_date - _DATE_TOLERANCE
            ):
                raise ValueError("fixed node dates are not chronological")
        else:
            alpha = float(parameter_values[layout.alpha_parameter_by_node_id[node_id]])
            node_date = parent_date + (alpha * (upper_date - parent_date))
            if node_date <= (parent_date + _DATE_TOLERANCE) or node_date >= (
                upper_date - _DATE_TOLERANCE
            ):
                raise ValueError("candidate dates are not chronological")
        node_dates[node_id] = node_date
    for tip_name, tip_date in tip_dates.items():
        node_dates[layout.tip_node_id_by_name[tip_name]] = tip_date
    return node_dates


def _solve_penalized_node_log_rates(
    tree: PhyloTree,
    *,
    layout: _PenalizedLikelihoodTreeLayout,
    node_dates: Mapping[str, float],
    smoothing_parameter: float,
) -> tuple[dict[str, float], float]:
    edge_count = sum(1 for _parent, _child in tree.iter_edges())
    node_count = len(layout.all_nodes)
    data_matrix = numpy.zeros((edge_count, node_count), dtype=float)
    penalty_matrix = numpy.zeros((edge_count, node_count), dtype=float)
    observed_log_rates = numpy.zeros(edge_count, dtype=float)

    for edge_index, (parent, child) in enumerate(tree.iter_edges()):
        parent_node_id = parent.node_id or ""
        child_node_id = child.node_id or ""
        duration = node_dates[child_node_id] - node_dates[parent_node_id]
        if duration <= _DATE_TOLERANCE:
            raise ValueError("candidate branch duration is not positive")
        observed_branch_length = float(child.branch_length or 0.0)
        observed_log_rates[edge_index] = math.log(observed_branch_length / duration)
        parent_index = layout.node_index[parent_node_id]
        child_index = layout.node_index[child_node_id]
        data_matrix[edge_index][parent_index] = 0.5
        data_matrix[edge_index][child_index] = 0.5
        penalty_matrix[edge_index][parent_index] = 1.0
        penalty_matrix[edge_index][child_index] = -1.0

    normal_matrix = (data_matrix.T @ data_matrix) + (
        smoothing_parameter * (penalty_matrix.T @ penalty_matrix)
    )
    right_hand_side = data_matrix.T @ observed_log_rates
    try:
        solution = numpy.linalg.solve(normal_matrix, right_hand_side)
    except numpy.linalg.LinAlgError as error:
        raise ValueError("penalized likelihood rate system is singular") from error
    node_log_rates = {
        node.node_id or "": float(solution[layout.node_index[node.node_id or ""]])
        for node in layout.all_nodes
    }
    return node_log_rates, float(numpy.linalg.cond(normal_matrix, numpy.inf))


def _build_penalized_node_rows(
    tree: PhyloTree,
    *,
    node_dates: Mapping[str, float],
    node_log_rates: Mapping[str, float],
    node_rates: Mapping[str, float],
    tip_dates: Mapping[str, float],
    maximum_tip_date: float,
) -> list[PenalizedLikelihoodDatingNodeRow]:
    rows: list[PenalizedLikelihoodDatingNodeRow] = []
    for node in tree.iter_nodes(order="preorder"):
        node_id = node.node_id or ""
        if node is tree.root:
            node_kind = "root"
        elif node.is_leaf():
            node_kind = "tip"
        else:
            node_kind = "internal"
        rows.append(
            PenalizedLikelihoodDatingNodeRow(
                node_id=node_id,
                node_kind=node_kind,
                node_label=node.name,
                descendant_taxa=node.descendant_taxa,
                estimated_date=node_dates[node_id],
                fixed_tip_date=node.is_leaf() and node.name in tip_dates,
                time_height=maximum_tip_date - node_dates[node_id],
                estimated_log_rate=node_log_rates[node_id],
                estimated_rate=node_rates[node_id],
            )
        )
    return rows


def _build_penalized_branch_rows(
    tree: PhyloTree,
    *,
    node_dates: Mapping[str, float],
    node_log_rates: Mapping[str, float],
    node_rates: Mapping[str, float],
    smoothing_parameter: float,
) -> list[PenalizedLikelihoodDatingBranchRow]:
    rows: list[PenalizedLikelihoodDatingBranchRow] = []
    for parent, child in tree.iter_edges():
        parent_node_id = parent.node_id or ""
        child_node_id = child.node_id or ""
        parent_date = node_dates[parent_node_id]
        child_date = node_dates[child_node_id]
        duration = child_date - parent_date
        observed_branch_length = float(child.branch_length or 0.0)
        observed_log_rate = math.log(observed_branch_length / duration)
        fitted_log_rate = (
            node_log_rates[parent_node_id] + node_log_rates[child_node_id]
        ) / 2.0
        estimated_branch_rate = math.exp(fitted_log_rate)
        fitted_branch_length = duration * estimated_branch_rate
        data_score_contribution = (fitted_log_rate - observed_log_rate) ** 2
        smoothing_penalty_contribution = (
            smoothing_parameter
            * (node_log_rates[parent_node_id] - node_log_rates[child_node_id]) ** 2
        )
        rows.append(
            PenalizedLikelihoodDatingBranchRow(
                branch_id=child_node_id,
                child_name=child.name,
                descendant_taxa=child.descendant_taxa,
                parent_date=parent_date,
                child_date=child_date,
                fitted_time_duration=duration,
                observed_branch_length=observed_branch_length,
                observed_log_rate=observed_log_rate,
                fitted_log_rate=fitted_log_rate,
                estimated_branch_rate=estimated_branch_rate,
                fitted_branch_length=fitted_branch_length,
                data_score_contribution=data_score_contribution,
                smoothing_penalty_contribution=smoothing_penalty_contribution,
            )
        )
    return rows


def _build_dated_tree(
    tree: PhyloTree,
    node_dates: Mapping[str, float],
) -> PhyloTree:
    dated_tree = tree.copy()
    for parent, child in dated_tree.iter_edges():
        parent_node_id = parent.node_id or ""
        child_node_id = child.node_id or ""
        duration = node_dates[child_node_id] - node_dates[parent_node_id]
        if duration < -_DATE_TOLERANCE:
            raise ValueError("candidate dated tree is not chronological")
        child.branch_length = max(duration, 0.0)
    return dated_tree


def write_penalized_likelihood_dating_summary_tsv(
    path: Path,
    report: PenalizedLikelihoodDatingReport,
) -> Path:
    """Write one summary row for one penalized likelihood dating run."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "tree_path",
        "metadata_path",
        "taxon_column",
        "date_column",
        "tip_count",
        "internal_node_count",
        "branch_count",
        "parameter_count",
        "minimum_tip_date",
        "maximum_tip_date",
        "root_date",
        "smoothing_parameter",
        "data_score",
        "penalty_score",
        "total_score",
        "condition_number",
        "optimizer_name",
        "optimization_pass_count",
        "function_evaluation_count",
        "converged",
    ]
    values = [
        report.tree_path or "",
        report.metadata_path or "",
        report.taxon_column,
        report.date_column,
        str(report.tip_count),
        str(report.internal_node_count),
        str(report.branch_count),
        str(report.parameter_count),
        format(report.minimum_tip_date, ".15g"),
        format(report.maximum_tip_date, ".15g"),
        format(report.root_date, ".15g"),
        format(report.smoothing_parameter, ".15g"),
        format(report.data_score, ".15g"),
        format(report.penalty_score, ".15g"),
        format(report.total_score, ".15g"),
        format(report.condition_number, ".15g"),
        report.optimizer_name,
        str(report.optimization_pass_count),
        str(report.function_evaluation_count),
        str(report.converged).lower(),
    ]
    path.write_text(
        "\n".join(["\t".join(columns), "\t".join(values)]) + "\n",
        encoding="utf-8",
    )
    return path


def write_penalized_likelihood_node_dates_tsv(
    path: Path,
    report: PenalizedLikelihoodDatingReport,
) -> Path:
    """Write one dated-node row per node in tree preorder."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "node_id",
        "node_kind",
        "node_label",
        "descendant_taxa",
        "estimated_date",
        "fixed_tip_date",
        "time_height",
        "estimated_log_rate",
        "estimated_rate",
    ]
    lines = ["\t".join(columns)]
    lines.extend(
        "\t".join(
            [
                row.node_id,
                row.node_kind,
                row.node_label or "",
                "|".join(row.descendant_taxa),
                format(row.estimated_date, ".15g"),
                str(row.fixed_tip_date).lower(),
                format(row.time_height, ".15g"),
                format(row.estimated_log_rate, ".15g"),
                format(row.estimated_rate, ".15g"),
            ]
        )
        for row in report.node_rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_penalized_likelihood_branch_rate_tsv(
    path: Path,
    report: PenalizedLikelihoodDatingReport,
) -> Path:
    """Write one branch rate row per edge in tree preorder."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "branch_id",
        "child_name",
        "descendant_taxa",
        "parent_date",
        "child_date",
        "fitted_time_duration",
        "observed_branch_length",
        "observed_log_rate",
        "fitted_log_rate",
        "estimated_branch_rate",
        "fitted_branch_length",
        "data_score_contribution",
        "smoothing_penalty_contribution",
    ]
    lines = ["\t".join(columns)]
    lines.extend(
        "\t".join(
            [
                row.branch_id,
                row.child_name or "",
                "|".join(row.descendant_taxa),
                format(row.parent_date, ".15g"),
                format(row.child_date, ".15g"),
                format(row.fitted_time_duration, ".15g"),
                format(row.observed_branch_length, ".15g"),
                format(row.observed_log_rate, ".15g"),
                format(row.fitted_log_rate, ".15g"),
                format(row.estimated_branch_rate, ".15g"),
                format(row.fitted_branch_length, ".15g"),
                format(row.data_score_contribution, ".15g"),
                format(row.smoothing_penalty_contribution, ".15g"),
            ]
        )
        for row in report.branch_rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_penalized_likelihood_dating_run_json(
    path: Path,
    report: PenalizedLikelihoodDatingReport,
) -> Path:
    """Write the full penalized likelihood dating report as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_penalized_likelihood_dating_artifacts(
    out_dir: Path,
    report: PenalizedLikelihoodDatingReport,
) -> dict[str, Path]:
    """Write governed artifact outputs for one penalized likelihood dating run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    dated_tree_path = write_newick(
        out_dir / "dated_tree.nwk",
        loads_newick(report.dated_tree_newick),
    )
    summary_path = write_penalized_likelihood_dating_summary_tsv(
        out_dir / "summary.tsv",
        report,
    )
    node_dates_path = write_penalized_likelihood_node_dates_tsv(
        out_dir / "node_dates.tsv",
        report,
    )
    branch_rates_path = write_penalized_likelihood_branch_rate_tsv(
        out_dir / "branch_rates.tsv",
        report,
    )
    run_json_path = write_penalized_likelihood_dating_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "dated_tree_path": dated_tree_path,
        "summary_path": summary_path,
        "node_dates_path": node_dates_path,
        "branch_rates_path": branch_rates_path,
        "run_json_path": run_json_path,
    }
