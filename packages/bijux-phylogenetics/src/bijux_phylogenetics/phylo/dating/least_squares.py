from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict
import json
import math
from pathlib import Path

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    matrix_condition_number,
    matrix_vector_multiply,
)
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.newick import dumps_newick, loads_newick, write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    MetadataJoinError,
    PhylogeneticsError,
    UnrootedTreeError,
)

from .inputs import load_tip_dates_for_tree, validate_tip_dates_against_tree
from .models import (
    LeastSquaresDatingBranchRow,
    LeastSquaresDatingNodeRow,
    LeastSquaresDatingReport,
)

_DATE_TOLERANCE = 1e-9


def fit_least_squares_dating(
    tree: PhyloTree,
    tip_dates: Mapping[str, float],
    *,
    tree_path: Path | None = None,
    metadata_path: Path | None = None,
    taxon_column: str = "taxon",
    date_column: str = "date",
) -> LeastSquaresDatingReport:
    """Fit one rooted tree to fixed tip dates by closed-form least squares."""
    if tree.rooted is not True:
        raise UnrootedTreeError(
            "least-squares dating requires one rooted tree",
            code="least_squares_dating_requires_rooted_tree",
        )
    if any(length is None for length in tree.branch_lengths()):
        raise InvalidBranchLengthError(
            "least-squares dating requires complete branch lengths"
        )
    validate_tip_dates_against_tree(tree, tip_dates)
    internal_nodes = list(tree.iter_internal_nodes(order="preorder"))
    if not internal_nodes:
        raise PhylogeneticsError(
            "least-squares dating requires at least one internal node",
            code="least_squares_dating_error",
        )
    if len(set(tip_dates.values())) < 2:
        raise PhylogeneticsError(
            "least-squares dating requires variation in tip dates",
            code="least_squares_dating_error",
        )

    coefficient_matrix, observations, internal_index = _build_linear_system(
        tree,
        tip_dates=tip_dates,
        internal_nodes=internal_nodes,
    )
    parameter_values, condition_number = _solve_linear_least_squares(
        coefficient_matrix,
        observations,
    )
    estimated_clock_rate = parameter_values[0]
    if estimated_clock_rate <= 0.0:
        raise PhylogeneticsError(
            "least-squares dating fit requires a positive clock rate",
            code="least_squares_dating_error",
        )

    node_dates = {
        node_id: parameter_values[index] / estimated_clock_rate
        for node_id, index in internal_index.items()
    }
    for tip_name, tip_date in tip_dates.items():
        node_dates[_tip_node_id(tree, tip_name)] = tip_date
    dated_tree = _build_dated_tree(tree, node_dates)
    branch_rows = _build_branch_rows(
        tree,
        node_dates=node_dates,
        estimated_clock_rate=estimated_clock_rate,
    )
    residual_sum_squares = sum(row.residual * row.residual for row in branch_rows)
    maximum_tip_date = max(tip_dates.values())
    node_rows = _build_node_rows(
        tree,
        node_dates=node_dates,
        tip_dates=tip_dates,
        maximum_tip_date=maximum_tip_date,
    )
    return LeastSquaresDatingReport(
        tree_newick=dumps_newick(tree),
        dated_tree_newick=dumps_newick(dated_tree),
        taxa=sorted(tree.tip_names),
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        branch_count=len(branch_rows),
        parameter_count=len(parameter_values),
        tree_path=None if tree_path is None else str(tree_path),
        metadata_path=None if metadata_path is None else str(metadata_path),
        taxon_column=taxon_column,
        date_column=date_column,
        minimum_tip_date=min(tip_dates.values()),
        maximum_tip_date=maximum_tip_date,
        root_date=node_dates[tree.root.node_id or ""],
        estimated_clock_rate=estimated_clock_rate,
        residual_sum_squares=residual_sum_squares,
        condition_number=condition_number,
        exact_fit=math.isclose(residual_sum_squares, 0.0, abs_tol=1e-12),
        optimizer_name="closed-form-linear-least-squares",
        converged=True,
        node_rows=node_rows,
        branch_rows=branch_rows,
    )


def fit_least_squares_dating_from_metadata(
    tree_path: Path,
    metadata_path: Path,
    *,
    taxon_column: str | None = None,
    date_column: str = "date",
) -> LeastSquaresDatingReport:
    """Fit least-squares node dates from one rooted substitution tree and one tip-date table."""
    validate_tree_path(tree_path, require_rooted=True)
    tree = load_tree(tree_path)
    tree.rooted = True
    if any(length is None for length in tree.branch_lengths()):
        raise InvalidBranchLengthError(
            "least-squares dating requires complete branch lengths"
        )
    tip_dates, resolved_taxon_column = load_tip_dates_for_tree(
        metadata_path,
        tree_taxa=tree.tip_names,
        taxon_column=taxon_column,
        date_column=date_column,
    )
    return fit_least_squares_dating(
        tree,
        tip_dates,
        tree_path=tree_path,
        metadata_path=metadata_path,
        taxon_column=resolved_taxon_column,
        date_column=date_column,
    )


def _build_linear_system(
    tree: PhyloTree,
    *,
    tip_dates: Mapping[str, float],
    internal_nodes: list[TreeNode],
) -> tuple[list[list[float]], list[float], dict[str, int]]:
    internal_index = {
        node.node_id or "": index for index, node in enumerate(internal_nodes, start=1)
    }
    branch_count = sum(1 for _parent, _child in tree.iter_edges())
    coefficient_matrix = [
        [0.0 for _ in range(len(internal_nodes) + 1)] for _ in range(branch_count)
    ]
    observations = [0.0 for _ in range(branch_count)]
    for row_index, (parent, child) in enumerate(tree.iter_edges()):
        observations[row_index] = float(child.branch_length or 0.0)
        coefficient_matrix[row_index][internal_index[parent.node_id or ""]] = -1.0
        if child.is_leaf():
            if child.name is None:
                raise MetadataJoinError("least-squares dating requires named tips")
            coefficient_matrix[row_index][0] = float(tip_dates[child.name])
        else:
            coefficient_matrix[row_index][internal_index[child.node_id or ""]] = 1.0
    return coefficient_matrix, observations, internal_index


def _solve_linear_least_squares(
    coefficient_matrix: list[list[float]],
    observations: list[float],
) -> tuple[list[float], float]:
    transposed = list(zip(*coefficient_matrix, strict=False))
    normal_matrix = [
        [
            sum(
                transposed[row_index][k] * coefficient_matrix[k][column_index]
                for k in range(len(coefficient_matrix))
            )
            for column_index in range(len(coefficient_matrix[0]))
        ]
        for row_index in range(len(coefficient_matrix[0]))
    ]
    normal_vector = [
        sum(
            transposed[row_index][k] * observations[k]
            for k in range(len(coefficient_matrix))
        )
        for row_index in range(len(coefficient_matrix[0]))
    ]
    try:
        inverse = invert_matrix(normal_matrix)
    except ValueError as error:
        raise PhylogeneticsError(
            "least-squares dating system is singular for the supplied tree and tip dates",
            code="least_squares_dating_error",
        ) from error
    parameter_values = matrix_vector_multiply(inverse, normal_vector)
    return parameter_values, matrix_condition_number(normal_matrix)


def _tip_node_id(tree: PhyloTree, tip_name: str) -> str:
    for node in tree.iter_leaves():
        if node.name == tip_name:
            if node.node_id is None:
                raise PhylogeneticsError(
                    "least-squares dating requires stable node ids",
                    code="least_squares_dating_error",
                )
            return node.node_id
    raise MetadataJoinError(f"tip-date taxon '{tip_name}' is absent from the tree")


def _build_dated_tree(tree: PhyloTree, node_dates: Mapping[str, float]) -> PhyloTree:
    dated_tree = tree.copy()
    for parent, child in dated_tree.iter_edges():
        parent_node_id = parent.node_id or ""
        child_node_id = child.node_id or ""
        duration = node_dates[child_node_id] - node_dates[parent_node_id]
        if duration < -_DATE_TOLERANCE:
            raise PhylogeneticsError(
                "least-squares dating fit produced non-chronological node ages",
                code="least_squares_dating_error",
            )
        child.branch_length = max(duration, 0.0)
    return dated_tree


def _build_node_rows(
    tree: PhyloTree,
    *,
    node_dates: Mapping[str, float],
    tip_dates: Mapping[str, float],
    maximum_tip_date: float,
) -> list[LeastSquaresDatingNodeRow]:
    rows: list[LeastSquaresDatingNodeRow] = []
    for node in tree.iter_nodes(order="preorder"):
        if node.node_id is None:
            raise PhylogeneticsError(
                "least-squares dating requires stable node ids",
                code="least_squares_dating_error",
            )
        if node is tree.root:
            node_kind = "root"
        elif node.is_leaf():
            node_kind = "tip"
        else:
            node_kind = "internal"
        estimated_date = node_dates[node.node_id]
        rows.append(
            LeastSquaresDatingNodeRow(
                node_id=node.node_id,
                node_kind=node_kind,
                node_label=node.name,
                descendant_taxa=node.descendant_taxa,
                estimated_date=estimated_date,
                fixed_tip_date=node.is_leaf() and node.name in tip_dates,
                time_height=maximum_tip_date - estimated_date,
            )
        )
    return rows


def _build_branch_rows(
    tree: PhyloTree,
    *,
    node_dates: Mapping[str, float],
    estimated_clock_rate: float,
) -> list[LeastSquaresDatingBranchRow]:
    rows: list[LeastSquaresDatingBranchRow] = []
    for parent, child in tree.iter_edges():
        if parent.node_id is None or child.node_id is None:
            raise PhylogeneticsError(
                "least-squares dating requires stable node ids",
                code="least_squares_dating_error",
            )
        parent_date = node_dates[parent.node_id]
        child_date = node_dates[child.node_id]
        fitted_time_duration = child_date - parent_date
        if fitted_time_duration < -_DATE_TOLERANCE:
            raise PhylogeneticsError(
                "least-squares dating fit produced non-chronological branch durations",
                code="least_squares_dating_error",
            )
        fitted_time_duration = max(fitted_time_duration, 0.0)
        observed_branch_length = float(child.branch_length or 0.0)
        fitted_branch_length = estimated_clock_rate * fitted_time_duration
        rows.append(
            LeastSquaresDatingBranchRow(
                branch_id=child.node_id,
                child_name=child.name,
                descendant_taxa=child.descendant_taxa,
                parent_date=parent_date,
                child_date=child_date,
                fitted_time_duration=fitted_time_duration,
                observed_branch_length=observed_branch_length,
                fitted_branch_length=fitted_branch_length,
                residual=observed_branch_length - fitted_branch_length,
            )
        )
    return rows


def write_least_squares_dating_summary_tsv(
    path: Path,
    report: LeastSquaresDatingReport,
) -> Path:
    """Write one summary row for one least-squares dating run."""
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
        "estimated_clock_rate",
        "residual_sum_squares",
        "condition_number",
        "exact_fit",
        "optimizer_name",
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
        format(report.estimated_clock_rate, ".15g"),
        format(report.residual_sum_squares, ".15g"),
        format(report.condition_number, ".15g"),
        str(report.exact_fit).lower(),
        report.optimizer_name,
        str(report.converged).lower(),
    ]
    path.write_text(
        "\n".join(["\t".join(columns), "\t".join(values)]) + "\n",
        encoding="utf-8",
    )
    return path


def write_least_squares_node_dates_tsv(
    path: Path,
    report: LeastSquaresDatingReport,
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
            ]
        )
        for row in report.node_rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_least_squares_branch_residuals_tsv(
    path: Path,
    report: LeastSquaresDatingReport,
) -> Path:
    """Write one branch residual row per edge in tree preorder."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "branch_id",
        "child_name",
        "descendant_taxa",
        "parent_date",
        "child_date",
        "fitted_time_duration",
        "observed_branch_length",
        "fitted_branch_length",
        "residual",
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
                format(row.fitted_branch_length, ".15g"),
                format(row.residual, ".15g"),
            ]
        )
        for row in report.branch_rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_least_squares_dating_run_json(
    path: Path,
    report: LeastSquaresDatingReport,
) -> Path:
    """Write the full least-squares dating report as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_least_squares_dating_artifacts(
    out_dir: Path,
    report: LeastSquaresDatingReport,
) -> dict[str, Path]:
    """Write governed artifact outputs for one least-squares dating run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    dated_tree_path = write_newick(
        out_dir / "dated_tree.nwk",
        loads_newick(report.dated_tree_newick),
    )
    summary_path = write_least_squares_dating_summary_tsv(
        out_dir / "summary.tsv",
        report,
    )
    node_dates_path = write_least_squares_node_dates_tsv(
        out_dir / "node_dates.tsv",
        report,
    )
    branch_residuals_path = write_least_squares_branch_residuals_tsv(
        out_dir / "branch_residuals.tsv",
        report,
    )
    run_json_path = write_least_squares_dating_run_json(out_dir / "run.json", report)
    return {
        "dated_tree_path": dated_tree_path,
        "summary_path": summary_path,
        "node_dates_path": node_dates_path,
        "branch_residuals_path": branch_residuals_path,
        "run_json_path": run_json_path,
    }
