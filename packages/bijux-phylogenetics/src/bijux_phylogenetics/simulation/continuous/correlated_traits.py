from __future__ import annotations

from math import sqrt
from pathlib import Path
import random

from bijux_phylogenetics.io.trees import load_tree


def _normalize_trait_names(trait_names: list[str] | tuple[str, ...]) -> list[str]:
    normalized = [name.strip() for name in trait_names if name.strip()]
    if len(normalized) < 2:
        raise ValueError(
            "trait_names must contain at least two non-empty trait identifiers"
        )
    if len(set(normalized)) != len(normalized):
        raise ValueError("trait_names must be unique")
    return normalized


def _normalize_root_states(
    trait_names: list[str],
    root_states: list[float] | tuple[float, ...] | None,
) -> list[float]:
    if root_states is None:
        return [0.0 for _ in trait_names]
    if len(root_states) != len(trait_names):
        raise ValueError(
            "root_states length must match the number of trait_names when provided"
        )
    return [float(value) for value in root_states]


def _normalize_square_matrix(
    *,
    matrix: list[list[float]] | tuple[tuple[float, ...], ...],
    size: int,
    matrix_name: str,
) -> list[list[float]]:
    if len(matrix) != size:
        raise ValueError(
            f"{matrix_name} must contain exactly {size} rows for the requested traits"
        )
    normalized: list[list[float]] = []
    for row_index, row in enumerate(matrix, start=1):
        if len(row) != size:
            raise ValueError(
                f"{matrix_name} row {row_index} must contain exactly {size} values"
            )
        normalized.append([float(value) for value in row])
    return normalized


def _covariance_matrix_from_correlation_matrix(
    *,
    correlation_matrix: list[list[float]],
    trait_standard_deviations: list[float],
) -> list[list[float]]:
    covariance_matrix: list[list[float]] = []
    for row_index, row in enumerate(correlation_matrix):
        covariance_row: list[float] = []
        for column_index, correlation in enumerate(row):
            covariance_row.append(
                correlation
                * trait_standard_deviations[row_index]
                * trait_standard_deviations[column_index]
            )
        covariance_matrix.append(covariance_row)
    return covariance_matrix


def _resolve_correlated_brownian_covariance_matrix(
    *,
    trait_names: list[str],
    evolutionary_covariance_matrix: list[list[float]]
    | tuple[tuple[float, ...], ...]
    | None,
    evolutionary_correlation_matrix: list[list[float]]
    | tuple[tuple[float, ...], ...]
    | None,
    trait_standard_deviations: list[float] | tuple[float, ...] | None,
) -> list[list[float]]:
    if evolutionary_covariance_matrix is not None:
        return _normalize_square_matrix(
            matrix=evolutionary_covariance_matrix,
            size=len(trait_names),
            matrix_name="evolutionary_covariance_matrix",
        )
    if evolutionary_correlation_matrix is None:
        return [
            [
                1.0 if row_index == column_index else 0.0
                for column_index in range(len(trait_names))
            ]
            for row_index in range(len(trait_names))
        ]
    correlation_matrix = _normalize_square_matrix(
        matrix=evolutionary_correlation_matrix,
        size=len(trait_names),
        matrix_name="evolutionary_correlation_matrix",
    )
    if trait_standard_deviations is None or len(trait_standard_deviations) != len(
        trait_names
    ):
        raise ValueError(
            "trait_standard_deviations must be provided with one value per trait when using evolutionary_correlation_matrix"
        )
    return _covariance_matrix_from_correlation_matrix(
        correlation_matrix=correlation_matrix,
        trait_standard_deviations=[float(value) for value in trait_standard_deviations],
    )


def _cholesky_lower(matrix: list[list[float]]) -> list[list[float]]:
    size = len(matrix)
    lower = [[0.0 for _ in range(size)] for _ in range(size)]
    for row in range(size):
        for column in range(row + 1):
            value = matrix[row][column] - sum(
                lower[row][inner] * lower[column][inner] for inner in range(column)
            )
            if row == column:
                if value <= 0.0:
                    raise ValueError(
                        "correlated Brownian simulation requires a positive-definite evolutionary covariance matrix"
                    )
                lower[row][column] = sqrt(value)
            else:
                lower[row][column] = value / lower[column][column]
    return lower


def _simulate_correlated_brownian_node_values(
    tree,
    *,
    root_states: list[float],
    covariance_lower: list[list[float]],
    rng: random.Random,
) -> dict[str, list[float]]:
    from bijux_phylogenetics.ancestral.common import node_signature

    node_values: dict[str, list[float]] = {}

    def propagate(state: list[float], branch_length: float) -> list[float]:
        if branch_length == 0.0:
            return list(state)
        standard_normal = [rng.gauss(0.0, 1.0) for _ in state]
        scaled_noise = [sqrt(branch_length) * value for value in standard_normal]
        increments = [
            sum(
                covariance_lower[row_index][column_index] * scaled_noise[column_index]
                for column_index in range(row_index + 1)
            )
            for row_index in range(len(state))
        ]
        return [state[index] + increments[index] for index in range(len(state))]

    def visit(node, state: list[float]) -> None:
        node_values[node_signature(node)] = list(state)
        if node.is_leaf():
            return
        for child in node.children:
            branch_length = max(child.branch_length or 0.0, 0.0)
            visit(child, propagate(state, branch_length))

    visit(tree.root, list(root_states))
    return node_values


def _build_correlated_brownian_collection_summary_rows(
    trait_names: list[str],
    root_states: list[float],
    covariance_matrix: list[list[float]],
    simulations,
):
    from .._statistics import (
        _mean,
        _median,
        _round_float,
        _sample_correlation,
        _sample_covariance,
        _sample_standard_deviation,
    )
    from ..contracts import ContinuousTraitSimulationSummaryRow

    if not simulations:
        return []
    rows: list[ContinuousTraitSimulationSummaryRow] = []
    for trait_name, root_state in zip(trait_names, root_states, strict=True):
        rows.append(
            ContinuousTraitSimulationSummaryRow(
                row_kind="root_state",
                label=trait_name,
                mean_value=_round_float(root_state),
            )
        )
    for left_index, left_trait in enumerate(trait_names):
        left_variance = covariance_matrix[left_index][left_index]
        for right_index, right_trait in enumerate(
            trait_names[left_index:], start=left_index
        ):
            right_variance = covariance_matrix[right_index][right_index]
            covariance = covariance_matrix[left_index][right_index]
            correlation = 0.0
            if left_variance > 0.0 and right_variance > 0.0:
                correlation = covariance / sqrt(left_variance * right_variance)
            rows.append(
                ContinuousTraitSimulationSummaryRow(
                    row_kind="evolutionary_covariance",
                    label=f"{left_trait}|{right_trait}",
                    covariance=_round_float(covariance),
                    correlation=_round_float(correlation),
                )
            )
    dimension_labels = [
        (taxon, trait_name)
        for taxon in [row.taxon for row in simulations[0].traits[:: len(trait_names)]]
        for trait_name in trait_names
    ]
    values_by_dimension = {
        dimension_label: [
            next(
                row.value
                for row in simulation.traits
                if row.taxon == dimension_label[0] and row.trait == dimension_label[1]
            )
            for simulation in simulations
        ]
        for dimension_label in dimension_labels
    }
    for dimension_label in dimension_labels:
        values = values_by_dimension[dimension_label]
        rows.append(
            ContinuousTraitSimulationSummaryRow(
                row_kind="tip_distribution",
                label=f"{dimension_label[0]}|{dimension_label[1]}",
                mean_value=_mean(values),
                standard_deviation=_sample_standard_deviation(values),
                minimum=_round_float(min(values)),
                median=_median(values),
                maximum=_round_float(max(values)),
            )
        )
    for left_index, left_label in enumerate(dimension_labels):
        left_values = values_by_dimension[left_label]
        for right_label in dimension_labels[left_index:]:
            right_values = values_by_dimension[right_label]
            rows.append(
                ContinuousTraitSimulationSummaryRow(
                    row_kind="tip_covariance",
                    label=(
                        f"{left_label[0]}|{left_label[1]}||"
                        f"{right_label[0]}|{right_label[1]}"
                    ),
                    covariance=_sample_covariance(left_values, right_values),
                    correlation=_sample_correlation(left_values, right_values),
                )
            )
    return sorted(rows, key=lambda row: (row.row_kind, row.label))


def simulate_correlated_brownian_traits(
    tree_path: Path,
    *,
    trait_names: list[str] | tuple[str, ...],
    evolutionary_covariance_matrix: list[list[float]]
    | tuple[tuple[float, ...], ...]
    | None = None,
    evolutionary_correlation_matrix: list[list[float]]
    | tuple[tuple[float, ...], ...]
    | None = None,
    trait_standard_deviations: list[float] | tuple[float, ...] | None = None,
    root_states: list[float] | tuple[float, ...] | None = None,
    seed: int = 1,
):
    from bijux_phylogenetics.ancestral.common import node_signature

    from .._statistics import _round_float
    from ..contracts import (
        CorrelatedContinuousTraitSimulationReport,
        SimulatedCorrelatedContinuousTrait,
    )

    normalized_trait_names = _normalize_trait_names(trait_names)
    normalized_root_states = _normalize_root_states(
        normalized_trait_names,
        root_states,
    )
    covariance_matrix = _resolve_correlated_brownian_covariance_matrix(
        trait_names=normalized_trait_names,
        evolutionary_covariance_matrix=evolutionary_covariance_matrix,
        evolutionary_correlation_matrix=evolutionary_correlation_matrix,
        trait_standard_deviations=trait_standard_deviations,
    )
    covariance_lower = _cholesky_lower(covariance_matrix)
    tree = load_tree(tree_path)
    node_values = _simulate_correlated_brownian_node_values(
        tree,
        root_states=normalized_root_states,
        covariance_lower=covariance_lower,
        rng=random.Random(seed),  # nosec B311
    )
    traits: list[SimulatedCorrelatedContinuousTrait] = []
    for node in tree.iter_leaves():
        if node.name is None:
            continue
        values = node_values[node_signature(node)]
        for trait_name, value in zip(normalized_trait_names, values, strict=True):
            traits.append(
                SimulatedCorrelatedContinuousTrait(
                    taxon=node.name,
                    trait=trait_name,
                    value=_round_float(value),
                )
            )
    return CorrelatedContinuousTraitSimulationReport(
        model="multivariate-brownian-motion",
        tree_path=tree_path,
        tip_count=tree.tip_count,
        trait_names=list(normalized_trait_names),
        seed=seed,
        root_states=list(normalized_root_states),
        evolutionary_covariance_matrix=[
            [_round_float(value) for value in row] for row in covariance_matrix
        ],
        traits=sorted(traits, key=lambda row: (row.taxon, row.trait)),
    )


def simulate_correlated_brownian_trait_collection(
    tree_path: Path,
    *,
    trait_names: list[str] | tuple[str, ...],
    evolutionary_covariance_matrix: list[list[float]]
    | tuple[tuple[float, ...], ...]
    | None = None,
    evolutionary_correlation_matrix: list[list[float]]
    | tuple[tuple[float, ...], ...]
    | None = None,
    trait_standard_deviations: list[float] | tuple[float, ...] | None = None,
    root_states: list[float] | tuple[float, ...] | None = None,
    replicates: int = 128,
    seed: int = 1,
):
    from .._statistics import (
        _round_float,
    )
    from ..contracts import CorrelatedContinuousTraitSimulationCollectionReport

    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    normalized_trait_names = _normalize_trait_names(trait_names)
    normalized_root_states = _normalize_root_states(
        normalized_trait_names,
        root_states,
    )
    covariance_matrix = _resolve_correlated_brownian_covariance_matrix(
        trait_names=normalized_trait_names,
        evolutionary_covariance_matrix=evolutionary_covariance_matrix,
        evolutionary_correlation_matrix=evolutionary_correlation_matrix,
        trait_standard_deviations=trait_standard_deviations,
    )
    tree = load_tree(tree_path)
    simulations = [
        simulate_correlated_brownian_traits(
            tree_path,
            trait_names=normalized_trait_names,
            evolutionary_covariance_matrix=covariance_matrix,
            root_states=normalized_root_states,
            seed=seed + index - 1,
        )
        for index in range(1, replicates + 1)
    ]
    return CorrelatedContinuousTraitSimulationCollectionReport(
        model="multivariate-brownian-motion",
        tree_path=tree_path,
        tip_count=tree.tip_count,
        branch_count=sum(1 for _ in tree.iter_edges()),
        trait_names=list(normalized_trait_names),
        replicate_count=replicates,
        seed=seed,
        root_states=list(normalized_root_states),
        evolutionary_covariance_matrix=[
            [_round_float(value) for value in row] for row in covariance_matrix
        ],
        simulations=simulations,
        rows=_build_correlated_brownian_collection_summary_rows(
            normalized_trait_names,
            normalized_root_states,
            covariance_matrix,
            simulations,
        ),
    )


def write_correlated_continuous_trait_table(path: Path, report) -> Path:
    from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

    return write_taxon_rows(
        path,
        columns=["taxon", "trait", "value"],
        rows=[
            {
                "taxon": row.taxon,
                "trait": row.trait,
                "value": format(row.value, ".15g"),
            }
            for row in report.traits
        ],
    )


def write_correlated_continuous_trait_collection_table(path: Path, report) -> Path:
    from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

    return write_taxon_rows(
        path,
        columns=["replicate_index", "taxon", "trait", "value"],
        rows=[
            {
                "replicate_index": str(replicate_index),
                "taxon": row.taxon,
                "trait": row.trait,
                "value": format(row.value, ".15g"),
            }
            for replicate_index, simulation in enumerate(report.simulations, start=1)
            for row in simulation.traits
        ],
    )


def write_correlated_continuous_trait_collection_summary_table(
    path: Path, report
) -> Path:
    from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

    return write_taxon_rows(
        path,
        columns=[
            "row_kind",
            "label",
            "mean_value",
            "standard_deviation",
            "minimum",
            "median",
            "maximum",
            "covariance",
            "correlation",
        ],
        rows=[
            {
                "row_kind": row.row_kind,
                "label": row.label,
                "mean_value": (
                    "" if row.mean_value is None else format(row.mean_value, ".15g")
                ),
                "standard_deviation": (
                    ""
                    if row.standard_deviation is None
                    else format(row.standard_deviation, ".15g")
                ),
                "minimum": "" if row.minimum is None else format(row.minimum, ".15g"),
                "median": "" if row.median is None else format(row.median, ".15g"),
                "maximum": "" if row.maximum is None else format(row.maximum, ".15g"),
                "covariance": (
                    "" if row.covariance is None else format(row.covariance, ".15g")
                ),
                "correlation": (
                    "" if row.correlation is None else format(row.correlation, ".15g")
                ),
            }
            for row in report.rows
        ],
    )
