from __future__ import annotations

from math import exp, sqrt
from pathlib import Path
import random

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


def simulate_brownian_traits(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float | None = None,
    sigma_squared: float | None = None,
    seed: int = 1,
):
    from .._state_propagation import (
        _resolve_brownian_sigma_parameters,
        _simulate_brownian_node_values,
    )

    sigma, sigma_squared = _resolve_brownian_sigma_parameters(
        sigma=sigma,
        sigma_squared=sigma_squared,
    )
    tree = load_tree(tree_path)
    rng = random.Random(seed)  # nosec B311
    node_values = _simulate_brownian_node_values(
        tree,
        root_state=root_state,
        sigma=sigma,
        rng=rng,
    )
    return _build_continuous_trait_simulation_report(
        tree=tree,
        tree_path=tree_path,
        model="brownian-motion",
        seed=seed,
        root_state=root_state,
        sigma=sigma,
        sigma_squared=sigma_squared,
        alpha=None,
        theta=None,
        rate_change=None,
        node_values=node_values,
    )


def simulate_ou_traits(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float = 1.0,
    alpha: float = 1.0,
    theta: float = 0.0,
    seed: int = 1,
):
    from .._state_propagation import _iter_node_trait_values

    if sigma < 0.0:
        raise ValueError(f"sigma must be nonnegative, got {sigma}")
    if alpha < 0.0:
        raise ValueError(f"alpha must be nonnegative, got {alpha}")
    tree = load_tree(tree_path)
    rng = random.Random(seed)  # nosec B311

    def propagate(state: float, branch_length: float) -> float:
        if branch_length == 0.0:
            return state
        if alpha == 0.0:
            return state + rng.gauss(0.0, sigma * sqrt(branch_length))
        mean = theta + (state - theta) * exp(-alpha * branch_length)
        variance = (
            (sigma**2) * (1.0 - exp(-2.0 * alpha * branch_length)) / (2.0 * alpha)
        )
        return mean + rng.gauss(0.0, sqrt(max(variance, 0.0)))

    node_values = _iter_node_trait_values(
        tree, root_state=root_state, propagate=propagate
    )
    return _build_continuous_trait_simulation_report(
        tree=tree,
        tree_path=tree_path,
        model="ornstein-uhlenbeck",
        seed=seed,
        root_state=root_state,
        sigma=sigma,
        sigma_squared=sigma * sigma,
        alpha=alpha,
        theta=theta,
        rate_change=None,
        node_values=node_values,
    )


def simulate_early_burst_traits(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float = 1.0,
    rate_change: float = 1.0,
    seed: int = 1,
):
    from bijux_phylogenetics.comparative.evolutionary_modes import (
        transform_tree_for_evolutionary_mode,
    )

    from .._state_propagation import _iter_node_trait_values

    if sigma < 0.0:
        raise ValueError(f"sigma must be nonnegative, got {sigma}")
    if rate_change < 0.0:
        raise ValueError(f"rate_change must be nonnegative, got {rate_change}")
    tree = load_tree(tree_path)
    transformed_tree = transform_tree_for_evolutionary_mode(
        tree,
        mode="early-burst",
        parameter_value=rate_change,
    )
    rng = random.Random(seed)  # nosec B311
    node_values = _iter_node_trait_values(
        transformed_tree,
        root_state=root_state,
        propagate=lambda state, branch_length: (
            state + rng.gauss(0.0, sigma * sqrt(branch_length))
        ),
    )
    return _build_continuous_trait_simulation_report(
        tree=transformed_tree,
        tree_path=tree_path,
        model="early-burst",
        seed=seed,
        root_state=root_state,
        sigma=sigma,
        sigma_squared=sigma * sigma,
        alpha=None,
        theta=None,
        rate_change=rate_change,
        node_values=node_values,
    )


def simulate_speciational_traits(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float | None = None,
    sigma_squared: float | None = None,
    seed: int = 1,
):
    from .._state_propagation import (
        _resolve_brownian_sigma_parameters,
        _simulate_brownian_node_values,
    )

    sigma, sigma_squared = _resolve_brownian_sigma_parameters(
        sigma=sigma,
        sigma_squared=sigma_squared,
    )
    tree = _build_speciational_tree(load_tree(tree_path))
    rng = random.Random(seed)  # nosec B311
    node_values = _simulate_brownian_node_values(
        tree,
        root_state=root_state,
        sigma=sigma,
        rng=rng,
    )
    return _build_continuous_trait_simulation_report(
        tree=tree,
        tree_path=tree_path,
        model="speciational",
        seed=seed,
        root_state=root_state,
        sigma=sigma,
        sigma_squared=sigma_squared,
        alpha=None,
        theta=None,
        rate_change=None,
        node_values=node_values,
    )


def _build_continuous_trait_simulation_report(
    *,
    tree: PhyloTree,
    tree_path: Path,
    model: str,
    seed: int,
    root_state: float,
    sigma: float,
    sigma_squared: float,
    alpha: float | None,
    theta: float | None,
    rate_change: float | None,
    node_values: dict[str, float],
):
    from bijux_phylogenetics.ancestral.common import (
        node_descendant_taxa,
        node_signature,
    )

    from .._state_propagation import _tip_values_from_node_map
    from ..contracts import (
        ContinuousTraitSimulationReport,
        SimulatedContinuousNode,
        SimulatedContinuousTrait,
    )

    values = _tip_values_from_node_map(tree, node_values)
    return ContinuousTraitSimulationReport(
        model=model,
        tree_path=tree_path,
        tip_count=tree.tip_count,
        seed=seed,
        root_state=root_state,
        sigma=sigma,
        sigma_squared=sigma_squared,
        alpha=alpha,
        theta=theta,
        rate_change=rate_change,
        traits=[
            SimulatedContinuousTrait(taxon=taxon, value=value)
            for taxon, value in sorted(values.items())
        ],
        node_values=[
            SimulatedContinuousNode(
                node=node_signature(node),
                node_name=node.name,
                is_tip=node.is_leaf(),
                descendant_taxa=node_descendant_taxa(node),
                value=float(format(node_values[node_signature(node)], ".15g")),
            )
            for node in tree.iter_nodes()
        ],
    )


def _build_speciational_tree(tree: PhyloTree) -> PhyloTree:
    transformed_tree = tree.copy()
    for _, child in transformed_tree.iter_edges():
        branch_length = float(child.branch_length or 0.0)
        if branch_length < 0.0:
            raise ValueError(
                "speciational simulation requires nonnegative branch lengths"
            )
        child.branch_length = 0.0 if branch_length == 0.0 else 1.0
    return transformed_tree


def _build_continuous_collection_summary_rows(simulations):
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
    taxa = [row.taxon for row in simulations[0].traits]
    values_by_taxon = {
        taxon: [simulation.traits[index].value for simulation in simulations]
        for index, taxon in enumerate(taxa)
    }
    rows: list[ContinuousTraitSimulationSummaryRow] = []
    for taxon in taxa:
        values = values_by_taxon[taxon]
        rows.append(
            ContinuousTraitSimulationSummaryRow(
                row_kind="tip_distribution",
                label=taxon,
                mean_value=_mean(values),
                standard_deviation=_sample_standard_deviation(values),
                minimum=_round_float(min(values)),
                median=_median(values),
                maximum=_round_float(max(values)),
            )
        )
    for left_index, left_taxon in enumerate(taxa):
        for right_taxon in taxa[left_index:]:
            left_values = values_by_taxon[left_taxon]
            right_values = values_by_taxon[right_taxon]
            rows.append(
                ContinuousTraitSimulationSummaryRow(
                    row_kind="tip_covariance",
                    label=f"{left_taxon}|{right_taxon}",
                    covariance=_sample_covariance(left_values, right_values),
                    correlation=_sample_correlation(left_values, right_values),
                )
            )
    return sorted(rows, key=lambda row: (row.row_kind, row.label))


def simulate_brownian_trait_collection(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float | None = None,
    sigma_squared: float | None = None,
    replicates: int = 128,
    seed: int = 1,
):
    from .._state_propagation import (
        _resolve_brownian_sigma_parameters,
        _simulate_brownian_node_values,
    )
    from ..contracts import ContinuousTraitSimulationCollectionReport

    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    sigma, sigma_squared = _resolve_brownian_sigma_parameters(
        sigma=sigma,
        sigma_squared=sigma_squared,
    )
    tree = load_tree(tree_path)
    simulations = [
        _build_continuous_trait_simulation_report(
            tree=tree,
            tree_path=tree_path,
            model="brownian-motion",
            seed=seed + index - 1,
            root_state=root_state,
            sigma=sigma,
            sigma_squared=sigma_squared,
            alpha=None,
            theta=None,
            rate_change=None,
            node_values=_simulate_brownian_node_values(
                tree,
                root_state=root_state,
                sigma=sigma,
                rng=random.Random(seed + index - 1),  # nosec B311
            ),
        )
        for index in range(1, replicates + 1)
    ]
    return ContinuousTraitSimulationCollectionReport(
        model="brownian-motion",
        tree_path=tree_path,
        tip_count=tree.tip_count,
        branch_count=sum(1 for _ in tree.iter_edges()),
        replicate_count=replicates,
        seed=seed,
        root_state=root_state,
        sigma=sigma,
        sigma_squared=sigma_squared,
        simulations=simulations,
        rows=_build_continuous_collection_summary_rows(simulations),
    )


def simulate_speciational_trait_collection(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float | None = None,
    sigma_squared: float | None = None,
    replicates: int = 128,
    seed: int = 1,
):
    from .._state_propagation import (
        _resolve_brownian_sigma_parameters,
        _simulate_brownian_node_values,
    )
    from ..contracts import ContinuousTraitSimulationCollectionReport

    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    sigma, sigma_squared = _resolve_brownian_sigma_parameters(
        sigma=sigma,
        sigma_squared=sigma_squared,
    )
    tree = _build_speciational_tree(load_tree(tree_path))
    simulations = [
        _build_continuous_trait_simulation_report(
            tree=tree,
            tree_path=tree_path,
            model="speciational",
            seed=seed + index - 1,
            root_state=root_state,
            sigma=sigma,
            sigma_squared=sigma_squared,
            alpha=None,
            theta=None,
            rate_change=None,
            node_values=_simulate_brownian_node_values(
                tree,
                root_state=root_state,
                sigma=sigma,
                rng=random.Random(seed + index - 1),  # nosec B311
            ),
        )
        for index in range(1, replicates + 1)
    ]
    return ContinuousTraitSimulationCollectionReport(
        model="speciational",
        tree_path=tree_path,
        tip_count=tree.tip_count,
        branch_count=sum(1 for _ in tree.iter_edges()),
        replicate_count=replicates,
        seed=seed,
        root_state=root_state,
        sigma=sigma,
        sigma_squared=sigma_squared,
        simulations=simulations,
        rows=_build_continuous_collection_summary_rows(simulations),
    )


def write_continuous_trait_table(path: Path, report) -> Path:
    from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

    return write_taxon_rows(
        path,
        columns=["taxon", "value"],
        rows=[
            {"taxon": row.taxon, "value": format(row.value, ".15g")}
            for row in report.traits
        ],
    )


def write_continuous_trait_collection_table(path: Path, report) -> Path:
    from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

    return write_taxon_rows(
        path,
        columns=["replicate_index", "taxon", "value"],
        rows=[
            {
                "replicate_index": str(replicate_index),
                "taxon": row.taxon,
                "value": format(row.value, ".15g"),
            }
            for replicate_index, simulation in enumerate(report.simulations, start=1)
            for row in simulation.traits
        ],
    )


def write_continuous_trait_collection_summary_table(path: Path, report) -> Path:
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
