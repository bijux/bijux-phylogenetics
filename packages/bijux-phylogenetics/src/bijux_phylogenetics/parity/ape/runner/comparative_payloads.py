from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative import compute_diversification_gamma_statistic
from bijux_phylogenetics.comparative.covariance import summarize_brownian_covariance
from bijux_phylogenetics.comparative.signal import (
    compute_phylogenetic_independent_contrasts,
)
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_tree_simulation_fixture,
)
from bijux_phylogenetics.phylo.branch_lengths.branching_times import (
    compute_tree_branching_times,
)
from bijux_phylogenetics.phylo.branch_lengths.node_depths import (
    compute_tree_node_depths,
)
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    assess_tree_ultrametricity,
)
from bijux_phylogenetics.simulation import (
    simulate_coalescent_trees,
    simulate_random_trees,
)


def _build_bijux_brownian_covariance_rows(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = summarize_brownian_covariance(input_fixture)
    return {
        "tip_count": len(report.taxa),
        "rooted": report.tree_is_rooted,
        "tip_labels": report.taxa,
        "pair_count": len(report.rows),
        "tree_is_ultrametric": report.tree_is_ultrametric,
        "minimum_root_to_tip_depth": report.minimum_root_to_tip_depth,
        "maximum_root_to_tip_depth": report.maximum_root_to_tip_depth,
        "minimum_branch_length": report.minimum_branch_length,
        "maximum_branch_length": report.maximum_branch_length,
        "matrix_dimension": report.matrix_dimension,
        "matrix_rank": report.matrix_rank,
        "singular": report.singular,
        "near_singular": report.near_singular,
        "positive_definite": report.positive_definite,
        "condition_number": (
            None if math.isinf(report.condition_number) else report.condition_number
        ),
        "raw_log_determinant": report.raw_log_determinant,
    }, [
        {
            "left_taxon": row.left_taxon,
            "right_taxon": row.right_taxon,
            "shared_ancestry_covariance": row.shared_ancestry_covariance,
        }
        for row in report.rows
    ]


def _build_bijux_independent_contrast_rows(
    input_fixture: Path,
    *,
    trait_table_path: Path,
    trait_name: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_phylogenetic_independent_contrasts(
        input_fixture,
        trait_table_path,
        trait=trait_name,
    )
    rows = sorted(
        [
            {
                "node_id": row.node_id,
                "node": row.node,
                "left_taxa": "|".join(row.left_taxa),
                "right_taxa": "|".join(row.right_taxa),
                "contrast": row.contrast,
                "expected_variance": row.expected_variance,
            }
            for row in report.contrasts
        ],
        key=lambda row: int(row["node_id"]),
    )
    return {
        "trait": report.trait,
        "taxon_count": report.taxon_count,
        "contrast_count": len(report.contrasts),
        "tree_is_ultrametric": report.input_audit.tree_is_ultrametric,
        "minimum_root_to_tip_depth": report.input_audit.minimum_root_to_tip_depth,
        "maximum_root_to_tip_depth": report.input_audit.maximum_root_to_tip_depth,
    }, rows


def _build_bijux_tree_node_depth_rows(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_tree_node_depths(input_fixture)
    return {
        "node_count": report.node_count,
        "tip_count": report.tip_count,
        "internal_node_count": report.internal_node_count,
        "rooted": report.rooted,
        "tip_labels": report.tip_labels,
        "tree_is_ultrametric": report.tree_is_ultrametric,
        "zero_branch_length_count": report.zero_branch_length_count,
        "minimum_tip_depth": report.minimum_tip_depth,
        "maximum_tip_depth": report.maximum_tip_depth,
        "minimum_internal_depth": report.minimum_internal_depth,
        "maximum_internal_depth": report.maximum_internal_depth,
    }, [
        {
            "node_id": row.node_id,
            "node_kind": row.node_kind,
            "node_label": row.node_label or "",
            "descendant_taxa": "|".join(row.descendant_taxa),
            "branch_length_depth": row.branch_length_depth,
            "branch_length": "" if row.branch_length is None else row.branch_length,
        }
        for row in report.rows
    ]


def _build_bijux_branching_time_rows(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_tree_branching_times(input_fixture)
    return {
        "internal_node_count": report.internal_node_count,
        "rooted": report.rooted,
        "tip_labels": report.tip_labels,
        "tree_is_ultrametric": report.tree_is_ultrametric,
        "root_age": report.root_age,
        "zero_branch_length_count": report.zero_branch_length_count,
        "minimum_tip_depth": report.minimum_tip_depth,
        "maximum_tip_depth": report.maximum_tip_depth,
        "max_tip_depth_deviation": report.max_tip_depth_deviation,
        "tolerance": report.tolerance,
    }, [
        {
            "node_id": row.node_id,
            "node_kind": row.node_kind,
            "node_label": row.node_label or "",
            "descendant_taxa": "|".join(row.descendant_taxa),
            "node_depth": row.node_depth,
            "branching_time": row.branching_time,
        }
        for row in report.rows
    ]


def _build_bijux_diversification_gamma_rows(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_diversification_gamma_statistic(input_fixture)
    return {
        "tip_count": report.tip_count,
        "rooted": report.rooted,
        "ultrametric": report.ultrametric,
        "bifurcating": report.bifurcating,
        "root_age": report.root_age,
        "branching_time_count": report.branching_time_count,
        "interval_count": report.interval_count,
        "minimum_branching_time": report.minimum_branching_time,
        "maximum_branching_time": report.maximum_branching_time,
        "gamma_statistic": report.gamma_statistic,
    }, [
        {
            "tip_count": report.tip_count,
            "rooted": report.rooted,
            "ultrametric": report.ultrametric,
            "bifurcating": report.bifurcating,
            "root_age": report.root_age,
            "branching_time_count": report.branching_time_count,
            "interval_count": report.interval_count,
            "minimum_branching_time": report.minimum_branching_time,
            "maximum_branching_time": report.maximum_branching_time,
            "gamma_statistic": report.gamma_statistic,
        }
    ]


def _build_bijux_tree_simulation_envelope_rows(
    fixture_id: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    fixture = get_shared_tree_simulation_fixture(fixture_id)
    if fixture.simulation_model == "random-tree":
        _trees, report = simulate_random_trees(
            tree_count=fixture.replicate_count,
            tip_count=fixture.tip_count,
            seed=fixture.seed,
            branch_length_model=fixture.branch_length_model or "uniform",
        )
    elif fixture.simulation_model == "coalescent":
        _trees, report = simulate_coalescent_trees(
            tree_count=fixture.replicate_count,
            tip_count=fixture.tip_count,
            population_size=fixture.population_size or 1.0,
            seed=fixture.seed,
        )
    else:
        raise ValueError(
            f"unsupported governed simulation model {fixture.simulation_model!r}"
        )
    return {
        "simulation_model": report.model,
        "reference_function": fixture.reference_function,
        "tree_count": report.tree_count,
        "tip_count": report.tip_count,
        "seed": report.seed,
        "branch_length_model": report.branch_length_model,
        "population_size": report.population_size,
        "rooted": report.rooted,
        "binary": report.binary,
        "pooled_branch_count": report.pooled_branch_count,
        "envelope_metric_count": len(report.envelope_metrics),
    }, [
        {
            "metric": row.metric,
            "sample_scope": row.sample_scope,
            "observation_count": row.observation_count,
            "mean": row.mean,
            "standard_deviation": row.standard_deviation,
            "minimum": row.minimum,
            "median": row.median,
            "maximum": row.maximum,
        }
        for row in report.envelope_metrics
    ]


def _build_bijux_tree_ultrametric_rows(
    input_fixture: Path,
    *,
    tolerance: float,
    option: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = assess_tree_ultrametricity(
        input_fixture,
        tolerance=tolerance,
        option=option,
    )
    return {
        "tip_count": len(report.tip_labels),
        "rooted": report.rooted,
        "tip_labels": report.tip_labels,
        "ultrametric": report.ultrametric,
        "criterion_name": report.criterion_name,
        "criterion_value": report.criterion_value,
        "tolerance": report.tolerance,
        "option": report.option,
        "minimum_tip_depth": report.minimum_tip_depth,
        "maximum_tip_depth": report.maximum_tip_depth,
        "mean_tip_depth": report.mean_tip_depth,
        "max_tip_depth_deviation": report.max_tip_depth_deviation,
        "root_age": report.root_age,
        "offending_taxa": report.offending_taxa,
    }, [
        {
            "node_id": row.node_id,
            "tip_label": row.tip_label,
            "root_to_tip_depth": row.root_to_tip_depth,
            "deviation_from_mean_depth": row.deviation_from_mean_depth,
            "deviation_from_min_depth": row.deviation_from_min_depth,
            "deviation_from_max_depth": row.deviation_from_max_depth,
            "is_offending_taxon": row.is_offending_taxon,
        }
        for row in report.rows
    ]
