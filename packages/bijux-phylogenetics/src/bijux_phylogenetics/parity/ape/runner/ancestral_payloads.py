from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    load_continuous_dataset,
    load_discrete_dataset,
    node_signature,
)
from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete import (
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.phylo.topology.node_identity import build_ape_internal_node_map

from .normalization import _coerce_table_cell


def _build_bijux_continuous_ancestral_rows(
    input_fixture: Path,
    *,
    trait_table_path: Path,
    trait_name: str,
    trait_taxon_column: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    dataset = load_continuous_dataset(
        input_fixture,
        trait_table_path,
        trait=trait_name,
        taxon_column=trait_taxon_column,
    )
    report = reconstruct_continuous_ancestral_states(
        input_fixture,
        trait_table_path,
        trait=trait_name,
        taxon_column=trait_taxon_column,
        model="brownian",
    )
    internal_node_map = {
        node_signature(node): node_id
        for node_id, node in build_ape_internal_node_map(dataset.tree).items()
    }
    rows = sorted(
        [
            {
                "node_id": internal_node_map[estimate.node],
                "node": estimate.node,
                "estimate": estimate.estimate,
                "standard_error": estimate.standard_error,
                "lower_95_interval": estimate.lower_95_interval,
                "upper_95_interval": estimate.upper_95_interval,
            }
            for estimate in report.estimates
            if not estimate.is_tip
        ],
        key=lambda row: int(row["node_id"]),
    )
    diagnostics = report.brownian_fit_diagnostics
    return {
        "trait": report.trait,
        "taxon_count": report.taxon_count,
        "excluded_taxon_count": len(report.dropped_missing_taxa)
        + len(report.dropped_non_numeric_taxa),
        "dropped_missing_taxa": report.dropped_missing_taxa,
        "dropped_non_numeric_taxa": report.dropped_non_numeric_taxa,
        "internal_node_count": len(rows),
        "method": "pic",
        "tree_is_ultrametric": (
            None if diagnostics is None else diagnostics.tree_is_ultrametric
        ),
        "minimum_root_to_tip_depth": (
            None if diagnostics is None else diagnostics.minimum_root_to_tip_depth
        ),
        "maximum_root_to_tip_depth": (
            None if diagnostics is None else diagnostics.maximum_root_to_tip_depth
        ),
    }, rows


def _build_bijux_discrete_ancestral_rows(
    input_fixture: Path,
    *,
    trait_table_path: Path,
    trait_name: str,
    trait_taxon_column: str,
    ancestral_model: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    dataset = load_discrete_dataset(
        input_fixture,
        trait_table_path,
        trait=trait_name,
        taxon_column=trait_taxon_column,
    )
    report = reconstruct_discrete_ancestral_states(
        input_fixture,
        trait_table_path,
        trait=trait_name,
        taxon_column=trait_taxon_column,
        model=ancestral_model,
    )
    internal_node_map = {
        node_signature(node): node_id
        for node_id, node in build_ape_internal_node_map(dataset.tree).items()
    }
    rows = sorted(
        [
            {
                "node_id": internal_node_map[estimate.node],
                "node": estimate.node,
                "state": _coerce_table_cell(state),
                "posterior_probability": probability,
                "most_likely_state": _coerce_table_cell(estimate.most_likely_state),
                "max_posterior_probability": estimate.confidence,
            }
            for estimate in report.estimates
            if not estimate.is_tip
            for state, probability in sorted(estimate.state_probabilities.items())
        ],
        key=lambda row: (int(row["node_id"]), str(row["state"])),
    )
    transition_rows = [
        {
            "source_state": row.source_state,
            "target_state": row.target_state,
            "transition_allowed": row.transition_allowed,
            "step_distance": row.step_distance,
            "rate": row.rate,
        }
        for row in report.transition_rate_rows
    ]
    return {
        "trait": report.trait,
        "taxon_count": report.taxon_count,
        "excluded_taxon_count": len(report.dropped_missing_taxa),
        "dropped_missing_taxa": report.dropped_missing_taxa,
        "internal_node_count": len(
            [estimate for estimate in report.estimates if not estimate.is_tip]
        ),
        "model": report.model,
        "state_count": len(report.observed_states),
        "state_labels": report.observed_states,
        "log_likelihood": report.log_likelihood,
        "parameter_count": report.parameter_count,
        "aic": report.aic,
        "overparameterized": report.overparameterized,
        "baseline_model": (
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.baseline_model
        ),
        "baseline_delta_aic": (
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.delta_aic
        ),
        "preferred_model_by_aic": (
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.preferred_model_by_aic
        ),
        "transition_rate_rows": transition_rows,
    }, rows
