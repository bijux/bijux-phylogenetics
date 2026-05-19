from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.signal import (
    compute_blombergs_k,
    compute_phylogenetic_signal_test,
    estimate_pagels_lambda,
)

from ..registry import PhytoolsParityCase


def build_signal_case_payload(
    case: PhytoolsParityCase,
    *,
    tree_path: Path,
    traits_path: Path | None,
) -> tuple[dict[str, object], list[dict[str, object]] | None] | None:
    if case.operation == "phylogenetic-signal-lambda":
        report = estimate_pagels_lambda(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
        )
        return (
            {
                "taxon_count": report.taxon_count,
                "trait_name": report.trait,
                "lambda_value": report.lambda_value,
                "log_likelihood": report.log_likelihood,
                "null_log_likelihood": report.null_log_likelihood,
                "brownian_log_likelihood": report.brownian_log_likelihood,
                "tree_is_ultrametric": report.input_audit.tree_is_ultrametric,
                "pruned_missing_value_taxa": list(
                    report.input_audit.pruned_missing_value_taxa
                ),
                "warning_count": len(report.input_audit.warnings),
            },
            None,
        )
    if case.operation == "phylogenetic-signal-k":
        signal_test = compute_phylogenetic_signal_test(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
            permutations=case.permutation_count or 199,
            seed=case.permutation_seed or 1,
        )
        report = compute_blombergs_k(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
        )
        return (
            {
                "taxon_count": report.taxon_count,
                "trait_name": report.trait,
                "k": report.k,
                "p_value": signal_test.p_value,
                "permutation_count": signal_test.permutations,
                "permutation_seed": signal_test.seed,
                "null_distribution_count": len(signal_test.permutation_rows),
                "simulated_k_minimum": signal_test.null_distribution_minimum,
                "simulated_k_mean": signal_test.null_distribution_mean,
                "simulated_k_maximum": signal_test.null_distribution_maximum,
                "generalized_mean": report.generalized_mean,
                "observed_mean_square": report.observed_mean_square,
                "phylogenetic_mean_square": report.phylogenetic_mean_square,
                "expected_mean_square_ratio": report.expected_mean_square_ratio,
                "tree_is_ultrametric": report.input_audit.tree_is_ultrametric,
                "pruned_missing_value_taxa": list(
                    report.input_audit.pruned_missing_value_taxa
                ),
                "warning_count": len(report.input_audit.warnings),
            },
            None,
        )
    return None
