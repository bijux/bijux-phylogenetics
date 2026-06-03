from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from ..registry import PhytoolsParityCase


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_rows_table(path: Path) -> list[dict[str, object]]:
    string_identity_fields = {
        "row_kind",
        "label",
        "source_state",
        "target_state",
        "node",
        "state",
    }
    boolean_fields = {
        "transition_allowed",
    }
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: list[dict[str, object]] = []
        for row in reader:
            parsed: dict[str, object] = {}
            for key, value in row.items():
                if value is None or value == "":
                    parsed[key] = ""
                    continue
                if key in string_identity_fields:
                    parsed[key] = value
                    continue
                if key in boolean_fields:
                    lowered = value.lower()
                    if lowered == "true":
                        parsed[key] = True
                        continue
                    if lowered == "false":
                        parsed[key] = False
                        continue
                try:
                    parsed[key] = int(value)
                    continue
                except ValueError:
                    pass
                try:
                    parsed[key] = float(value)
                    continue
                except ValueError:
                    parsed[key] = value
            rows.append(parsed)
        return rows


def reference_rows_filename(case: PhytoolsParityCase) -> str | None:
    if not case.compare_rows:
        return None
    if case.operation == "discrete-fit-mk":
        return "fitmk-rate-matrix.tsv"
    if case.operation in {
        "discrete-stochastic-map",
        "discrete-stochastic-map-count",
        "discrete-stochastic-map-description",
        "discrete-stochastic-map-density",
        "simulate-discrete-history",
    }:
        return "stochastic-map-summary-rows.tsv"
    if case.operation == "simulate-continuous-brownian":
        return "fastbm-summary-rows.tsv"
    if case.operation == "simulate-continuous-correlated-brownian":
        return "simcorrs-summary-rows.tsv"
    if case.operation == "comparative-pgls-brownian":
        return "pgls-summary-rows.tsv"
    if case.operation == "phylogenetic-residuals":
        return "phyl-resid-summary-rows.tsv"
    if case.operation == "phylogenetic-anova":
        return "phyl-anova-summary-rows.tsv"
    if case.operation == "discrete-ancestral-rerooting":
        return "rerooting-method-node-probabilities.tsv"
    if case.operation == "continuous-ancestral-fast-anc":
        return "fast-anc-node-estimates.tsv"
    if case.operation == "continuous-ancestral-anc-ml":
        return "anc-ml-node-estimates.tsv"
    return None


def _isclose(left: object, right: object, *, tolerance: float) -> bool:
    if isinstance(left, bool) or isinstance(right, bool):
        return left == right
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return math.isclose(
            float(left),
            float(right),
            rel_tol=tolerance,
            abs_tol=tolerance,
        )
    return left == right


def _field_tolerance(case: PhytoolsParityCase, key: str) -> float:
    if case.field_tolerances and key in case.field_tolerances:
        return case.field_tolerances[key]
    return case.tolerance


def _row_field_tolerance(case: PhytoolsParityCase, key: str) -> float:
    if case.row_field_tolerances and key in case.row_field_tolerances:
        return case.row_field_tolerances[key]
    return case.tolerance


def mismatch_reason(
    case: PhytoolsParityCase,
    *,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
) -> str | None:
    if reference_summary is None or bijux_summary is None:
        return "summary_missing"
    if case.operation == "phylogenetic-signal-lambda":
        compare_keys = ("taxon_count", "trait_name", "lambda_value", "log_likelihood")
    elif case.operation == "phylogenetic-signal-k":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "k",
            "p_value",
            "permutation_count",
            "permutation_seed",
            "simulated_k_minimum",
            "simulated_k_mean",
        )
    elif case.operation == "discrete-fit-mk":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "excluded_taxon_count",
            "excluded_taxa",
            "model",
            "state_count",
            "parameter_count",
            "log_likelihood",
            "aic",
            "aicc",
            "overparameterized",
            "baseline_model",
            "preferred_model_by_aic",
        )
    elif case.operation in {
        "discrete-stochastic-map",
        "discrete-stochastic-map-count",
        "discrete-stochastic-map-description",
    }:
        compare_keys = (
            "taxon_count",
            "trait_name",
            "excluded_taxon_count",
            "excluded_taxa",
            "model",
            "state_count",
            "parameter_count",
            "log_likelihood",
            "aic",
            "aicc",
            "overparameterized",
            "baseline_model",
            "preferred_model_by_aic",
            "requested_replicate_count",
            "successful_replicate_count",
            "simulation_failure_count",
            "seed",
            "mean_total_transition_count",
            "lower_95_total_transition_count",
            "upper_95_total_transition_count",
        )
        if case.operation == "discrete-stochastic-map":
            compare_keys = (
                compare_keys[:16]
                + ("conditioned_on_node_estimates",)
                + compare_keys[16:]
            )
        elif case.operation == "discrete-stochastic-map-description":
            compare_keys = compare_keys + ("branch_count",)
    elif case.operation == "discrete-stochastic-map-density":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "excluded_taxon_count",
            "excluded_taxa",
            "model",
            "state_count",
            "parameter_count",
            "log_likelihood",
            "aic",
            "aicc",
            "overparameterized",
            "baseline_model",
            "preferred_model_by_aic",
            "requested_replicate_count",
            "successful_replicate_count",
            "simulation_failure_count",
            "seed",
            "branch_count",
            "focal_state",
            "baseline_state",
            "resolution",
            "total_tree_depth",
        )
    elif case.operation == "simulate-discrete-history":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "branch_count",
            "state_count",
            "requested_replicate_count",
            "successful_replicate_count",
            "fixed_root_state",
            "root_prior_probabilities",
            "seed",
            "mean_total_transition_count",
            "lower_95_total_transition_count",
            "upper_95_total_transition_count",
        )
    elif case.operation == "simulate-continuous-brownian":
        compare_keys = (
            "taxon_count",
            "branch_count",
            "requested_replicate_count",
            "successful_replicate_count",
            "seed",
            "root_state",
            "sigma_squared",
        )
    elif case.operation == "simulate-continuous-correlated-brownian":
        compare_keys = (
            "taxon_count",
            "branch_count",
            "trait_count",
            "requested_replicate_count",
            "successful_replicate_count",
            "seed",
        )
    elif case.operation == "comparative-pgls-brownian":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "formula",
            "analysis_taxon_count",
            "coefficient_count",
            "model_matrix_row_count",
            "model_matrix_column_count",
            "categorical_predictor_count",
            "interaction_term_count",
            "lambda_value",
            "lambda_estimation_mode",
            "log_likelihood",
            "aic",
        )
    elif case.operation == "phylogenetic-residuals":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "predictor_name",
            "method",
            "excluded_taxon_count",
            "excluded_taxa",
        )
        if case.comparative_lambda_value is None or not math.isclose(
            case.comparative_lambda_value, 1.0, abs_tol=1e-12
        ):
            compare_keys = compare_keys + ("lambda_value", "log_likelihood")
    elif case.operation == "phylogenetic-anova":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "group_column",
            "excluded_taxon_count",
            "excluded_taxa",
            "group_count",
            "simulation_count",
            "seed",
            "pairwise_adjustment_method",
            "brownian_sigma_squared",
            "sum_of_squares_between",
            "sum_of_squares_within",
            "mean_square_between",
            "mean_square_within",
            "f_statistic",
            "p_value",
            "low_sample_group_count",
        )
    elif case.operation == "discrete-ancestral-rerooting":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "excluded_taxon_count",
            "excluded_taxa",
            "model",
            "state_count",
            "internal_node_count",
            "root_prior_mode",
            "phytools_rerooting_method_comparable",
        )
    elif case.operation == "continuous-ancestral-fast-anc":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "internal_node_count",
            "excluded_taxon_count",
            "excluded_taxa",
            "tree_is_ultrametric",
        )
    elif case.operation == "continuous-ancestral-anc-ml":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "internal_node_count",
            "excluded_taxon_count",
            "excluded_taxa",
            "tree_is_ultrametric",
            "sigma_squared",
            "log_likelihood",
        )
    else:
        return "unsupported_operation"
    for key in compare_keys:
        if key not in reference_summary or key not in bijux_summary:
            return f"summary_field_missing:{key}"
        if not _isclose(
            reference_summary[key],
            bijux_summary[key],
            tolerance=_field_tolerance(case, key),
        ):
            return f"summary_mismatch:{key}"
    return None


def row_mismatch_reason(
    case: PhytoolsParityCase,
    *,
    reference_rows: list[dict[str, object]] | None,
    bijux_rows: list[dict[str, object]] | None,
) -> str | None:
    if not case.compare_rows:
        return None
    if case.operation not in {
        "discrete-fit-mk",
        "discrete-stochastic-map",
        "discrete-stochastic-map-count",
        "discrete-stochastic-map-description",
        "simulate-discrete-history",
        "simulate-continuous-brownian",
        "simulate-continuous-correlated-brownian",
        "comparative-pgls-brownian",
        "phylogenetic-residuals",
        "phylogenetic-anova",
        "discrete-ancestral-rerooting",
        "continuous-ancestral-fast-anc",
        "continuous-ancestral-anc-ml",
    }:
        return None
    if reference_rows is None or bijux_rows is None:
        return "rows_missing"
    if case.operation == "discrete-fit-mk":
        reference_rows = sorted(
            reference_rows,
            key=lambda row: (
                str(row.get("source_state", "")),
                str(row.get("target_state", "")),
            ),
        )
        bijux_rows = sorted(
            bijux_rows,
            key=lambda row: (
                str(row.get("source_state", "")),
                str(row.get("target_state", "")),
            ),
        )
    elif case.operation == "discrete-ancestral-rerooting":
        reference_rows = sorted(
            reference_rows,
            key=lambda row: (
                str(row.get("node", "")),
                str(row.get("state", "")),
            ),
        )
        bijux_rows = sorted(
            bijux_rows,
            key=lambda row: (
                str(row.get("node", "")),
                str(row.get("state", "")),
            ),
        )
    elif case.operation in {
        "discrete-stochastic-map",
        "discrete-stochastic-map-count",
        "discrete-stochastic-map-description",
        "simulate-discrete-history",
        "simulate-continuous-brownian",
        "phylogenetic-residuals",
        "phylogenetic-anova",
    }:
        reference_rows = sorted(
            reference_rows,
            key=lambda row: (
                str(row.get("row_kind", "")),
                str(row.get("label", "")),
            ),
        )
        bijux_rows = sorted(
            bijux_rows,
            key=lambda row: (
                str(row.get("row_kind", "")),
                str(row.get("label", "")),
            ),
        )
    else:
        reference_rows = sorted(
            reference_rows, key=lambda row: str(row.get("node", ""))
        )
        bijux_rows = sorted(bijux_rows, key=lambda row: str(row.get("node", "")))
    if len(reference_rows) != len(bijux_rows):
        return "row_count_mismatch"
    if case.operation == "discrete-fit-mk":
        compare_keys = (
            "source_state",
            "target_state",
            "transition_allowed",
            "step_distance",
            "rate",
        )
    elif case.operation == "discrete-ancestral-rerooting":
        compare_keys = (
            "node",
            "state",
            "probability",
        )
    elif case.operation in {
        "discrete-stochastic-map",
        "discrete-stochastic-map-count",
        "discrete-stochastic-map-description",
        "simulate-discrete-history",
    }:
        compare_keys = (
            "row_kind",
            "label",
            "mean_value",
            "lower_95_interval",
            "upper_95_interval",
            "presence_fraction",
        )
    elif (
        case.operation == "simulate-continuous-brownian"
        or case.operation == "simulate-continuous-correlated-brownian"
    ):
        compare_keys = (
            "row_kind",
            "label",
            "mean_value",
            "standard_deviation",
            "minimum",
            "median",
            "maximum",
            "covariance",
            "correlation",
        )
    elif case.operation == "comparative-pgls-brownian":
        compare_keys = (
            "row_kind",
            "label",
            "value",
        )
    else:
        compare_keys = (
            ("node", "estimate", "standard_error")
            if case.operation == "continuous-ancestral-fast-anc"
            else (
                "node",
                "estimate",
                "standard_error",
                "lower_95_interval",
                "upper_95_interval",
            )
        )
    for reference_row, bijux_row in zip(reference_rows, bijux_rows, strict=True):
        if case.operation == "phylogenetic-residuals":
            row_kind = str(reference_row.get("row_kind", ""))
            if row_kind != str(bijux_row.get("row_kind", "")):
                return "row_mismatch:row_kind"
            compare_keys = (
                ("row_kind", "label", "value")
                if row_kind == "coefficient_estimate"
                else (
                    "row_kind",
                    "label",
                    "observed_value",
                    "fitted_value",
                    "residual",
                )
            )
        elif case.operation == "phylogenetic-anova":
            row_kind = str(reference_row.get("row_kind", ""))
            if row_kind != str(bijux_row.get("row_kind", "")):
                return "row_mismatch:row_kind"
            compare_keys = (
                (
                    "row_kind",
                    "label",
                    "taxon_count",
                    "taxa",
                    "mean",
                    "variance",
                    "minimum",
                    "maximum",
                )
                if row_kind == "group_summary"
                else (
                    "row_kind",
                    "label",
                    "left_taxon_count",
                    "right_taxon_count",
                    "observed_t_statistic",
                    "uncorrected_p_value",
                    "adjusted_p_value",
                )
            )
        for key in compare_keys:
            if key not in reference_row or key not in bijux_row:
                return f"row_field_missing:{key}"
            if not _isclose(
                reference_row[key],
                bijux_row[key],
                tolerance=_row_field_tolerance(case, key),
            ):
                return f"row_mismatch:{key}"
    return None
