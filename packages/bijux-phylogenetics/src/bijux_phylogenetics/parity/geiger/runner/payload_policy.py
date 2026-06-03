from __future__ import annotations

from bijux_phylogenetics.parity.geiger.registry import GeigerParityCase


def comparison_modes(
    candidate_model_names: tuple[str, ...] | None,
) -> tuple[str, ...] | None:
    if candidate_model_names is None:
        return None
    mode_names = {
        "BM": "brownian",
        "white": "white-noise",
        "lambda": "pagel-lambda",
        "kappa": "pagel-kappa",
        "delta": "pagel-delta",
        "OU": "ornstein-uhlenbeck",
        "EB": "early-burst",
    }
    return tuple(
        mode_names.get(model_name, model_name) for model_name in candidate_model_names
    )


def standard_error_policy() -> str:
    return "fitcontinuous-standard-error-explicitly-excluded-this-round"


def discrete_missing_value_policy() -> str:
    return "prune-overlapping-missing-values"


def missing_value_policy() -> str:
    return "prune-tree-tip-overlap-with-missing-or-nonnumeric-trait-values"


def parameter_bound_policy(case: GeigerParityCase) -> str:
    if case.python_mode in {"brownian", "white-noise"}:
        return "closed-form-without-parameter-bounds"
    return "governed-bounded-grid-search"


def bijux_optimizer_result(
    case: GeigerParityCase,
    report,
) -> dict[str, object]:
    if case.operation == "fit-discrete-mk":
        diagnostics = report.optimizer_diagnostics
        return {
            "optimizer_name": diagnostics.optimizer_name,
            "parameter_count": diagnostics.parameter_count,
            "initial_candidate_count": diagnostics.initial_candidate_count,
            "best_initial_scale": diagnostics.best_initial_scale,
            "converged": diagnostics.converged,
            "iteration_count": diagnostics.iteration_count,
            "function_evaluation_count": diagnostics.function_evaluation_count,
            "simplex_shrink_count": diagnostics.simplex_shrink_count,
            "hit_lower_parameter_bound": diagnostics.hit_lower_parameter_bound,
            "hit_upper_parameter_bound": diagnostics.hit_upper_parameter_bound,
        }
    if report.optimizer_diagnostics is not None:
        diagnostics = report.optimizer_diagnostics
        return {
            "optimizer_name": diagnostics.optimizer_name,
            "parameter_search_strategy": diagnostics.parameter_search_strategy,
            "converged": diagnostics.converged,
            "lower_bound": diagnostics.lower_bound,
            "upper_bound": diagnostics.upper_bound,
            "starting_parameter_policy": diagnostics.starting_parameter_policy,
            "starting_parameter_value": diagnostics.starting_parameter_value,
            "starting_parameter_log_likelihood": (
                diagnostics.starting_parameter_log_likelihood
            ),
            "coarse_grid_point_count": diagnostics.coarse_grid_point_count,
            "fine_grid_point_count": diagnostics.fine_grid_point_count,
            "function_evaluation_count": diagnostics.function_evaluation_count,
            "coarse_best_parameter": diagnostics.coarse_best_parameter,
            "coarse_best_log_likelihood": diagnostics.coarse_best_log_likelihood,
            "fine_search_start": diagnostics.fine_search_start,
            "fine_search_stop": diagnostics.fine_search_stop,
            "hit_lower_boundary": diagnostics.hit_lower_boundary,
            "hit_upper_boundary": diagnostics.hit_upper_boundary,
            "profile_rows": (
                None
                if report.optimizer_profile_rows is None
                else [
                    {
                        "parameter_name": report.parameter_name,
                        "parameter_value": row.parameter_value,
                        "log_likelihood": row.log_likelihood,
                    }
                    for row in report.optimizer_profile_rows
                ]
            ),
        }
    if case.python_mode in {"brownian", "white-noise"}:
        return {
            "optimizer_name": "closed-form-profile-solution",
            "parameter_search": "none",
            "converged": True,
            "parameter_count": 2,
        }
    if case.python_mode == "ornstein-uhlenbeck":
        return {
            "optimizer_name": "governed-two-stage-grid-search",
            "parameter_search": "bounded-grid-search",
            "converged": True,
            "parameter_count": 3,
            "coarse_grid_point_count": 81,
            "fine_grid_point_count": 81,
        }
    return {
        "optimizer_name": "governed-two-stage-grid-search",
        "parameter_search": "bounded-grid-search",
        "converged": True,
        "parameter_count": 3,
        "coarse_grid_point_count": 81,
        "fine_grid_point_count": 81,
    }
