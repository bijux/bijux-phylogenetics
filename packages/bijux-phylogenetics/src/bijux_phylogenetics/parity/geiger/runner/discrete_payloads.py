from __future__ import annotations

import math

from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
from bijux_phylogenetics.parity.geiger.registry import GeigerParityCase

from .payload_policy import bijux_optimizer_result, discrete_missing_value_policy


def discrete_rate_rows(report) -> list[dict[str, object]]:
    return [
        {
            "source_state": row.source_state,
            "target_state": row.target_state,
            "transition_allowed": row.transition_allowed,
            "step_distance": row.step_distance,
            "rate": row.rate,
        }
        for row in report.transition_rate_rows
    ]


def build_bijux_discrete_case_payload(
    case: GeigerParityCase,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    tree_path, traits_path = case.input_fixtures
    report = fit_discrete_mk_model(
        tree_path,
        traits_path,
        trait=case.trait_name,
        taxon_column=case.taxon_column,
        model=case.python_mode,
        transform=case.discrete_transform_name,
        lambda_bounds=(0.0, 1.0) if case.lambda_bounds is None else case.lambda_bounds,
        kappa_bounds=(0.0, 1.0) if case.kappa_bounds is None else case.kappa_bounds,
        delta_bounds=(math.exp(-5.0), 3.0)
        if case.delta_bounds is None
        else case.delta_bounds,
        early_burst_bounds=(-10.0, 10.0)
        if case.early_burst_bounds is None
        else case.early_burst_bounds,
    )
    input_audit = report.input_audit
    missing_value_taxa = list(input_audit.pruned_missing_value_taxa)
    excluded_taxa = sorted(
        set(input_audit.missing_from_traits) | set(missing_value_taxa)
    )
    diagnostics = report.optimizer_diagnostics
    transform_fit = report.transform_fit
    summary = {
        "taxon_count": report.taxon_count,
        "trait_name": report.trait,
        "model_name": case.model_name,
        "transform_name": (
            None
            if transform_fit is None
            else {
                "lambda": "pagel-lambda",
                "kappa": "pagel-kappa",
                "delta": "pagel-delta",
                "early-burst": "early-burst",
            }.get(transform_fit.transform_name, transform_fit.transform_name)
        ),
        "observed_state_count": len(input_audit.observed_states),
        "state_order": list(report.state_order),
        "excluded_taxon_count": len(excluded_taxa),
        "excluded_taxa": excluded_taxa,
        "missing_value_taxa": missing_value_taxa,
        "missing_from_traits": list(input_audit.missing_from_traits),
        "extra_trait_taxa": list(input_audit.extra_trait_taxa),
        "missing_value_policy": discrete_missing_value_policy(),
        "log_likelihood": report.log_likelihood,
        "parameter_count": report.parameter_count,
        "aic": report.aic,
        "aicc": report.aicc,
        "parameter_name": None
        if transform_fit is None
        else transform_fit.parameter_name,
        "parameter_value": (
            None if transform_fit is None else transform_fit.parameter_value
        ),
        "hit_lower_parameter_boundary": (
            diagnostics.hit_lower_parameter_bound
            if transform_fit is None
            else transform_fit.hit_lower_parameter_boundary
        ),
        "hit_upper_parameter_boundary": (
            diagnostics.hit_upper_parameter_bound
            if transform_fit is None
            else transform_fit.hit_upper_parameter_boundary
        ),
        "optimizer_settings": case.optimizer_settings,
        "optimizer_result": {
            **bijux_optimizer_result(case, report),
            "profile_rows": (
                None
                if transform_fit is None
                else [
                    {
                        "parameter_name": transform_fit.parameter_name,
                        "parameter_value": row.transform_parameter_value,
                        "log_likelihood": row.log_likelihood,
                    }
                    for row in transform_fit.profile_rows
                ]
            ),
        },
    }
    return summary, discrete_rate_rows(report)
