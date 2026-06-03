from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Never

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    log_determinant,
    stable_covariance,
    student_t_quantile,
    student_t_two_sided_p_value,
)
from bijux_phylogenetics.comparative.common import (
    build_brownian_covariance_matrix,
    build_ou_covariance_matrix,
    load_comparative_dataset,
    tip_root_depths,
)
from bijux_phylogenetics.comparative.continuous.model_fitting import _alpha_grid
from bijux_phylogenetics.comparative.pgls.design import inspect_pgls_inputs
from bijux_phylogenetics.comparative.pgls.fitting import (
    _build_pgls_diagnostics,
    _fit_gls,
    _gls_log_likelihood,
    _quadratic_form,
)
from bijux_phylogenetics.comparative.pgls.models import (
    ComparativeFormulaSpecification,
    PGLSCoefficient,
    PGLSDiagnosticsReport,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.runtime.errors import ComparativeMethodError


@dataclass(slots=True)
class OUCovarianceRow:
    """One pairwise covariance observation for an OU-style regression surface."""

    left_taxon: str
    right_taxon: str
    is_diagonal: bool
    covariance_value: float
    shared_path_length: float
    left_root_depth: float
    right_root_depth: float


@dataclass(slots=True)
class OUPGLSAlphaProfileRow:
    """One likelihood-profile row across candidate OU alpha values."""

    alpha: float
    log_likelihood: float
    delta_log_likelihood: float
    within_95_confidence_interval: bool


@dataclass(slots=True)
class OUCovarianceModelFit:
    """Reviewer-facing comparative regression fit under one OU covariance."""

    response: str
    formula: ComparativeFormulaSpecification
    predictors: list[str]
    interaction_terms: list[str]
    encoded_columns: list[str]
    taxon_count: int
    alpha: float
    alpha_estimation_mode: str
    log_likelihood: float
    aic: float
    residual_variance: float
    r_squared: float
    coefficients: list[PGLSCoefficient]
    fitted_values: list[float]
    residuals: list[float]
    taxa: list[str]
    diagnostics: PGLSDiagnosticsReport


@dataclass(slots=True)
class OUCovariancePGLSReport:
    """Owned OU-covariance comparative regression workflow."""

    tree_path: Path
    traits_path: Path
    response: str
    formula: ComparativeFormulaSpecification
    taxon_count: int
    tree_is_ultrametric: bool
    minimum_root_to_tip_depth: float
    maximum_root_to_tip_depth: float
    minimum_branch_length: float
    maximum_branch_length: float
    raw_log_determinant: float
    positive_definite_before_stabilization: bool
    alpha: float
    alpha_estimation_mode: str
    lower_95_confidence_interval: float | None
    upper_95_confidence_interval: float | None
    alpha_profile_rows: list[OUPGLSAlphaProfileRow]
    rows: list[OUCovarianceRow]
    model: OUCovarianceModelFit


def summarize_ou_covariance_pgls(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    alpha: float | str = "estimate",
) -> OUCovariancePGLSReport:
    """Fit one comparative regression under stationary-root OU covariance."""
    try:
        input_report = inspect_pgls_inputs(
            tree_path,
            traits_path,
            response=response,
            predictors=predictors,
            formula=formula,
            taxon_column=taxon_column,
        )
        if not input_report.ready:
            raise ComparativeMethodError("; ".join(input_report.blockers))
        dataset = load_comparative_dataset(
            tree_path,
            traits_path,
            trait=input_report.formula_audit.response_column,
            taxon_column=taxon_column,
            minimum_taxa=len(input_report.encoded_columns) + 1,
            require_rooted=True,
            require_binary=False,
        )
        taxa = list(input_report.analysis_taxa)
        design_matrix = [
            [
                row.encoded_values[column]
                for column in input_report.model_matrix.encoded_columns
            ]
            for row in input_report.model_matrix.rows
        ]
        response_values = [row.response_value for row in input_report.model_matrix.rows]
        resolved_alpha, alpha_mode, profile_rows, interval_bounds = _resolve_alpha_fit(
            dataset=dataset,
            taxa=taxa,
            design_matrix=design_matrix,
            response_values=response_values,
            alpha=alpha,
        )
        raw_covariance = build_ou_covariance_matrix(
            dataset.tree, taxa, alpha=resolved_alpha
        )
        root_depths = tip_root_depths(dataset.tree, taxa)
        ultrametric_summary = summarize_ultrametric_tip_depths(
            root_depths,
            tolerance=1e-12,
        )
        shared_paths = build_brownian_covariance_matrix(dataset.tree, taxa)
        minimum_branch_length, maximum_branch_length = _branch_length_range(tree_path)
        raw_log_determinant = _validate_raw_ou_covariance(
            tree_path=tree_path,
            taxa=taxa,
            covariance_matrix=raw_covariance,
            minimum_branch_length=minimum_branch_length,
        )
        stabilized_covariance = stable_covariance(raw_covariance)
        fit = _fit_ou_covariance_model(
            response=input_report.response,
            formula=input_report.formula,
            predictors=list(input_report.formula.predictors),
            interaction_terms=list(input_report.formula.interaction_terms),
            encoded_columns=list(input_report.model_matrix.encoded_columns),
            taxa=taxa,
            response_values=response_values,
            design_matrix=design_matrix,
            covariance=stabilized_covariance,
            alpha=resolved_alpha,
            alpha_mode=alpha_mode,
        )
        return OUCovariancePGLSReport(
            tree_path=tree_path,
            traits_path=traits_path,
            response=input_report.response,
            formula=input_report.formula,
            taxon_count=len(taxa),
            tree_is_ultrametric=ultrametric_summary.ultrametric,
            minimum_root_to_tip_depth=ultrametric_summary.minimum_tip_depth,
            maximum_root_to_tip_depth=ultrametric_summary.maximum_tip_depth,
            minimum_branch_length=minimum_branch_length,
            maximum_branch_length=maximum_branch_length,
            raw_log_determinant=raw_log_determinant,
            positive_definite_before_stabilization=True,
            alpha=resolved_alpha,
            alpha_estimation_mode=alpha_mode,
            lower_95_confidence_interval=interval_bounds[0],
            upper_95_confidence_interval=interval_bounds[1],
            alpha_profile_rows=profile_rows,
            rows=_build_covariance_rows(
                taxa=taxa,
                covariance_matrix=raw_covariance,
                shared_paths=shared_paths,
                root_depths=root_depths,
            ),
            model=fit,
        )
    except ComparativeMethodError as error:
        _reraise_ou_input_error(tree_path, error)


def write_ou_covariance_table(path: Path, report: OUCovariancePGLSReport) -> Path:
    """Write the pairwise OU covariance ledger as TSV or CSV."""
    return write_taxon_rows(
        path,
        columns=[
            "left_taxon",
            "right_taxon",
            "is_diagonal",
            "covariance_value",
            "shared_path_length",
            "left_root_depth",
            "right_root_depth",
            "alpha",
            "alpha_estimation_mode",
            "tree_is_ultrametric",
            "minimum_root_to_tip_depth",
            "maximum_root_to_tip_depth",
            "minimum_branch_length",
            "maximum_branch_length",
            "raw_log_determinant",
            "positive_definite_before_stabilization",
        ],
        rows=[
            {
                "left_taxon": row.left_taxon,
                "right_taxon": row.right_taxon,
                "is_diagonal": "true" if row.is_diagonal else "false",
                "covariance_value": format(row.covariance_value, ".15g"),
                "shared_path_length": format(row.shared_path_length, ".15g"),
                "left_root_depth": format(row.left_root_depth, ".15g"),
                "right_root_depth": format(row.right_root_depth, ".15g"),
                "alpha": format(report.alpha, ".15g"),
                "alpha_estimation_mode": report.alpha_estimation_mode,
                "tree_is_ultrametric": (
                    "true" if report.tree_is_ultrametric else "false"
                ),
                "minimum_root_to_tip_depth": format(
                    report.minimum_root_to_tip_depth, ".15g"
                ),
                "maximum_root_to_tip_depth": format(
                    report.maximum_root_to_tip_depth, ".15g"
                ),
                "minimum_branch_length": format(report.minimum_branch_length, ".15g"),
                "maximum_branch_length": format(report.maximum_branch_length, ".15g"),
                "raw_log_determinant": format(report.raw_log_determinant, ".15g"),
                "positive_definite_before_stabilization": (
                    "true" if report.positive_definite_before_stabilization else "false"
                ),
            }
            for row in report.rows
        ],
    )


def write_ou_alpha_profile_table(path: Path, report: OUCovariancePGLSReport) -> Path:
    """Write the fitted OU alpha likelihood profile as TSV or CSV."""
    return write_taxon_rows(
        path,
        columns=[
            "alpha_estimation_mode",
            "alpha",
            "log_likelihood",
            "delta_log_likelihood",
            "within_95_confidence_interval",
            "profile_lower_95_confidence_interval",
            "profile_upper_95_confidence_interval",
        ],
        rows=[
            {
                "alpha_estimation_mode": report.alpha_estimation_mode,
                "alpha": format(row.alpha, ".15g"),
                "log_likelihood": format(row.log_likelihood, ".15g"),
                "delta_log_likelihood": format(row.delta_log_likelihood, ".15g"),
                "within_95_confidence_interval": (
                    "true" if row.within_95_confidence_interval else "false"
                ),
                "profile_lower_95_confidence_interval": (
                    ""
                    if report.lower_95_confidence_interval is None
                    else format(report.lower_95_confidence_interval, ".15g")
                ),
                "profile_upper_95_confidence_interval": (
                    ""
                    if report.upper_95_confidence_interval is None
                    else format(report.upper_95_confidence_interval, ".15g")
                ),
            }
            for row in report.alpha_profile_rows
        ],
    )


def _resolve_alpha_fit(
    *,
    dataset,
    taxa: list[str],
    design_matrix: list[list[float]],
    response_values: list[float],
    alpha: float | str,
) -> tuple[float, str, list[OUPGLSAlphaProfileRow], tuple[float | None, float | None]]:
    if alpha != "estimate":
        fixed_alpha = float(alpha)
        if fixed_alpha <= 0.0:
            raise ComparativeMethodError("OU alpha must be positive")
        log_likelihood = _ou_gls_log_likelihood(
            dataset=dataset,
            taxa=taxa,
            design_matrix=design_matrix,
            response_values=response_values,
            alpha=fixed_alpha,
        )
        return (
            fixed_alpha,
            "fixed",
            [
                OUPGLSAlphaProfileRow(
                    alpha=fixed_alpha,
                    log_likelihood=log_likelihood,
                    delta_log_likelihood=0.0,
                    within_95_confidence_interval=True,
                )
            ],
            (None, None),
        )
    profile: list[tuple[float, float]] = []
    for candidate in _alpha_grid(dataset):
        log_likelihood = _ou_gls_log_likelihood(
            dataset=dataset,
            taxa=taxa,
            design_matrix=design_matrix,
            response_values=response_values,
            alpha=candidate,
        )
        profile.append((candidate, log_likelihood))
    best_alpha, best_log_likelihood = max(profile, key=lambda item: item[1])
    supported = [
        candidate
        for candidate, log_likelihood in profile
        if log_likelihood >= best_log_likelihood - 1.92
    ]
    return (
        best_alpha,
        "estimated",
        [
            OUPGLSAlphaProfileRow(
                alpha=candidate,
                log_likelihood=log_likelihood,
                delta_log_likelihood=best_log_likelihood - log_likelihood,
                within_95_confidence_interval=candidate in supported,
            )
            for candidate, log_likelihood in profile
        ],
        (min(supported), max(supported)),
    )


def _ou_gls_log_likelihood(
    *,
    dataset,
    taxa: list[str],
    design_matrix: list[list[float]],
    response_values: list[float],
    alpha: float,
) -> float:
    covariance = stable_covariance(
        build_ou_covariance_matrix(dataset.tree, taxa, alpha=alpha)
    )
    inverse_covariance = invert_matrix(covariance)
    _, _, fitted_values = _fit_gls(design_matrix, response_values, inverse_covariance)
    residuals = [
        observed - fitted
        for observed, fitted in zip(response_values, fitted_values, strict=True)
    ]
    return _gls_log_likelihood(
        response_values,
        residuals,
        inverse_covariance,
        covariance,
    )


def _fit_ou_covariance_model(
    *,
    response: str,
    formula: ComparativeFormulaSpecification,
    predictors: list[str],
    interaction_terms: list[str],
    encoded_columns: list[str],
    taxa: list[str],
    response_values: list[float],
    design_matrix: list[list[float]],
    covariance: list[list[float]],
    alpha: float,
    alpha_mode: str,
) -> OUCovarianceModelFit:
    inverse_covariance = invert_matrix(covariance)
    coefficients, covariance_of_betas, fitted_values = _fit_gls(
        design_matrix, response_values, inverse_covariance
    )
    residuals = [
        observed - fitted
        for observed, fitted in zip(response_values, fitted_values, strict=True)
    ]
    degrees_of_freedom = len(response_values) - len(coefficients)
    residual_variance = (
        _quadratic_form(residuals, inverse_covariance) / degrees_of_freedom
    )
    critical_value = student_t_quantile(0.975, degrees_of_freedom)
    coefficient_rows: list[PGLSCoefficient] = []
    for index, name in enumerate(encoded_columns):
        standard_error = math.sqrt(
            max(covariance_of_betas[index][index] * residual_variance, 0.0)
        )
        test_statistic = coefficients[index] / standard_error if standard_error else 0.0
        p_value = student_t_two_sided_p_value(test_statistic, degrees_of_freedom)
        interval_radius = critical_value * standard_error
        coefficient_rows.append(
            PGLSCoefficient(
                name=name,
                estimate=coefficients[index],
                standard_error=standard_error,
                test_statistic=test_statistic,
                p_value=p_value,
                lower_95_confidence_interval=coefficients[index] - interval_radius,
                upper_95_confidence_interval=coefficients[index] + interval_radius,
                degrees_of_freedom=degrees_of_freedom,
                inference_distribution="student-t",
            )
        )
    mean_response = sum(response_values) / len(response_values)
    total_sum_of_squares = sum(
        (value - mean_response) ** 2 for value in response_values
    )
    residual_sum_of_squares = sum(value * value for value in residuals)
    r_squared = (
        1.0 - (residual_sum_of_squares / total_sum_of_squares)
        if total_sum_of_squares
        else 1.0
    )
    log_likelihood = _gls_log_likelihood(
        response_values,
        residuals,
        inverse_covariance,
        covariance,
    )
    parameter_count = len(coefficients) + (1 if alpha_mode == "estimated" else 0)
    return OUCovarianceModelFit(
        response=response,
        formula=formula,
        predictors=predictors,
        interaction_terms=interaction_terms,
        encoded_columns=encoded_columns,
        taxon_count=len(taxa),
        alpha=alpha,
        alpha_estimation_mode=alpha_mode,
        log_likelihood=log_likelihood,
        aic=(2.0 * parameter_count) - (2.0 * log_likelihood),
        residual_variance=residual_variance,
        r_squared=r_squared,
        coefficients=coefficient_rows,
        fitted_values=fitted_values,
        residuals=residuals,
        taxa=taxa,
        diagnostics=_build_pgls_diagnostics(
            taxa,
            response_values,
            fitted_values,
            residuals,
            residual_variance,
            design_matrix,
            inverse_covariance,
        ),
    )


def _reraise_ou_input_error(tree_path: Path, error: ComparativeMethodError) -> Never:
    details = error.details or {}
    failure_reason = details.get("failure_reason")
    evidence = details.get("evidence", {})
    if failure_reason == "comparative_negative_branch_lengths":
        raise ComparativeMethodError(
            "OU covariance is invalid: tree contains negative branch lengths",
            details={
                "failure_reason": "ou_covariance_negative_branch_lengths",
                "scientific_explanation": (
                    "OU covariance is invalid on a tree with negative branch lengths because stationary variance depends on non-negative evolutionary distance."
                ),
                "likely_causes": [
                    "the tree file contains one or more negative branch lengths",
                ],
                "actionable_fixes": [
                    "repair or re-estimate the tree so every non-root branch length is non-negative",
                    "inspect the tree for scaling or export errors that introduced negative lengths",
                ],
                "evidence": {
                    "tree_path": str(tree_path),
                    "minimum_branch_length": evidence.get("minimum_branch_length"),
                },
            },
        ) from error
    if failure_reason == "comparative_branch_lengths_incomplete":
        raise ComparativeMethodError(
            "OU covariance requires complete branch lengths",
            details={
                "failure_reason": "ou_covariance_branch_lengths_incomplete",
                "scientific_explanation": (
                    "OU covariance needs complete numeric branch lengths because the stationary process depends on branch-length-scaled attraction toward the optimum."
                ),
                "likely_causes": [
                    "the tree was exported without complete branch lengths",
                    "one or more branches have blank or missing lengths",
                ],
                "actionable_fixes": [
                    "rerun tree inference or export with branch lengths preserved",
                    "inspect the tree file for missing branch-length fields",
                ],
                "evidence": {"tree_path": str(tree_path)},
            },
        ) from error
    raise error


def _validate_raw_ou_covariance(
    *,
    tree_path: Path,
    taxa: list[str],
    covariance_matrix: list[list[float]],
    minimum_branch_length: float,
) -> float:
    if minimum_branch_length < 0.0:
        raise ComparativeMethodError(
            "OU covariance is invalid: tree contains negative branch lengths",
            details={
                "failure_reason": "ou_covariance_negative_branch_lengths",
                "scientific_explanation": (
                    "OU covariance is invalid on a tree with negative branch lengths because stationary variance depends on non-negative evolutionary distance."
                ),
                "likely_causes": [
                    "the tree file contains one or more negative branch lengths",
                ],
                "actionable_fixes": [
                    "repair or re-estimate the tree so every non-root branch length is non-negative",
                    "inspect the tree for scaling or export errors that introduced negative lengths",
                ],
                "evidence": {
                    "tree_path": str(tree_path),
                    "minimum_branch_length": minimum_branch_length,
                },
            },
        )
    diagonal = [covariance_matrix[index][index] for index in range(len(taxa))]
    non_positive_taxa = [
        taxon for taxon, value in zip(taxa, diagonal, strict=True) if value <= 0.0
    ]
    if non_positive_taxa:
        raise ComparativeMethodError(
            "OU covariance is invalid: non-positive stationary variance for "
            + ", ".join(non_positive_taxa)
        )
    try:
        invert_matrix(covariance_matrix)
        return log_determinant(covariance_matrix)
    except ValueError as error:
        raise ComparativeMethodError(
            f"OU covariance is invalid before stabilization: {tree_path.name}: {error}"
        ) from error


def _branch_length_range(tree_path: Path) -> tuple[float, float]:
    from bijux_phylogenetics.io.trees import load_tree

    tree = load_tree(tree_path)
    branch_lengths = [
        node.branch_length
        for node in tree.iter_nodes()
        if node is not tree.root and node.branch_length is not None
    ]
    if not branch_lengths:
        raise ComparativeMethodError(
            "OU covariance requires complete branch lengths",
            details={
                "failure_reason": "ou_covariance_branch_lengths_incomplete",
                "scientific_explanation": (
                    "OU covariance needs complete numeric branch lengths because the stationary process depends on branch-length-scaled attraction toward the optimum."
                ),
                "likely_causes": [
                    "the tree was exported without complete branch lengths",
                    "one or more branches have blank or missing lengths",
                ],
                "actionable_fixes": [
                    "rerun tree inference or export with branch lengths preserved",
                    "inspect the tree file for missing branch-length fields",
                ],
                "evidence": {"tree_path": str(tree_path)},
            },
        )
    return min(branch_lengths), max(branch_lengths)


def _build_covariance_rows(
    *,
    taxa: list[str],
    covariance_matrix: list[list[float]],
    shared_paths: list[list[float]],
    root_depths: dict[str, float],
) -> list[OUCovarianceRow]:
    rows: list[OUCovarianceRow] = []
    for row_index, left_taxon in enumerate(taxa):
        for column_index, right_taxon in enumerate(taxa):
            rows.append(
                OUCovarianceRow(
                    left_taxon=left_taxon,
                    right_taxon=right_taxon,
                    is_diagonal=row_index == column_index,
                    covariance_value=covariance_matrix[row_index][column_index],
                    shared_path_length=shared_paths[row_index][column_index],
                    left_root_depth=root_depths[left_taxon],
                    right_root_depth=root_depths[right_taxon],
                )
            )
    return rows
