from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
import tempfile

from bijux_phylogenetics.comparative.continuous import (
    BrownianMotionFitReport,
    ComparativeParameterInterval,
    ComparativeResidualSummary,
    OUIdentifiabilityWarning,
    OUTraitModelReport,
    compare_brownian_and_ou_models,
    fit_brownian_motion_model,
    fit_ornstein_uhlenbeck_model,
)
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    _alpha_grid,
    _estimate_lambda_for_values,
)
from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    load_comparative_dataset,
)
from bijux_phylogenetics.comparative.pgls import inspect_pgls_inputs, run_pgls
from bijux_phylogenetics.comparative.signal import (
    compute_blombergs_k,
    compute_phylogenetic_independent_contrasts,
    estimate_pagels_lambda,
)

@dataclass(slots=True)
class ComparativeReferenceObservation:
    """One deterministic comparative-model validation example."""

    case: str
    model: str
    trait: str
    source: str
    expected_parameters: dict[str, float]
    observed_parameters: dict[str, float]
    tolerance: float
    passed: bool


@dataclass(slots=True)
class ComparativeReferenceValidationReport:
    """Validation of standalone comparative models against durable reference examples."""

    observations: list[ComparativeReferenceObservation]
    all_passed: bool


@dataclass(slots=True)
class ComparativeParameterIntervalAuditRow:
    """Audit row showing whether an external reference estimate is covered by one interval."""

    model: str
    parameter: str
    estimate: float
    lower_95: float
    upper_95: float
    reference_estimate: float
    interval_method: str
    contains_reference: bool
    reaches_search_boundary: bool
    boundary_note: str | None


@dataclass(slots=True)
class ComparativeParameterUncertaintyAudit:
    """Audit of BM/OU interval behavior against trusted reference estimates."""

    rows: list[ComparativeParameterIntervalAuditRow]
    warnings: list[str]
    all_reference_estimates_covered: bool


@dataclass(slots=True)
class ComparativeOUIdentifiabilityCase:
    """One fixture-level OU identifiability audit case."""

    case: str
    alpha: float
    warning_kinds: list[str]


@dataclass(slots=True)
class ComparativeOUIdentifiabilityAudit:
    """Whether the expected OU failure modes are detected on reference fixtures."""

    cases: list[ComparativeOUIdentifiabilityCase]
    detected_warning_kinds: list[str]
    expected_warning_kinds: list[str]
    all_expected_warning_kinds_detected: bool


@dataclass(slots=True)
class ComparativeResidualDiagnosticSurface:
    """Reviewer-facing residual diagnostics for one comparative analysis surface."""

    analysis: str
    residual_variance: float
    max_abs_standardized_residual: float
    phylogenetic_residual_lambda: float | None
    max_leverage: float | None
    outlier_taxa: list[str]
    high_leverage_taxa: list[str]
    warnings: list[str]


@dataclass(slots=True)
class ComparativeSensitivitySummary:
    """Compact summary of leave-one-taxon-out comparative sensitivity."""

    model: str
    influential_taxa: list[str]
    max_abs_delta_log_likelihood: float
    max_abs_delta_primary_parameter: float


@dataclass(slots=True)
class ComparativeMethodMaturityReport:
    """Integrated comparative audit over one user-supplied response trait workflow."""

    tree_path: Path
    traits_path: Path
    trait: str
    selected_model: str
    reference_validation_passed: bool
    residual_diagnostics: list[ComparativeResidualDiagnosticSurface]
    sensitivity: ComparativeSensitivitySummary
    warnings: list[str]


@dataclass(slots=True)
class LeaveOneTaxonOutRow:
    """Influence of dropping one taxon from a continuous-trait fit."""

    dropped_taxon: str
    taxon_count: int
    primary_parameter: float
    delta_primary_parameter: float
    log_likelihood: float
    delta_log_likelihood: float


@dataclass(slots=True)
class ComparativeSensitivityReport:
    """Leave-one-taxon-out sensitivity over one comparative trait model."""

    tree_path: Path
    traits_path: Path
    trait: str
    model: str
    baseline_primary_parameter: float
    baseline_log_likelihood: float
    rows: list[LeaveOneTaxonOutRow]
    most_influential_taxa: list[str]




def validate_comparative_reference_examples() -> ComparativeReferenceValidationReport:
    """Validate comparative-model outputs against checked-in external reference expectations."""
    fixtures = _load_comparative_reference_fixture()
    root = Path(__file__).resolve().parents[3] / "tests/fixtures"
    tree = root / "trees/example_tree.nwk"
    traits = root / "metadata/example_traits_comparative.tsv"

    contrasts = compute_phylogenetic_independent_contrasts(
        tree, traits, trait="response"
    )
    contrast_lookup = {row.node: row.contrast for row in contrasts.contrasts}
    pgls = run_pgls(
        tree,
        traits,
        response="response",
        predictors=["predictor_one"],
        lambda_value=1.0,
    )
    brownian = fit_brownian_motion_model(tree, traits, trait="response")
    ou = fit_ornstein_uhlenbeck_model(tree, traits, trait="response")
    blomberg = compute_blombergs_k(tree, traits, trait="response")
    pagel_lambda = estimate_pagels_lambda(tree, traits, trait="response")

    observed_by_case = {
        "brownian-example-tree": {
            "root_state": brownian.root_state,
            "rate": brownian.rate,
            "log_likelihood": brownian.log_likelihood,
        },
        "ou-example-tree-grid": {
            "alpha": ou.alpha,
            "theta": ou.theta,
            "log_likelihood": ou.log_likelihood,
        },
        "pic-example-tree": {
            node: contrast_lookup[node] for node in ("A|B", "C|D", "A|B|C|D")
        },
        "pgls-example-tree-brownian": {
            "intercept": pgls.coefficients[0].estimate,
            "predictor_one": pgls.coefficients[1].estimate,
            "log_likelihood": pgls.log_likelihood,
            **{
                f"residual_{taxon}": residual
                for taxon, residual in zip(pgls.taxa, pgls.residuals, strict=True)
            },
        },
        "blomberg-k-example-tree": {
            "k": blomberg.k,
        },
        "pagel-lambda-example-tree": {
            "lambda_value": pagel_lambda.lambda_value,
            "log_likelihood": pagel_lambda.log_likelihood,
        },
    }

    observations: list[ComparativeReferenceObservation] = []
    for example in fixtures["observations"]:
        observed = observed_by_case[example["case"]]
        tolerance = float(example["tolerance"])
        passed = all(
            math.isclose(
                observed[name], expected_value, rel_tol=tolerance, abs_tol=tolerance
            )
            for name, expected_value in example["expected_parameters"].items()
        )
        observations.append(
            ComparativeReferenceObservation(
                case=example["case"],
                model=example["model"],
                trait=example["trait"],
                source=example["source"],
                expected_parameters={
                    name: float(value)
                    for name, value in example["expected_parameters"].items()
                },
                observed_parameters={
                    name: float(observed[name])
                    for name in example["expected_parameters"]
                },
                tolerance=tolerance,
                passed=passed,
            )
        )
    return ComparativeReferenceValidationReport(
        observations=observations,
        all_passed=all(observation.passed for observation in observations),
    )


def audit_comparative_parameter_uncertainty() -> ComparativeParameterUncertaintyAudit:
    """Audit BM/OU parameter intervals against external reference estimates."""
    reference = {
        observation.case: observation
        for observation in validate_comparative_reference_examples().observations
    }
    root = Path(__file__).resolve().parents[3] / "tests/fixtures"
    tree = root / "trees/example_tree.nwk"
    traits = root / "metadata/example_traits_comparative.tsv"
    brownian = fit_brownian_motion_model(tree, traits, trait="response")
    ou = fit_ornstein_uhlenbeck_model(tree, traits, trait="response")

    rows: list[ComparativeParameterIntervalAuditRow] = []
    for interval in brownian.confidence_intervals:
        reference_estimate = reference["brownian-example-tree"].expected_parameters[
            interval.name
        ]
        rows.append(
            ComparativeParameterIntervalAuditRow(
                model="brownian",
                parameter=interval.name,
                estimate=interval.estimate,
                lower_95=interval.lower_95,
                upper_95=interval.upper_95,
                reference_estimate=reference_estimate,
                interval_method=interval.method,
                contains_reference=interval.lower_95
                <= reference_estimate
                <= interval.upper_95,
                reaches_search_boundary=False,
                boundary_note=None,
            )
        )
    for interval in ou.confidence_intervals:
        if interval.name not in {"alpha", "theta"}:
            continue
        reference_estimate = reference["ou-example-tree-grid"].expected_parameters[
            interval.name
        ]
        reaches_boundary = interval.name == "alpha" and math.isclose(
            interval.upper_95, interval.estimate, abs_tol=1e-9
        )
        boundary_note = (
            "supported alpha interval reaches the upper grid boundary and should be interpreted with the boundary warning"
            if reaches_boundary
            else None
        )
        rows.append(
            ComparativeParameterIntervalAuditRow(
                model="ou",
                parameter=interval.name,
                estimate=interval.estimate,
                lower_95=interval.lower_95,
                upper_95=interval.upper_95,
                reference_estimate=reference_estimate,
                interval_method=interval.method,
                contains_reference=interval.lower_95
                <= reference_estimate
                <= interval.upper_95,
                reaches_search_boundary=reaches_boundary,
                boundary_note=boundary_note,
            )
        )
    warnings = [row.boundary_note for row in rows if row.boundary_note]
    return ComparativeParameterUncertaintyAudit(
        rows=rows,
        warnings=warnings,
        all_reference_estimates_covered=all(row.contains_reference for row in rows),
    )


def audit_ou_identifiability_reference_examples() -> ComparativeOUIdentifiabilityAudit:
    """Verify that all expected OU warning modes are triggered on built-in fixtures."""
    root = Path(__file__).resolve().parents[3] / "tests/fixtures"
    cases = [
        (
            "example-tree-small-n",
            root / "trees/example_tree.nwk",
            root / "metadata/example_traits_comparative.tsv",
            "response",
        ),
        (
            "example-tree-weak-pull",
            root / "trees/example_tree_six_taxa.nwk",
            root / "metadata/example_traits_comparative_multiple.tsv",
            "response_growth",
        ),
    ]
    observations: list[ComparativeOUIdentifiabilityCase] = []
    detected: set[str] = set()
    for case, tree, traits, trait in cases:
        report = fit_ornstein_uhlenbeck_model(tree, traits, trait=trait)
        warning_kinds = [warning.kind for warning in report.identifiability_warnings]
        detected.update(warning_kinds)
        observations.append(
            ComparativeOUIdentifiabilityCase(
                case=case,
                alpha=report.alpha,
                warning_kinds=warning_kinds,
            )
        )
    expected = [
        "small_sample_size",
        "boundary_alpha",
        "flat_likelihood",
        "weak_pull_to_optimum",
    ]
    return ComparativeOUIdentifiabilityAudit(
        cases=observations,
        detected_warning_kinds=sorted(detected),
        expected_warning_kinds=expected,
        all_expected_warning_kinds_detected=all(kind in detected for kind in expected),
    )


def assess_comparative_method_maturity(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativeMethodMaturityReport:
    """Summarize residual and sensitivity trust signals for one comparative workflow."""
    pgls_inputs = inspect_pgls_inputs(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
    )
    brownian = fit_brownian_motion_model(
        tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column
    )
    ou = fit_ornstein_uhlenbeck_model(
        tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column
    )
    model_comparison = compare_brownian_and_ou_models(
        tree_path,
        traits_path,
        trait=pgls_inputs.response,
        taxon_column=taxon_column,
    )
    pgls = run_pgls(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    pgls_dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=pgls_inputs.response,
        taxon_column=taxon_column,
        minimum_taxa=max(3, len(pgls.taxa)),
        require_rooted=True,
        require_binary=False,
    )
    residual_lambda = _estimate_lambda_for_values(
        ComparativeDataset(
            tree_path=pgls_dataset.tree_path,
            traits_path=pgls_dataset.traits_path,
            tree=pgls_dataset.tree,
            taxon_column=pgls_dataset.taxon_column,
            trait=pgls_dataset.trait,
            taxa=pgls.taxa,
            trait_values=[0.0] * len(pgls.taxa),
            covariance_matrix=_subset_covariance(
                pgls_dataset.covariance_matrix, pgls_dataset.taxa, pgls.taxa
            ),
            readiness=pgls_dataset.readiness,
        ),
        pgls.residuals,
    )
    pgls_leverage_cutoff = (2.0 * len(pgls.encoded_columns)) / len(pgls.taxa)
    pgls_high_leverage = [
        row.taxon
        for row in pgls.diagnostics.leverage_rows
        if row.leverage >= pgls_leverage_cutoff
    ]
    residual_surfaces = [
        ComparativeResidualDiagnosticSurface(
            analysis="brownian",
            residual_variance=brownian.residual_diagnostics.residual_variance,
            max_abs_standardized_residual=brownian.residual_diagnostics.max_abs_standardized_residual,
            phylogenetic_residual_lambda=brownian.residual_diagnostics.phylogenetic_residual_lambda,
            max_leverage=None,
            outlier_taxa=[
                row.taxon for row in brownian.residual_diagnostics.outlier_taxa
            ],
            high_leverage_taxa=[],
            warnings=list(brownian.residual_diagnostics.warnings),
        ),
        ComparativeResidualDiagnosticSurface(
            analysis="ou",
            residual_variance=ou.residual_diagnostics.residual_variance,
            max_abs_standardized_residual=ou.residual_diagnostics.max_abs_standardized_residual,
            phylogenetic_residual_lambda=ou.residual_diagnostics.phylogenetic_residual_lambda,
            max_leverage=None,
            outlier_taxa=[row.taxon for row in ou.residual_diagnostics.outlier_taxa],
            high_leverage_taxa=[],
            warnings=[
                *ou.residual_diagnostics.warnings,
                *[warning.message for warning in ou.identifiability_warnings],
            ],
        ),
        ComparativeResidualDiagnosticSurface(
            analysis="pgls",
            residual_variance=pgls.residual_variance,
            max_abs_standardized_residual=max(
                abs(row.standardized_residual) for row in pgls.diagnostics.leverage_rows
            ),
            phylogenetic_residual_lambda=residual_lambda,
            max_leverage=max(row.leverage for row in pgls.diagnostics.leverage_rows),
            outlier_taxa=[row.taxon for row in pgls.diagnostics.outlier_taxa],
            high_leverage_taxa=pgls_high_leverage,
            warnings=_build_pgls_residual_warnings(
                residual_lambda,
                outlier_taxa=[row.taxon for row in pgls.diagnostics.outlier_taxa],
                high_leverage_taxa=pgls_high_leverage,
            ),
        ),
    ]
    sensitivity = run_comparative_sensitivity_analysis(
        tree_path,
        traits_path,
        trait=pgls_inputs.response,
        model=model_comparison.better_model,
        taxon_column=taxon_column,
    )
    sensitivity_summary = ComparativeSensitivitySummary(
        model=sensitivity.model,
        influential_taxa=sensitivity.most_influential_taxa,
        max_abs_delta_log_likelihood=max(
            abs(row.delta_log_likelihood) for row in sensitivity.rows
        ),
        max_abs_delta_primary_parameter=max(
            abs(row.delta_primary_parameter) for row in sensitivity.rows
        ),
    )
    warnings = sorted(
        {
            *brownian.residual_diagnostics.warnings,
            *ou.residual_diagnostics.warnings,
            *[warning.message for warning in ou.identifiability_warnings],
            *pgls_inputs.warnings,
            *[warning for surface in residual_surfaces for warning in surface.warnings],
        }
    )
    return ComparativeMethodMaturityReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=pgls_inputs.response,
        selected_model=model_comparison.better_model,
        reference_validation_passed=validate_comparative_reference_examples().all_passed,
        residual_diagnostics=residual_surfaces,
        sensitivity=sensitivity_summary,
        warnings=warnings,
    )


def run_comparative_sensitivity_analysis(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    model: str,
    taxon_column: str | None = None,
) -> ComparativeSensitivityReport:
    """Quantify leave-one-taxon-out sensitivity for a standalone BM or OU fit."""
    if model not in {"brownian", "ou"}:
        raise ValueError(f"unsupported comparative sensitivity model: {model}")
    baseline_bm = fit_brownian_motion_model(
        tree_path, traits_path, trait=trait, taxon_column=taxon_column
    )
    baseline_ou = (
        None
        if model == "brownian"
        else fit_ornstein_uhlenbeck_model(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
        )
    )
    baseline_primary = baseline_bm.rate if model == "brownian" else baseline_ou.alpha
    baseline_log_likelihood = (
        baseline_bm.log_likelihood
        if model == "brownian"
        else baseline_ou.log_likelihood
    )
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=4,
        require_rooted=True,
        require_binary=False,
    )
    rows: list[LeaveOneTaxonOutRow] = []
    root = Path(tree_path)
    table = Path(traits_path)
    for dropped_taxon in dataset.taxa:
        reduced_tree, _ = _write_reduced_comparative_inputs(
            root, dataset.taxa, dropped_taxon
        )
        reduced_traits = _write_reduced_trait_table(
            table, dataset.taxa, dropped_taxon, taxon_column=dataset.taxon_column
        )
        try:
            if model == "brownian":
                reduced = fit_brownian_motion_model(
                    reduced_tree,
                    reduced_traits,
                    trait=trait,
                    taxon_column=dataset.taxon_column,
                )
                primary = reduced.rate
            else:
                reduced = fit_ornstein_uhlenbeck_model(
                    reduced_tree,
                    reduced_traits,
                    trait=trait,
                    taxon_column=dataset.taxon_column,
                )
                primary = reduced.alpha
            rows.append(
                LeaveOneTaxonOutRow(
                    dropped_taxon=dropped_taxon,
                    taxon_count=reduced.taxon_count,
                    primary_parameter=primary,
                    delta_primary_parameter=primary - baseline_primary,
                    log_likelihood=reduced.log_likelihood,
                    delta_log_likelihood=reduced.log_likelihood
                    - baseline_log_likelihood,
                )
            )
        finally:
            reduced_tree.unlink(missing_ok=True)
            reduced_traits.unlink(missing_ok=True)
    most_influential_taxa = [
        row.dropped_taxon
        for row in sorted(
            rows,
            key=lambda row: (
                abs(row.delta_log_likelihood),
                abs(row.delta_primary_parameter),
                row.dropped_taxon,
            ),
            reverse=True,
        )[:3]
    ]
    return ComparativeSensitivityReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        model=model,
        baseline_primary_parameter=baseline_primary,
        baseline_log_likelihood=baseline_log_likelihood,
        rows=sorted(rows, key=lambda row: row.dropped_taxon),
        most_influential_taxa=most_influential_taxa,
    )




def _load_comparative_reference_fixture() -> dict[str, object]:
    fixture_path = (
        Path(__file__).resolve().parents[3]
        / "tests/fixtures/expected/comparative_reference_validation.json"
    )
    return json.loads(fixture_path.read_text(encoding="utf-8"))


def _subset_covariance(
    covariance_matrix: list[list[float]],
    source_taxa: list[str],
    target_taxa: list[str],
) -> list[list[float]]:
    positions = {taxon: index for index, taxon in enumerate(source_taxa)}
    return [
        [
            covariance_matrix[positions[left_taxon]][positions[right_taxon]]
            for right_taxon in target_taxa
        ]
        for left_taxon in target_taxa
    ]


def _build_pgls_residual_warnings(
    residual_lambda: float,
    *,
    outlier_taxa: list[str],
    high_leverage_taxa: list[str],
) -> list[str]:
    warnings: list[str] = []
    if residual_lambda > 0.5:
        warnings.append("PGLS residuals retain moderate phylogenetic structure")
    if outlier_taxa:
        warnings.append("PGLS residual diagnostics identify one or more outlier taxa")
    if high_leverage_taxa:
        warnings.append(
            "PGLS leverage diagnostics identify one or more high-leverage taxa"
        )
    return warnings


def _write_reduced_comparative_inputs(
    tree_path: Path, taxa: list[str], dropped_taxon: str
) -> tuple[Path, list[str]]:
    from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
    from bijux_phylogenetics.io.newick import write_newick

    kept_taxa = [taxon for taxon in taxa if taxon != dropped_taxon]
    reduced_tree, _ = prune_tree_to_requested_taxa(tree_path, kept_taxa)
    out_path = Path(
        tempfile.mkstemp(prefix=f"bijux-comparative-{dropped_taxon}-", suffix=".nwk")[1]
    )
    write_newick(out_path, reduced_tree)
    return out_path, kept_taxa


def _write_reduced_trait_table(
    traits_path: Path,
    taxa: list[str],
    dropped_taxon: str,
    *,
    taxon_column: str,
) -> Path:
    from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows

    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    kept_taxa = {taxon for taxon in taxa if taxon != dropped_taxon}
    rows = [row for row in table.rows if row[table.taxon_column] in kept_taxa]
    out_path = Path(
        tempfile.mkstemp(
            prefix=f"bijux-comparative-{dropped_taxon}-", suffix=traits_path.suffix
        )[1]
    )
    write_taxon_rows(out_path, columns=table.columns, rows=rows)
    return out_path
