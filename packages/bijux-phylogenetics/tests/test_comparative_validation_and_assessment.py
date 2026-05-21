from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.assessment import (
    assess_comparative_method_maturity,
    run_comparative_sensitivity_analysis,
)
from bijux_phylogenetics.comparative.continuous import (
    compare_brownian_and_ou_models,
    fit_brownian_motion_model,
    fit_ornstein_uhlenbeck_model,
)
from bijux_phylogenetics.comparative.validation import (
    audit_comparative_parameter_uncertainty,
    audit_ou_identifiability_reference_examples,
    validate_comparative_reference_examples,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_fit_brownian_motion_model_reports_root_state_rate_and_intervals() -> None:
    report = fit_brownian_motion_model(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    assert report.taxon_count == 4
    assert math.isclose(report.root_state, 2.8055555543209874)
    assert math.isclose(report.rate, 4.774305191647407)
    assert [interval.name for interval in report.confidence_intervals] == [
        "root_state",
        "rate",
    ]
    assert report.residual_diagnostics.outlier_taxa == []


def test_fit_ornstein_uhlenbeck_model_reports_identifiability_warnings() -> None:
    report = fit_ornstein_uhlenbeck_model(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    assert math.isclose(report.alpha, 33.333333)
    assert math.isclose(report.theta, 2.7503175510416416)
    assert [warning.kind for warning in report.identifiability_warnings] == [
        "small_sample_size",
        "boundary_alpha",
        "flat_likelihood",
    ]
    assert report.convergence_status == "grid-search-converged"


def test_compare_brownian_and_ou_models_prefers_brownian_on_example_dataset() -> None:
    report = compare_brownian_and_ou_models(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    assert report.better_model == "brownian"
    assert [row.model for row in report.rows] == ["brownian", "ou"]
    assert report.rows[0].selected is True


def test_validate_comparative_reference_examples_passes() -> None:
    report = validate_comparative_reference_examples()
    assert report.all_passed is True
    assert len(report.observations) == 6
    assert {row.model for row in report.observations} == {
        "blombergs-k",
        "brownian",
        "independent-contrasts",
        "ou",
        "pagels-lambda",
        "pgls",
    }


def test_audit_comparative_parameter_uncertainty_covers_reference_estimates() -> None:
    audit = audit_comparative_parameter_uncertainty()
    assert audit.all_reference_estimates_covered is True
    alpha_row = next(
        row for row in audit.rows if row.model == "ou" and row.parameter == "alpha"
    )
    assert alpha_row.reaches_search_boundary is True
    assert alpha_row.boundary_note is not None


def test_audit_ou_identifiability_reference_examples_detects_all_expected_modes() -> (
    None
):
    audit = audit_ou_identifiability_reference_examples()
    assert audit.all_expected_warning_kinds_detected is True
    assert audit.detected_warning_kinds == [
        "boundary_alpha",
        "flat_likelihood",
        "small_sample_size",
        "weak_pull_to_optimum",
    ]


def test_run_comparative_sensitivity_analysis_reports_influential_taxa() -> None:
    report = run_comparative_sensitivity_analysis(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        model="brownian",
    )
    assert report.most_influential_taxa == ["A", "D", "B"]
    rows = {row.dropped_taxon: row for row in report.rows}
    assert rows["A"].delta_log_likelihood > 2.0
    assert rows["C"].delta_primary_parameter > 0.8


def test_assess_comparative_method_maturity_reports_residual_surfaces_and_sensitivity() -> (
    None
):
    report = assess_comparative_method_maturity(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value=1.0,
    )
    assert report.reference_validation_passed is True
    assert report.selected_model == "brownian"
    assert [surface.analysis for surface in report.residual_diagnostics] == [
        "brownian",
        "ou",
        "pgls",
    ]
    pgls_surface = next(
        surface for surface in report.residual_diagnostics if surface.analysis == "pgls"
    )
    assert pgls_surface.max_leverage is not None
    assert report.sensitivity.influential_taxa == ["A", "D", "B"]
