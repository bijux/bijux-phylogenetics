from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.parity import run_geiger_parity_cases
from tests.support.fake_geiger_parity import fake_geiger_rscript
from tests.support.geiger_fitcontinuous_brownian_reference import (
    GEIGER_FITCONTINUOUS_BROWNIAN_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitcontinuous_early_burst_reference import (
    GEIGER_FITCONTINUOUS_EARLY_BURST_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitcontinuous_lambda_reference import (
    GEIGER_FITCONTINUOUS_LAMBDA_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitcontinuous_model_comparison_reference import (
    GEIGER_FITCONTINUOUS_MODEL_COMPARISON_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitcontinuous_ou_reference import (
    GEIGER_FITCONTINUOUS_OU_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitcontinuous_white_reference import (
    GEIGER_FITCONTINUOUS_WHITE_REFERENCE_PAYLOADS,
)
from tests.support.geiger_optimizer_triage_reference import (
    geiger_optimizer_triage_reference_payloads,
)

pytestmark = pytest.mark.slow


def test_run_geiger_parity_cases_triages_same_likelihood_different_parameters(
    tmp_path: Path,
) -> None:
    reference_payloads = geiger_optimizer_triage_reference_payloads()
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads={
            "fitcontinuous-lambda-strong-signal-review": reference_payloads[
                "same_likelihood_different_parameters"
            ]
        },
    )

    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-lambda-strong-signal-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    triage_row = report.optimizer_triage_rows[0]
    assert triage_row.mismatch_type == "same_likelihood_different_parameters"
    assert triage_row.same_likelihood_within_tolerance is True
    assert triage_row.same_parameter_surface_within_tolerance is False
    assert triage_row.parameter_surface_comparable is True


def test_run_geiger_parity_cases_triages_different_likelihood_same_parameters(
    tmp_path: Path,
) -> None:
    reference_payloads = geiger_optimizer_triage_reference_payloads()
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads={
            "fitcontinuous-lambda-strong-signal-review": reference_payloads[
                "different_likelihood_same_parameters"
            ]
        },
    )

    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-lambda-strong-signal-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    triage_row = report.optimizer_triage_rows[0]
    assert triage_row.mismatch_type == "different_likelihood_same_parameters"
    assert triage_row.same_likelihood_within_tolerance is False
    assert triage_row.same_parameter_surface_within_tolerance is True
    assert triage_row.parameter_surface_comparable is True


def test_run_geiger_parity_cases_triages_boundary_solution_review(
    tmp_path: Path,
) -> None:
    reference_payloads = geiger_optimizer_triage_reference_payloads()
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads={
            "fitcontinuous-lambda-strong-signal-review": reference_payloads[
                "boundary_solution_review"
            ]
        },
    )

    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-lambda-strong-signal-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    triage_row = report.optimizer_triage_rows[0]
    assert triage_row.mismatch_type == "boundary_solution_review"
    assert triage_row.boundary_solution_detected is True
    assert triage_row.reference_boundary_detected is True


def test_run_geiger_parity_cases_builds_boundary_warning_rows(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=(
            GEIGER_FITCONTINUOUS_LAMBDA_REFERENCE_PAYLOADS
            | GEIGER_FITCONTINUOUS_OU_REFERENCE_PAYLOADS
            | GEIGER_FITCONTINUOUS_EARLY_BURST_REFERENCE_PAYLOADS
        ),
    )
    report = run_geiger_parity_cases(
        case_ids=[
            "fitcontinuous-lambda-weak-signal-review",
            "fitcontinuous-ou-lower-boundary-review",
            "fitcontinuous-eb-lower-boundary-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    rows_by_case = {row.case_id: row for row in report.boundary_warning_rows}
    lambda_row = rows_by_case["fitcontinuous-lambda-weak-signal-review"]
    assert lambda_row.affected_parameter == "lambda"
    assert lambda_row.reference_hit_lower_boundary is True
    assert "weak_phylogenetic_signal" in lambda_row.reference_boundary_warning_kinds
    assert lambda_row.reference_stable_conclusion_supported is False
    assert lambda_row.stable_conclusion_supported_match is True
    ou_row = rows_by_case["fitcontinuous-ou-lower-boundary-review"]
    assert ou_row.affected_parameter == "alpha"
    assert ou_row.reference_hit_lower_boundary is True
    assert "weak_pull_to_optimum" in ou_row.reference_boundary_warning_kinds
    assert ou_row.reference_stable_conclusion_supported is False
    eb_row = rows_by_case["fitcontinuous-eb-lower-boundary-review"]
    assert eb_row.affected_parameter == "rate_change"
    assert eb_row.reference_hit_lower_boundary is True
    assert "brownian_like_rate_change" in eb_row.reference_boundary_warning_kinds


def test_run_geiger_parity_cases_derives_flat_boundary_review_from_fake_profiles(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_LAMBDA_REFERENCE_PAYLOADS,
    )
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-lambda-weak-signal-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    row = report.boundary_warning_rows[0]
    assert row.reference_flat_likelihood_near_boundary is False
    assert row.bijux_flat_likelihood_near_boundary is True
    assert "flat_likelihood" not in row.reference_boundary_warning_kinds
    assert "flat_likelihood" in row.bijux_boundary_warning_kinds


def test_run_geiger_parity_cases_builds_single_fit_likelihood_policy_rows(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=(
            GEIGER_FITCONTINUOUS_BROWNIAN_REFERENCE_PAYLOADS
            | GEIGER_FITCONTINUOUS_WHITE_REFERENCE_PAYLOADS
        ),
    )
    report = run_geiger_parity_cases(
        case_ids=[
            "fitcontinuous-bm-example-tree",
            "fitcontinuous-white-weak-signal-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    rows_by_case = {row.case_id: row for row in report.likelihood_policy_rows}
    brownian = rows_by_case["fitcontinuous-bm-example-tree"]
    assert brownian.reference_likelihood_constant_policy == (
        "full-gaussian-loglikelihood-includes-normalizing-constant"
    )
    assert brownian.bijux_likelihood_constant_policy == (
        "full-gaussian-loglikelihood-includes-normalizing-constant"
    )
    assert brownian.case_level_raw_log_likelihood_comparable is True
    assert brownian.raw_log_likelihood_match_within_tolerance is True
    assert brownian.reference_aic_matches_raw_log_likelihood is True
    assert brownian.bijux_aic_matches_raw_log_likelihood is True
    assert brownian.ranking_permitted is False
    white = rows_by_case["fitcontinuous-white-weak-signal-review"]
    assert white.reference_aic_matches_raw_log_likelihood is True
    assert white.bijux_aic_matches_raw_log_likelihood is True


def test_run_geiger_parity_cases_builds_model_comparison_likelihood_policy_rows(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_MODEL_COMPARISON_REFERENCE_PAYLOADS,
    )
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-model-comparison-brownian-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    row = report.likelihood_policy_rows[0]
    assert row.case_id == "fitcontinuous-model-comparison-brownian-review"
    assert row.case_level_raw_log_likelihood_comparable is False
    assert row.raw_log_likelihood_match_within_tolerance is None
    assert row.reference_aic_matches_raw_log_likelihood is True
    assert row.bijux_aic_matches_raw_log_likelihood is True
    assert row.relative_aic_comparable is True
    assert row.ranking_permitted is True
    assert row.ranking_guard_outcome == "shared-fitcontinuous-policy-ranking-permitted"


def test_run_geiger_parity_cases_builds_model_confidence_rows(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_MODEL_COMPARISON_REFERENCE_PAYLOADS,
    )
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-model-comparison-brownian-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    rows_by_model = {row.candidate_model: row for row in report.model_confidence_rows}
    brownian_row = rows_by_model["brownian"]
    assert brownian_row.reference_best_model == "brownian"
    assert brownian_row.bijux_best_model == "brownian"
    assert brownian_row.reference_selected is True
    assert brownian_row.bijux_selected is True
    assert brownian_row.reference_akaike_weight is not None
    assert brownian_row.bijux_akaike_weight is not None
    assert brownian_row.akaike_weight_match_within_tolerance is True
    assert brownian_row.reference_within_delta_aic_threshold is True
    assert brownian_row.within_delta_aic_threshold_match is True
    runner_up_row = next(
        row for row in report.model_confidence_rows if row.reference_rank == 2
    )
    assert runner_up_row.candidate_model == "ornstein-uhlenbeck"
    assert runner_up_row.reference_within_delta_aic_threshold is True
    assert runner_up_row.reference_within_delta_aicc_threshold is False
    white_row = rows_by_model["white-noise"]
    assert white_row.reference_within_delta_aic_threshold is False
    assert white_row.bijux_within_delta_aic_threshold is False
    assert white_row.reference_uncertainty_class in {"limited", "moderate", "strong"}
    assert white_row.bijux_uncertainty_class in {"limited", "moderate", "strong"}


def test_run_geiger_parity_cases_builds_direct_parameterization_registry_rows(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_LAMBDA_REFERENCE_PAYLOADS,
    )
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-lambda-weak-signal-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    row = report.parameterization_registry_rows[0]
    assert row.case_id == "fitcontinuous-lambda-weak-signal-review"
    assert row.reference_parameter_name == "lambda"
    assert row.bijux_parameter_name == "lambda"
    assert row.canonical_parameter_name == "lambda"
    assert row.parameter_match_after_conversion is True
    assert row.parameter_bounds_match_after_conversion is True
    assert row.expected_divergence is False
    assert row.root_state_match_within_tolerance is True
    assert row.variance_match_within_tolerance is True
    assert row.log_likelihood_match_within_tolerance is True


def test_run_geiger_parity_cases_marks_expected_early_burst_parameterization_divergence(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_EARLY_BURST_REFERENCE_PAYLOADS,
    )
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-eb-early-burst-rate-recovery"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    row = report.parameterization_registry_rows[0]
    assert row.case_id == "fitcontinuous-eb-early-burst-rate-recovery"
    assert row.reference_parameter_name == "a"
    assert row.bijux_parameter_name == "rate_change"
    assert row.canonical_parameter_name == "rate_change"
    assert row.reference_parameter_value is not None
    assert row.converted_reference_parameter_value is not None
    assert row.bijux_parameter_value is not None
    assert row.parameter_match_after_conversion is True
    assert row.expected_divergence is True
    assert (
        row.expected_divergence_kind
        == "continuous-early-burst-sign-and-bound-convention"
    )
    assert row.expected_divergence_evidence is not None
    assert row.parameter_bounds_match_after_conversion is False
