from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.parity import run_geiger_parity_cases
from tests.support.fake_geiger_parity import fake_geiger_rscript
from tests.support.geiger_fitcontinuous_brownian_reference import (
    GEIGER_FITCONTINUOUS_BROWNIAN_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitcontinuous_delta_reference import (
    GEIGER_FITCONTINUOUS_DELTA_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitcontinuous_early_burst_reference import (
    GEIGER_FITCONTINUOUS_EARLY_BURST_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitcontinuous_kappa_reference import (
    GEIGER_FITCONTINUOUS_KAPPA_REFERENCE_PAYLOADS,
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
from tests.support.geiger_fitdiscrete_er_reference import (
    GEIGER_FITDISCRETE_ER_REFERENCE_PAYLOADS,
)

pytestmark = pytest.mark.slow


def test_run_geiger_parity_cases_passes_representative_shards_against_fake_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=(
            GEIGER_FITCONTINUOUS_BROWNIAN_REFERENCE_PAYLOADS
            | GEIGER_FITCONTINUOUS_OU_REFERENCE_PAYLOADS
            | GEIGER_FITDISCRETE_ER_REFERENCE_PAYLOADS
        ),
    )

    report = run_geiger_parity_cases(
        shard_ids=[
            "geiger-fitcontinuous-model-bm",
            "geiger-fitcontinuous-model-ou",
            "geiger-fitdiscrete-model-er",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 10
    assert report.passed_case_count == 10
    assert report.failed_case_count == 0
    assert report.skipped_case_count == 0
    assert report.all_passed is True
    assert [row.function_name for row in report.summary_rows] == [
        "geiger::fitContinuous(model='BM')",
        "geiger::fitContinuous(model='OU')",
        "geiger::fitDiscrete(model='ER')",
    ]
    observation = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-ou-lower-boundary-review"
    )
    assert observation.model_name == "OU"
    assert observation.r_version == "4.4.0"
    assert observation.geiger_version == "2.0.11"
    assert observation.optimizer_settings is not None
    assert observation.reference_rows is not None


def test_run_geiger_parity_cases_counts_skips_for_selected_shards_when_geiger_is_unavailable(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        geiger_available=False,
    )

    report = run_geiger_parity_cases(
        shard_ids=[
            "geiger-fitcontinuous-model-bm",
            "geiger-fitdiscrete-model-er",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 6
    assert report.passed_case_count == 0
    assert report.failed_case_count == 0
    assert report.skipped_case_count == 6
    assert report.all_passed is False
    assert all(
        item.mismatch_reason == "geiger_package_unavailable"
        for item in report.observations
    )


def test_run_geiger_parity_cases_governs_brownian_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_BROWNIAN_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitcontinuous-bm-example-tree",
            "fitcontinuous-bm-brownian-sigma-recovery",
            "fitcontinuous-bm-missing-values-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 3
    assert report.passed_case_count == 3
    missing_values = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-bm-missing-values-review"
    )
    assert missing_values.reference_summary is not None
    assert missing_values.reference_summary["excluded_taxa"] == ["Phy10", "Phy14"]
    assert missing_values.reference_summary["missing_value_taxa"] == ["Phy10"]
    assert missing_values.reference_summary["non_numeric_taxa"] == ["Phy14"]
    assert missing_values.bijux_summary is not None
    assert missing_values.bijux_summary["missing_value_policy"] == (
        "prune-tree-tip-overlap-with-missing-or-nonnumeric-trait-values"
    )
    assert missing_values.bijux_summary["standard_error_policy"] == (
        "fitcontinuous-standard-error-explicitly-excluded-this-round"
    )


def test_run_geiger_parity_cases_governs_ou_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_OU_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitcontinuous-ou-ou-parameter-recovery",
            "fitcontinuous-ou-missing-values-review",
            "fitcontinuous-ou-lower-boundary-review",
            "fitcontinuous-ou-bounded-control-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 4
    assert report.passed_case_count == 4
    missing_values = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-ou-missing-values-review"
    )
    assert missing_values.reference_summary is not None
    assert (
        missing_values.reference_summary["aicc"]
        > missing_values.reference_summary["aic"]
    )
    assert missing_values.reference_summary["missing_value_taxa"] == ["Phy10"]
    assert missing_values.bijux_summary is not None
    assert missing_values.bijux_summary["identifiability_warning_kinds"] == [
        "flat_likelihood",
        "weak_pull_to_optimum",
    ]
    lower_boundary = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-ou-lower-boundary-review"
    )
    assert lower_boundary.reference_summary is not None
    assert lower_boundary.reference_summary["hit_lower_parameter_boundary"] is True
    assert lower_boundary.bijux_summary is not None
    assert lower_boundary.bijux_summary["hit_lower_parameter_boundary"] is True
    assert lower_boundary.bijux_summary["identifiability_warning_kinds"] == [
        "boundary_alpha",
        "flat_likelihood",
        "weak_pull_to_optimum",
    ]
    bounded = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-ou-bounded-control-review"
    )
    assert bounded.reference_summary is not None
    assert bounded.reference_summary["hit_upper_parameter_boundary"] is True
    assert bounded.reference_summary["optimizer_settings"] == {
        "bijux_coarse_grid_point_count": 61,
        "bijux_fine_grid_point_count": 41,
        "bijux_initial_parameter_value": 0.45,
        "bijux_optimizer_name": "governed-two-stage-grid-search",
        "bijux_parameter_bounds": {"lower": 0.2, "upper": 1.0},
        "reference_control_policy": "fitcontinuous-explicit-control",
        "reference_control_settings": {
            "CI": 0.95,
            "FAIL": 1e200,
            "hessian": False,
            "method": "subplex",
            "niter": 8,
        },
    }
    assert bounded.bijux_summary is not None
    assert bounded.bijux_summary["optimizer_result"]["coarse_grid_point_count"] == 61
    assert bounded.bijux_summary["optimizer_result"]["fine_grid_point_count"] == 41
    assert bounded.bijux_summary["optimizer_result"]["starting_parameter_value"] == 0.45


def test_run_geiger_parity_cases_governs_lambda_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_LAMBDA_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitcontinuous-lambda-strong-signal-review",
            "fitcontinuous-lambda-weak-signal-review",
            "fitcontinuous-lambda-missing-values-review",
            "fitcontinuous-lambda-bounded-control-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 4
    assert report.passed_case_count == 4
    strong = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-lambda-strong-signal-review"
    )
    assert strong.reference_summary is not None
    assert strong.reference_summary["hit_upper_parameter_boundary"] is True
    assert strong.bijux_summary is not None
    assert strong.bijux_summary["identifiability_warning_kinds"] == [
        "boundary_lambda",
        "flat_likelihood",
        "brownian_limit",
    ]
    weak = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-lambda-weak-signal-review"
    )
    assert weak.reference_summary is not None
    assert weak.reference_summary["hit_lower_parameter_boundary"] is True
    assert weak.bijux_summary is not None
    assert weak.bijux_summary["parameter_value"] == 0.0
    assert weak.bijux_summary["identifiability_warning_kinds"] == [
        "boundary_lambda",
        "flat_likelihood",
        "weak_phylogenetic_signal",
    ]
    missing_values = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-lambda-missing-values-review"
    )
    assert missing_values.reference_summary is not None
    assert missing_values.reference_summary["missing_value_taxa"] == ["Phy10"]
    assert missing_values.reference_summary["non_numeric_taxa"] == ["Phy14"]
    assert missing_values.bijux_summary is not None
    assert missing_values.bijux_summary["excluded_taxa"] == ["Phy10", "Phy14"]
    bounded = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-lambda-bounded-control-review"
    )
    assert bounded.reference_summary is not None
    assert bounded.reference_summary["hit_upper_parameter_boundary"] is True
    assert bounded.reference_summary["optimizer_settings"] == {
        "bijux_coarse_grid_point_count": 61,
        "bijux_fine_grid_point_count": 41,
        "bijux_initial_parameter_value": 0.35,
        "bijux_optimizer_name": "governed-two-stage-grid-search",
        "bijux_parameter_bounds": {"lower": 0.2, "upper": 0.6},
        "reference_control_policy": "fitcontinuous-explicit-control",
        "reference_control_settings": {
            "CI": 0.95,
            "FAIL": 1e200,
            "hessian": False,
            "method": "subplex",
            "niter": 8,
        },
    }
    assert bounded.bijux_summary is not None
    assert (
        bounded.bijux_summary["optimizer_result"]["starting_parameter_policy"]
        == "user-provided-first-evaluation"
    )
    assert bounded.bijux_summary["optimizer_result"]["starting_parameter_value"] == 0.35


def test_run_geiger_parity_cases_governs_white_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_WHITE_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitcontinuous-white-strong-signal-review",
            "fitcontinuous-white-weak-signal-review",
            "fitcontinuous-white-missing-values-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 3
    assert report.passed_case_count == 3
    strong = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-white-strong-signal-review"
    )
    assert strong.reference_summary is not None
    assert strong.reference_summary["parameter_bound_policy"] == (
        "reference-default-without-explicit-bounds"
    )
    assert strong.bijux_summary is not None
    assert strong.bijux_summary["parameter_bound_policy"] == (
        "closed-form-without-parameter-bounds"
    )
    assert strong.bijux_summary["identifiability_warning_kinds"] == [
        "no_phylogenetic_correlation"
    ]
    weak = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-white-weak-signal-review"
    )
    assert weak.reference_summary is not None
    assert weak.reference_summary["rate"] > 1.0
    assert weak.bijux_summary is not None
    assert weak.bijux_summary["parameter_name"] is None
    missing_values = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-white-missing-values-review"
    )
    assert missing_values.reference_summary is not None
    assert missing_values.reference_summary["missing_value_taxa"] == ["Phy10"]
    assert missing_values.reference_summary["non_numeric_taxa"] == ["Phy14"]
    assert missing_values.bijux_summary is not None
    assert missing_values.bijux_summary["excluded_taxa"] == ["Phy10", "Phy14"]


def test_run_geiger_parity_cases_governs_kappa_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_KAPPA_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitcontinuous-kappa-strong-signal-review",
            "fitcontinuous-kappa-weak-signal-review",
            "fitcontinuous-kappa-missing-values-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 3
    assert report.passed_case_count == 3
    strong = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-kappa-strong-signal-review"
    )
    assert strong.reference_summary is not None
    assert strong.reference_summary["parameter_value"] > 1.0
    assert strong.bijux_summary is not None
    assert strong.bijux_summary["identifiability_warning_kinds"] == ["flat_likelihood"]
    weak = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-kappa-weak-signal-review"
    )
    assert weak.reference_summary is not None
    assert weak.reference_summary["hit_lower_parameter_boundary"] is True
    assert weak.bijux_summary is not None
    assert weak.bijux_summary["parameter_value"] == 0.0
    assert weak.bijux_summary["identifiability_warning_kinds"] == [
        "boundary_kappa",
        "flat_likelihood",
        "punctuational_limit",
    ]
    missing_values = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-kappa-missing-values-review"
    )
    assert missing_values.reference_summary is not None
    assert missing_values.reference_summary["missing_value_taxa"] == ["Phy10"]
    assert missing_values.reference_summary["non_numeric_taxa"] == ["Phy14"]
    assert missing_values.bijux_summary is not None
    assert missing_values.bijux_summary["excluded_taxa"] == ["Phy10", "Phy14"]


def test_run_geiger_parity_cases_governs_delta_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_DELTA_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitcontinuous-delta-strong-signal-review",
            "fitcontinuous-delta-weak-signal-review",
            "fitcontinuous-delta-missing-values-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 3
    assert report.passed_case_count == 3
    strong = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-delta-strong-signal-review"
    )
    assert strong.reference_summary is not None
    assert 1.0 < strong.reference_summary["parameter_value"] < 2.0
    assert strong.bijux_summary is not None
    assert strong.bijux_summary["identifiability_warning_kinds"] == ["flat_likelihood"]
    weak = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-delta-weak-signal-review"
    )
    assert weak.reference_summary is not None
    assert weak.reference_summary["hit_upper_parameter_boundary"] is True
    assert weak.bijux_summary is not None
    assert weak.bijux_summary["parameter_value"] == 3.0
    assert weak.bijux_summary["identifiability_warning_kinds"] == [
        "boundary_delta",
        "flat_likelihood",
        "late_change_limit",
    ]
    missing_values = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-delta-missing-values-review"
    )
    assert missing_values.reference_summary is not None
    assert missing_values.reference_summary["missing_value_taxa"] == ["Phy10"]
    assert missing_values.reference_summary["non_numeric_taxa"] == ["Phy14"]
    assert missing_values.bijux_summary is not None
    assert missing_values.bijux_summary["excluded_taxa"] == ["Phy10", "Phy14"]


def test_run_geiger_parity_cases_governs_early_burst_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_EARLY_BURST_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitcontinuous-eb-early-burst-rate-recovery",
            "fitcontinuous-eb-lower-boundary-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 2
    assert report.passed_case_count == 2
    recovery = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-eb-early-burst-rate-recovery"
    )
    assert recovery.reference_summary is not None
    assert recovery.reference_summary["aicc"] < recovery.reference_summary["aic"] + 2.0
    assert recovery.bijux_summary is not None
    assert recovery.bijux_summary["parameter_value"] > 4.0
    assert recovery.bijux_summary["identifiability_warning_kinds"] == [
        "flat_likelihood_profile"
    ]
    lower_boundary = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-eb-lower-boundary-review"
    )
    assert lower_boundary.reference_summary is not None
    assert lower_boundary.reference_summary["hit_lower_parameter_boundary"] is True
    assert lower_boundary.bijux_summary is not None
    assert lower_boundary.bijux_summary["hit_lower_parameter_boundary"] is True
    assert lower_boundary.bijux_summary["identifiability_warning_kinds"] == [
        "boundary_rate_change",
        "flat_likelihood_profile",
        "brownian_like_rate_change",
    ]


def test_run_geiger_parity_cases_governs_model_comparison_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_MODEL_COMPARISON_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitcontinuous-model-comparison-brownian-review",
            "fitcontinuous-model-comparison-ou-review",
            "fitcontinuous-model-comparison-early-burst-review",
            "fitcontinuous-model-comparison-white-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 4
    assert report.passed_case_count == 4
    brownian = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-model-comparison-brownian-review"
    )
    assert brownian.reference_summary is not None
    assert brownian.reference_summary["selected_model"] == "brownian"
    assert brownian.reference_summary["runner_up_model"] == "ornstein-uhlenbeck"
    assert brownian.reference_summary["warning_count"] == 1
    assert brownian.reference_rows is not None
    assert brownian.reference_rows[0]["model"] == "brownian"
    assert brownian.reference_rows[0]["selected"] is True
    assert brownian.reference_rows[-1]["model"] == "white-noise"
    assert brownian.reference_rows[-1]["comparable"] is True
    assert brownian.bijux_summary is not None
    assert brownian.bijux_summary["model_ranking"][:3] == [
        "brownian",
        "ornstein-uhlenbeck",
        "pagel-delta",
    ]
    white = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-model-comparison-white-review"
    )
    assert white.reference_summary is not None
    assert white.reference_summary["selected_model"] == "white-noise"
    assert white.reference_summary["runner_up_model"] == "pagel-lambda"
    assert white.bijux_summary is not None
    assert white.bijux_summary["noncomparable_model_count"] == 0


def test_run_geiger_parity_cases_persists_failure_artifacts_for_mismatches(
    tmp_path: Path,
) -> None:
    payloads = dict(GEIGER_FITCONTINUOUS_BROWNIAN_REFERENCE_PAYLOADS)
    mismatched_payload = dict(payloads["fitcontinuous-bm-example-tree"])
    mismatched_summary = dict(mismatched_payload["summary"])
    mismatched_summary["root_state"] = 999.0
    mismatched_payload["summary"] = mismatched_summary
    payloads["fitcontinuous-bm-example-tree"] = mismatched_payload

    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=payloads,
    )
    failure_root = tmp_path / "geiger-parity-failures"

    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-bm-example-tree"],
        rscript_executable=str(rscript),
        failure_root=failure_root,
    )

    assert report.case_count == 1
    assert report.failed_case_count == 1
    observation = report.observations[0]
    assert observation.status == "failed"
    assert observation.mismatch_reason == "summary_field_mismatch:root_state"
    assert observation.reproducible_artifact_root is not None
    artifact_root = observation.reproducible_artifact_root
    assert artifact_root.is_dir()
    assert (artifact_root / "case.json").is_file()
    assert (artifact_root / "reference-summary.json").is_file()
    assert (artifact_root / "bijux-summary.json").is_file()
    assert (artifact_root / "mismatch-reason.txt").read_text(
        encoding="utf-8"
    ) == "summary_field_mismatch:root_state"
    stored_summary = json.loads(
        (artifact_root / "reference-summary.json").read_text(encoding="utf-8")
    )
    assert stored_summary["root_state"] == 999.0
