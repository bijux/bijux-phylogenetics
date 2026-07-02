from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from bijux_phylogenetics.parity import (
    list_geiger_parity_cases,
    run_geiger_parity_cases,
    write_geiger_boundary_warning_table,
    write_geiger_likelihood_policy_table,
    write_geiger_model_confidence_table,
    write_geiger_optimizer_triage_table,
    write_geiger_parameterization_registry_table,
    write_geiger_parity_observation_table,
    write_geiger_parity_summary_table,
)
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
from tests.support.geiger_fitdiscrete_ard_reference import (
    GEIGER_FITDISCRETE_ARD_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitdiscrete_delta_reference import (
    GEIGER_FITDISCRETE_DELTA_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitdiscrete_early_burst_reference import (
    GEIGER_FITDISCRETE_EARLY_BURST_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitdiscrete_er_reference import (
    GEIGER_FITDISCRETE_ER_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitdiscrete_kappa_reference import (
    GEIGER_FITDISCRETE_KAPPA_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitdiscrete_lambda_reference import (
    GEIGER_FITDISCRETE_LAMBDA_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitdiscrete_sym_reference import (
    GEIGER_FITDISCRETE_SYM_REFERENCE_PAYLOADS,
)
from tests.support.geiger_optimizer_triage_reference import (
    geiger_optimizer_triage_reference_payloads,
)

pytestmark = pytest.mark.slow


def test_run_geiger_parity_cases_governs_fitdiscrete_er_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITDISCRETE_ER_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitdiscrete-er-binary-twenty-four-taxa",
            "fitdiscrete-er-multistate-missing-twenty-four-taxa",
            "fitdiscrete-er-multistate-tip-intersection-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 3
    assert report.passed_case_count == 3
    binary = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-er-binary-twenty-four-taxa"
    )
    assert binary.reference_summary is not None
    assert binary.reference_summary["log_likelihood"] == -9.078105640476831
    assert binary.reference_rows is not None
    assert binary.reference_rows == [
        {
            "source_state": "0",
            "target_state": "1",
            "transition_allowed": True,
            "step_distance": 1,
            "rate": 0.393523166730309,
        },
        {
            "source_state": "1",
            "target_state": "0",
            "transition_allowed": True,
            "step_distance": 1,
            "rate": 0.393523166730309,
        },
    ]
    assert binary.bijux_summary is not None
    assert binary.bijux_summary["missing_value_policy"] == (
        "prune-overlapping-missing-values"
    )
    missing = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-er-multistate-missing-twenty-four-taxa"
    )
    assert missing.reference_summary is not None
    assert missing.reference_summary["missing_value_taxa"] == ["Phy14"]
    assert missing.bijux_summary is not None
    assert missing.bijux_summary["excluded_taxa"] == ["Phy14"]
    mismatch = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-er-multistate-tip-intersection-review"
    )
    assert mismatch.reference_summary is not None
    assert mismatch.reference_summary["missing_from_traits"] == ["Phy9"]
    assert mismatch.reference_summary["extra_trait_taxa"] == ["PhyExtra"]
    assert mismatch.bijux_summary is not None
    assert mismatch.bijux_summary["missing_value_taxa"] == []


def test_run_geiger_parity_cases_governs_fitdiscrete_sym_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITDISCRETE_SYM_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitdiscrete-sym-three-state-twenty-four-taxa",
            "fitdiscrete-sym-four-state-twenty-four-taxa",
            "fitdiscrete-sym-multistate-missing-twenty-four-taxa",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 3
    assert report.passed_case_count == 3
    three_state = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-sym-three-state-twenty-four-taxa"
    )
    assert three_state.reference_summary is not None
    assert three_state.reference_summary["log_likelihood"] == -15.16327924942547
    assert three_state.reference_rows is not None
    assert three_state.reference_rows[0]["source_state"] == "central"
    assert three_state.reference_rows[0]["target_state"] == "north"
    assert three_state.reference_rows[0]["rate"] == 0.727456725972813
    four_state = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-sym-four-state-twenty-four-taxa"
    )
    assert four_state.reference_summary is not None
    assert four_state.reference_summary["observed_state_count"] == 4
    assert four_state.reference_rows is not None
    assert any(
        row["source_state"] == "south"
        and row["target_state"] == "west"
        and row["rate"] == 4.57831457247917
        for row in four_state.reference_rows
    )
    missing = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-sym-multistate-missing-twenty-four-taxa"
    )
    assert missing.reference_summary is not None
    assert missing.reference_summary["missing_value_taxa"] == ["Phy14"]
    assert missing.bijux_summary is not None
    assert missing.bijux_summary["excluded_taxa"] == ["Phy14"]


def test_run_geiger_parity_cases_governs_fitdiscrete_lambda_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITDISCRETE_LAMBDA_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitdiscrete-lambda-strong-signal-review",
            "fitdiscrete-lambda-weak-signal-review",
            "fitdiscrete-lambda-missing-values-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 3
    assert report.passed_case_count == 3
    strong = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-lambda-strong-signal-review"
    )
    assert strong.reference_summary is not None
    assert strong.reference_summary["transform_name"] == "pagel-lambda"
    assert strong.reference_summary["parameter_value"] == 1
    assert strong.reference_rows is not None
    assert strong.reference_rows[0]["rate"] == 0.393523199371205
    assert strong.bijux_summary is not None
    assert strong.bijux_summary["transform_name"] == "pagel-lambda"
    strong_case = next(
        item
        for item in list_geiger_parity_cases()
        if item.case_id == "fitdiscrete-lambda-strong-signal-review"
    )
    assert "parameter_value" in strong_case.comparison_fields
    weak = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-lambda-weak-signal-review"
    )
    assert weak.reference_summary is not None
    assert weak.reference_summary["trait_name"] == "er_binary_transform_weak_signal"
    assert weak.reference_summary["parameter_value"] < 1e-40
    assert weak.reference_rows is not None
    assert weak.reference_rows[0]["rate"] == 31.5169313693995
    assert weak.bijux_summary is not None
    assert weak.bijux_summary["parameter_value"] == 0.0
    weak_case = next(
        item
        for item in list_geiger_parity_cases()
        if item.case_id == "fitdiscrete-lambda-weak-signal-review"
    )
    assert weak_case.row_comparison_policy == "summary-only"
    assert "parameter_value" not in weak_case.comparison_fields
    missing = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-lambda-missing-values-review"
    )
    assert missing.reference_summary is not None
    assert missing.reference_summary["missing_value_taxa"] == ["Phy10"]
    assert missing.bijux_summary is not None
    assert missing.bijux_summary["excluded_taxa"] == ["Phy10"]


def test_run_geiger_parity_cases_accepts_flat_lambda_plateau_review_surface(
    tmp_path: Path,
) -> None:
    payloads = dict(GEIGER_FITDISCRETE_LAMBDA_REFERENCE_PAYLOADS)
    weak_payload = dict(payloads["fitdiscrete-lambda-weak-signal-review"])
    weak_summary = dict(weak_payload["summary"])
    weak_summary["parameter_value"] = 0.7215146376631235
    weak_payload["summary"] = weak_summary
    payloads["fitdiscrete-lambda-weak-signal-review"] = weak_payload

    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=payloads,
    )

    report = run_geiger_parity_cases(
        case_ids=["fitdiscrete-lambda-weak-signal-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 1
    assert report.passed_case_count == 1
    observation = report.observations[0]
    assert observation.passed is True
    assert observation.mismatch_reason is None
    assert observation.reference_summary is not None
    assert observation.reference_summary["parameter_value"] == 0.7215146376631235
    assert observation.bijux_summary is not None
    assert observation.bijux_summary["parameter_value"] == 0.0


def test_run_geiger_parity_cases_governs_fitdiscrete_kappa_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITDISCRETE_KAPPA_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitdiscrete-kappa-strong-signal-review",
            "fitdiscrete-kappa-weak-signal-review",
            "fitdiscrete-kappa-missing-values-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 3
    assert report.passed_case_count == 3
    strong = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-kappa-strong-signal-review"
    )
    assert strong.reference_summary is not None
    assert strong.reference_summary["transform_name"] == "pagel-kappa"
    assert strong.reference_summary["parameter_value"] == 0.901187058855041
    assert strong.reference_rows is not None
    assert strong.reference_rows[0]["rate"] == 0.348699760584072
    assert strong.bijux_summary is not None
    assert strong.bijux_summary["transform_name"] == "pagel-kappa"
    weak = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-kappa-weak-signal-review"
    )
    assert weak.reference_summary is not None
    assert weak.reference_summary["trait_name"] == "er_binary_transform_weak_signal"
    assert weak.reference_summary["parameter_value"] < 0.01
    assert weak.reference_rows is not None
    assert weak.reference_rows[0]["rate"] == 61.5648011375806
    assert weak.bijux_summary is not None
    assert weak.bijux_summary["parameter_value"] == 0.0
    missing = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-kappa-missing-values-review"
    )
    assert missing.reference_summary is not None
    assert missing.reference_summary["missing_value_taxa"] == ["Phy14"]
    assert missing.bijux_summary is not None
    assert missing.bijux_summary["excluded_taxa"] == ["Phy14"]


def test_run_geiger_parity_cases_governs_fitdiscrete_delta_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITDISCRETE_DELTA_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitdiscrete-delta-late-change-boundary-review",
            "fitdiscrete-delta-earliest-change-review",
            "fitdiscrete-delta-missing-values-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 3
    assert report.passed_case_count == 3
    boundary = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-delta-late-change-boundary-review"
    )
    assert boundary.reference_summary is not None
    assert boundary.reference_summary["transform_name"] == "pagel-delta"
    assert boundary.reference_summary["parameter_value"] == 2.999999
    assert boundary.reference_rows is not None
    assert boundary.reference_rows[0]["rate"] == 0.2282873
    assert boundary.bijux_summary is not None
    assert boundary.bijux_summary["transform_name"] == "pagel-delta"
    earliest = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-delta-earliest-change-review"
    )
    assert earliest.reference_summary is not None
    assert earliest.reference_summary["trait_name"] == "sym_three_state_truth"
    assert earliest.reference_summary["parameter_value"] == 0.006737947
    assert earliest.bijux_summary is not None
    assert abs(earliest.bijux_summary["parameter_value"] - 0.006737947) < 1e-12
    missing = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-delta-missing-values-review"
    )
    assert missing.reference_summary is not None
    assert missing.reference_summary["missing_value_taxa"] == ["Phy14"]
    assert missing.bijux_summary is not None
    assert missing.bijux_summary["excluded_taxa"] == ["Phy14"]


def test_run_geiger_parity_cases_governs_fitdiscrete_early_burst_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITDISCRETE_EARLY_BURST_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitdiscrete-early-burst-early-change-review",
            "fitdiscrete-early-burst-weak-signal-review",
            "fitdiscrete-early-burst-late-change-review",
            "fitdiscrete-early-burst-missing-values-review",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 4
    assert report.passed_case_count == 4
    early = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-early-burst-early-change-review"
    )
    assert early.reference_summary is not None
    assert early.reference_summary["transform_name"] == "early-burst"
    assert early.reference_summary["parameter_value"] == 2.22884
    assert early.reference_rows is not None
    assert early.reference_rows[0]["rate"] == 0.01886852
    assert early.bijux_summary is not None
    assert early.bijux_summary["transform_name"] == "early-burst"
    weak = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-early-burst-weak-signal-review"
    )
    assert weak.reference_summary is not None
    assert weak.reference_summary["trait_name"] == "er_binary_transform_weak_signal"
    assert weak.reference_summary["parameter_value"] > 9.0
    assert weak.bijux_summary is not None
    assert weak.bijux_summary["parameter_value"] == 0.0
    late = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-early-burst-late-change-review"
    )
    assert late.reference_summary is not None
    assert late.reference_summary["parameter_value"] < 0.0
    assert late.bijux_summary is not None
    assert late.bijux_summary["parameter_value"] < 0.0
    missing = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-early-burst-missing-values-review"
    )
    assert missing.reference_summary is not None
    assert missing.reference_summary["missing_value_taxa"] == ["Phy10"]
    assert missing.bijux_summary is not None
    assert missing.bijux_summary["excluded_taxa"] == ["Phy10"]


def test_run_geiger_parity_cases_governs_fitdiscrete_ard_reference_payloads(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITDISCRETE_ARD_REFERENCE_PAYLOADS,
    )

    report = run_geiger_parity_cases(
        case_ids=[
            "fitdiscrete-ard-binary-twenty-four-taxa",
            "fitdiscrete-ard-four-state-weak-identification-review",
            "fitdiscrete-ard-multistate-missing-twenty-four-taxa",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 3
    assert report.passed_case_count == 3
    binary = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-ard-binary-twenty-four-taxa"
    )
    assert binary.reference_summary is not None
    assert binary.reference_summary["log_likelihood"] == -10.750446663416824
    assert binary.reference_rows is not None
    assert binary.reference_rows[0]["source_state"] == "0"
    assert binary.reference_rows[0]["rate"] == 1.31428822804799
    assert binary.bijux_summary is not None
    assert binary.bijux_summary["optimizer_result"]["converged"] is True
    four_state = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-ard-four-state-weak-identification-review"
    )
    assert four_state.reference_summary is not None
    assert four_state.reference_summary["parameter_count"] == 12
    assert four_state.reference_rows is not None
    assert any(
        row["source_state"] == "west"
        and row["target_state"] == "south"
        and row["rate"] == 8.15930951548937
        for row in four_state.reference_rows
    )
    assert four_state.bijux_summary is not None
    assert four_state.bijux_summary["optimizer_result"]["converged"] is False
    assert (
        four_state.bijux_summary["optimizer_result"]["hit_lower_parameter_bound"]
        is True
    )
    missing = next(
        item
        for item in report.observations
        if item.case_id == "fitdiscrete-ard-multistate-missing-twenty-four-taxa"
    )
    assert missing.reference_summary is not None
    assert missing.reference_summary["missing_value_taxa"] == ["Phy14"]
    assert missing.reference_rows is not None
    assert missing.reference_rows[1]["rate"] == 2.37295469511749
    assert missing.bijux_summary is not None
    assert missing.bijux_summary["excluded_taxa"] == ["Phy14"]


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


def test_write_geiger_optimizer_triage_table_writes_rows(tmp_path: Path) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-lambda-weak-signal-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )
    triage_path = tmp_path / "geiger-optimizer-triage.tsv"

    write_geiger_optimizer_triage_table(triage_path, report)

    with triage_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["case_id"] == "fitcontinuous-lambda-weak-signal-review"
    assert rows[0]["mismatch_type"] == "no_algorithm_mismatch"


def test_run_geiger_parity_cases_builds_boundary_warning_rows(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=(
            GEIGER_FITCONTINUOUS_LAMBDA_REFERENCE_PAYLOADS
            | GEIGER_FITCONTINUOUS_OU_REFERENCE_PAYLOADS
            | GEIGER_FITCONTINUOUS_EARLY_BURST_REFERENCE_PAYLOADS
            | GEIGER_FITCONTINUOUS_KAPPA_REFERENCE_PAYLOADS
            | GEIGER_FITCONTINUOUS_DELTA_REFERENCE_PAYLOADS
        ),
    )
    report = run_geiger_parity_cases(
        case_ids=[
            "fitcontinuous-lambda-weak-signal-review",
            "fitcontinuous-ou-lower-boundary-review",
            "fitcontinuous-eb-lower-boundary-review",
            "fitcontinuous-kappa-weak-signal-review",
            "fitcontinuous-delta-weak-signal-review",
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
    kappa_row = rows_by_case["fitcontinuous-kappa-weak-signal-review"]
    assert kappa_row.affected_parameter == "kappa"
    assert kappa_row.reference_hit_lower_boundary is True
    assert "punctuational_limit" in kappa_row.reference_boundary_warning_kinds
    delta_row = rows_by_case["fitcontinuous-delta-weak-signal-review"]
    assert delta_row.affected_parameter == "delta"
    assert delta_row.reference_hit_upper_boundary is True
    assert "late_change_limit" in delta_row.reference_boundary_warning_kinds


def test_run_geiger_parity_cases_derives_flat_boundary_review_from_fake_profiles(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-lambda-weak-signal-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    row = report.boundary_warning_rows[0]
    assert row.reference_flat_likelihood_near_boundary is True
    assert row.bijux_flat_likelihood_near_boundary is True
    assert "flat_likelihood" in row.reference_boundary_warning_kinds
    assert "flat_likelihood" in row.bijux_boundary_warning_kinds


def test_write_geiger_boundary_warning_table_writes_rows(tmp_path: Path) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-lambda-weak-signal-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )
    boundary_path = tmp_path / "geiger-boundary-warning.tsv"

    write_geiger_boundary_warning_table(boundary_path, report)

    with boundary_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["case_id"] == "fitcontinuous-lambda-weak-signal-review"
    assert rows[0]["reference_hit_lower_boundary"] == "True"
    assert "weak_phylogenetic_signal" in rows[0]["bijux_boundary_warning_kinds"]


def test_run_geiger_parity_cases_builds_single_fit_likelihood_policy_rows(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")
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


def test_write_geiger_likelihood_policy_table_writes_rows(tmp_path: Path) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-bm-example-tree"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )
    policy_path = tmp_path / "geiger-likelihood-policy.tsv"

    write_geiger_likelihood_policy_table(policy_path, report)

    with policy_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["case_id"] == "fitcontinuous-bm-example-tree"
    assert rows[0]["case_level_raw_log_likelihood_comparable"] == "True"


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


def test_write_geiger_model_confidence_table_writes_rows(tmp_path: Path) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_MODEL_COMPARISON_REFERENCE_PAYLOADS,
    )
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-model-comparison-brownian-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )
    confidence_path = tmp_path / "geiger-model-confidence.tsv"

    write_geiger_model_confidence_table(confidence_path, report)

    with confidence_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 7
    assert rows[0]["case_id"] == "fitcontinuous-model-comparison-brownian-review"
    assert rows[0]["weight_basis"] == "AICc"


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


def test_write_geiger_parameterization_registry_table_writes_rows(
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
    registry_path = tmp_path / "geiger-parameterization-registry.tsv"

    write_geiger_parameterization_registry_table(registry_path, report)

    with registry_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["case_id"] == "fitcontinuous-lambda-weak-signal-review"
    assert rows[0]["reference_parameter_name"] == "lambda"


def test_write_geiger_parity_tables_writes_summary_and_observations(
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
    summary_path = tmp_path / "geiger-parity-summary.tsv"
    observation_path = tmp_path / "geiger-parity-observations.tsv"
    triage_path = tmp_path / "geiger-optimizer-triage.tsv"

    write_geiger_parity_summary_table(summary_path, report)
    write_geiger_parity_observation_table(observation_path, report)
    write_geiger_optimizer_triage_table(triage_path, report)

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("function_name\tcase_count")
    assert any("geiger::fitContinuous(model='OU')" in row for row in summary_rows[1:])
    with observation_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 10
    assert rows[0]["model_name"] in {
        "BM",
        "OU",
        "ER",
    }
    optimizer_settings = json.loads(rows[0]["optimizer_settings"])
    assert "reference_control_policy" in optimizer_settings
    with triage_path.open(encoding="utf-8", newline="") as handle:
        triage_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(triage_rows) == len(rows)
    assert triage_rows[0]["mismatch_type"] in {
        "no_algorithm_mismatch",
        "parameter_surface_not_applicable",
    }
