from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.parity import list_geiger_parity_cases, run_geiger_parity_cases
from tests.support.fake_geiger_parity import fake_geiger_rscript
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
