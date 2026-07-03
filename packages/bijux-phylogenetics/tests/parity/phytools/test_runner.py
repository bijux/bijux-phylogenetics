from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.parity import run_phytools_parity_cases
from tests.support.fake_phytools_parity import fake_phytools_rscript


@pytest.mark.slow
def test_run_phytools_parity_cases_passes_representative_shards_against_fake_reference_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_phytools_rscript(tmp_path / "fake-phytools-rscript")

    report = run_phytools_parity_cases(
        shard_ids=[
            "phytools-densitymap",
            "phytools-make-simmap-model-er",
        ],
        rscript_executable=str(rscript),
    )

    assert report.all_passed is True
    assert report.case_count == 5
    assert report.failed_case_count == 0
    assert report.skipped_case_count == 0
    assert [row.function_name for row in report.summary_rows] == [
        "phytools::densityMap",
        "phytools::make.simmap(model='ER')",
    ]
    first = report.observations[0]
    assert first.phytools_version == "2.5.2"
    assert first.r_version == "4.6.0"
    assert first.reference_summary is not None
    assert {observation.function_name for observation in report.observations} == {
        "phytools::densityMap",
        "phytools::make.simmap(model='ER')",
    }


@pytest.mark.slow
def test_run_phytools_parity_cases_passes_density_map_cases_against_fake_reference_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_phytools_rscript(tmp_path / "fake-phytools-rscript")

    report = run_phytools_parity_cases(
        case_ids=[
            "density-map-er-binary-twenty-four-taxa",
            "density-map-er-binary-missing-twenty-four-taxa",
        ],
        rscript_executable=str(rscript),
    )

    assert report.all_passed is True
    assert report.case_count == 2
    assert report.failed_case_count == 0
    assert [row.function_name for row in report.summary_rows] == [
        "phytools::densityMap"
    ]
    assert all(
        observation.function_name == "phytools::densityMap"
        for observation in report.observations
    )


def test_run_phytools_parity_cases_records_failure_artifacts(
    tmp_path: Path,
) -> None:
    rscript = fake_phytools_rscript(
        tmp_path / "fake-phytools-rscript",
        summary_overrides={
            "phylosig-k-strong-signal-twenty-four-taxa": {
                "taxon_count": 24,
                "trait_name": "signal_strong",
                "k": 0.5,
                "p_value": 0.005025125628140704,
                "permutation_count": 199,
                "permutation_seed": 17,
                "simulated_k_minimum": 0.004227447597570447,
                "simulated_k_mean": 0.029614826253456215,
            }
        },
    )

    report = run_phytools_parity_cases(
        case_ids=["phylosig-k-strong-signal-twenty-four-taxa"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "failures",
    )

    assert report.all_passed is False
    assert report.failed_case_count == 1
    observation = report.observations[0]
    assert observation.status == "failed"
    assert observation.mismatch_reason == "summary_mismatch:k"
    assert observation.reproducible_artifact_root is not None
    assert (observation.reproducible_artifact_root / "reference-summary.json").exists()
    assert (observation.reproducible_artifact_root / "bijux-summary.json").exists()


def test_run_phytools_parity_cases_marks_missing_rscript_as_skipped(
    tmp_path: Path,
) -> None:
    report = run_phytools_parity_cases(
        case_ids=["phylosig-lambda-non-ultrametric-strong-signal-twenty-four-taxa"],
        rscript_executable=str(tmp_path / "missing-rscript"),
        failure_root=tmp_path / "failures",
    )

    assert report.all_passed is False
    assert report.skipped_case_count == 1
    observation = report.observations[0]
    assert observation.status == "skipped"
    assert observation.mismatch_reason == "rscript_unavailable"
    assert observation.reproducible_artifact_root is not None
