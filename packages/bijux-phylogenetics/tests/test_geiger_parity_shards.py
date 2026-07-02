from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.parity import (
    list_geiger_parity_shards,
    run_geiger_parity_cases,
)
from tests.support.fake_geiger_parity import fake_geiger_rscript


def test_list_geiger_parity_shards_groups_cases_by_function_and_model() -> None:
    shards = list_geiger_parity_shards()

    lambda_shard = next(
        shard
        for shard in shards
        if "fitcontinuous-lambda-weak-signal-review" in shard.case_ids
    )

    assert lambda_shard.shard_id == "geiger-fitcontinuous-model-lambda"
    assert lambda_shard.function_name == "geiger::fitContinuous(model='lambda')"
    assert lambda_shard.model_name == "lambda"
    assert lambda_shard.case_count >= 1


def test_run_geiger_parity_cases_supports_shard_selection(tmp_path: Path) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")
    lambda_shard = next(
        shard
        for shard in list_geiger_parity_shards()
        if "fitcontinuous-lambda-weak-signal-review" in shard.case_ids
    )

    report = run_geiger_parity_cases(
        shard_ids=[lambda_shard.shard_id],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-shard-failures",
    )

    assert report.case_count == lambda_shard.case_count
    assert {observation.case_id for observation in report.observations} == set(
        lambda_shard.case_ids
    )
    assert {observation.model_name for observation in report.observations} == {
        "lambda"
    }
    assert all(
        observation.function_name == "geiger::fitContinuous(model='lambda')"
        for observation in report.observations
    )
