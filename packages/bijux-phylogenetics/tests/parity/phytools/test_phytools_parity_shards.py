from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.parity import (
    list_phytools_parity_shards,
    run_phytools_parity_cases,
)
from tests.support.fake_phytools_parity import fake_phytools_rscript

pytestmark = pytest.mark.usefixtures("enable_live_parity_access")


def test_list_phytools_parity_shards_groups_cases_by_reference_call() -> None:
    shards = list_phytools_parity_shards()

    simmap_shard = next(
        shard
        for shard in shards
        if "simmap-er-binary-twenty-four-taxa" in shard.case_ids
    )

    assert simmap_shard.shard_id == "phytools-make-simmap-model-er"
    assert simmap_shard.function_name == "phytools::make.simmap(model='ER')"
    assert simmap_shard.case_count >= 1


def test_run_phytools_parity_cases_supports_shard_selection(tmp_path: Path) -> None:
    rscript = fake_phytools_rscript(tmp_path / "fake-phytools-rscript")
    simmap_shard = next(
        shard
        for shard in list_phytools_parity_shards()
        if "simmap-er-binary-twenty-four-taxa" in shard.case_ids
    )

    report = run_phytools_parity_cases(
        shard_ids=[simmap_shard.shard_id],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "phytools-shard-failures",
    )

    assert report.all_passed is True
    assert report.case_count == simmap_shard.case_count
    assert {observation.case_id for observation in report.observations} == set(
        simmap_shard.case_ids
    )
    assert {observation.function_name for observation in report.observations} == {
        "phytools::make.simmap(model='ER')"
    }
