from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.parity import (
    list_ape_parity_shards,
    run_ape_parity_cases,
)
from tests.support.fake_reference_parity import fake_ape_rscript

pytestmark = pytest.mark.usefixtures("enable_live_parity_access")


def test_list_ape_parity_shards_groups_cases_by_owned_function_surface() -> None:
    shards = list_ape_parity_shards()

    mrca_shard = next(
        shard for shard in shards if "get-mrca-balanced-two-tip" in shard.case_ids
    )

    assert mrca_shard.shard_id == "find_tree_mrca"
    assert mrca_shard.function_name == "ape::getMRCA"
    assert mrca_shard.case_count >= 1


def test_run_ape_parity_cases_supports_shard_selection(tmp_path: Path) -> None:
    rscript = fake_ape_rscript(tmp_path / "fake-ape-rscript")
    mrca_shard = next(
        shard
        for shard in list_ape_parity_shards()
        if "get-mrca-balanced-two-tip" in shard.case_ids
    )

    report = run_ape_parity_cases(
        shard_ids=[mrca_shard.shard_id],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-shard-failures",
    )

    assert report.all_passed is True
    assert report.case_count == mrca_shard.case_count
    assert {observation.case_id for observation in report.observations} == set(
        mrca_shard.case_ids
    )
    assert {
        observation.python_function_name for observation in report.observations
    } == {"find_tree_mrca"}
