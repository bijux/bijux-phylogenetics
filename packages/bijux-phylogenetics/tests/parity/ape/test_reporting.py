from __future__ import annotations

import csv
from pathlib import Path

import pytest

from bijux_phylogenetics.parity import (
    run_ape_parity_cases,
    write_ape_parity_observation_table,
    write_ape_parity_summary_table,
)
from tests.support.fake_reference_parity import fake_ape_rscript


@pytest.mark.slow
def test_write_ape_parity_tables_writes_summary_and_observations(
    tmp_path: Path,
) -> None:
    rscript = fake_ape_rscript(tmp_path / "fake-ape-rscript")
    report = run_ape_parity_cases(
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-parity-failures",
    )
    summary_path = tmp_path / "ape-parity-summary.tsv"
    observation_path = tmp_path / "ape-parity-observations.tsv"

    write_ape_parity_summary_table(summary_path, report)
    write_ape_parity_observation_table(observation_path, report)

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0] == (
        "function_name\tcase_count\tpassed_case_count\tfailed_case_count\tskipped_case_count"
    )
    with observation_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 180
    assert rows[0]["function_name"] == "ape::read.tree"
    assert rows[0]["fixture_kind"] == "tree"
    assert rows[0]["fixture_id"]
    assert rows[0]["status"] == "passed"
    assert rows[0]["bijux_version"]
