from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from bijux_phylogenetics.parity import (
    run_phytools_parity_cases,
    write_phytools_parity_observation_table,
    write_phytools_parity_summary_table,
)
from tests.support.fake_phytools_parity import fake_phytools_rscript


@pytest.mark.slow
def test_write_phytools_parity_tables_writes_summary_and_observations(
    tmp_path: Path,
) -> None:
    rscript = fake_phytools_rscript(tmp_path / "fake-phytools-rscript")
    report = run_phytools_parity_cases(
        shard_ids=[
            "phytools-make-simmap-model-er",
            "phytools-sim-history",
        ],
        rscript_executable=str(rscript),
    )
    summary_path = tmp_path / "phytools-parity-summary.tsv"
    observation_path = tmp_path / "phytools-parity-observations.tsv"

    write_phytools_parity_summary_table(summary_path, report)
    write_phytools_parity_observation_table(observation_path, report)

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("function_name\tcase_count")
    with observation_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    assert len(rows) == 8
    input_fixtures = json.loads(rows[0]["input_fixtures"])
    assert isinstance(input_fixtures, list)
    assert rows[0]["phytools_version"] == "2.5.2"
