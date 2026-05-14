from __future__ import annotations

import csv
import json
from pathlib import Path

from bijux_phylogenetics.reference_parity import (
    validate_reference_parity_examples,
    write_reference_parity_observation_table,
    write_reference_parity_summary_table,
)


def test_validate_reference_parity_examples_passes() -> None:
    report = validate_reference_parity_examples()
    assert report.all_passed is True
    assert report.case_count == 10
    assert report.failed_case_count == 0
    assert report.covered_methods == [
        "blombergs-k",
        "branch-score-distance",
        "brownian-trait-model",
        "consensus-tree-generation",
        "ornstein-uhlenbeck-trait-model",
        "pagels-lambda",
        "pgls",
        "phylogenetic-independent-contrasts",
        "posterior-clade-frequencies",
        "robinson-foulds-distance",
    ]
    assert report.reference_tools["DendroPy"] == "5.0.8"
    assert report.reference_tools["phytools"] == "2.5.2"


def test_validate_reference_parity_examples_records_failure_modes_and_inputs() -> None:
    report = validate_reference_parity_examples()
    pgls = next(item for item in report.observations if item.method == "pgls")
    rf = next(
        item
        for item in report.observations
        if item.method == "robinson-foulds-distance"
    )
    consensus = next(
        item
        for item in report.observations
        if item.method == "consensus-tree-generation"
    )
    assert pgls.input_fixtures[0].name == "example_tree.nwk"
    assert pgls.reference_tool == "ape+nlme"
    assert "floating-point linear algebra" in pgls.tolerance_reason
    assert pgls.expected_failure_mode == "model_assumption"
    assert pgls.taxon_overlap_policy is None
    assert rf.expected_output["robinson_foulds_distance"] == 2
    assert rf.expected_output["normalized_robinson_foulds"] == 1.0
    assert rf.expected_failure_mode == "topology"
    assert rf.shared_taxa == ["A", "B", "C", "D"]
    assert rf.left_only_taxa == []
    assert rf.right_only_taxa == []
    assert consensus.observed_output["unrooted_robinson_foulds"] == 0
    assert consensus.observed_output["consensus_splits"] == [
        "C|D||A|B|E|F",
        "E|F||A|B|C|D",
    ]


def test_write_reference_parity_tables_writes_summary_and_observations(
    tmp_path: Path,
) -> None:
    report = validate_reference_parity_examples()
    summary_path = tmp_path / "reference-parity-summary.tsv"
    observation_path = tmp_path / "reference-parity-observations.tsv"
    write_reference_parity_summary_table(summary_path, report)
    write_reference_parity_observation_table(observation_path, report)
    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("suite\tmethod\tcase_count")
    assert any("robinson-foulds-distance" in row for row in summary_rows[1:])
    with observation_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    assert rows[0]["expected_failure_mode"] in {
        "topology",
        "branch_length",
        "model_assumption",
        "numerical_tolerance",
    }
    parsed_payload = json.loads(rows[0]["observed_output"])
    assert isinstance(parsed_payload, dict)
