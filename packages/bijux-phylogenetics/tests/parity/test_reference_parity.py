from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from bijux_phylogenetics.parity import (
    validate_reference_parity_examples,
    write_reference_parity_observation_table,
    write_reference_parity_summary_table,
)


def test_validate_reference_parity_examples_passes() -> None:
    report = validate_reference_parity_examples()
    assert report.all_passed is True
    assert report.case_count == 14
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
    pgls = next(
        item
        for item in report.observations
        if item.case == "pgls-example-tree-brownian"
    )
    pgls_categorical = next(
        item
        for item in report.observations
        if item.case == "pgls-example-tree-brownian-categorical"
    )
    pgls_interaction = next(
        item
        for item in report.observations
        if item.case == "pgls-example-tree-brownian-interaction"
    )
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
    overlap_rf = next(
        item
        for item in report.observations
        if item.case == "robinson-foulds-prune-to-shared-taxa"
    )
    overlap_branch_score = next(
        item
        for item in report.observations
        if item.case == "branch-score-distance-prune-to-shared-taxa"
    )
    assert pgls.input_fixtures[0].name == "example_tree.nwk"
    assert pgls.reference_tool == "ape+nlme"
    assert "floating-point linear algebra" in pgls.tolerance_reason
    assert pgls.expected_failure_mode == "model_assumption"
    assert pgls.taxon_overlap_policy is None
    assert math.isclose(
        pgls.observed_output["aic"],
        pgls.expected_output["aic"],
        rel_tol=1e-6,
        abs_tol=1e-6,
    )
    assert math.isclose(
        pgls.observed_output["coefficient.predictor_one.p_value"],
        pgls.expected_output["coefficient.predictor_one.p_value"],
        rel_tol=1e-6,
        abs_tol=1e-6,
    )
    assert math.isclose(
        pgls_categorical.observed_output["coefficient.habitat[tundra].estimate"],
        pgls_categorical.expected_output["coefficient.habitat[tundra].estimate"],
        rel_tol=1e-6,
        abs_tol=1e-6,
    )
    assert pgls_categorical.observed_output["model_matrix.G.diet[herbivore]"] == 1.0
    assert (
        "stats::model.matrix treatment-coded categorical predictors"
        in pgls_categorical.reference_source
    )
    assert math.isclose(
        pgls_interaction.observed_output[
            "coefficient.habitat[tundra]:diet[herbivore].estimate"
        ],
        pgls_interaction.expected_output[
            "coefficient.habitat[tundra]:diet[herbivore].estimate"
        ],
        rel_tol=1e-6,
        abs_tol=1e-6,
    )
    assert (
        pgls_interaction.observed_output[
            "model_matrix.H.habitat[tundra]:diet[herbivore]"
        ]
        == 1.0
    )
    assert rf.expected_output["robinson_foulds_distance"] == 2
    assert rf.expected_output["normalized_robinson_foulds"] == 1.0
    assert rf.expected_failure_mode == "topology"
    assert rf.shared_taxa == ["A", "B", "C", "D"]
    assert rf.left_only_taxa == []
    assert rf.right_only_taxa == []
    assert overlap_rf.expected_failure_mode == "missing_taxa_policy"
    assert overlap_rf.taxon_overlap_policy == "prune-to-shared"
    assert overlap_rf.shared_taxa == ["A", "B", "C"]
    assert overlap_rf.left_only_taxa == ["D"]
    assert overlap_rf.right_only_taxa == ["E"]
    assert overlap_rf.observed_output["robinson_foulds_distance"] == 0
    assert overlap_branch_score.expected_failure_mode == "missing_taxa_policy"
    assert overlap_branch_score.observed_output["branch_score_distance"] == 0.0
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
    assert any(row["case"] == "pgls-example-tree-brownian-interaction" for row in rows)
