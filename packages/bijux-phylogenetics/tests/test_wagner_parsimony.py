from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    load_parsimony_character_matrix,
    score_wagner,
    write_wagner_artifacts,
)
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_wagner_surface() -> None:
    assert (
        parsimony_api.load_parsimony_character_matrix is load_parsimony_character_matrix
    )
    assert parsimony_api.score_wagner is score_wagner
    assert parsimony_api.write_wagner_artifacts is write_wagner_artifacts


def test_score_wagner_matches_hand_computed_ordinal_fixture() -> None:
    report = score_wagner(
        fixture("fitch_tree.nwk"),
        fixture("wagner_ordinal_matrix.tsv"),
    )

    assert report.total_cost == 6
    assert [
        (row.character_id, row.weighted_step_count, row.optimal_root_states)
        for row in report.step_rows
    ] == [
        ("char01_gradient", 3, ["1", "2"]),
        ("char02_distance_gap", 3, ["3"]),
    ]
    costs = {
        (row.character_id, row.node, row.state): (row.cost, row.is_optimal_state)
        for row in report.node_cost_rows
    }
    assert costs[("char01_gradient", "A|B", "0")] == (1, True)
    assert costs[("char01_gradient", "A|B", "1")] == (1, True)
    assert costs[("char01_gradient", "A|B", "2")] == (3, False)
    assert costs[("char01_gradient", "A|B", "3")] == (5, False)
    assert costs[("char01_gradient", "C|D", "0")] == (5, False)
    assert costs[("char01_gradient", "C|D", "1")] == (3, False)
    assert costs[("char01_gradient", "C|D", "2")] == (1, True)
    assert costs[("char01_gradient", "C|D", "3")] == (1, True)
    assert costs[("char01_gradient", "A|B|C|D", "0")] == (4, False)
    assert costs[("char01_gradient", "A|B|C|D", "1")] == (3, True)
    assert costs[("char01_gradient", "A|B|C|D", "2")] == (3, True)
    assert costs[("char01_gradient", "A|B|C|D", "3")] == (4, False)
    assert costs[("char02_distance_gap", "A|B", "0")] == (3, True)
    assert costs[("char02_distance_gap", "A|B", "3")] == (3, True)
    assert costs[("char02_distance_gap", "C|D", "0")] == (6, False)
    assert costs[("char02_distance_gap", "C|D", "3")] == (0, True)
    assert costs[("char02_distance_gap", "A|B|C|D", "0")] == (6, False)
    assert costs[("char02_distance_gap", "A|B|C|D", "3")] == (3, True)


def test_score_wagner_accepts_named_states_when_explicit_order_is_supplied() -> None:
    report = score_wagner(
        fixture("fitch_tree.nwk"),
        fixture("wagner_named_state_matrix.tsv"),
        state_order=["low", "medium", "high", "very_high"],
    )

    assert report.total_cost == 6
    assert [row.state_order for row in report.step_rows] == [
        ["low", "medium", "high", "very_high"],
        ["low", "medium", "high", "very_high"],
    ]


def test_score_wagner_rejects_nonordinal_labels_without_explicit_order() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        score_wagner(
            fixture("fitch_tree.nwk"),
            fixture("wagner_named_state_matrix.tsv"),
        )

    assert error_info.value.code == "parsimony_state_order_required"
    assert error_info.value.details["character_id"] == "char01_gradient"
    assert error_info.value.details["observed_states"] == [
        "high",
        "low",
        "medium",
        "very_high",
    ]


def test_write_wagner_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = score_wagner(
        fixture("fitch_tree.nwk"),
        fixture("wagner_ordinal_matrix.tsv"),
    )

    outputs = write_wagner_artifacts(tmp_path / "wagner-run", report)

    assert set(outputs) == {"steps_path", "node_costs_path", "run_json_path"}
    assert (
        outputs["steps_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "character_id\tweighted_step_count\tobserved_states\tstate_order\toptimal_root_states\tcharacter_weight\tweighted_score\n"
        )
    )
    assert (
        outputs["node_costs_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "character_id\tnode\tnode_name\tdescendant_taxa\tstate\tcost\tis_optimal_state\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "ordered-wagner"
    assert payload["total_cost"] == 6
    assert payload["total_weighted_score"] == 6.0
