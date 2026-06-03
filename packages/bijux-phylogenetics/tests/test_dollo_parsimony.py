from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import score_dollo, write_dollo_artifacts
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_dollo_surface() -> None:
    assert parsimony_api.score_dollo is score_dollo
    assert parsimony_api.write_dollo_artifacts is write_dollo_artifacts


def test_score_dollo_matches_hand_computed_binary_fixture() -> None:
    report = score_dollo(
        fixture("dollo_tree_5_taxa.nwk"),
        fixture("dollo_binary_matrix.tsv"),
    )

    assert report.total_gains == 2
    assert report.total_losses == 3
    assert [
        (
            row.character_id,
            row.derived_taxon_count,
            row.gain_node,
            row.total_losses,
            row.impossible_state_warning,
        )
        for row in report.step_rows
    ] == [
        ("char01_gain_only", 2, "D|E", 0, None),
        ("char02_root_gain_with_losses", 2, "A|B|C|D|E", 3, None),
        ("char03_absent", 0, None, 0, None),
    ]
    assert [
        (row.character_id, row.change_kind, row.node, row.descendant_taxa)
        for row in report.branch_change_rows
    ] == [
        ("char01_gain_only", "gain", "D|E", ["D", "E"]),
        (
            "char02_root_gain_with_losses",
            "gain",
            "A|B|C|D|E",
            ["A", "B", "C", "D", "E"],
        ),
        ("char02_root_gain_with_losses", "loss", "A", ["A"]),
        ("char02_root_gain_with_losses", "loss", "C", ["C"]),
        ("char02_root_gain_with_losses", "loss", "E", ["E"]),
    ]


def test_score_dollo_rejects_multistate_traits_without_binarization() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        score_dollo(
            fixture("dollo_tree_5_taxa.nwk"),
            fixture("dollo_multistate_matrix.tsv"),
        )

    assert error_info.value.code == "parsimony_matrix_multistate_not_binarized"
    assert error_info.value.details["character_id"] == "char01_needs_binarization"
    assert error_info.value.details["invalid_states"] == ["2"]


def test_write_dollo_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = score_dollo(
        fixture("dollo_tree_5_taxa.nwk"),
        fixture("dollo_binary_matrix.tsv"),
    )

    outputs = write_dollo_artifacts(tmp_path / "dollo-run", report)

    assert set(outputs) == {"steps_path", "branch_changes_path", "run_json_path"}
    assert (
        outputs["steps_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "character_id\tderived_taxon_count\tgain_node\tgain_node_name\tgain_descendant_taxa\ttotal_losses\timpossible_state_warning\tstep_count\tcharacter_weight\tweighted_score\n"
        )
    )
    assert (
        outputs["branch_changes_path"]
        .read_text(encoding="utf-8")
        .startswith("character_id\tchange_kind\tnode\tnode_name\tdescendant_taxa\n")
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "dollo"
    assert payload["total_gains"] == 2
    assert payload["total_losses"] == 3
    assert payload["total_weighted_score"] == 5.0
