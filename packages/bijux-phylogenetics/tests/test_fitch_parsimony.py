from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    load_fitch_character_matrix,
    score_fitch,
    write_fitch_artifacts,
)
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_root_exports_parsimony_gateway() -> None:
    assert parsimony_api.score_fitch is score_fitch
    assert parsimony_api.load_fitch_character_matrix is load_fitch_character_matrix
    assert parsimony_api.write_fitch_artifacts is write_fitch_artifacts


def test_score_fitch_matches_hand_computed_binary_fixture() -> None:
    report = score_fitch(
        fixture("fitch_tree.nwk"),
        fixture("fitch_binary_matrix.tsv"),
    )

    assert report.total_steps == 2
    assert [(row.character_id, row.step_count) for row in report.step_rows] == [
        ("char01_split", 1),
        ("char02_left_conflict", 1),
    ]
    rows = {
        (row.character_id, row.node): row.state_set for row in report.node_state_rows
    }
    assert rows[("char01_split", "A")] == ["0"]
    assert rows[("char01_split", "B")] == ["0"]
    assert rows[("char01_split", "C")] == ["1"]
    assert rows[("char01_split", "D")] == ["1"]
    assert rows[("char01_split", "A|B")] == ["0"]
    assert rows[("char01_split", "C|D")] == ["1"]
    assert rows[("char01_split", "A|B|C|D")] == ["0", "1"]
    assert rows[("char02_left_conflict", "A")] == ["0"]
    assert rows[("char02_left_conflict", "B")] == ["1"]
    assert rows[("char02_left_conflict", "C")] == ["1"]
    assert rows[("char02_left_conflict", "D")] == ["1"]
    assert rows[("char02_left_conflict", "A|B")] == ["0", "1"]
    assert rows[("char02_left_conflict", "C|D")] == ["1"]
    assert rows[("char02_left_conflict", "A|B|C|D")] == ["1"]


def test_score_fitch_matches_hand_computed_four_state_fixture() -> None:
    report = score_fitch(
        fixture("fitch_tree.nwk"),
        fixture("fitch_four_state_matrix.tsv"),
    )

    assert report.total_steps == 5
    assert [(row.character_id, row.step_count) for row in report.step_rows] == [
        ("char03_four_way", 3),
        ("char04_mixed_overlap", 2),
    ]
    rows = {
        (row.character_id, row.node): row.state_set for row in report.node_state_rows
    }
    assert rows[("char03_four_way", "A|B")] == ["blue", "red"]
    assert rows[("char03_four_way", "C|D")] == ["green", "yellow"]
    assert rows[("char03_four_way", "A|B|C|D")] == [
        "blue",
        "green",
        "red",
        "yellow",
    ]
    assert rows[("char04_mixed_overlap", "A|B")] == ["red"]
    assert rows[("char04_mixed_overlap", "C|D")] == ["blue", "green"]
    assert rows[("char04_mixed_overlap", "A|B|C|D")] == ["blue", "green", "red"]


def test_load_fitch_character_matrix_rejects_empty_character_matrices() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        load_fitch_character_matrix(fixture("fitch_empty_matrix.tsv"))

    assert error_info.value.code == "parsimony_matrix_empty"
    assert error_info.value.details["character_count"] == 0


def test_load_fitch_character_matrix_honors_explicit_taxon_column() -> None:
    matrix = load_fitch_character_matrix(
        fixture("fitch_binary_species_matrix.tsv"),
        taxon_column="species",
    )

    assert matrix.taxon_column == "species"
    assert matrix.character_ids == ["char01_split"]
    assert matrix.states_by_taxon["A"]["char01_split"] == "0"


def test_load_fitch_character_matrix_rejects_unknown_states() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        load_fitch_character_matrix(fixture("fitch_unknown_state_matrix.tsv"))

    assert error_info.value.code == "parsimony_matrix_unknown_state"
    assert error_info.value.details["taxon"] == "C"
    assert error_info.value.details["character_id"] == "char01_split"


def test_score_fitch_rejects_missing_tree_taxa() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        score_fitch(
            fixture("fitch_tree.nwk"),
            fixture("fitch_missing_taxon_matrix.tsv"),
        )

    assert error_info.value.code == "parsimony_matrix_missing_taxa"
    assert error_info.value.details["missing_taxa"] == ["D"]


def test_write_fitch_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = score_fitch(
        fixture("fitch_tree.nwk"),
        fixture("fitch_binary_matrix.tsv"),
    )

    outputs = write_fitch_artifacts(tmp_path / "fitch-run", report)

    assert set(outputs) == {"steps_path", "node_state_sets_path", "run_json_path"}
    assert (
        outputs["steps_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "character_id\tstep_count\tobserved_states\tcharacter_weight\tweighted_score\n"
        )
    )
    assert (
        outputs["node_state_sets_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "character_id\tnode\tnode_name\tdescendant_taxa\tis_tip\tobserved_state\tstate_set\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "unordered-fitch"
    assert payload["total_steps"] == 2
    assert payload["total_weighted_score"] == 2.0
