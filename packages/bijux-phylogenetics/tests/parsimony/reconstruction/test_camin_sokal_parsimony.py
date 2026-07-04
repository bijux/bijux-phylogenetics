from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    score_camin_sokal,
    score_fitch,
    write_camin_sokal_artifacts,
)
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_camin_sokal_surface() -> None:
    assert parsimony_api.score_camin_sokal is score_camin_sokal
    assert parsimony_api.write_camin_sokal_artifacts is write_camin_sokal_artifacts


def test_score_camin_sokal_matches_hand_computed_irreversible_fixture() -> None:
    report = score_camin_sokal(
        fixture("camin_sokal_tree_5_taxa.nwk"),
        fixture("camin_sokal_binary_matrix.tsv"),
    )

    assert report.root_state == "0"
    assert report.total_gains == 3
    assert [
        (
            row.character_id,
            row.derived_taxon_count,
            row.gain_count,
            row.root_state,
        )
        for row in report.step_rows
    ] == [
        ("char01_single_gain_clade", 2, 1, "0"),
        ("char02_repeated_gains", 2, 2, "0"),
        ("char03_absent", 0, 0, "0"),
    ]
    assert [
        (row.character_id, row.change_kind, row.node, row.descendant_taxa)
        for row in report.branch_change_rows
    ] == [
        ("char01_single_gain_clade", "gain", "D|E", ["D", "E"]),
        ("char02_repeated_gains", "gain", "B", ["B"]),
        ("char02_repeated_gains", "gain", "D", ["D"]),
    ]


def test_unordered_fitch_does_not_provide_irreversible_change_placement() -> None:
    fitch_report = score_fitch(
        fixture("camin_sokal_tree_5_taxa.nwk"),
        fixture("camin_sokal_unordered_gap_matrix.tsv"),
    )
    camin_sokal_report = score_camin_sokal(
        fixture("camin_sokal_tree_5_taxa.nwk"),
        fixture("camin_sokal_unordered_gap_matrix.tsv"),
    )

    root_rows = [
        row
        for row in fitch_report.node_state_rows
        if row.character_id == "char01_disjoint_presences" and row.node == "A|B|C|D|E"
    ]
    assert [row.state_set for row in root_rows] == [["0", "1"]]
    assert [
        (row.character_id, row.node, row.descendant_taxa)
        for row in camin_sokal_report.branch_change_rows
    ] == [
        ("char01_disjoint_presences", "A", ["A"]),
        ("char01_disjoint_presences", "C", ["C"]),
    ]


def test_score_camin_sokal_rejects_multistate_traits_without_binarization() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        score_camin_sokal(
            fixture("camin_sokal_tree_5_taxa.nwk"),
            fixture("camin_sokal_multistate_matrix.tsv"),
        )

    assert error_info.value.code == "parsimony_matrix_multistate_not_binarized"
    assert error_info.value.details["character_id"] == "char01_needs_binarization"
    assert error_info.value.details["invalid_states"] == ["2"]


def test_write_camin_sokal_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = score_camin_sokal(
        fixture("camin_sokal_tree_5_taxa.nwk"),
        fixture("camin_sokal_binary_matrix.tsv"),
    )

    outputs = write_camin_sokal_artifacts(tmp_path / "camin-sokal-run", report)

    assert set(outputs) == {"steps_path", "branch_changes_path", "run_json_path"}
    assert (
        outputs["steps_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "character_id\tderived_taxon_count\tgain_count\troot_state\tcharacter_weight\tweighted_score\n"
        )
    )
    assert (
        outputs["branch_changes_path"]
        .read_text(encoding="utf-8")
        .startswith("character_id\tchange_kind\tnode\tnode_name\tdescendant_taxa\n")
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "camin-sokal"
    assert payload["root_state"] == "0"
    assert payload["total_gains"] == 3
    assert payload["total_weighted_score"] == 3.0
