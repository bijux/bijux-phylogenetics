from __future__ import annotations

import json
from pathlib import Path

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    reconstruct_acctran,
    write_parsimony_reconstruction_artifacts,
)
from bijux_phylogenetics.parsimony.reconstruction import _reconstruct_deltran

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_acctran_surface() -> None:
    assert parsimony_api.reconstruct_acctran is reconstruct_acctran
    assert (
        parsimony_api.write_parsimony_reconstruction_artifacts
        is write_parsimony_reconstruction_artifacts
    )


def test_reconstruct_acctran_resolves_ambiguous_changes_toward_earlier_branches() -> (
    None
):
    acctran_report = reconstruct_acctran(
        fixture("acctran_tree_5_taxa.nwk"),
        fixture("acctran_ambiguous_matrix.tsv"),
    )
    deltran_report = _reconstruct_deltran(
        fixture("acctran_tree_5_taxa.nwk"),
        fixture("acctran_ambiguous_matrix.tsv"),
    )

    assert acctran_report.total_steps == 2
    assert deltran_report.total_steps == 2
    assert [
        (row.character_id, row.step_count, row.root_state)
        for row in acctran_report.step_rows
    ] == [
        ("char01_ambiguous", 2, "0"),
    ]
    assert [
        (
            row.branch_id,
            row.parent_node,
            row.child_node,
            row.character_id,
            row.change_from,
            row.change_to,
            row.ambiguous,
        )
        for row in acctran_report.branch_change_rows
    ] == [
        ("C|D|E", "A|B|C|D|E", "C|D|E", "char01_ambiguous", "0", "1", False),
        ("D", "D|E", "D", "char01_ambiguous", "1", "0", False),
    ]
    assert [
        (
            row.branch_id,
            row.parent_node,
            row.child_node,
            row.character_id,
            row.change_from,
            row.change_to,
            row.ambiguous,
        )
        for row in deltran_report.branch_change_rows
    ] == [
        ("C", "C|D|E", "C", "char01_ambiguous", "0", "1", False),
        ("E", "D|E", "E", "char01_ambiguous", "0", "1", False),
    ]
    resolved_states = {
        (row.character_id, row.node): row.resolved_state
        for row in acctran_report.node_state_rows
    }
    assert resolved_states[("char01_ambiguous", "A|B|C|D|E")] == "0"
    assert resolved_states[("char01_ambiguous", "C|D|E")] == "1"
    assert resolved_states[("char01_ambiguous", "D|E")] == "1"
    assert resolved_states[("char01_ambiguous", "D")] == "0"
    assert resolved_states[("char01_ambiguous", "E")] == "1"
    ancestral_rows = {
        (row.character_id, row.node_id): (
            row.clade_id,
            row.possible_states,
            row.chosen_state,
            row.method,
            row.ambiguous,
        )
        for row in acctran_report.ancestral_state_rows
    }
    assert ancestral_rows[("char01_ambiguous", "A|B|C|D|E")] == (
        "A|B|C|D|E",
        ["0", "1"],
        "0",
        "acctran",
        True,
    )
    assert ancestral_rows[("char01_ambiguous", "A|B")] == (
        "A|B",
        ["0"],
        "0",
        "acctran",
        False,
    )
    assert ancestral_rows[("char01_ambiguous", "C|D|E")] == (
        "C|D|E",
        ["1"],
        "1",
        "acctran",
        False,
    )
    assert ancestral_rows[("char01_ambiguous", "D|E")] == (
        "D|E",
        ["0", "1"],
        "1",
        "acctran",
        True,
    )


def test_write_acctran_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = reconstruct_acctran(
        fixture("acctran_tree_5_taxa.nwk"),
        fixture("acctran_ambiguous_matrix.tsv"),
    )

    outputs = write_parsimony_reconstruction_artifacts(tmp_path / "acctran-run", report)

    assert set(outputs) == {
        "steps_path",
        "node_states_path",
        "ancestral_states_path",
        "branch_changes_path",
        "run_json_path",
    }
    assert (
        outputs["steps_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "character_id\tstep_count\tobserved_states\troot_state\tcharacter_weight\tweighted_score\n"
        )
    )
    assert (
        outputs["node_states_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "character_id\tnode\tnode_name\tdescendant_taxa\tresolved_state\tis_tip\tobserved_state\n"
        )
    )
    assert (
        outputs["ancestral_states_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "node_id\tclade_id\tcharacter_id\tpossible_states\tchosen_state\tmethod\tambiguous\n"
        )
    )
    assert (
        outputs["branch_changes_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "branch_id\tcharacter_id\tparent_node\tchild_node\tchild_node_name\tchild_descendant_taxa\tchange_from\tchange_to\tambiguous\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "acctran"
    assert payload["total_steps"] == 2
    assert payload["total_weighted_score"] == 2.0
    assert payload["branch_change_rows"] == [
        {
            "ambiguous": False,
            "branch_id": "C|D|E",
            "character_id": "char01_ambiguous",
            "change_from": "0",
            "change_to": "1",
            "child_descendant_taxa": ["C", "D", "E"],
            "child_node": "C|D|E",
            "child_node_name": None,
            "parent_node": "A|B|C|D|E",
            "parent_state": "0",
        },
        {
            "ambiguous": False,
            "branch_id": "D",
            "character_id": "char01_ambiguous",
            "change_from": "1",
            "change_to": "0",
            "child_descendant_taxa": ["D"],
            "child_node": "D",
            "child_node_name": "D",
            "parent_node": "D|E",
            "parent_state": "1",
        },
    ]
    assert payload["ancestral_state_rows"] == [
        {
            "ambiguous": True,
            "character_id": "char01_ambiguous",
            "chosen_state": "0",
            "clade_id": "A|B|C|D|E",
            "method": "acctran",
            "node_id": "A|B|C|D|E",
            "possible_states": ["0", "1"],
        },
        {
            "ambiguous": False,
            "character_id": "char01_ambiguous",
            "chosen_state": "0",
            "clade_id": "A|B",
            "method": "acctran",
            "node_id": "A|B",
            "possible_states": ["0"],
        },
        {
            "ambiguous": False,
            "character_id": "char01_ambiguous",
            "chosen_state": "1",
            "clade_id": "C|D|E",
            "method": "acctran",
            "node_id": "C|D|E",
            "possible_states": ["1"],
        },
        {
            "ambiguous": True,
            "character_id": "char01_ambiguous",
            "chosen_state": "1",
            "clade_id": "D|E",
            "method": "acctran",
            "node_id": "D|E",
            "possible_states": ["0", "1"],
        },
    ]


def test_reconstruct_acctran_marks_ambiguous_branch_change_mapping() -> None:
    report = reconstruct_acctran(
        fixture("branch_change_mapping_tree_4_taxa.nwk"),
        fixture("branch_change_mapping_multistate_matrix.tsv"),
    )

    assert [
        (
            row.branch_id,
            row.parent_node,
            row.child_node,
            row.change_from,
            row.change_to,
            row.ambiguous,
        )
        for row in report.branch_change_rows
    ] == [
        ("C|D", "A|B|C|D", "C|D", "0", "1", True),
        ("D", "C|D", "D", "1", "2", False),
    ]
