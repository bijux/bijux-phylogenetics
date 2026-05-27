from __future__ import annotations

import json
from pathlib import Path

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    reconstruct_acctran,
    reconstruct_deltran,
    write_parsimony_reconstruction_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_deltran_surface() -> None:
    assert parsimony_api.reconstruct_deltran is reconstruct_deltran
    assert (
        parsimony_api.write_parsimony_reconstruction_artifacts
        is write_parsimony_reconstruction_artifacts
    )


def test_reconstruct_deltran_delays_ambiguous_changes_toward_terminal_branches() -> None:
    deltran_report = reconstruct_deltran(
        fixture("acctran_tree_5_taxa.nwk"),
        fixture("acctran_ambiguous_matrix.tsv"),
    )
    acctran_report = reconstruct_acctran(
        fixture("acctran_tree_5_taxa.nwk"),
        fixture("acctran_ambiguous_matrix.tsv"),
    )

    assert deltran_report.total_steps == 2
    assert acctran_report.total_steps == 2
    assert [(row.character_id, row.step_count, row.root_state) for row in deltran_report.step_rows] == [
        ("char01_ambiguous", 2, "0"),
    ]
    assert [
        (row.character_id, row.node, row.change_from, row.change_to)
        for row in deltran_report.branch_change_rows
    ] == [
        ("char01_ambiguous", "C", "0", "1"),
        ("char01_ambiguous", "E", "0", "1"),
    ]
    assert [
        (row.character_id, row.node, row.change_from, row.change_to)
        for row in acctran_report.branch_change_rows
    ] == [
        ("char01_ambiguous", "C|D|E", "0", "1"),
        ("char01_ambiguous", "D", "1", "0"),
    ]
    resolved_states = {
        (row.character_id, row.node): row.resolved_state
        for row in deltran_report.node_state_rows
    }
    assert resolved_states[("char01_ambiguous", "A|B|C|D|E")] == "0"
    assert resolved_states[("char01_ambiguous", "C|D|E")] == "0"
    assert resolved_states[("char01_ambiguous", "D|E")] == "0"
    assert resolved_states[("char01_ambiguous", "C")] == "1"
    assert resolved_states[("char01_ambiguous", "E")] == "1"


def test_write_deltran_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = reconstruct_deltran(
        fixture("acctran_tree_5_taxa.nwk"),
        fixture("acctran_ambiguous_matrix.tsv"),
    )

    outputs = write_parsimony_reconstruction_artifacts(tmp_path / "deltran-run", report)

    assert set(outputs) == {
        "steps_path",
        "node_states_path",
        "branch_changes_path",
        "run_json_path",
    }
    assert outputs["steps_path"].read_text(encoding="utf-8").startswith(
        "character_id\tstep_count\tobserved_states\troot_state\n"
    )
    assert outputs["node_states_path"].read_text(encoding="utf-8").startswith(
        "character_id\tnode\tnode_name\tdescendant_taxa\tresolved_state\tis_tip\tobserved_state\n"
    )
    assert outputs["branch_changes_path"].read_text(encoding="utf-8").startswith(
        "character_id\tparent_node\tparent_state\tnode\tnode_name\tdescendant_taxa\tchange_from\tchange_to\n"
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "deltran"
    assert payload["total_steps"] == 2
