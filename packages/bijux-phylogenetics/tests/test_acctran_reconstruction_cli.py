from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_acctran_cli_writes_governed_artifacts(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "acctran-cli"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "acctran",
            str(fixture("acctran_tree_5_taxa.nwk")),
            str(fixture("acctran_ambiguous_matrix.tsv")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["algorithm"] == "acctran"
    assert payload["metrics"]["taxon_count"] == 5
    assert payload["metrics"]["character_count"] == 1
    assert payload["metrics"]["total_steps"] == 2
    assert (out_dir / "steps.tsv").is_file()
    assert (out_dir / "ancestral_states.tsv").is_file()
    assert (out_dir / "resolved_states.tsv").is_file()
    assert (out_dir / "branch_changes.tsv").is_file()
    assert (out_dir / "run.json").is_file()
    assert (
        (out_dir / "branch_changes.tsv")
        .read_text(encoding="utf-8")
        .startswith(
            "branch_id\tcharacter_id\tparent_node\tchild_node\tchild_node_name\tchild_descendant_taxa\tchange_from\tchange_to\tambiguous\n"
        )
    )


def test_phylo_parsimony_acctran_cli_preserves_ambiguous_branch_change_flag(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "acctran-branch-map-cli"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "acctran",
            str(fixture("branch_change_mapping_tree_4_taxa.nwk")),
            str(fixture("branch_change_mapping_multistate_matrix.tsv")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert "C|D\tchar01_branch_mapping\tA|B|C|D\tC|D" in (
        out_dir / "branch_changes.tsv"
    ).read_text(encoding="utf-8")
    assert "\t0\t1\ttrue\n" in (out_dir / "branch_changes.tsv").read_text(
        encoding="utf-8"
    )
