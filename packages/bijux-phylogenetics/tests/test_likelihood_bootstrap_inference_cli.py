from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"
pytestmark = pytest.mark.slow


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_phylo_likelihood_bootstrap_tree_cli_writes_governed_seeded_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    from bijux_phylogenetics.command_line import main

    left_out_dir = tmp_path / "likelihood-bootstrap-left"
    right_out_dir = tmp_path / "likelihood-bootstrap-right"
    command = [
        "phylo",
        "likelihood",
        "bootstrap-tree",
        str(fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta")),
        "--model",
        "jc69",
        "--search-method",
        "nni",
        "--start-tree-count",
        "2",
        "--start-tree-seed",
        "5",
        "--replicate-count",
        "4",
        "--bootstrap-seed",
        "9",
        "--upper-branch-length-bound",
        "1.0",
        "--out-dir",
    ]

    left_exit_code = main([*command, str(left_out_dir), "--json"])
    left_payload = json.loads(capsys.readouterr().out)
    right_exit_code = main([*command, str(right_out_dir), "--json"])
    right_payload = json.loads(capsys.readouterr().out)

    assert left_exit_code == 0
    assert right_exit_code == 0
    assert left_payload["status"] == "ok"
    assert right_payload["status"] == "ok"
    assert left_payload["metrics"]["selected_reference_model_name"] == "JC69"
    assert left_payload["metrics"]["taxon_count"] == 4
    assert left_payload["metrics"]["site_count"] == 12
    assert left_payload["metrics"]["pattern_count"] == 2
    assert left_payload["metrics"]["replicate_count"] == 4
    assert left_payload["metrics"]["support_row_count"] == 2
    assert left_payload["metrics"]["search_method"] == "nni"
    assert math.isclose(
        left_payload["metrics"]["reference_log_likelihood"],
        -34.13524969797671,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert right_payload["metrics"] == left_payload["metrics"]
    assert (left_out_dir / "final_tree.nwk").is_file()
    assert (left_out_dir / "replicate_trees.nwk").is_file()
    assert (left_out_dir / "replicate_draws.tsv").is_file()
    assert (left_out_dir / "clade_support.tsv").is_file()
    assert (left_out_dir / "consensus_tree.nwk").is_file()
    assert (left_out_dir / "clade_frequencies.tsv").is_file()
    assert (left_out_dir / "run.json").is_file()
    assert (left_out_dir / "replicate_trees.nwk").read_text(encoding="utf-8") == (
        right_out_dir / "replicate_trees.nwk"
    ).read_text(encoding="utf-8")
    assert (left_out_dir / "clade_support.tsv").read_text(encoding="utf-8") == (
        right_out_dir / "clade_support.tsv"
    ).read_text(encoding="utf-8")
