from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main

pytestmark = pytest.mark.slow

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_phylo_likelihood_infer_tree_cli_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "likelihood-tree-inference"

    exit_code = main(
        [
            "phylo",
            "likelihood",
            "infer-tree",
            str(fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta")),
            "--start-tree-count",
            "3",
            "--start-tree-seed",
            "17",
            "--upper-branch-length-bound",
            "1.0",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["taxon_count"] == 4
    assert payload["metrics"]["site_count"] == 12
    assert payload["metrics"]["pattern_count"] == 2
    assert payload["metrics"]["start_tree_count"] == 3
    assert payload["metrics"]["best_run_source_label"] == "stepwise-addition-tree"
    assert payload["metrics"]["search_method"] == "nni"
    assert payload["metrics"]["selected_model_name"] in {
        "JC69",
        "K80",
        "F81",
        "HKY85",
        "GTR",
    }
    assert math.isclose(
        payload["metrics"]["best_final_log_likelihood"],
        -27.213282538844922,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert (out_dir / "stepwise_start_tree.nwk").is_file()
    assert (out_dir / "start_trees.nwk").is_file()
    assert (out_dir / "final_tree.nwk").is_file()
    assert (out_dir / "best_trees.nwk").is_file()
    assert (out_dir / "likelihood_table.tsv").is_file()
    assert (out_dir / "model_table.tsv").is_file()
    assert (out_dir / "search_trace.tsv").is_file()
    assert (out_dir / "run.json").is_file()
