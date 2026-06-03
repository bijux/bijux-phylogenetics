from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.io.newick import load_newick_tree_set

pytestmark = pytest.mark.slow

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_phylo_likelihood_placement_cli_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "likelihood-placement"

    exit_code = main(
        [
            "phylo",
            "likelihood",
            "placement",
            str(fixture("trees", "likelihood_placement_reference_tree_4_taxa.nwk")),
            str(
                fixture(
                    "alignments",
                    "likelihood_placement_reference_alignment_4_taxa.fasta",
                )
            ),
            str(
                fixture(
                    "alignments",
                    "likelihood_placement_query_alignment_2_taxa.fasta",
                )
            ),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["model_name"] == "JC69"
    assert payload["metrics"]["reference_taxon_count"] == 4
    assert payload["metrics"]["edge_count"] == 6
    assert payload["metrics"]["query_count"] == 2
    assert payload["metrics"]["site_count"] == 12
    assert payload["metrics"]["placement_count"] == 12
    assert payload["metrics"]["total_function_evaluation_count"] > 0
    assert len(payload["outputs"]) == 4
    assert (out_dir / "summary.tsv").is_file()
    assert (out_dir / "alternative_placements.tsv").is_file()
    assert (out_dir / "best_placements.nwk").is_file()
    assert (out_dir / "run.json").is_file()
    assert len(load_newick_tree_set(out_dir / "best_placements.nwk")) == 2
