from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_phylo_likelihood_strict_clock_cli_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "strict-clock-likelihood"

    exit_code = main(
        [
            "phylo",
            "likelihood",
            "strict-clock",
            str(fixture("trees", "strict_clock_time_tree_4_taxa.nwk")),
            str(
                fixture("alignments", "strict_clock_likelihood_alignment_4_taxa.fasta")
            ),
            "--model",
            "jc69",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["model_name"] == "JC69"
    assert payload["metrics"]["taxon_count"] == 4
    assert payload["metrics"]["site_count"] == 12
    assert payload["metrics"]["pattern_count"] == 2
    assert payload["metrics"]["branch_count"] == 6
    assert math.isclose(
        payload["metrics"]["optimized_clock_rate"],
        0.5217383669434831,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        payload["metrics"]["optimized_log_likelihood"],
        -64.38812241070909,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert payload["metrics"]["function_evaluation_count"] == 51
    assert payload["metrics"]["converged"] is True
    assert (out_dir / "scaled_tree.nwk").is_file()
    assert (out_dir / "branch_rates.tsv").is_file()
    assert (out_dir / "branch_likelihood_diagnostics.tsv").is_file()
    assert (out_dir / "site_log_likelihoods.tsv").is_file()
    assert (out_dir / "run.json").is_file()
