from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_phylo_likelihood_local_clock_cli_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "local-clock-likelihood"

    exit_code = main(
        [
            "phylo",
            "likelihood",
            "local-clock",
            str(fixture("trees", "strict_clock_time_tree_4_taxa.nwk")),
            str(fixture("alignments", "local_clock_likelihood_alignment_4_taxa.fasta")),
            str(fixture("metadata", "local_clock_regimes_4_taxa.tsv")),
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
    assert payload["metrics"]["model_name"] == "JC69 local-clock"
    assert payload["metrics"]["taxon_count"] == 4
    assert payload["metrics"]["site_count"] == 60
    assert payload["metrics"]["pattern_count"] == 4
    assert payload["metrics"]["branch_count"] == 6
    assert payload["metrics"]["regime_count"] == 3
    assert payload["metrics"]["preferred_model_by_aic"] == "local-clock"
    assert (
        payload["metrics"]["strict_clock_aic"] > payload["metrics"]["local_clock_aic"]
    )
    assert payload["metrics"]["converged"] is True
    assert (out_dir / "scaled_tree.nwk").is_file()
    assert (out_dir / "branch_rates.tsv").is_file()
    assert (out_dir / "regimes.tsv").is_file()
    assert (out_dir / "branch_likelihood_diagnostics.tsv").is_file()
    assert (out_dir / "site_log_likelihoods.tsv").is_file()
    assert (out_dir / "run.json").is_file()
