from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_adapter_beast_log_cli_writes_posterior_decomposition_bundle(
    tmp_path: Path,
    capsys,
) -> None:
    log_path = tmp_path / "classified-beast.log"
    summary_path = tmp_path / "classified-beast-summary.tsv"
    decomposition_path = tmp_path / "classified-beast-decomposition.tsv"
    log_path.write_text(
        "# posterior decomposition fixture\n"
        "state\tposterior\tlikelihood\tprior\tclockRate\n"
        "0\t-510.0\t-490.0\t-20.0\t0.0010\n"
        "1000\t-505.0\t-486.0\t-19.0\t0.0011\n"
        "2000\t-500.0\t-482.0\t-18.0\t0.0012\n"
        "3000\t-497.0\t-479.0\t-18.0\t0.0013\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "adapter",
            "beast-log",
            str(log_path),
            "--burnin-fraction",
            "0.25",
            "--summary-out",
            str(summary_path),
            "--decomposition-out",
            str(decomposition_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["decomposition_verified"] is True
    assert payload["metrics"]["decomposition_maximum_absolute_delta"] == 0.0
    assert payload["data"]["decomposition"]["prior_term_source"] == "logged"
    assert summary_path.exists()
    assert decomposition_path.exists()


def test_adapter_mrbayes_parameters_cli_writes_posterior_decomposition_bundle(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "mrbayes-parameters.tsv"
    decomposition_path = tmp_path / "mrbayes-decomposition.tsv"

    exit_code = main(
        [
            "adapter",
            "mrbayes-parameters",
            str(fixture("mrbayes", "partitioned-analysis.run1.p")),
            "--summary-out",
            str(summary_path),
            "--decomposition-out",
            str(decomposition_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["decomposition_verified"] is True
    assert payload["data"]["decomposition"]["posterior_term_source"] == (
        "derived_from_likelihood_and_prior"
    )
    assert summary_path.exists()
    assert decomposition_path.exists()
