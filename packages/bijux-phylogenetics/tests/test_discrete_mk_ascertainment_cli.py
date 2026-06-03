from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_comparative_discrete_mk_cli_reports_lewis_ascertainment_metrics(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "discrete-mk",
            str(fixture("trees", "example_tree.nwk")),
            str(
                fixture(
                    "metadata", "example_traits_discrete_mk_variable_only_four_taxa.tsv"
                )
            ),
            "--trait",
            "state",
            "--taxon-column",
            "taxon",
            "--ascertainment",
            "lewis-variable-only",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["ascertainment_policy"] == "lewis-variable-only"
    assert payload["metrics"]["ascertainment_conditioning_log_probability"] is not None
    assert payload["metrics"]["invariant_pattern_log_probability"] is not None
