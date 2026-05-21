from __future__ import annotations

import json

from bijux_phylogenetics.command_line import main


def test_validate_sim_char_reference_cli_reports_governed_cases(capsys) -> None:
    exit_code = main(
        [
            "simulate",
            "validate-sim-char-reference",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["case_count"] == 3
    assert payload["metrics"]["all_passed"] is True
    assert [row["case_id"] for row in payload["data"]["observations"]] == [
        "brownian-internal-long-branch-envelope",
        "speciational-internal-long-branch-envelope",
        "discrete-rate-matrix-internal-long-branch-envelope",
    ]
