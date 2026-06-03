from __future__ import annotations

import json

from bijux_phylogenetics.command_line import main


def test_topology_distance_reference_cli_reports_passing_cases(capsys) -> None:
    exit_code = main(["topology", "distance-reference", "--json"])
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["case_count"] == 27
    assert payload["metrics"]["external_case_count"] == 25
    assert payload["metrics"]["policy_case_count"] == 2
    assert payload["metrics"]["all_passed"] is True
