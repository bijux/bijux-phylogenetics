from __future__ import annotations

import json

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.parity import (
    run_ape_parity_cases,
    run_geiger_parity_cases,
    run_phytools_parity_cases,
)
from bijux_phylogenetics.parity.live_access import RUN_LIVE_PARITY_ENV_VAR
from bijux_phylogenetics.runtime.errors import ExternalEngineAccessError


def test_checked_fixture_parity_cli_runs_without_live_opt_in(
    capsys,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv(RUN_LIVE_PARITY_ENV_VAR, raising=False)

    exit_code = main(
        [
            "parity",
            "--reference-method",
            "pgls",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["reference_source"] == "checked-fixture"
    assert payload["metrics"]["method_count"] == 1


@pytest.mark.parametrize(
    ("reference_source", "selection_flag", "selection_value"),
    [
        ("ape-live", "--ape-case", "get-mrca-balanced-two-tip"),
        ("phytools-live", "--phytools-case", "simmap-er-binary-twenty-four-taxa"),
        (
            "geiger-live",
            "--geiger-case",
            "fitcontinuous-lambda-weak-signal-review",
        ),
    ],
)
def test_live_parity_cli_requires_explicit_opt_in(
    reference_source: str,
    selection_flag: str,
    selection_value: str,
    capsys,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv(RUN_LIVE_PARITY_ENV_VAR, raising=False)

    exit_code = main(
        [
            "parity",
            "--reference-source",
            reference_source,
            selection_flag,
            selection_value,
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "external_engine_access_forbidden"
    assert payload["errors"][0]["details"]["surface"] == (
        "parity command live reference dispatch"
    )
    assert payload["errors"][0]["details"]["reference_source"] == reference_source
    assert (
        payload["errors"][0]["details"]["environment_variable"]
        == RUN_LIVE_PARITY_ENV_VAR
    )


def test_live_parity_cli_reports_invalid_environment_values(
    capsys,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv(RUN_LIVE_PARITY_ENV_VAR, "sometimes")

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "ape-live",
            "--ape-case",
            "get-mrca-balanced-two-tip",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "live_parity_opt_in_invalid"
    assert payload["errors"][0]["details"]["received_value"] == "sometimes"


@pytest.mark.parametrize(
    "runner",
    [
        run_ape_parity_cases,
        run_phytools_parity_cases,
        run_geiger_parity_cases,
    ],
)
def test_live_parity_runners_require_explicit_opt_in(
    runner,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv(RUN_LIVE_PARITY_ENV_VAR, raising=False)

    with pytest.raises(ExternalEngineAccessError, match="RUN_LIVE_PARITY=1"):
        runner()
