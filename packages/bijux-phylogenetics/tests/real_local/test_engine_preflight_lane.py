from __future__ import annotations

import json

import pytest

from bijux_phylogenetics.command_line import main

from ..support.external_engines import (
    real_beast_executable,
    require_alignment_validation_matrix_executables,
)

pytestmark = [pytest.mark.real_local, pytest.mark.engine_real]


def test_phylo_preflight_cli_reports_real_alignment_workflow_readiness(capsys) -> None:
    executables = require_alignment_validation_matrix_executables()

    exit_code = main(
        [
            "phylo",
            "preflight",
            "--workflow",
            "fasta-to-tree",
            "--mafft-executable",
            executables["mafft"],
            "--trimal-executable",
            executables["trimal"],
            "--iqtree-executable",
            executables["iqtree"],
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["selected_workflow"] == "fasta-to-tree"
    assert payload["metrics"]["selected_workflow_status"] in {"ready", "caution"}

    selected = next(
        workflow
        for workflow in payload["data"]["workflows"]
        if workflow["workflow_id"] == "fasta-to-tree"
    )
    assert selected["runnable"] is True
    assert selected["readiness_status"] in {"ready", "caution"}
    assert selected["blocking_engines"] == []

    for engine_id in ("mafft", "trimal", "iqtree"):
        engine = next(
            item
            for item in payload["data"]["engines"]
            if item["engine_id"] == engine_id
        )
        assert engine["available"] is True
        assert engine["executable_path"] == executables[engine_id]
        assert engine["support_status"] in {"tested", "untested"}
        assert engine["version_text"]


def test_phylo_preflight_cli_reports_real_beast_workflow_gate(capsys) -> None:
    executable = real_beast_executable()
    arguments = [
        "phylo",
        "preflight",
        "--workflow",
        "beast-posterior",
        "--json",
    ]
    if executable is not None:
        arguments.extend(["--beast-executable", str(executable)])

    exit_code = main(arguments)
    payload = json.loads(capsys.readouterr().out)

    if executable is None:
        assert exit_code == 2
        assert payload["status"] == "error"
        assert payload["errors"][0]["code"] == "engine_preflight_workflow_blocked"
        assert payload["errors"][0]["details"]["workflow_id"] == "beast-posterior"
        assert payload["errors"][0]["details"]["blocking_engines"] == ["BEAST"]
        return

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["selected_workflow"] == "beast-posterior"
    assert payload["metrics"]["selected_workflow_status"] in {"ready", "caution"}

    beast_engine = next(
        item for item in payload["data"]["engines"] if item["engine_id"] == "beast"
    )
    assert beast_engine["available"] is True
    assert beast_engine["executable_path"] == str(executable)
    assert beast_engine["support_status"] in {"tested", "untested"}
    assert beast_engine["version_text"]
