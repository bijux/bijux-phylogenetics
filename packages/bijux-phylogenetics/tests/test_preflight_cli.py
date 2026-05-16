from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def _write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _fake_mafft(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys

if "--version" in sys.argv[1:]:
    print("MAFFT v7.520")
    raise SystemExit(0)

raise SystemExit(1)
""",
    )


def _fake_trimal(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys

if "--version" in sys.argv[1:]:
    print("trimAl v1.5.0")
    raise SystemExit(0)

raise SystemExit(1)
""",
    )


def _fake_iqtree(path: Path, *, version: str = "2.4.0") -> Path:
    return _write_executable(
        path,
        f"""#!/usr/bin/env python3
import sys

if "--version" in sys.argv[1:]:
    print("IQ-TREE multicore version {version}")
    raise SystemExit(0)

raise SystemExit(1)
""",
    )


def test_phylo_preflight_cli_reports_ready_selected_workflow(
    tmp_path: Path,
    capsys,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")

    exit_code = main(
        [
            "phylo",
            "preflight",
            "--workflow",
            "fasta-to-tree",
            "--mafft-executable",
            str(mafft),
            "--trimal-executable",
            str(trimal),
            "--iqtree-executable",
            str(iqtree),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["selected_workflow"] == "fasta-to-tree"
    assert payload["metrics"]["selected_workflow_status"] == "ready"
    assert payload["data"]["selected_workflow"] == "fasta-to-tree"
    selected = next(
        workflow
        for workflow in payload["data"]["workflows"]
        if workflow["workflow_id"] == "fasta-to-tree"
    )
    assert selected["runnable"] is True
    assert selected["blocking_engines"] == []
    mafft_status = next(
        engine
        for engine in payload["data"]["engines"]
        if engine["engine_id"] == "mafft"
    )
    assert mafft_status["executable_path"] == str(mafft)
    assert mafft_status["support_status"] == "tested"


def test_phylo_preflight_cli_reports_caution_for_untested_versions(
    tmp_path: Path,
    capsys,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture", version="4.0.0")

    exit_code = main(
        [
            "phylo",
            "preflight",
            "--workflow",
            "fasta-to-tree",
            "--mafft-executable",
            str(mafft),
            "--trimal-executable",
            str(trimal),
            "--iqtree-executable",
            str(iqtree),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["selected_workflow_status"] == "caution"
    selected = next(
        workflow
        for workflow in payload["data"]["workflows"]
        if workflow["workflow_id"] == "fasta-to-tree"
    )
    assert selected["runnable"] is True
    assert selected["caution_engines"] == ["IQ-TREE"]


def test_phylo_preflight_cli_fails_early_for_blocked_selected_workflow(
    tmp_path: Path,
    capsys,
) -> None:
    missing_beast = tmp_path / "missing-beast"

    exit_code = main(
        [
            "phylo",
            "preflight",
            "--workflow",
            "beast-posterior",
            "--beast-executable",
            str(missing_beast),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "engine_preflight_workflow_blocked"
    assert payload["errors"][0]["details"]["workflow_id"] == "beast-posterior"
    assert payload["errors"][0]["details"]["blocking_engines"] == ["BEAST"]
