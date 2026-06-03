from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.engines import run_model_selection

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def _write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _fake_iqtree(path: Path, *, version: str = "2.9.9") -> Path:
    return _write_executable(
        path,
        f"""#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    print("IQ-TREE multicore version {version}")
    raise SystemExit(0)

prefix = Path(args[args.index("-pre") + 1])
prefix.parent.mkdir(parents=True, exist_ok=True)
prefix.with_suffix(".iqtree").write_text(
    " No. Model         -LnL         df  AIC          AICc         BIC\\n"
    "  1  GTR+G         123.456      12  270.912      330.912      272.912\\n"
    "  2  HKY+G         124.000      10  268.000      320.000      269.000\\n"
    "  3  JC            130.500      5   271.000      300.000      271.500\\n"
    "Akaike Information Criterion:           HKY+G\\n"
    "Corrected Akaike Information Criterion: JC\\n"
    "Bayesian Information Criterion:         GTR+G\\n"
    "Best-fit model according to BIC: GTR+G\\n"
    "Log-likelihood of the tree: -123.456\\n",
    encoding="utf-8",
)
prefix.with_suffix(".log").write_text(
    "IQ-TREE fixture model-selection log\\nBEST SCORE FOUND : -123.456\\n",
    encoding="utf-8",
)
prefix.with_suffix(".model").write_text(
    "Best-fit model: GTR+G\\n",
    encoding="utf-8",
)
raise SystemExit(0)
""",
    )


def test_phylo_replay_cli_reports_equivalent_outputs(tmp_path: Path, capsys) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")
    workflow = run_model_selection(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "model-selection",
        executable=executable,
        prefix="example",
    )

    exit_code = main(
        [
            "phylo",
            "replay",
            str(workflow.manifest_path),
            "--out-dir",
            str(tmp_path / "model-selection-replay"),
            "--iqtree-executable",
            str(executable),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["workflow"] == "model-selection"
    assert payload["metrics"]["outputs_equivalent"] is True
    assert payload["data"]["outputs_equivalent"] is True


def test_phylo_replay_cli_fails_for_changed_inputs(tmp_path: Path, capsys) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = tmp_path / "input.fasta"
    input_path.write_text(
        fixture("alignments/example_alignment.fasta").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    workflow = run_model_selection(
        input_path,
        out_dir=tmp_path / "model-selection",
        executable=executable,
        prefix="example",
    )
    input_path.write_text(
        input_path.read_text(encoding="utf-8").replace("ACTGACTG", "ACTGACTA", 1),
        encoding="utf-8",
    )

    exit_code = main(
        [
            "phylo",
            "replay",
            str(workflow.manifest_path),
            "--out-dir",
            str(tmp_path / "model-selection-replay"),
            "--iqtree-executable",
            str(executable),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "manifest_replay_input_changed"


def test_phylo_replay_cli_reports_engine_version_drift(tmp_path: Path, capsys) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture", version="2.9.9")
    drifted_executable = _fake_iqtree(tmp_path / "iqtree-drifted", version="3.0.0")
    workflow = run_model_selection(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "model-selection",
        executable=executable,
        prefix="example",
    )

    exit_code = main(
        [
            "phylo",
            "replay",
            str(workflow.manifest_path),
            "--out-dir",
            str(tmp_path / "model-selection-replay"),
            "--iqtree-executable",
            str(drifted_executable),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["engine_version_drift_count"] == 1
    assert payload["metrics"]["outputs_equivalent"] is True
    assert payload["data"]["engine_version_drift_detected"] is True
    assert payload["data"]["engine_version_drift"][0]["label"] == "iqtree"
    assert payload["data"]["engine_version_drift"][0]["matched"] is False
