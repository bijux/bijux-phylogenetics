from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
from tests.support.fake_external_engines import fake_iqtree, fake_mafft, fake_trimal
from tests.test_workflow_config import _write_config


def _sentinel_engine(path: Path, sentinel_path: Path, version_text: str) -> Path:
    path.write_text(
        (
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "from pathlib import Path\n"
            f"sentinel = Path({str(sentinel_path)!r})\n"
            "sentinel.write_text('invoked\\n', encoding='utf-8')\n"
            "if '--version' in sys.argv:\n"
            f"    print({version_text!r})\n"
            "    raise SystemExit(0)\n"
            "raise SystemExit(0)\n"
        ),
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def test_phylo_run_cli_executes_workflow_config_and_validates_bundle(
    tmp_path: Path, capsys
) -> None:
    input_path = tmp_path / "input.fasta"
    input_path.write_text(
        ">A\nACTGACTG\n>B\nACTGACTA\n>C\nACTGACTC\n>D\nACTGACTT\n",
        encoding="utf-8",
    )
    metadata_path = tmp_path / "metadata.tsv"
    metadata_path.write_text(
        "taxon\tregion\nA\twest\nB\teast\nC\tnorth\nD\tsouth\n",
        encoding="utf-8",
    )
    traits_path = tmp_path / "traits.tsv"
    traits_path.write_text(
        "taxon\tbody_mass\nA\t1.2\nB\t1.4\nC\t1.1\nD\t1.7\n",
        encoding="utf-8",
    )
    config_path = _write_config(
        tmp_path / "workflow-config.yaml",
        input_path=input_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        mafft_executable=fake_mafft(tmp_path / "mafft-fixture"),
        trimal_executable=fake_trimal(tmp_path / "trimal-fixture"),
        iqtree_executable=fake_iqtree(tmp_path / "iqtree-fixture"),
    )

    exit_code = main(["phylo", "run", str(config_path), "--json"])
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["workflow"] == "fasta-to-tree"
    assert payload["metrics"]["bundle_validation_passed"] is True
    assert payload["metrics"]["metadata_present"] is True
    assert payload["metrics"]["traits_present"] is True


def test_phylo_run_cli_rejects_invalid_config_before_engine_preflight(
    tmp_path: Path, capsys
) -> None:
    input_path = tmp_path / "input.fasta"
    input_path.write_text(">A\nACTGACTG\n>B\nACTGACTA\n", encoding="utf-8")
    sentinel_path = tmp_path / "engine-invoked.txt"
    config_path = _write_config(
        tmp_path / "workflow-config.yaml",
        input_path=input_path,
        mafft_executable=_sentinel_engine(
            tmp_path / "mafft-sentinel",
            sentinel_path,
            "mafft v7.999",
        ),
        trimal_executable=_sentinel_engine(
            tmp_path / "trimal-sentinel",
            sentinel_path,
            "trimAl v2.0",
        ),
        iqtree_executable=_sentinel_engine(
            tmp_path / "iqtree-sentinel",
            sentinel_path,
            "IQ-TREE multicore version 2.9.9",
        ),
        alignment_mode="not-a-real-mode",
    )

    exit_code = main(["phylo", "run", str(config_path), "--json"])
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "workflow_config_invalid"
    assert sentinel_path.exists() is False
