from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.engines import run_fasta_to_tree_workflow
from tests.support.fake_external_engines import fake_iqtree, fake_mafft, fake_trimal


def _run_example_fasta_to_tree(tmp_path: Path):
    input_path = tmp_path / "input.fasta"
    input_path.write_text(
        ">A\nACTGACTG\n>B\nACTGACTA\n>C\nACTGACTC\n>D\nACTGACTT\n",
        encoding="utf-8",
    )
    return run_fasta_to_tree_workflow(
        input_path,
        out_dir=tmp_path / "workflow",
        prefix="example",
        mafft_executable=fake_mafft(tmp_path / "mafft-fixture"),
        trimal_executable=fake_trimal(tmp_path / "trimal-fixture"),
        iqtree_executable=fake_iqtree(tmp_path / "iqtree-fixture"),
        bootstrap_replicates=1000,
    )


def test_phylo_bundle_cli_exports_and_validates_workflow_result_bundle(
    tmp_path: Path, capsys
) -> None:
    workflow = _run_example_fasta_to_tree(tmp_path)
    bundle_root = tmp_path / "bundle"

    exit_code = main(
        [
            "phylo",
            "bundle",
            str(workflow.manifest_path),
            "--out-dir",
            str(bundle_root),
            "--json",
        ]
    )
    bundle_payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert bundle_payload["status"] == "ok"
    assert bundle_payload["metrics"]["workflow"] == "fasta-to-tree"
    assert bundle_payload["metrics"]["validation_passed"] is True
    assert bundle_payload["metrics"]["copied_step_manifest_count"] == 5

    exit_code = main(["phylo", "validate-bundle", str(bundle_root), "--json"])
    validate_payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert validate_payload["status"] == "ok"
    assert validate_payload["data"]["valid"] is True


def test_phylo_validate_bundle_cli_fails_for_missing_required_workflow_report(
    tmp_path: Path, capsys
) -> None:
    workflow = _run_example_fasta_to_tree(tmp_path)
    bundle_root = tmp_path / "bundle"

    exit_code = main(
        [
            "phylo",
            "bundle",
            str(workflow.manifest_path),
            "--out-dir",
            str(bundle_root),
            "--json",
        ]
    )
    assert exit_code == 0
    _ = capsys.readouterr()
    (bundle_root / "reports" / "workflow-report.html").unlink()

    exit_code = main(["phylo", "validate-bundle", str(bundle_root), "--json"])
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "workflow_bundle_validation_failed"
