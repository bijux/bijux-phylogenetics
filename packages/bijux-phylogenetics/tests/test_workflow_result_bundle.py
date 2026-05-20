from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines import run_fasta_to_tree_workflow
from bijux_phylogenetics.engines.bundles import (
    export_workflow_result_bundle,
    validate_workflow_result_bundle,
)
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


def test_export_workflow_result_bundle_copies_fasta_to_tree_handoff_surface(
    tmp_path: Path,
) -> None:
    workflow = _run_example_fasta_to_tree(tmp_path)

    bundle = export_workflow_result_bundle(
        workflow.manifest_path,
        bundle_root=tmp_path / "bundle",
    )
    validation = validate_workflow_result_bundle(bundle.bundle_root)

    assert bundle.workflow == "fasta-to-tree"
    assert bundle.includes_input_files is True
    assert bundle.copied_output_count >= 7
    assert bundle.copied_step_manifest_count == 5
    assert bundle.copied_step_output_count >= 5
    assert bundle.copied_report_count == 6
    assert bundle.report_path.exists()
    assert bundle.config_path.exists()
    assert bundle.rerun_path.exists()
    assert bundle.workflow_manifest_path.exists()
    assert validation.valid is True
    assert validation.issues == []


def test_validate_workflow_result_bundle_detects_missing_required_workflow_report(
    tmp_path: Path,
) -> None:
    workflow = _run_example_fasta_to_tree(tmp_path)
    bundle = export_workflow_result_bundle(
        workflow.manifest_path,
        bundle_root=tmp_path / "bundle",
    )
    bundle.report_path.unlink()

    validation = validate_workflow_result_bundle(bundle.bundle_root)

    assert validation.valid is False
    assert any(issue.label == "workflow_report" for issue in validation.issues)


def test_validate_workflow_result_bundle_detects_missing_required_final_tree(
    tmp_path: Path,
) -> None:
    workflow = _run_example_fasta_to_tree(tmp_path)
    bundle = export_workflow_result_bundle(
        workflow.manifest_path,
        bundle_root=tmp_path / "bundle",
    )
    tree_output = next(
        file.relative_path
        for file in bundle.files
        if file.role == "workflow_output" and file.label == "tree"
    )
    (bundle.bundle_root / tree_output).unlink()

    validation = validate_workflow_result_bundle(bundle.bundle_root)

    assert validation.valid is False
    assert any(issue.label == "tree" for issue in validation.issues)
