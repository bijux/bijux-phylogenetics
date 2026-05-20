from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.engines.inference import (
    load_phylo_workflow_config,
    run_phylo_workflow_config,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError
from tests.support.fake_external_engines import fake_iqtree, fake_mafft, fake_trimal


def _write_config(
    path: Path,
    *,
    input_path: Path,
    metadata_path: Path | None = None,
    traits_path: Path | None = None,
    mafft_executable: Path,
    trimal_executable: Path,
    iqtree_executable: Path,
    out_dir: str = "artifacts/run",
    prefix: str = "example-study",
    alignment_mode: str = "auto",
) -> Path:
    metadata_block = ""
    if metadata_path is not None:
        metadata_block = (
            f"  metadata: {metadata_path.name}\n  metadata_taxon_column: taxon\n"
        )
    traits_block = ""
    if traits_path is not None:
        traits_block = f"  traits: {traits_path.name}\n  traits_taxon_column: taxon\n"
    path.write_text(
        (
            "workflow: fasta-to-tree\n"
            "inputs:\n"
            f"  fasta: {input_path.name}\n"
            f"{metadata_block}"
            f"{traits_block}"
            "engines:\n"
            f"  mafft_executable: {mafft_executable}\n"
            f"  trimal_executable: {trimal_executable}\n"
            f"  iqtree_executable: {iqtree_executable}\n"
            "alignment:\n"
            f"  mode: {alignment_mode}\n"
            "  sequence_type: dna\n"
            "trimming:\n"
            "  mode: gap-threshold\n"
            "  gap_threshold: 0.2\n"
            "inference:\n"
            "  bootstrap_replicates: 1000\n"
            "  seed: 17\n"
            "  threads: 2\n"
            "output:\n"
            f"  out_dir: {out_dir}\n"
            f"  prefix: {prefix}\n"
            "resources:\n"
            "  timeout_seconds: 45\n"
            "  resume: false\n"
            "  incomplete_run_policy: reject\n"
        ),
        encoding="utf-8",
    )
    return path


def test_load_phylo_workflow_config_resolves_relative_inputs_and_outputs(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "input.fasta"
    input_path.write_text(">A\nACTGACTG\n>B\nACTGACTA\n", encoding="utf-8")
    metadata_path = tmp_path / "metadata.tsv"
    metadata_path.write_text("taxon\tregion\nA\twest\nB\teast\n", encoding="utf-8")
    traits_path = tmp_path / "traits.tsv"
    traits_path.write_text("taxon\tbody_mass\nA\t1.2\nB\t1.4\n", encoding="utf-8")
    config_path = _write_config(
        tmp_path / "workflow-config.yaml",
        input_path=input_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        mafft_executable=tmp_path / "mafft",
        trimal_executable=tmp_path / "trimal",
        iqtree_executable=tmp_path / "iqtree2",
    )

    config = load_phylo_workflow_config(config_path)

    assert config.workflow == "fasta-to-tree"
    assert config.input_fasta_path == input_path.resolve()
    assert config.metadata_path == metadata_path.resolve()
    assert config.traits_path == traits_path.resolve()
    assert config.out_dir == (tmp_path / "artifacts/run").resolve()
    assert (
        config.bundle_root
        == (tmp_path / "artifacts/run/example-study.result-bundle").resolve()
    )
    assert config.iqtree_seed == 17
    assert config.iqtree_threads == 2
    assert config.timeout_seconds == 45.0
    assert config.resolved_payload["workflow"] == "fasta-to-tree"


def test_load_phylo_workflow_config_rejects_invalid_modes_before_execution(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "input.fasta"
    input_path.write_text(">A\nACTGACTG\n>B\nACTGACTA\n", encoding="utf-8")
    config_path = _write_config(
        tmp_path / "workflow-config.yaml",
        input_path=input_path,
        mafft_executable=tmp_path / "mafft",
        trimal_executable=tmp_path / "trimal",
        iqtree_executable=tmp_path / "iqtree2",
        alignment_mode="bad-mode",
    )

    with pytest.raises(PhylogeneticsError) as error_info:
        load_phylo_workflow_config(config_path)

    assert error_info.value.code == "workflow_config_invalid"
    assert error_info.value.details["field"] == "alignment.mode"


def test_run_phylo_workflow_config_exports_valid_bundle_with_auxiliary_inputs(
    tmp_path: Path,
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
    mafft_executable = fake_mafft(tmp_path / "mafft-fixture")
    trimal_executable = fake_trimal(tmp_path / "trimal-fixture")
    iqtree_executable = fake_iqtree(tmp_path / "iqtree-fixture")
    config_path = _write_config(
        tmp_path / "workflow-config.yaml",
        input_path=input_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        mafft_executable=mafft_executable,
        trimal_executable=trimal_executable,
        iqtree_executable=iqtree_executable,
    )

    report = run_phylo_workflow_config(config_path)

    assert report.bundle_validation.valid is True
    assert report.output_paths["bundle_root"] == report.bundle_report.bundle_root
    bundle_manifest = json.loads(
        report.bundle_report.bundle_manifest_path.read_text(encoding="utf-8")
    )
    input_labels = {entry["label"] for entry in bundle_manifest["input_files"]}
    assert "workflow-config-source.yaml" in input_labels
    assert f"metadata-{metadata_path.name}" in input_labels
    assert f"traits-{traits_path.name}" in input_labels
    bundle_config = json.loads(
        report.bundle_report.config_path.read_text(encoding="utf-8")
    )
    assert bundle_config["inputs"]["metadata"] == str(metadata_path.resolve())
    assert bundle_config["inputs"]["traits"] == str(traits_path.resolve())
    assert bundle_config["inference"]["seed"] == 17
