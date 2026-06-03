from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines.inference import run_fasta_to_tree_workflow
from bijux_phylogenetics.io.artifact_schema import (
    validate_artifact_schema,
    validate_fasta_to_tree_manifest_schema,
    validate_fasta_to_tree_model_table_schema,
    validate_fasta_to_tree_support_table_schema,
    validate_run_manifest_schema,
)
from tests.support.fake_external_engines import fake_iqtree, fake_mafft, fake_trimal

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def _run_schema_fixture_workflow(tmp_path: Path):
    engine_root = tmp_path / "engines"
    engine_root.mkdir(parents=True, exist_ok=True)
    mafft = fake_mafft(engine_root / "mafft")
    trimal = fake_trimal(engine_root / "trimal")
    iqtree = fake_iqtree(engine_root / "iqtree2")
    return run_fasta_to_tree_workflow(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "fasta-to-tree",
        prefix="example",
        mafft_executable=mafft,
        trimal_executable=trimal,
        iqtree_executable=iqtree,
        bootstrap_replicates=1000,
    )


def test_fasta_to_tree_outputs_match_stable_artifact_schemas(tmp_path: Path) -> None:
    report = _run_schema_fixture_workflow(tmp_path)

    assert report.method_tier.tier == "supported"
    assert validate_fasta_to_tree_model_table_schema(
        report.output_paths["model_table"]
    ).valid
    assert validate_fasta_to_tree_support_table_schema(
        report.output_paths["support_table"]
    ).valid
    assert validate_fasta_to_tree_manifest_schema(report.output_paths["manifest"]).valid
    assert validate_run_manifest_schema(report.run_manifest_path).valid


def test_named_artifact_schema_validation_accepts_fasta_to_tree_outputs(
    tmp_path: Path,
) -> None:
    report = _run_schema_fixture_workflow(tmp_path)

    assert validate_artifact_schema(
        report.output_paths["model_table"],
        "fasta_to_tree_model_tsv",
    ).valid
    assert validate_artifact_schema(
        report.output_paths["support_table"],
        "fasta_to_tree_support_tsv",
    ).valid
    assert validate_artifact_schema(
        report.output_paths["manifest"],
        "fasta_to_tree_manifest_json",
    ).valid
    assert validate_artifact_schema(
        report.run_manifest_path,
        "run_manifest_json",
    ).valid
