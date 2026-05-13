from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.cli import main
from bijux_phylogenetics.engines.fasta_to_tree import run_fasta_to_tree_workflow
from tests.support.external_engines import require_alignment_engine_executables

REPO_ROOT = Path(__file__).resolve().parents[3]
REAL_INPUT_ROOT = Path(
    "packages/bijux-phylogenetics/tests/fixtures/fasta_to_tree/real/inputs"
)
EXPECTED_ROOT = REPO_ROOT / "packages/bijux-phylogenetics/tests/fixtures/expected/fasta_to_tree"
USER_FACING_SUFFIXES = (
    ".aln",
    ".trimmed.aln",
    ".tree",
    ".log",
    ".model.tsv",
    ".support.tsv",
)

def _assert_matches_expected_bundle(actual_root: Path, slug: str) -> None:
    expected_root = EXPECTED_ROOT / slug
    for suffix in USER_FACING_SUFFIXES:
        expected_path = expected_root / f"{slug}{suffix}"
        actual_path = actual_root / f"{slug}{suffix}"
        assert actual_path.exists(), actual_path
        assert actual_path.read_bytes() == expected_path.read_bytes(), actual_path.name


@pytest.mark.parametrize(
    ("slug", "sequence_type"),
    (
        ("gnathostome-ortholog-proteins", "protein"),
        ("gnathostome-ortholog-coding-sequences", "dna"),
        ("strnog-enog411bqtj-proteins", "protein"),
    ),
)
def test_run_fasta_to_tree_workflow_matches_real_output_golden(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    slug: str,
    sequence_type: str,
) -> None:
    executables = require_alignment_engine_executables()
    monkeypatch.chdir(REPO_ROOT)
    input_path = REAL_INPUT_ROOT / f"{slug}.fasta"
    out_dir = tmp_path / slug

    run_fasta_to_tree_workflow(
        input_path,
        out_dir=out_dir,
        prefix=slug,
        sequence_type=sequence_type,
        mafft_executable=executables["mafft"],
        trimal_executable=executables["trimal"],
        iqtree_executable=executables["iqtree2"],
        iqtree_seed=1,
        iqtree_threads=1,
        bootstrap_replicates=1000,
    )

    _assert_matches_expected_bundle(out_dir, slug)


def test_adapter_fasta_to_tree_cli_matches_real_output_golden(
    tmp_path: Path,
    capsys,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    executables = require_alignment_engine_executables()
    monkeypatch.chdir(REPO_ROOT)
    slug = "gnathostome-ortholog-proteins"
    input_path = REAL_INPUT_ROOT / f"{slug}.fasta"
    out_dir = tmp_path / slug

    exit_code = main(
        [
            "adapter",
            "fasta-to-tree",
            str(input_path),
            "--out-dir",
            str(out_dir),
            "--prefix",
            slug,
            "--sequence-type",
            "protein",
            "--mafft-executable",
            executables["mafft"],
            "--trimal-executable",
            executables["trimal"],
            "--iqtree-executable",
            executables["iqtree2"],
            "--iqtree-seed",
            "1",
            "--iqtree-threads",
            "1",
            "--bootstrap-replicates",
            "1000",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["selected_model"] == "Q.insect+I"
    assert payload["metrics"]["sequence_type"] == "protein"
    _assert_matches_expected_bundle(out_dir, slug)
