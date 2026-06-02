from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.engines.inference import run_fasta_to_tree_workflow
from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree

from .support.external_engines import require_alignment_engine_executables

pytestmark = [
    pytest.mark.real_local,
    pytest.mark.evaluation,
    pytest.mark.scientific_validation,
    pytest.mark.slow,
]

REPO_ROOT = Path(__file__).resolve().parents[3]
REAL_INPUT_ROOT = Path(
    "packages/bijux-phylogenetics/tests/fixtures/fasta_to_tree/real/inputs"
)
EXPECTED_ROOT = (
    REPO_ROOT / "packages/bijux-phylogenetics/tests/fixtures/expected/fasta_to_tree"
)
USER_FACING_SUFFIXES = (
    ".aln",
    ".trimmed.aln",
    ".tree",
    ".log",
    ".methods-summary.md",
    ".model.tsv",
    ".support.tsv",
    ".manifest.json",
    ".run.json",
)

REAL_DATASETS = (
    (
        "gnathostome-ortholog-proteins",
        "protein",
        REAL_INPUT_ROOT / "gnathostome-ortholog-proteins.fasta",
    ),
    (
        "gnathostome-ortholog-coding-sequences",
        "dna",
        REAL_INPUT_ROOT / "gnathostome-ortholog-coding-sequences.fasta",
    ),
    (
        "strnog-enog411bqtj-proteins",
        "protein",
        REAL_INPUT_ROOT / "strnog-enog411bqtj-proteins.fasta",
    ),
    (
        "vertebrate-homolog-proteins",
        "protein",
        REAL_INPUT_ROOT / "vertebrate-homolog-proteins.fasta",
    ),
    (
        "influenza-a-ha-reference-panel",
        "dna",
        Path(
            "packages/bijux-phylogenetics/src/bijux_phylogenetics/resources/datasets/viruses/influenza_a_ha_reference_panel/sequences.fasta"
        ),
    ),
    (
        "pleistocene-bear-cytb-fragments",
        "dna",
        Path(
            "packages/bijux-phylogenetics/src/bijux_phylogenetics/resources/datasets/ancient_dna/pleistocene_bear_cytb_fragments/sequences.fasta"
        ),
    ),
)


def _read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(encoding="utf-8"), delimiter="\t"))


def _alignment_signature(path: Path) -> tuple[tuple[str, ...], int]:
    records = load_fasta_alignment(path)
    identifiers = tuple(record.identifier for record in records)
    alignment_length = len(records[0].sequence)
    return identifiers, alignment_length


def _tree_signature(path: Path) -> tuple[tuple[str, ...], int]:
    tree = load_tree(path)
    tip_names = tuple(sorted(tree.tip_names))
    internal_node_count = tree.internal_node_count
    return tip_names, internal_node_count


def _support_signature(path: Path) -> tuple[tuple[str, float], ...]:
    rows = _read_tsv_rows(path)
    return tuple(
        sorted(
            (
                row["descendant_taxa"],
                round(float(row["support"]), 6),
            )
            for row in rows
        )
    )


def _model_signature(path: Path) -> dict[str, str]:
    rows = _read_tsv_rows(path)
    assert len(rows) == 1
    return rows[0]


def _normalize_paths(values: list[str]) -> tuple[str, ...]:
    return tuple(sorted(Path(value).name for value in values))


def _normalize_checksum_keys(checksums: dict[str, str]) -> tuple[str, ...]:
    return tuple(sorted(Path(path).name for path in checksums))


def _normalize_checksum_map(checksums: dict[str, str]) -> dict[str, str]:
    return {Path(path).name: value for path, value in checksums.items()}


def _parse_run_arguments(arguments: list[str]) -> dict[str, str]:
    parsed: dict[str, str] = {}
    index = 1
    while index < len(arguments):
        flag = arguments[index]
        if index + 1 >= len(arguments):
            break
        parsed[flag] = arguments[index + 1]
        index += 2
    return parsed


def _assert_matches_expected_bundle(
    actual_root: Path,
    slug: str,
    *,
    sequence_type: str,
    input_path: Path,
) -> None:
    expected_root = EXPECTED_ROOT / slug
    for suffix in USER_FACING_SUFFIXES:
        expected_path = expected_root / f"{slug}{suffix}"
        actual_path = actual_root / f"{slug}{suffix}"
        assert actual_path.exists(), actual_path
        assert expected_path.exists(), expected_path

    actual_alignment = actual_root / f"{slug}.aln"
    expected_alignment = expected_root / f"{slug}.aln"
    actual_trimmed = actual_root / f"{slug}.trimmed.aln"
    expected_trimmed = expected_root / f"{slug}.trimmed.aln"
    actual_tree = actual_root / f"{slug}.tree"
    expected_tree = expected_root / f"{slug}.tree"
    actual_log = actual_root / f"{slug}.log"
    expected_log = expected_root / f"{slug}.log"
    actual_model = actual_root / f"{slug}.model.tsv"
    expected_model = expected_root / f"{slug}.model.tsv"
    actual_support = actual_root / f"{slug}.support.tsv"
    expected_support = expected_root / f"{slug}.support.tsv"
    actual_manifest = actual_root / f"{slug}.manifest.json"
    expected_manifest = expected_root / f"{slug}.manifest.json"
    actual_run_manifest = actual_root / f"{slug}.run.json"
    expected_run_manifest = expected_root / f"{slug}.run.json"

    assert _alignment_signature(actual_alignment) == _alignment_signature(
        expected_alignment
    )
    assert _alignment_signature(actual_trimmed) == _alignment_signature(
        expected_trimmed
    )
    assert _tree_signature(actual_tree) == _tree_signature(expected_tree)
    assert _support_signature(actual_support) == _support_signature(expected_support)

    actual_model_row = _model_signature(actual_model)
    expected_model_row = _model_signature(expected_model)
    for key in (
        "workflow",
        "engine_name",
        "sequence_type",
        "selected_model",
        "report_selected_model",
        "artifact_selected_model",
        "model_consistent",
    ):
        assert actual_model_row[key] == expected_model_row[key]
    assert actual_model_row["sequence_type"] == sequence_type

    actual_log_text = actual_log.read_text(encoding="utf-8")
    expected_log_text = expected_log.read_text(encoding="utf-8")
    for stable_line in (
        "workflow: fasta-to-tree",
        f"sequence_type: {sequence_type}",
        f"selected_model: {expected_model_row['selected_model']}",
        f"run_manifest_path: {slug}.run.json",
    ):
        assert stable_line in actual_log_text
        assert stable_line in expected_log_text
    for section in (
        "[alignment]",
        "[trimming]",
        "[model_selection]",
        "[maximum_likelihood]",
        "[bootstrap_support]",
    ):
        assert section in actual_log_text

    actual_manifest_payload = _read_json(actual_manifest)
    expected_manifest_payload = _read_json(expected_manifest)
    for key in (
        "sequence_type",
        "selected_model",
        "alignment_mode",
        "trimming_mode",
        "trim_gap_threshold",
        "iqtree_seed",
        "iqtree_threads",
        "bootstrap_replicates",
    ):
        assert actual_manifest_payload[key] == expected_manifest_payload[key]
    assert actual_manifest_payload["sequence_type"] == sequence_type
    assert (
        actual_manifest_payload["commands"].keys()
        == expected_manifest_payload["commands"].keys()
    )
    assert (
        actual_manifest_payload["engine_versions"].keys()
        == expected_manifest_payload["engine_versions"].keys()
    )
    assert (
        actual_manifest_payload["step_manifests"].keys()
        == expected_manifest_payload["step_manifests"].keys()
    )
    assert _normalize_checksum_map(actual_manifest_payload["input_checksums"]) == (
        _normalize_checksum_map(expected_manifest_payload["input_checksums"])
    )
    assert _normalize_paths(list(actual_manifest_payload["output_paths"].values())) == (
        _normalize_paths(list(expected_manifest_payload["output_paths"].values()))
    )
    assert _normalize_checksum_keys(actual_manifest_payload["output_checksums"]) == (
        _normalize_checksum_keys(expected_manifest_payload["output_checksums"])
    )
    assert all(
        len(value) == 64
        for value in actual_manifest_payload["input_checksums"].values()
    )
    assert all(
        len(value) == 64
        for value in actual_manifest_payload["output_checksums"].values()
    )

    actual_run_payload = _read_json(actual_run_manifest)
    expected_run_payload = _read_json(expected_run_manifest)
    assert (
        actual_run_payload["command"]
        == expected_run_payload["command"]
        == "run_fasta_to_tree_workflow"
    )
    actual_arguments = _parse_run_arguments(actual_run_payload["arguments"])
    expected_arguments = _parse_run_arguments(expected_run_payload["arguments"])
    for key in (
        "--prefix",
        "--sequence-type",
        "--alignment-mode",
        "--trimming-mode",
        "--trim-gap-threshold",
        "--iqtree-seed",
        "--iqtree-threads",
        "--bootstrap-replicates",
        "--resume",
        "--incomplete-run-policy",
    ):
        assert actual_arguments[key] == expected_arguments[key]
    assert actual_arguments["--prefix"] == slug
    assert actual_arguments["--sequence-type"] == sequence_type
    assert _normalize_checksum_map(actual_run_payload["input_checksums"]) == (
        _normalize_checksum_map(expected_run_payload["input_checksums"])
    )
    assert str(input_path) in actual_run_payload["input_paths"]
    assert _normalize_paths(actual_run_payload["output_paths"]) == _normalize_paths(
        expected_run_payload["output_paths"]
    )
    assert _normalize_checksum_keys(actual_run_payload["output_checksums"]) == (
        _normalize_checksum_keys(expected_run_payload["output_checksums"])
    )
    assert all(
        len(value) == 64 for value in actual_run_payload["input_checksums"].values()
    )
    assert all(
        len(value) == 64 for value in actual_run_payload["output_checksums"].values()
    )


@pytest.mark.parametrize(
    ("slug", "sequence_type", "input_path"),
    REAL_DATASETS,
)
def test_run_fasta_to_tree_workflow_matches_real_output_golden(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    slug: str,
    sequence_type: str,
    input_path: Path,
) -> None:
    executables = require_alignment_engine_executables()
    monkeypatch.chdir(REPO_ROOT)
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

    _assert_matches_expected_bundle(
        out_dir,
        slug,
        sequence_type=sequence_type,
        input_path=input_path,
    )


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
    _assert_matches_expected_bundle(
        out_dir,
        slug,
        sequence_type="protein",
        input_path=input_path,
    )
