from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.cli import main
from bijux_phylogenetics.core.alignment import AlignmentRecord
from bijux_phylogenetics.io.fasta import write_fasta_alignment

FIXTURES = Path(__file__).parent / "fixtures" / "alignments"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_alignment_coding_cli_reports_genetic_code_and_invalid_codon_metrics(
    capsys,
) -> None:
    exit_code = main(
        [
            "alignment",
            "coding",
            str(fixture("example_alignment_coding.fasta")),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["genetic_code_id"] == 1
    assert payload["metrics"]["frameshift_like_sequence_count"] == 1
    assert payload["metrics"]["invalid_codon_count"] == 0
    assert payload["metrics"]["stop_codon_count"] == 2


def test_alignment_translate_cli_reports_genetic_code_and_invalid_codon_metrics(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "translated.fasta"
    exit_code = main(
        [
            "alignment",
            "translate",
            str(fixture("example_alignment_coding.fasta")),
            "--out",
            str(output_path),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        ">A\nME*\n>B\nMEW\n>C\nMXW\n>D\nM*W\n"
    )
    assert payload["metrics"]["genetic_code_id"] == 1
    assert payload["metrics"]["translated_sequence_count"] == 4
    assert payload["metrics"]["invalid_codon_count"] == 0
    assert payload["metrics"]["stop_codon_count"] == 2


def test_alignment_translate_cli_honors_configurable_genetic_code(
    tmp_path: Path, capsys
) -> None:
    input_path = tmp_path / "coding-mito.fasta"
    output_path = tmp_path / "translated.fasta"
    write_fasta_alignment(
        input_path,
        [
            AlignmentRecord(identifier="mito_triplet", sequence="ATGTGAGGG"),
        ],
    )

    exit_code = main(
        [
            "alignment",
            "translate",
            str(input_path),
            "--out",
            str(output_path),
            "--genetic-code",
            "2",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == ">mito_triplet\nMWG\n"
    assert payload["metrics"]["genetic_code_id"] == 2
    assert payload["metrics"]["stop_codon_count"] == 0
