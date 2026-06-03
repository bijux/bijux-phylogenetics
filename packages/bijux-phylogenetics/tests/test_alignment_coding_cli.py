from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.io.fasta import write_fasta_alignment
from bijux_phylogenetics.phylo.alignment import AlignmentRecord

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
    codon_table_path = tmp_path / "codon-validation.tsv"
    exclusion_table_path = tmp_path / "excluded-sequences.tsv"
    exit_code = main(
        [
            "alignment",
            "translate",
            str(fixture("example_alignment_coding.fasta")),
            "--out",
            str(output_path),
            "--codon-validation-out",
            str(codon_table_path),
            "--excluded-sequences-out",
            str(exclusion_table_path),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        ">A\nME*\n>B\nMEW\n>C\nMXW\n>D\nM*W\n"
    )
    assert codon_table_path.read_text(encoding="utf-8").splitlines()[:3] == [
        "identifier\tcodon_index\tnucleotide_start\tcodon\tamino_acid\ttranslation_status",
        "A\t1\t1\tATG\tM\ttranslated",
        "A\t2\t4\tGAA\tE\ttranslated",
    ]
    assert (
        exclusion_table_path.read_text(encoding="utf-8") == "identifier\treason\tnote\n"
    )
    assert payload["metrics"]["genetic_code_id"] == 1
    assert payload["metrics"]["translated_sequence_count"] == 4
    assert payload["metrics"]["invalid_codon_count"] == 1
    assert payload["metrics"]["stop_codon_count"] == 2
    assert payload["metrics"]["internal_stop_sequence_count"] == 1
    assert payload["metrics"]["terminal_stop_sequence_count"] == 1
    assert payload["metrics"]["trailing_partial_codon_sequence_count"] == 0
    assert payload["metrics"]["excluded_sequence_count"] == 0


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


def test_alignment_translate_cli_reports_trailing_partial_codon_warning(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "translated.fasta"

    exit_code = main(
        [
            "alignment",
            "translate",
            str(fixture("example_alignment_coding_frame_error.fasta")),
            "--out",
            str(output_path),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == ">frame_error\nME\n"
    assert payload["metrics"]["translated_alignment_length"] == 2
    assert payload["metrics"]["dropped_trailing_nucleotide_count"] == 2
    assert payload["metrics"]["trailing_partial_codon_sequence_count"] == 1
    assert payload["warnings"] == [
        "sequence length not a multiple of 3: 2 nucleotides dropped"
    ]
