from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def test_alignment_concatenate_cli_writes_supermatrix_outputs_and_manifest(
    tmp_path: Path,
    capsys,
) -> None:
    input_dir = fixture("concatenation/mixed-locus-supermatrix/inputs")
    alignment_out = tmp_path / "mixed-locus-supermatrix.aln.fasta"
    partitions_out = tmp_path / "mixed-locus-supermatrix.partitions.txt"
    matrix_out = tmp_path / "mixed-locus-supermatrix.matrix.tsv"
    taxa_out = tmp_path / "mixed-locus-supermatrix.taxa.tsv"
    loci_out = tmp_path / "mixed-locus-supermatrix.loci.tsv"
    manifest_path = tmp_path / "mixed-locus-supermatrix.manifest.json"

    exit_code = main(
        [
            "alignment",
            "concatenate",
            str(input_dir / "alpha-dna.fasta"),
            str(input_dir / "beta-protein.fasta"),
            str(input_dir / "gamma-dna.fasta"),
            "--data-type",
            "DNA",
            "--data-type",
            "PROTEIN",
            "--data-type",
            "DNA",
            "--out",
            str(alignment_out),
            "--partitions-out",
            str(partitions_out),
            "--matrix-out",
            str(matrix_out),
            "--taxa-out",
            str(taxa_out),
            "--loci-out",
            str(loci_out),
            "--manifest",
            str(manifest_path),
            "--json",
        ]
    )

    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"] == {
        "alignment_length": 11,
        "locus_count": 3,
        "taxon_count": 5,
    }
    assert alignment_out.read_text(encoding="utf-8") == fixture(
        "expected/concatenation/mixed-locus-supermatrix.aln.fasta"
    ).read_text(encoding="utf-8")
    assert partitions_out.read_text(encoding="utf-8") == fixture(
        "expected/concatenation/mixed-locus-supermatrix.partitions.txt"
    ).read_text(encoding="utf-8")
    assert matrix_out.read_text(encoding="utf-8") == fixture(
        "expected/concatenation/mixed-locus-supermatrix.matrix.tsv"
    ).read_text(encoding="utf-8")
    assert taxa_out.read_text(encoding="utf-8").startswith(
        "taxon\tcovered_locus_count\t"
    )
    assert loci_out.read_text(encoding="utf-8").startswith(
        "locus_name\tcovered_taxon_count\t"
    )
    assert manifest_path.exists()


def test_alignment_concatenate_cli_honors_explicit_locus_names(
    tmp_path: Path,
) -> None:
    input_dir = fixture("concatenation/mixed-locus-supermatrix/inputs")
    partitions_out = tmp_path / "renamed.partitions.txt"

    exit_code = main(
        [
            "alignment",
            "concatenate",
            str(input_dir / "alpha-dna.fasta"),
            str(input_dir / "gamma-dna.fasta"),
            "--locus-name",
            "first-marker",
            "--locus-name",
            "second-marker",
            "--out",
            str(tmp_path / "renamed.aln.fasta"),
            "--partitions-out",
            str(partitions_out),
            "--matrix-out",
            str(tmp_path / "renamed.matrix.tsv"),
            "--json",
        ]
    )

    assert exit_code == 0
    assert partitions_out.read_text(encoding="utf-8") == (
        "DNA,first-marker = 1-4\nDNA,second-marker = 5-7\n"
    )
