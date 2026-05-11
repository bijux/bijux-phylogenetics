from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines.fasta_to_tree import (
    FastaToTreeModelRow,
    FastaToTreeSupportRow,
    infer_unaligned_sequence_type,
    write_fasta_to_tree_model_table,
    write_fasta_to_tree_support_table,
)


def test_infer_unaligned_sequence_type_detects_dna_and_protein_inputs() -> None:
    assert (
        infer_unaligned_sequence_type(
            [("A", "ACTGACTG"), ("B", "ACTGACTA"), ("C", "ACTGACGG")]
        )
        == "dna"
    )
    assert (
        infer_unaligned_sequence_type(
            [("P1", "MKTWFLIM"), ("P2", "MKTWYLIM"), ("P3", "MRTWFLVM")]
        )
        == "protein"
    )


def test_write_fasta_to_tree_tables_emits_expected_tsv(tmp_path: Path) -> None:
    model_path = tmp_path / "example.model.tsv"
    support_path = tmp_path / "example.support.tsv"

    write_fasta_to_tree_model_table(
        model_path,
        [
            FastaToTreeModelRow(
                workflow="model-selection",
                engine_name="iqtree",
                sequence_type="dna",
                selected_model="GTR+G",
                report_selected_model="GTR+G",
                artifact_selected_model="GTR+G",
                model_consistent=True,
                alignment_path=Path("example.aln"),
                trimmed_alignment_path=Path("example.trimmed.aln"),
                manifest_path=Path("engine-artifacts/example/model.manifest.json"),
            )
        ],
    )
    write_fasta_to_tree_support_table(
        support_path,
        [
            FastaToTreeSupportRow(
                node="A|B",
                descendant_taxa=("A", "B"),
                support=95.0,
                support_fraction=0.95,
                is_backbone=True,
            )
        ],
    )

    assert model_path.read_text(encoding="utf-8") == (
        "workflow\tengine_name\tsequence_type\tselected_model\t"
        "report_selected_model\tartifact_selected_model\tmodel_consistent\t"
        "alignment_path\ttrimmed_alignment_path\tmanifest_path\n"
        "model-selection\tiqtree\tdna\tGTR+G\tGTR+G\tGTR+G\ttrue\t"
        "example.aln\texample.trimmed.aln\t"
        "engine-artifacts/example/model.manifest.json\n"
    )
    assert support_path.read_text(encoding="utf-8") == (
        "node\tdescendant_taxa\tsupport\tsupport_fraction\tis_backbone\n"
        "A|B\tA,B\t95\t0.95\ttrue\n"
    )
