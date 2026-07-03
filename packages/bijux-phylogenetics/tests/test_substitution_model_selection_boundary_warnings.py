from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    compare_nucleotide_substitution_models_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_substitution_model_selection_rows_emit_structured_boundary_warnings() -> None:
    report = compare_nucleotide_substitution_models_from_alignment(
        fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta"),
        candidate_models=("JC69+G", "K80+I"),
        max_coordinate_passes=2,
    )

    rows = {row.model_name: row for row in report.rows}

    jc69_gamma_row = rows["JC69+G"]
    assert jc69_gamma_row.warnings == ["alpha hit upper search boundary"]
    assert [
        warning.warning_kind for warning in jc69_gamma_row.boundary_warning_rows
    ] == ["flat-alpha"]
    assert [
        warning.boundary_side for warning in jc69_gamma_row.boundary_warning_rows
    ] == ["upper"]

    k80_invariant_row = rows["K80+I"]
    assert k80_invariant_row.warnings == [
        "invariant_proportion hit lower search boundary"
    ]
    assert [
        warning.warning_kind for warning in k80_invariant_row.boundary_warning_rows
    ] == ["invariant-boundary"]
    assert [
        warning.boundary_side for warning in k80_invariant_row.boundary_warning_rows
    ] == ["lower"]
