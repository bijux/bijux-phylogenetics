from __future__ import annotations

from pathlib import Path

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.parsimony import (
    FitchCharacterMatrix,
    ParsimonyBremerSupportReport,
    ParsimonyBremerSupportRow,
    compute_parsimony_bremer_support,
)

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_bremer_support_surface() -> None:
    assert parsimony_api.ParsimonyBremerSupportRow is ParsimonyBremerSupportRow
    assert parsimony_api.ParsimonyBremerSupportReport is ParsimonyBremerSupportReport
    assert (
        parsimony_api.compute_parsimony_bremer_support
        is compute_parsimony_bremer_support
    )


def test_parsimony_bremer_support_reports_positive_decay_index() -> None:
    matrix = FitchCharacterMatrix(
        matrix_path=None,
        taxon_column="taxon",
        character_ids=["char01_terminal_a", "char02_clade_bd"],
        states_by_taxon={
            "A": {"char01_terminal_a": "1", "char02_clade_bd": "0"},
            "B": {"char01_terminal_a": "0", "char02_clade_bd": "1"},
            "C": {"char01_terminal_a": "0", "char02_clade_bd": "0"},
            "D": {"char01_terminal_a": "0", "char02_clade_bd": "1"},
        },
    )

    report = compute_parsimony_bremer_support(
        loads_newick("((A,(B,D)),C);"),
        matrix,
        method="dollo",
    )

    assert report.algorithm == "parsimony-bremer-support"
    assert report.method == "dollo"
    assert report.candidate_tree_count == 15
    assert report.reference_tree_newick == "((A,(B,D)),C);"
    assert report.reference_tree_score == 2.0
    assert report.optimal_score == 2.0
    assert report.optimal_tree_count == 3
    assert report.optimal_tree_newick == "((A,(B,D)),C);"
    assert report.reference_tree_score_delta_from_optimal == 0.0
    assert report.reference_tree_is_optimal is True
    assert [
        (
            row.branch_id,
            row.descendant_taxa,
            row.shortest_lacking_score,
            row.decay_index,
            row.shortest_lacking_tree_count,
            row.shortest_lacking_tree_newick,
        )
        for row in report.bremer_rows
    ] == [
        ("B|D", ["B", "D"], 3.0, 1.0, 6, "(((A,B),D),C);"),
        ("A|B|D", ["A", "B", "D"], 2.0, 0.0, 2, "((A,C),(B,D));"),
    ]
