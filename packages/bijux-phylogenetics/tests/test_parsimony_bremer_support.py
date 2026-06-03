from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import loads_newick
import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    FitchCharacterMatrix,
    ParsimonyBremerSupportReport,
    ParsimonyBremerSupportRow,
    compute_parsimony_bremer_support,
    write_parsimony_bremer_support_artifacts,
)
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

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
    assert (
        parsimony_api.write_parsimony_bremer_support_artifacts
        is write_parsimony_bremer_support_artifacts
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


def test_write_parsimony_bremer_support_artifacts_materializes_outputs(
    tmp_path: Path,
) -> None:
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

    outputs = write_parsimony_bremer_support_artifacts(
        tmp_path / "bremer-support",
        report,
    )

    assert set(outputs) == {
        "reference_tree_path",
        "optimal_tree_path",
        "bremer_support_path",
        "run_json_path",
    }
    assert (
        outputs["bremer_support_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "branch_id\tnode_name\tdescendant_taxa\tshortest_lacking_score\tdecay_index\tshortest_lacking_tree_count\tshortest_lacking_tree_newick\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "parsimony-bremer-support"
    assert payload["optimal_score"] == 2.0
    assert payload["bremer_rows"][0]["branch_id"] == "B|D"
    assert payload["bremer_rows"][0]["decay_index"] == 1.0


def test_parsimony_bremer_support_keeps_zero_decay_when_excluding_tree_is_still_optimal() -> (
    None
):
    report = compute_parsimony_bremer_support(
        loads_newick("(((A,B),C),D);"),
        fixture("bootstrap_matrix.tsv"),
        method="fitch",
    )

    assert report.reference_tree_score == 5.0
    assert report.optimal_score == 5.0
    assert report.optimal_tree_count == 5
    assert [
        (
            row.branch_id,
            row.shortest_lacking_score,
            row.decay_index,
            row.shortest_lacking_tree_count,
            row.shortest_lacking_tree_newick,
        )
        for row in report.bremer_rows
    ] == [
        ("A|B", 5.0, 0.0, 2, "((A,(C,D)),B);"),
        ("A|B|C", 5.0, 0.0, 4, "(((A,B),D),C);"),
    ]


def test_parsimony_bremer_support_rejects_nonbinary_reference_tree() -> None:
    with pytest.raises(ParsimonyAnalysisError, match="rooted binary tree"):
        compute_parsimony_bremer_support(
            loads_newick("(A,B,C,D);"),
            fixture("bootstrap_matrix.tsv"),
            method="fitch",
        )
