from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    ParsimonyPlacementAlternativeRow,
    ParsimonyPlacementQuerySummary,
    ParsimonyPlacementReport,
    place_parsimony_queries,
    write_parsimony_placement_artifacts,
)
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_placement_surface() -> None:
    assert (
        parsimony_api.ParsimonyPlacementAlternativeRow
        is ParsimonyPlacementAlternativeRow
    )
    assert (
        parsimony_api.ParsimonyPlacementQuerySummary is ParsimonyPlacementQuerySummary
    )
    assert parsimony_api.ParsimonyPlacementReport is ParsimonyPlacementReport
    assert parsimony_api.place_parsimony_queries is place_parsimony_queries
    assert (
        parsimony_api.write_parsimony_placement_artifacts
        is write_parsimony_placement_artifacts
    )


def test_parsimony_placement_reports_equal_best_and_unique_best_queries() -> None:
    report = place_parsimony_queries(
        fixture("placement_reference_tree_4_taxa.nwk"),
        fixture("placement_reference_matrix.tsv"),
        fixture("placement_query_matrix.tsv"),
    )

    assert report.algorithm == "parsimony-placement"
    assert report.method == "unordered-fitch"
    assert report.reference_taxon_count == 4
    assert report.character_count == 2
    assert report.edge_count == 6
    assert report.query_count == 2
    assert report.reference_total_steps == 2
    assert report.reference_total_weighted_score == 2.0

    summaries_by_query = {row.query_id: row for row in report.query_summaries}
    assert set(summaries_by_query) == {"Q_C", "Q_TIE"}

    q_c = summaries_by_query["Q_C"]
    assert q_c.best_child_name == "C"
    assert q_c.best_descendant_taxa == ["C"]
    assert q_c.best_additional_steps == 0
    assert q_c.equally_best_placement_count == 1
    assert "Q_C" in q_c.selected_best_tree_newick

    q_tie = summaries_by_query["Q_TIE"]
    assert q_tie.best_child_name is None
    assert q_tie.best_descendant_taxa == ["A", "B"]
    assert q_tie.best_additional_steps == 0
    assert q_tie.equally_best_placement_count == 4
    assert "Q_TIE" in q_tie.selected_best_tree_newick

    rows_by_query: dict[str, list[ParsimonyPlacementAlternativeRow]] = {}
    for row in report.alternative_rows:
        rows_by_query.setdefault(row.query_id, []).append(row)

    assert [row.placement_rank for row in rows_by_query["Q_C"]] == [1, 2, 3, 4, 5, 6]
    assert [row.placement_rank for row in rows_by_query["Q_TIE"]] == [1, 2, 3, 4, 5, 6]
    assert [
        (row.child_name, row.descendant_taxa)
        for row in rows_by_query["Q_TIE"]
        if row.is_equally_best
    ] == [
        (None, ["A", "B"]),
        ("A", ["A"]),
        ("B", ["B"]),
        (None, ["C", "D"]),
    ]


def test_parsimony_placement_rejects_overlapping_query_taxa(tmp_path: Path) -> None:
    query_path = tmp_path / "overlapping_query.tsv"
    query_path.write_text(
        "taxon\tchar01_split\tchar02_cd_conflict\nA\t0\t0\n",
        encoding="utf-8",
    )

    with pytest.raises(ParsimonyAnalysisError) as error:
        place_parsimony_queries(
            fixture("placement_reference_tree_4_taxa.nwk"),
            fixture("placement_reference_matrix.tsv"),
            query_path,
        )

    assert error.value.code == "parsimony_placement_query_taxon_overlap"
    assert error.value.details["overlapping_taxa"] == ["A"]


def test_parsimony_placement_requires_matching_query_characters(
    tmp_path: Path,
) -> None:
    query_path = tmp_path / "mismatched_query.tsv"
    query_path.write_text(
        "taxon\tchar01_split\tchar99_other\nQ_BAD\t0\t0\n",
        encoding="utf-8",
    )

    with pytest.raises(ParsimonyAnalysisError) as error:
        place_parsimony_queries(
            fixture("placement_reference_tree_4_taxa.nwk"),
            fixture("placement_reference_matrix.tsv"),
            query_path,
        )

    assert error.value.code == "parsimony_placement_query_character_mismatch"
    assert error.value.details["reference_character_ids"] == [
        "char01_split",
        "char02_cd_conflict",
    ]
