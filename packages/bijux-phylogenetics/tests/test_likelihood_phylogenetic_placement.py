from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.io.fasta.core import write_fasta_alignment
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    place_queries_by_likelihood_from_alignment,
)
from bijux_phylogenetics.runtime.errors import (
    AlignmentTaxonMismatchError,
    InvalidAlignmentError,
)

pytestmark = pytest.mark.slow

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_package_likelihood_gateway_exports_placement_surface() -> None:
    assert (
        likelihood_api.place_queries_by_likelihood_from_alignment
        is place_queries_by_likelihood_from_alignment
    )


def test_likelihood_placement_reports_best_edges_weights_and_alternative_rows() -> None:
    report = place_queries_by_likelihood_from_alignment(
        fixture("trees", "likelihood_placement_reference_tree_4_taxa.nwk"),
        fixture("alignments", "likelihood_placement_reference_alignment_4_taxa.fasta"),
        fixture("alignments", "likelihood_placement_query_alignment_2_taxa.fasta"),
    )

    assert report.model_name == "JC69"
    assert report.reference_taxa == ["A", "B", "C", "D"]
    assert report.edge_count == 6
    assert report.query_count == 2
    assert report.site_count == 12
    assert report.max_coordinate_passes == 12

    summaries_by_query = {row.query_id: row for row in report.query_summaries}
    assert set(summaries_by_query) == {"Q_A", "Q_C"}
    assert summaries_by_query["Q_A"].best_child_name == "A"
    assert summaries_by_query["Q_A"].best_descendant_taxa == ["A"]
    assert summaries_by_query["Q_C"].best_child_name == "C"
    assert summaries_by_query["Q_C"].best_descendant_taxa == ["C"]
    assert summaries_by_query["Q_A"].candidate_placement_count == 6
    assert summaries_by_query["Q_C"].candidate_placement_count == 6
    assert "Q_A" in summaries_by_query["Q_A"].best_tree_newick
    assert "Q_C" in summaries_by_query["Q_C"].best_tree_newick

    rows_by_query: dict[str, list[object]] = {}
    for row in report.alternative_placements:
        rows_by_query.setdefault(row.query_id, []).append(row)
    assert set(rows_by_query) == {"Q_A", "Q_C"}

    for query_id, rows in rows_by_query.items():
        assert len(rows) == 6
        assert [row.placement_rank for row in rows] == [1, 2, 3, 4, 5, 6]
        assert math.isclose(
            sum(row.likelihood_weight_ratio for row in rows),
            1.0,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        assert rows[0].log_likelihood >= rows[1].log_likelihood
        assert math.isclose(
            rows[0].optimized_proximal_length + rows[0].optimized_distal_length,
            rows[0].original_branch_length,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        assert rows[0].optimized_pendant_length >= report.lower_pendant_length_bound
        assert rows[0].optimized_pendant_length <= report.upper_pendant_length_bound
        assert rows[0].child_name == summaries_by_query[query_id].best_child_name


def test_likelihood_placement_rejects_overlapping_query_identifiers(
    tmp_path: Path,
) -> None:
    query_path = write_fasta_alignment(
        tmp_path / "overlapping-query.fasta",
        [AlignmentRecord(identifier="A", sequence="AAAAAAAAAAAA")],
    )

    with pytest.raises(AlignmentTaxonMismatchError) as error:
        place_queries_by_likelihood_from_alignment(
            fixture("trees", "likelihood_placement_reference_tree_4_taxa.nwk"),
            fixture(
                "alignments",
                "likelihood_placement_reference_alignment_4_taxa.fasta",
            ),
            query_path,
        )

    assert "do not overlap the reference tree taxa" in str(error.value)


def test_likelihood_placement_requires_matching_query_alignment_length(
    tmp_path: Path,
) -> None:
    query_path = write_fasta_alignment(
        tmp_path / "short-query.fasta",
        [AlignmentRecord(identifier="Q_SHORT", sequence="AAAA")],
    )

    with pytest.raises(InvalidAlignmentError) as error:
        place_queries_by_likelihood_from_alignment(
            fixture("trees", "likelihood_placement_reference_tree_4_taxa.nwk"),
            fixture(
                "alignments",
                "likelihood_placement_reference_alignment_4_taxa.fasta",
            ),
            query_path,
        )

    assert "same aligned length" in str(error.value)
