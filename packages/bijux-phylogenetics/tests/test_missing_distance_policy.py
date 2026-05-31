from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.distance import (
    apply_missing_distance_policy,
    build_distance_tree,
    build_tree_from_imported_distance_matrix,
)
from bijux_phylogenetics.runtime.errors import (
    InvalidAlignmentError,
    InvalidDistanceMatrixError,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_missing_distance_policy_reports_distinct_resolution_results() -> None:
    identifiers = ["A", "B", "C", "D"]
    pair_distances = {
        ("A", "B"): 2.0,
        ("A", "D"): 5.0,
        ("B", "C"): 3.0,
        ("B", "D"): 11.0,
        ("C", "D"): 7.0,
    }

    with pytest.raises(
        InvalidDistanceMatrixError,
        match="missing-distance policy 'reject' blocks incomplete distance pairs",
    ):
        apply_missing_distance_policy(
            identifiers,
            pair_distances,
            policy="reject",
        )

    _, mean_report = apply_missing_distance_policy(
        identifiers,
        pair_distances,
        policy="mean-impute",
    )
    _, nearest_report = apply_missing_distance_policy(
        identifiers,
        pair_distances,
        policy="nearest-valid",
    )
    _, triangle_report = apply_missing_distance_policy(
        identifiers,
        pair_distances,
        policy="triangle-bound",
    )

    assert mean_report.imputed_rows[0].imputed_distance == 5.6
    assert nearest_report.imputed_rows[0].imputed_distance == 2.0
    assert triangle_report.imputed_rows[0].imputed_distance == 5.0


def test_build_tree_from_imported_distance_matrix_applies_missing_distance_policy() -> (
    None
):
    with pytest.raises(
        InvalidDistanceMatrixError,
        match="missing-distance policy 'reject' blocks incomplete distance pairs",
    ):
        build_tree_from_imported_distance_matrix(
            fixture("example_distance_matrix_missing_pair_four_taxon.tsv"),
            method="neighbor-joining",
            missing_distance_policy="reject",
        )

    _, report = build_tree_from_imported_distance_matrix(
        fixture("example_distance_matrix_missing_pair_four_taxon.tsv"),
        method="neighbor-joining",
        missing_distance_policy="triangle-bound",
    )

    assert report.missing_distance_policy_report.policy == "triangle-bound"
    assert report.missing_distance_policy_report.missing_pairs == ["A/C"]
    assert report.missing_distance_policy_report.imputed_rows[0].imputed_distance == 5.0


def test_build_distance_tree_applies_missing_distance_policy_to_alignment_pairs() -> (
    None
):
    with pytest.raises(
        InvalidAlignmentError,
        match="missing-distance policy 'reject' blocks incomplete distance pairs",
    ):
        build_distance_tree(
            fixture("example_alignment_distance_missing_pair.fasta"),
            method="neighbor-joining",
            model="p-distance",
            missing_distance_policy="reject",
        )

    _, mean_report = build_distance_tree(
        fixture("example_alignment_distance_missing_pair.fasta"),
        method="neighbor-joining",
        model="p-distance",
        missing_distance_policy="mean-impute",
    )
    _, nearest_report = build_distance_tree(
        fixture("example_alignment_distance_missing_pair.fasta"),
        method="neighbor-joining",
        model="p-distance",
        missing_distance_policy="nearest-valid",
    )
    _, triangle_report = build_distance_tree(
        fixture("example_alignment_distance_missing_pair.fasta"),
        method="neighbor-joining",
        model="p-distance",
        missing_distance_policy="triangle-bound",
    )

    assert (
        mean_report.missing_distance_policy_report.imputed_rows[0].imputed_distance
        == 0.75
    )
    assert (
        nearest_report.missing_distance_policy_report.imputed_rows[0].imputed_distance
        == 0.5
    )
    assert (
        triangle_report.missing_distance_policy_report.imputed_rows[0].imputed_distance
        == 1.5
    )
