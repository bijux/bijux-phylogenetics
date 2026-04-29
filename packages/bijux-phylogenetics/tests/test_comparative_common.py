from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.comparative.common import (
    build_brownian_covariance_matrix,
    load_comparative_dataset,
    summarize_numeric_trait,
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.errors import ComparativeMethodError
from bijux_phylogenetics.io.trees import load_tree


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


def test_comparative_readiness_reports_numeric_trait_pruning() -> None:
    report = summarize_numeric_trait_readiness(
        fixture("example_tree.nwk"),
        fixture("example_traits_validate.tsv"),
        trait="height_cm",
    )
    assert report.ready is True
    assert report.analysis_taxa == ["A", "B", "C", "D"]
    assert report.blockers == []
    assert report.warnings == []


def test_numeric_trait_summary_uses_phylogenetic_overlap_taxa() -> None:
    summary = summarize_numeric_trait(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    assert summary.taxa == ["A", "B", "C", "D"]
    assert math.isclose(summary.mean, 2.75)
    assert math.isclose(summary.variance, 1.0833333333333333)
    assert math.isclose(summary.minimum, 1.5)
    assert math.isclose(summary.maximum, 4.0)


def test_brownian_covariance_matrix_matches_shared_path_lengths() -> None:
    tree = load_tree(fixture("example_tree.nwk"))
    matrix = build_brownian_covariance_matrix(tree, ["A", "B", "C", "D"])
    expected = [
        [0.3, 0.2, 0.0, 0.0],
        [0.2, 0.3, 0.0, 0.0],
        [0.0, 0.0, 0.3, 0.1],
        [0.0, 0.0, 0.1, 0.3],
    ]
    for row, expected_row in zip(matrix, expected, strict=True):
        for value, expected_value in zip(row, expected_row, strict=True):
            assert math.isclose(value, expected_value)


def test_load_comparative_dataset_requires_binary_tree_for_contrasts() -> None:
    with pytest.raises(ComparativeMethodError):
        load_comparative_dataset(
            fixture("example_tree_star.nwk"),
            fixture("example_traits_comparative.tsv"),
            trait="response",
            require_binary=True,
        )
