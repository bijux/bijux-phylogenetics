from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative import (
    rescale_tree_early_burst,
    rescale_tree_pagel_delta,
    rescale_tree_pagel_kappa,
    rescale_tree_pagel_lambda,
    rescale_tree_white_noise,
)
from bijux_phylogenetics.comparative.common import build_brownian_covariance_matrix
from bijux_phylogenetics.io.newick import loads_newick
from tests.support.geiger_rescale_reference import (
    GEIGER_RESCALE_REFERENCE_PAYLOADS,
)

FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures"
EXAMPLE_TREE = FIXTURE_ROOT / "trees" / "example_tree.nwk"


def _covariance_from_newick(newick: str) -> list[list[float]]:
    tree = loads_newick(newick)
    return build_brownian_covariance_matrix(tree, tree.tip_names)


def _sorted_branch_lengths_from_newick(newick: str) -> list[float]:
    tree = loads_newick(newick)
    branch_lengths: list[float] = []

    def visit(node) -> None:
        for child in node.children:
            branch_lengths.append(float(child.branch_length or 0.0))
            visit(child)

    visit(tree.root)
    return sorted(branch_lengths)


def _assert_matrix_close(
    observed: list[list[float]],
    expected: list[list[float]],
) -> None:
    assert len(observed) == len(expected)
    for observed_row, expected_row in zip(observed, expected, strict=True):
        assert len(observed_row) == len(expected_row)
        for observed_value, expected_value in zip(
            observed_row,
            expected_row,
            strict=True,
        ):
            assert math.isclose(
                observed_value,
                expected_value,
                rel_tol=0.0,
                abs_tol=1e-9,
            )


def test_rescale_tree_pagel_lambda_matches_governed_geiger_reference() -> None:
    report = rescale_tree_pagel_lambda(EXAMPLE_TREE, lambda_value=0.5)
    reference = GEIGER_RESCALE_REFERENCE_PAYLOADS["pagel_lambda_half"]

    for observed, expected in zip(
        _sorted_branch_lengths_from_newick(report.transformed_tree_newick),
        sorted(reference["edge_lengths"]),
        strict=True,
    ):
        assert math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-9)
    _assert_matrix_close(
        _covariance_from_newick(report.transformed_tree_newick),
        reference["covariance"],
    )


def test_rescale_tree_pagel_kappa_matches_governed_geiger_reference() -> None:
    report = rescale_tree_pagel_kappa(EXAMPLE_TREE, kappa=0.5)
    reference = GEIGER_RESCALE_REFERENCE_PAYLOADS["pagel_kappa_half"]

    for observed, expected in zip(
        _sorted_branch_lengths_from_newick(report.transformed_tree_newick),
        sorted(reference["edge_lengths"]),
        strict=True,
    ):
        assert math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-9)
    _assert_matrix_close(
        _covariance_from_newick(report.transformed_tree_newick),
        reference["covariance"],
    )


def test_rescale_tree_pagel_delta_matches_governed_geiger_reference() -> None:
    report = rescale_tree_pagel_delta(EXAMPLE_TREE, delta=0.5)
    reference = GEIGER_RESCALE_REFERENCE_PAYLOADS["pagel_delta_half"]

    for observed, expected in zip(
        _sorted_branch_lengths_from_newick(report.transformed_tree_newick),
        sorted(reference["edge_lengths"]),
        strict=True,
    ):
        assert math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-9)
    _assert_matrix_close(
        _covariance_from_newick(report.transformed_tree_newick),
        reference["covariance"],
    )


def test_rescale_tree_early_burst_matches_governed_geiger_reference() -> None:
    report = rescale_tree_early_burst(EXAMPLE_TREE, rate_change=2.0)
    reference = GEIGER_RESCALE_REFERENCE_PAYLOADS["early_burst_negative_two"]

    for observed, expected in zip(
        _sorted_branch_lengths_from_newick(report.transformed_tree_newick),
        sorted(reference["edge_lengths"]),
        strict=True,
    ):
        assert math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-9)
    _assert_matrix_close(
        _covariance_from_newick(report.transformed_tree_newick),
        reference["covariance"],
    )


def test_rescale_tree_white_noise_matches_governed_geiger_reference() -> None:
    report = rescale_tree_white_noise(EXAMPLE_TREE, sigsq=2.5)
    reference = GEIGER_RESCALE_REFERENCE_PAYLOADS["white_noise_two_point_five"]

    for observed, expected in zip(
        _sorted_branch_lengths_from_newick(report.transformed_tree_newick),
        sorted(reference["edge_lengths"]),
        strict=True,
    ):
        assert math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-9)
    _assert_matrix_close(
        _covariance_from_newick(report.transformed_tree_newick),
        reference["covariance"],
    )
