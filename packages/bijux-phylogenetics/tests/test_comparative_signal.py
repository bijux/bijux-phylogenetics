from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.signal import (
    compute_blombergs_k,
    compute_phylogenetic_independent_contrasts,
    compute_phylogenetic_signal_test,
    estimate_pagels_lambda,
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


def test_independent_contrasts_return_expected_internal_node_values() -> None:
    report = compute_phylogenetic_independent_contrasts(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    assert report.taxon_count == 4
    assert len(report.contrasts) == 3
    assert report.contrasts[0].node == "A|B"
    assert math.isclose(report.contrasts[0].contrast, -3.3541019662496847)
    assert report.contrasts[1].node == "C|D"
    assert math.isclose(report.contrasts[1].contrast, -2.3717082451262845)
    assert report.contrasts[2].node == "A|B|C|D"
    assert math.isclose(report.root_estimate, 2.8055555555555554)


def test_blombergs_k_and_pagels_lambda_return_stable_positive_signal() -> None:
    k_report = compute_blombergs_k(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    lambda_report = estimate_pagels_lambda(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    assert k_report.k > 0.0
    assert 0.0 <= lambda_report.lambda_value <= 1.0
    assert lambda_report.log_likelihood >= lambda_report.null_log_likelihood


def test_phylogenetic_signal_test_returns_permutation_p_value() -> None:
    report = compute_phylogenetic_signal_test(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        permutations=19,
        seed=7,
    )
    assert report.permutations == 19
    assert 0.0 < report.p_value <= 1.0
    assert 0.0 <= report.estimated_lambda <= 1.0
