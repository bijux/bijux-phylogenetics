from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.comparative.signal import (
    compute_blombergs_k,
    compute_phylogenetic_independent_contrasts,
    compute_phylogenetic_signal_test,
    estimate_pagels_lambda,
)
from bijux_phylogenetics.errors import ComparativeMethodError

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


def test_phylogenetic_signal_test_reuses_seeded_permutation_path() -> None:
    left = compute_phylogenetic_signal_test(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        trait="response_growth",
        permutations=11,
        seed=17,
    )
    right = compute_phylogenetic_signal_test(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        trait="response_growth",
        permutations=11,
        seed=17,
    )
    different_seed = compute_phylogenetic_signal_test(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        trait="response_growth",
        permutations=11,
        seed=18,
    )
    assert left.seed == 17
    assert right.seed == 17
    assert left.permutation_rows == right.permutation_rows
    assert left.p_value == right.p_value
    assert left.permutation_rows != different_seed.permutation_rows


def test_phylogenetic_signal_reports_pruned_missing_values_explicitly() -> None:
    report = compute_blombergs_k(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_brownian_missing.tsv"),
        trait="response_growth",
    )
    assert report.taxon_count == 4
    assert report.input_audit.missing_value_policy == "prune-overlapping-missing-values"
    assert report.input_audit.pruned_missing_value_taxa == ["B"]
    assert (
        "one or more overlapping taxa have missing trait values and will be pruned"
        in report.input_audit.warnings
    )


def test_phylogenetic_signal_accepts_rooted_non_ultrametric_tree_by_policy() -> None:
    report = estimate_pagels_lambda(
        fixture("example_tree_internal_long_branch.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    assert report.input_audit.tree_is_ultrametric is False
    assert report.input_audit.ultrametric_policy == (
        "accept-rooted-trees-and-report-ultrametricity"
    )
    assert report.input_audit.minimum_root_to_tip_depth == 0.2
    assert report.input_audit.maximum_root_to_tip_depth == 1.1


def test_phylogenetic_signal_rejects_constant_trait_values(tmp_path: Path) -> None:
    traits_path = tmp_path / "constant-traits.tsv"
    traits_path.write_text(
        "taxon\tresponse\nA\t2.0\nB\t2.0\nC\t2.0\nD\t2.0\n",
        encoding="utf-8",
    )
    with pytest.raises(ComparativeMethodError) as error:
        compute_phylogenetic_signal_test(
            fixture("example_tree.nwk"),
            traits_path,
            trait="response",
            permutations=7,
            seed=3,
        )
    assert str(error.value) == (
        "phylogenetic signal requires at least two distinct numeric trait values after pruning"
    )
