from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.comparative.common import load_comparative_dataset
from bijux_phylogenetics.comparative.signal import (
    compute_blombergs_k,
    compute_phylogenetic_independent_contrasts,
    compute_phylogenetic_signal_test,
    estimate_pagels_lambda,
    evaluate_pagels_lambda_likelihood,
    evaluate_pagels_lambda_likelihood_from_dataset,
)
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

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
    assert report.input_audit.missing_value_policy == "prune-overlapping-missing-values"
    assert report.input_audit.pruned_missing_value_taxa == []
    assert report.input_audit.tree_is_ultrametric is True
    assert len(report.contrasts) == 3
    assert report.contrasts[0].node_id == 6
    assert report.contrasts[0].node == "A|B"
    assert math.isclose(report.contrasts[0].contrast, -3.3541019662496847)
    assert report.contrasts[1].node_id == 7
    assert report.contrasts[1].node == "C|D"
    assert math.isclose(report.contrasts[1].contrast, -2.3717082451262845)
    assert report.contrasts[2].node_id == 5
    assert report.contrasts[2].node == "A|B|C|D"
    assert math.isclose(report.root_estimate, 2.8055555555555554)


def test_independent_contrasts_accept_pectinate_non_ultrametric_tree_by_policy() -> (
    None
):
    report = compute_phylogenetic_independent_contrasts(
        fixture("example_tree_ladderized.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )

    assert report.input_audit.tree_is_ultrametric is False
    assert report.input_audit.ultrametric_policy == (
        "accept-rooted-trees-and-report-ultrametricity"
    )
    assert len(report.contrasts) == 3
    assert [row.node_id for row in report.contrasts] == [7, 6, 5]
    assert math.isclose(report.contrasts[1].contrast, -0.5)
    assert math.isclose(report.contrasts[2].expected_variance, 0.26)
    assert math.isclose(report.root_estimate, 3.3846153846153846)


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
    assert (
        report.null_distribution_minimum
        <= report.null_distribution_mean
        <= report.null_distribution_maximum
    )


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
    assert left.null_distribution_mean == right.null_distribution_mean
    assert left.permutation_rows != different_seed.permutation_rows


@pytest.mark.slow
def test_phylogenetic_signal_test_distinguishes_strong_and_weak_signal_permutation_summaries() -> (
    None
):
    strong = compute_phylogenetic_signal_test(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="signal_strong",
        permutations=199,
        seed=17,
    )
    weak = compute_phylogenetic_signal_test(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="signal_weak",
        permutations=199,
        seed=17,
    )

    assert strong.observed_k > weak.observed_k
    assert strong.p_value < weak.p_value
    assert strong.null_distribution_mean > weak.null_distribution_mean


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


@pytest.mark.parametrize(
    ("runner", "extra_kwargs"),
    [
        (
            compute_blombergs_k,
            {},
        ),
        (
            estimate_pagels_lambda,
            {},
        ),
        (
            compute_phylogenetic_signal_test,
            {"permutations": 9, "seed": 17},
        ),
    ],
)
def test_signal_methods_report_pruned_missing_values_explicitly(
    runner,
    extra_kwargs: dict[str, int],
) -> None:
    report = runner(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_brownian_missing.tsv"),
        trait="response_growth",
        **extra_kwargs,
    )

    assert report.taxon_count == 4
    assert report.input_audit.missing_value_policy == (
        "prune-overlapping-missing-values"
    )
    assert report.input_audit.pruned_missing_value_taxa == ["B"]
    assert (
        "one or more overlapping taxa have missing trait values and will be pruned"
        in report.input_audit.warnings
    )


def test_independent_contrasts_reports_pruned_missing_values_explicitly() -> None:
    report = compute_phylogenetic_independent_contrasts(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_brownian_missing.tsv"),
        trait="response_growth",
    )

    assert report.taxon_count == 4
    assert report.input_audit.missing_value_policy == "prune-overlapping-missing-values"
    assert report.input_audit.pruned_missing_value_taxa == ["B"]
    assert report.input_audit.warnings == [
        "trait table contains taxa absent from the tree",
        "one or more overlapping taxa have missing trait values and will be pruned",
        "one or more overlapping taxa have non-numeric trait values and will be pruned",
    ]


def test_independent_contrasts_reject_negative_branch_lengths() -> None:
    with pytest.raises(ComparativeMethodError) as error:
        compute_phylogenetic_independent_contrasts(
            fixture("example_tree_negative_length.nwk"),
            fixture("example_traits_three_taxa.tsv"),
            trait="response",
        )

    assert (
        error.value.details["failure_reason"] == "comparative_negative_branch_lengths"
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


@pytest.mark.parametrize(
    ("runner", "extra_kwargs"),
    [
        (
            compute_blombergs_k,
            {},
        ),
        (
            compute_phylogenetic_signal_test,
            {"permutations": 9, "seed": 5},
        ),
    ],
)
def test_signal_methods_accept_rooted_non_ultrametric_tree_by_policy(
    runner,
    extra_kwargs: dict[str, int],
) -> None:
    report = runner(
        fixture("example_tree_internal_long_branch.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        **extra_kwargs,
    )

    assert report.input_audit.tree_is_ultrametric is False
    assert report.input_audit.ultrametric_policy == (
        "accept-rooted-trees-and-report-ultrametricity"
    )
    assert report.input_audit.minimum_root_to_tip_depth == 0.2
    assert report.input_audit.maximum_root_to_tip_depth == 1.1


def test_pagels_lambda_fixed_likelihood_surface_matches_boundary_reports() -> None:
    tree = fixture("example_tree_phytools_non_ultrametric_twenty_four_taxa.nwk")
    traits = fixture(
        "example_traits_phytools_signal_non_ultrametric_twenty_four_taxa.tsv"
    )

    estimate = estimate_pagels_lambda(tree, traits, trait="signal_strong")
    null_report = evaluate_pagels_lambda_likelihood(
        tree,
        traits,
        trait="signal_strong",
        lambda_value=0.0,
    )
    brownian_report = evaluate_pagels_lambda_likelihood(
        tree,
        traits,
        trait="signal_strong",
        lambda_value=1.0,
    )

    assert math.isclose(null_report.log_likelihood, estimate.null_log_likelihood)
    assert math.isclose(
        brownian_report.log_likelihood, estimate.brownian_log_likelihood
    )
    assert math.isclose(
        estimate.likelihood_ratio_statistic,
        2.0 * (estimate.log_likelihood - estimate.null_log_likelihood),
    )
    assert 0.0 <= estimate.likelihood_ratio_p_value <= 1.0


def test_pagels_lambda_fixed_likelihood_from_dataset_matches_path_surface() -> None:
    tree = fixture("example_tree_phytools_non_ultrametric_twenty_four_taxa.nwk")
    traits = fixture(
        "example_traits_phytools_signal_non_ultrametric_twenty_four_taxa.tsv"
    )
    dataset = load_comparative_dataset(tree, traits, trait="signal_strong")

    path_report = evaluate_pagels_lambda_likelihood(
        tree,
        traits,
        trait="signal_strong",
        lambda_value=0.35,
    )
    dataset_report = evaluate_pagels_lambda_likelihood_from_dataset(
        dataset,
        lambda_value=0.35,
    )

    assert math.isclose(dataset_report.log_likelihood, path_report.log_likelihood)
    assert dataset_report.lambda_value == path_report.lambda_value
    assert dataset_report.taxon_count == path_report.taxon_count
    assert dataset_report.trait == path_report.trait


def test_pagels_lambda_reports_optimizer_diagnostics_for_strong_and_weak_signal() -> (
    None
):
    strong_report = estimate_pagels_lambda(
        fixture("example_tree_phytools_non_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_non_ultrametric_twenty_four_taxa.tsv"),
        trait="signal_strong",
    )
    weak_report = estimate_pagels_lambda(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="signal_weak",
    )

    strong_diagnostics = strong_report.optimizer_diagnostics
    weak_diagnostics = weak_report.optimizer_diagnostics

    assert strong_diagnostics.optimizer_name == "two-stage-grid-search"
    assert strong_diagnostics.coarse_grid_point_count == 21
    assert strong_diagnostics.fine_grid_point_count == 11
    assert (
        strong_diagnostics.function_evaluation_count
        == strong_diagnostics.coarse_grid_point_count
        + strong_diagnostics.fine_grid_point_count
    )
    assert strong_diagnostics.hit_upper_boundary is True
    assert weak_diagnostics.hit_lower_boundary is True
    assert len(strong_report.profile_rows) == strong_diagnostics.fine_grid_point_count
    assert any(row.within_95_confidence_interval for row in weak_report.profile_rows)
    assert all(row.delta_log_likelihood >= 0.0 for row in strong_report.profile_rows)


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


@pytest.mark.parametrize(
    ("runner", "extra_kwargs"),
    [
        (
            compute_blombergs_k,
            {},
        ),
        (
            estimate_pagels_lambda,
            {},
        ),
    ],
)
def test_signal_estimates_reject_constant_trait_values(
    tmp_path: Path,
    runner,
    extra_kwargs: dict[str, int],
) -> None:
    traits_path = tmp_path / "constant-traits.tsv"
    traits_path.write_text(
        "taxon\tresponse\nA\t2.0\nB\t2.0\nC\t2.0\nD\t2.0\n",
        encoding="utf-8",
    )
    with pytest.raises(ComparativeMethodError) as error:
        runner(
            fixture("example_tree.nwk"),
            traits_path,
            trait="response",
            **extra_kwargs,
        )

    assert str(error.value) == (
        "phylogenetic signal requires at least two distinct numeric trait values after pruning"
    )
