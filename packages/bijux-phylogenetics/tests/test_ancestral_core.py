from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.ancestral.common import (
    load_continuous_dataset,
    load_discrete_dataset,
)
from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.errors import AncestralReconstructionError

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


def test_load_continuous_dataset_prunes_missing_taxa_and_warns_on_small_sample() -> (
    None
):
    dataset = load_continuous_dataset(
        fixture("example_tree.nwk"),
        fixture("example_traits.tsv"),
        trait="value",
    )
    assert dataset.taxa == ["A", "B", "C"]
    assert dataset.dropped_missing_taxa == ["D"]
    assert dataset.dropped_non_numeric_taxa == []
    assert dataset.warnings == [
        "continuous trait reconstruction is using only 3 taxa; results may be unstable",
        "one or more taxa were excluded because the continuous trait value was missing",
    ]


def test_continuous_reconstruction_reports_internal_estimates_and_intervals() -> None:
    report = reconstruct_continuous_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        model="brownian",
    )
    internal_estimates = {
        estimate.node: estimate for estimate in report.estimates if not estimate.is_tip
    }
    assert report.model == "brownian"
    assert report.taxon_count == 4
    assert math.isclose(internal_estimates["A|B"].estimate, 2.25)
    assert math.isclose(internal_estimates["C|D"].estimate, 3.25)
    assert math.isclose(internal_estimates["A|B|C|D"].estimate, 2.8055555555555554)
    assert (
        internal_estimates["A|B|C|D"].lower_95_interval
        < internal_estimates["A|B|C|D"].estimate
    )
    assert (
        internal_estimates["A|B|C|D"].upper_95_interval
        > internal_estimates["A|B|C|D"].estimate
    )
    assert report.brownian_fit_diagnostics is not None
    assert report.brownian_fit_diagnostics.covariance_model == "brownian-shared-path"
    assert report.brownian_fit_diagnostics.tree_is_ultrametric is True
    assert report.brownian_fit_diagnostics.covariance_matrix_dimension == 4
    assert report.brownian_fit_diagnostics.covariance_matrix_rank == 4
    assert report.brownian_fit_diagnostics.covariance_singular is False
    assert report.brownian_fit_diagnostics.covariance_positive_definite is True
    assert report.brownian_fit_diagnostics.covariance_condition_number > 0.0
    assert math.isfinite(report.brownian_fit_diagnostics.log_likelihood)
    assert report.brownian_fit_diagnostics.residual_sigma_squared > 0.0
    assert internal_estimates["A|B|C|D"].interpretation in {
        "moderate uncertainty",
        "narrow uncertainty",
        "broad uncertainty",
        "unstable node estimate",
    }
    assert 0.0 <= internal_estimates["A|B|C|D"].confidence <= 1.0
    assert isinstance(internal_estimates["A|B|C|D"].downstream_risks, list)


def test_continuous_reconstruction_supports_ou_model() -> None:
    report = reconstruct_continuous_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        model="ou",
        alpha=1.5,
    )
    root = next(estimate for estimate in report.estimates if estimate.node == "A|B|C|D")
    assert report.model == "ou"
    assert report.alpha == 1.5
    assert report.brownian_fit_diagnostics is None
    assert 2.5 < root.estimate < 3.1


def test_continuous_reconstruction_rejects_missing_branch_lengths() -> None:
    with pytest.raises(AncestralReconstructionError):
        reconstruct_continuous_ancestral_states(
            fixture("example_tree_no_lengths.nwk"),
            fixture("example_traits_comparative.tsv"),
            trait="response",
        )


def test_load_discrete_dataset_reports_sparse_states() -> None:
    dataset = load_discrete_dataset(
        fixture("example_tree.nwk"),
        fixture("example_traits_ancestral_sparse.tsv"),
        trait="habitat",
    )
    assert dataset.sparse_states == ["desert"]
    assert dataset.warnings == [
        "one or more discrete states are represented by fewer than two taxa and should be interpreted cautiously"
    ]


def test_discrete_reconstruction_reports_ambiguous_root_probabilities() -> None:
    report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="habitat",
    )
    estimates = {
        estimate.node: estimate for estimate in report.estimates if not estimate.is_tip
    }
    assert estimates["A|B"].most_likely_state == "forest"
    assert estimates["C|D"].most_likely_state == "tundra"
    assert estimates["A|B|C|D"].ambiguous is True
    assert estimates["A|B|C|D"].state_probabilities == {"forest": 0.5, "tundra": 0.5}
    assert estimates["A|B|C|D"].unstable is True
    assert estimates["A|B|C|D"].interpretation == "unstable node state"
    assert estimates["A|B|C|D"].downstream_risks
    assert report.weak_support_nodes == ["A|B|C|D"]


def test_discrete_reconstruction_supports_likelihood_models() -> None:
    report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="symmetric",
    )
    assert report.model == "symmetric"
    assert report.estimates[0].confidence > 0.0
    assert isinstance(report.unstable_nodes, list)


def test_discrete_reconstruction_supports_ordered_state_models() -> None:
    report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
        state_ordering="ordered",
        ordered_states=["north", "south", "island"],
    )
    assert report.state_ordering == "ordered"
    assert report.ordered_states == ["north", "south", "island"]


def test_discrete_reconstruction_rejects_single_observed_state() -> None:
    with pytest.raises(AncestralReconstructionError):
        reconstruct_discrete_ancestral_states(
            fixture("example_tree.nwk"),
            fixture("example_traits_ancestral_single_state.tsv"),
            trait="habitat",
        )
