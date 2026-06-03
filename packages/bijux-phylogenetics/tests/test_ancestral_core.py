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
    reconstruct_continuous_ancestral_states_from_dataset,
)
from bijux_phylogenetics.ancestral.discrete import (
    reconstruct_discrete_ancestral_states,
    reconstruct_discrete_ancestral_states_from_dataset,
)
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

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
    assert dataset.dropped_missing_taxa == []
    assert dataset.dropped_non_numeric_taxa == []
    assert dataset.warnings == [
        "continuous trait reconstruction is using only 3 taxa; results may be unstable",
        "one or more tree taxa were excluded because they were absent from the trait table",
        "one or more trait rows were excluded because their taxa were absent from the tree",
    ]
    assert dataset.alignment_report.dropped_tree_taxa == ["D"]
    assert dataset.alignment_report.dropped_trait_taxa == ["E"]
    assert dataset.alignment_report.dropped_missing_value_taxa == []


def test_load_continuous_dataset_prunes_nonnumeric_taxa_explicitly() -> None:
    dataset = load_continuous_dataset(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_brownian_missing.tsv"),
        trait="response_growth",
    )
    assert dataset.taxa == ["A", "D", "E", "F"]
    assert dataset.dropped_missing_taxa == ["B"]
    assert dataset.dropped_non_numeric_taxa == ["C"]
    assert dataset.warnings == [
        "one or more trait rows were excluded because their taxa were absent from the tree",
        "one or more taxa were excluded because the continuous trait value was missing",
        "one or more taxa were excluded because the continuous trait value was not numeric",
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


def test_continuous_reconstruction_from_dataset_matches_path_surface() -> None:
    dataset = load_continuous_dataset(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )

    direct_report = reconstruct_continuous_ancestral_states_from_dataset(
        dataset,
        model="brownian",
    )
    path_report = reconstruct_continuous_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        model="brownian",
    )

    assert direct_report.tree_path == path_report.tree_path
    assert direct_report.traits_path == path_report.traits_path
    assert direct_report.analysis_tree_newick == path_report.analysis_tree_newick
    assert direct_report.warnings == path_report.warnings
    assert direct_report.estimates == path_report.estimates
    assert (
        direct_report.brownian_fit_diagnostics == path_report.brownian_fit_diagnostics
    )


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


def test_continuous_reconstruction_supports_fast_anc_estimator() -> None:
    report = reconstruct_continuous_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        model="brownian",
        estimator="fast-anc",
    )
    internal_estimates = {
        estimate.node: estimate for estimate in report.estimates if not estimate.is_tip
    }

    assert report.estimator == "fast-anc"
    assert math.isclose(internal_estimates["A|B"].estimate, 2.361111111111111)
    assert math.isclose(internal_estimates["C|D"].estimate, 3.027777777777778)
    assert math.isclose(internal_estimates["A|B|C|D"].estimate, 2.8055555555555554)
    assert math.isclose(
        internal_estimates["A|B"].standard_error,
        0.5319039487535212,
    )
    assert math.isclose(
        internal_estimates["A|B|C|D"].standard_error,
        0.8410139872493032,
    )
    assert report.brownian_fit_diagnostics is not None


def test_continuous_fast_anc_reconstruction_from_dataset_matches_path_surface() -> None:
    dataset = load_continuous_dataset(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )

    direct_report = reconstruct_continuous_ancestral_states_from_dataset(
        dataset,
        model="brownian",
        estimator="fast-anc",
    )
    path_report = reconstruct_continuous_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        model="brownian",
        estimator="fast-anc",
    )

    assert direct_report.estimator == "fast-anc"
    assert direct_report.analysis_tree_newick == path_report.analysis_tree_newick
    assert direct_report.warnings == path_report.warnings
    assert direct_report.estimates == path_report.estimates
    assert (
        direct_report.brownian_fit_diagnostics == path_report.brownian_fit_diagnostics
    )


def test_continuous_reconstruction_supports_anc_ml_estimator() -> None:
    report = reconstruct_continuous_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        model="brownian",
        estimator="anc-ml",
    )
    internal_estimates = {
        estimate.node: estimate for estimate in report.estimates if not estimate.is_tip
    }

    assert report.estimator == "anc-ml"
    assert report.optimizer_diagnostics is not None
    assert report.optimizer_diagnostics.optimizer_name == "closed-form-profile-solution"
    assert report.optimizer_diagnostics.converged is True
    assert report.optimizer_diagnostics.iteration_count == 0
    assert math.isclose(internal_estimates["A|B|C|D"].estimate, 2.8055555555555554)
    assert math.isclose(internal_estimates["A|B"].estimate, 2.3611111111111112)
    assert math.isclose(internal_estimates["C|D"].estimate, 3.0277777777777781)
    assert math.isclose(
        internal_estimates["A|B|C|D"].standard_error,
        0.5946866992509994,
        rel_tol=0.0,
        abs_tol=1e-8,
    )
    assert math.isclose(
        internal_estimates["A|B"].standard_error,
        0.37611289276824167,
        rel_tol=0.0,
        abs_tol=1e-8,
    )
    assert math.isclose(
        internal_estimates["C|D"].standard_error,
        0.49755058957489323,
        rel_tol=0.0,
        abs_tol=1e-8,
    )
    assert math.isclose(
        internal_estimates["A|B|C|D"].lower_95_interval,
        1.6399696250235967,
        rel_tol=0.0,
        abs_tol=5e-8,
    )
    assert math.isclose(
        internal_estimates["A|B|C|D"].upper_95_interval,
        3.971141486087514,
        rel_tol=0.0,
        abs_tol=5e-8,
    )
    assert report.brownian_fit_diagnostics is not None
    assert math.isclose(
        report.brownian_fit_diagnostics.log_likelihood,
        -6.1189469566554999,
    )
    assert math.isclose(
        report.brownian_fit_diagnostics.residual_sigma_squared,
        3.1828704323963097,
        rel_tol=0.0,
        abs_tol=1e-6,
    )


def test_continuous_anc_ml_reconstruction_from_dataset_matches_path_surface() -> None:
    dataset = load_continuous_dataset(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )

    direct_report = reconstruct_continuous_ancestral_states_from_dataset(
        dataset,
        model="brownian",
        estimator="anc-ml",
    )
    path_report = reconstruct_continuous_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        model="brownian",
        estimator="anc-ml",
    )

    assert direct_report.estimator == "anc-ml"
    assert direct_report.analysis_tree_newick == path_report.analysis_tree_newick
    assert direct_report.warnings == path_report.warnings
    assert direct_report.estimates == path_report.estimates
    assert (
        direct_report.brownian_fit_diagnostics == path_report.brownian_fit_diagnostics
    )
    assert direct_report.optimizer_diagnostics == path_report.optimizer_diagnostics


def test_continuous_reconstruction_rejects_incompatible_fast_anc_estimator() -> None:
    with pytest.raises(ValueError):
        reconstruct_continuous_ancestral_states(
            fixture("example_tree.nwk"),
            fixture("example_traits_comparative.tsv"),
            trait="response",
            model="ou",
            estimator="fast-anc",
        )


def test_continuous_reconstruction_rejects_incompatible_anc_ml_estimator() -> None:
    with pytest.raises(ValueError):
        reconstruct_continuous_ancestral_states(
            fixture("example_tree.nwk"),
            fixture("example_traits_comparative.tsv"),
            trait="response",
            model="ou",
            estimator="anc-ml",
        )


def test_continuous_reconstruction_reports_singular_anc_ml_covariance() -> None:
    with pytest.raises(
        ValueError,
        match="requires a nonsingular full Brownian covariance matrix",
    ):
        reconstruct_continuous_ancestral_states(
            fixture("example_tree_ultrametric_zero_internal.nwk"),
            fixture("example_traits_comparative.tsv"),
            trait="response",
            model="brownian",
            estimator="anc-ml",
        )


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


def test_discrete_reconstruction_from_dataset_matches_path_surface() -> None:
    dataset = load_discrete_dataset(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="habitat",
    )

    direct_report = reconstruct_discrete_ancestral_states_from_dataset(
        dataset,
        model="equal-rates",
    )
    path_report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="habitat",
        model="equal-rates",
    )

    assert direct_report.tree_path == path_report.tree_path
    assert direct_report.traits_path == path_report.traits_path
    assert direct_report.analysis_tree_newick == path_report.analysis_tree_newick
    assert direct_report.warnings == path_report.warnings
    assert direct_report.estimates == path_report.estimates
    assert direct_report.transition_rate_rows == path_report.transition_rate_rows
    assert direct_report.optimizer_diagnostics == path_report.optimizer_diagnostics


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
    assert report.rerooting_method_compatibility.comparable is True
    assert report.rerooting_method_compatibility.reference_model == "SYM"
    assert report.rerooting_method_compatibility.reference_root_prior_mode == "equal"
    assert report.estimates[0].confidence > 0.0
    assert isinstance(report.unstable_nodes, list)


@pytest.mark.slow
def test_discrete_reconstruction_matches_governed_ard_probability_fixture() -> None:
    report = reconstruct_discrete_ancestral_states(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_geography_biased.tsv"),
        trait="region",
        taxon_column="taxon",
        model="all-rates-different",
    )
    internal_estimates = {
        estimate.node: estimate for estimate in report.estimates if not estimate.is_tip
    }

    expected_probabilities = {
        "A|B|C|D|E|F": {
            "island": 0.532202345485066,
            "north": 0.0563498359339132,
            "south": 0.41144781858102,
        },
        "A|B|C|D": {
            "island": 0.134323458950047,
            "north": 0.806889166520943,
            "south": 0.05878737452901,
        },
        "A|B": {
            "island": 0.0640940564112071,
            "north": 0.914130423465446,
            "south": 0.021775520123347,
        },
        "C|D": {
            "island": 0.0640940564112071,
            "north": 0.914130423465446,
            "south": 0.021775520123347,
        },
        "E|F": {
            "island": 0.389284284222312,
            "north": 0.00013599404628342,
            "south": 0.610579721731405,
        },
    }
    for node, expected_node_probabilities in expected_probabilities.items():
        observed_probabilities = internal_estimates[node].state_probabilities
        for state, expected_probability in expected_node_probabilities.items():
            assert math.isclose(
                observed_probabilities[state],
                expected_probability,
                rel_tol=0.0,
                abs_tol=1e-15,
            )


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
    assert report.rerooting_method_compatibility.comparable is False
    assert (
        "phytools::rerootingMethod does not provide a governed ordered-transition parity surface in this repository"
        in report.rerooting_method_compatibility.notes
    )


def test_discrete_reconstruction_rejects_meristic_parity_claim() -> None:
    with pytest.raises(ValueError) as excinfo:
        reconstruct_discrete_ancestral_states(
            fixture("example_tree.nwk"),
            fixture("example_traits_geography.tsv"),
            trait="region",
            model="meristic",
            state_ordering="ordered",
            ordered_states=["north", "south", "island"],
        )

    assert "explicitly excluded this round" in str(excinfo.value)
    assert "integer-state meristic contract" in str(excinfo.value)


@pytest.mark.slow
def test_discrete_reconstruction_reports_rerooting_method_assumption_limits() -> None:
    fixed_root_report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
        root_prior_mode="fixed",
        fixed_root_state="north",
    )
    ard_report = reconstruct_discrete_ancestral_states(
        fixture("example_tree_ladderized.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="all-rates-different",
    )

    assert fixed_root_report.rerooting_method_compatibility.comparable is False
    assert (
        "phytools::rerootingMethod inherits fitMk's default equal root prior; empirical or fixed root-prior runs remain Bijux sensitivity scenarios without direct rerootingMethod parity"
        in fixed_root_report.warnings
    )
    assert ard_report.rerooting_method_compatibility.comparable is False
    assert (
        "phytools::rerootingMethod is invalid for non-symmetric Q matrices such as all-rates-different models in phytools 2.5.2"
        in ard_report.warnings
    )


def test_discrete_reconstruction_rejects_single_observed_state() -> None:
    with pytest.raises(AncestralReconstructionError):
        reconstruct_discrete_ancestral_states(
            fixture("example_tree.nwk"),
            fixture("example_traits_ancestral_single_state.tsv"),
            trait="habitat",
        )
