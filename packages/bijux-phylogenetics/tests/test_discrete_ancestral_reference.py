from __future__ import annotations

import pytest

from bijux_phylogenetics.ancestral.discrete.review import (
    validate_discrete_ancestral_reference_examples,
)

pytestmark = pytest.mark.slow


def test_validate_discrete_ancestral_reference_examples_reports_passing_cases() -> None:
    report = validate_discrete_ancestral_reference_examples()

    assert report.case_count == 10
    assert report.external_case_count == 6
    assert report.all_passed is True
    assert {observation.category for observation in report.observations} == {
        "ambiguity-policy",
        "external-marginal-probability",
        "irreversible-constraint-policy",
        "ordered-constraint-policy",
        "root-prior-policy",
    }


def test_validate_discrete_ancestral_reference_examples_tracks_expanded_ace_cases() -> (
    None
):
    report = validate_discrete_ancestral_reference_examples()
    observation = next(
        row
        for row in report.observations
        if row.case_id == "geography_biased_all_rates_different"
    )

    assert observation.model == "all-rates-different"
    assert observation.passed is True
    assert observation.observed_metrics["root_state"] == "island"
    assert observation.observed_metrics["root_ambiguous"] is False
    assert (
        observation.observed_metrics["max_probability_delta"] <= observation.tolerance
    )
    assert len(observation.probability_rows) == 15


def test_validate_discrete_ancestral_reference_examples_covers_policy_surfaces() -> (
    None
):
    report = validate_discrete_ancestral_reference_examples()
    observations = {row.case_id: row for row in report.observations}

    root_prior = observations["root_prior_sensitivity_policy"]
    ordered = observations["ordered_transition_constraints"]
    irreversible = observations["irreversible_loss_constraints"]
    ambiguous = observations["ambiguous_internal_states"]

    assert root_prior.observed_metrics["state_changed_node_count"] == 2
    assert root_prior.observed_metrics["support_changed_node_count"] == 1
    assert root_prior.observed_metrics["fixed_root_state"] == "island"
    assert ordered.observed_metrics["preferred_ordering"] == "ordered"
    assert ordered.observed_metrics["north_to_island_allowed"] is False
    assert irreversible.observed_metrics["preferred_constraint"] == "constrained"
    assert irreversible.observed_metrics["absent_to_present_allowed"] is False
    assert ambiguous.observed_metrics["ambiguous_nodes"] == ["A|B|C|D"]
    assert ambiguous.observed_metrics["root_probabilities"] == {
        "forest": 0.5,
        "tundra": 0.5,
    }
