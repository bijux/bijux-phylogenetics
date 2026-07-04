from __future__ import annotations

import csv

import pytest

from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
    summarize_continuous_ancestral_report,
)
from bijux_phylogenetics.ancestral.discrete import (
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.comparative.covariance import (
    summarize_comparative_covariance_audit,
)
from bijux_phylogenetics.comparative.pgls import inspect_pgls_inputs
from bijux_phylogenetics.comparative.regression import (
    summarize_phylogenetic_anova,
)
from bijux_phylogenetics.comparative.signal import (
    compute_phylogenetic_independent_contrasts,
)
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_trait_table_fixture,
    get_shared_tree_fixture,
    list_shared_trait_table_fixtures,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError


def _read_taxa(fixture_id: str) -> list[str]:
    fixture = get_shared_trait_table_fixture(fixture_id)
    with fixture.path.open(encoding="utf-8", newline="") as handle:
        return [
            row[fixture.taxon_column] for row in csv.DictReader(handle, delimiter="\t")
        ]


def test_shared_trait_table_fixture_catalog_covers_required_goal_cases() -> None:
    fixtures = list_shared_trait_table_fixtures()
    feature_tags = {tag for fixture in fixtures for tag in fixture.feature_tags}

    assert {
        "continuous-trait",
        "continuous-ancestral-states",
        "binary-discrete-trait",
        "multistate-discrete-trait",
        "missing-trait-values",
        "extra-taxa-not-in-tree",
        "tree-taxa-missing-from-trait-table",
        "duplicated-taxon-row-negative-case",
        "constant-trait-negative-case",
        "categorical-predictor-case",
        "misordered-taxon-rows",
        "phylogenetic-independent-contrasts",
    } <= feature_tags


def test_shared_trait_table_fixture_lookup_preserves_durable_ids() -> None:
    fixture = get_shared_trait_table_fixture("binary_discrete_match")

    assert fixture.relative_path == "metadata/example_traits_binary.tsv"
    assert fixture.tree_fixture_id == "balanced_rooted_ultrametric"
    assert fixture.primary_trait_columns == ("presence",)
    assert "binary-discrete-trait" in fixture.feature_tags
    assert fixture.path.is_file()


def test_shared_trait_table_fixture_catalog_supports_governed_phytools_panels() -> None:
    small_fixture = get_shared_trait_table_fixture(
        "phytools_signal_panel_twenty_four_taxa"
    )
    large_fixture = get_shared_trait_table_fixture(
        "phytools_signal_panel_one_hundred_twenty_eight_taxa"
    )
    missing_fixture = get_shared_trait_table_fixture(
        "phytools_signal_missing_twenty_four_taxa"
    )
    mismatch_fixture = get_shared_trait_table_fixture(
        "phytools_signal_mismatch_twenty_four_taxa"
    )

    assert small_fixture.tree_fixture_id == "phytools_ultrametric_twenty_four_taxa"
    assert (
        large_fixture.tree_fixture_id
        == "phytools_ultrametric_one_hundred_twenty_eight_taxa"
    )
    assert small_fixture.row_count == 24
    assert large_fixture.row_count == 128
    assert {
        "signal_strong",
        "signal_weak",
        "binary_state",
        "region_state",
    } <= set(small_fixture.primary_trait_columns)
    assert "missing-trait-values" in missing_fixture.feature_tags
    assert "tree-trait-taxon-mismatch" in mismatch_fixture.feature_tags


def test_shared_trait_table_fixture_catalog_tracks_governed_standard_error_columns() -> (
    None
):
    fixture = get_shared_trait_table_fixture(
        "geiger_continuous_model_panel_twenty_four_taxa"
    )

    assert fixture.standard_error_columns == ("ou_truth_standard_error",)
    assert "ou_truth" in fixture.primary_trait_columns


def test_shared_trait_table_fixture_catalog_tracks_tree_trait_mismatch_surface() -> (
    None
):
    fixture = get_shared_trait_table_fixture("continuous_tree_mismatch")
    tree_fixture = get_shared_tree_fixture(fixture.tree_fixture_id)

    table_taxa = set(_read_taxa("continuous_tree_mismatch"))
    tree_taxa = set(load_tree(tree_fixture.path).tip_names)

    report = summarize_continuous_ancestral_report(
        reconstruct_continuous_ancestral_states(
            tree_fixture.path,
            fixture.path,
            trait="value",
        )
    )

    assert table_taxa - tree_taxa == {"E"}
    assert tree_taxa - table_taxa == {"D"}
    assert report.analyzed_taxon_count == 3
    assert report.missing_tip_taxon_count == 1


def test_shared_trait_table_fixture_catalog_supports_binary_and_multistate_ancestral_cases() -> (
    None
):
    tree = get_shared_tree_fixture("balanced_rooted_ultrametric").path
    binary_fixture = get_shared_trait_table_fixture("binary_discrete_match")
    multistate_fixture = get_shared_trait_table_fixture("multistate_discrete_match")
    constant_fixture = get_shared_trait_table_fixture("constant_trait_negative")

    binary_report = reconstruct_discrete_ancestral_states(
        tree,
        binary_fixture.path,
        trait="presence",
        model="equal-rates",
    )
    multistate_report = reconstruct_discrete_ancestral_states(
        tree,
        multistate_fixture.path,
        trait="region",
        model="equal-rates",
    )
    binary_root = next(
        estimate for estimate in binary_report.estimates if estimate.node == "A|B|C|D"
    )
    multistate_root = next(
        estimate
        for estimate in multistate_report.estimates
        if estimate.node == "A|B|C|D"
    )
    assert set(binary_root.state_probabilities) == {"0", "1"}
    assert set(multistate_root.state_probabilities) == {"north", "south", "island"}
    with pytest.raises(AncestralReconstructionError):
        reconstruct_discrete_ancestral_states(
            tree,
            constant_fixture.path,
            trait="habitat",
            model="equal-rates",
        )


def test_shared_trait_table_fixture_catalog_supports_comparative_mismatch_and_duplicate_cases() -> (
    None
):
    tree = get_shared_tree_fixture("balanced_rooted_ultrametric").path
    duplicate_fixture = get_shared_trait_table_fixture("duplicate_taxon_negative")
    missing_fixture = get_shared_trait_table_fixture("missing_trait_values")

    duplicate_report = summarize_comparative_covariance_audit(
        tree,
        duplicate_fixture.path,
        analysis="pgls",
        response="response",
        predictors=["predictor_one"],
        lambda_value=1.0,
    )
    missing_report = inspect_pgls_inputs(
        tree,
        missing_fixture.path,
        response="response",
        predictors=["predictor_one", "habitat"],
    )

    assert duplicate_report.duplicate_trait_taxa == ["A"]
    assert "trait table contains duplicate taxon keys" in duplicate_report.blockers
    assert missing_report.ready is False
    assert [row.taxon for row in missing_report.formula_audit.excluded_taxa] == ["D"]
    assert (
        "PGLS overfit guard requires at least one residual degree of freedom after predictor encoding"
        in missing_report.blockers
    )


def test_shared_trait_table_fixture_catalog_preserves_categorical_order_independence() -> (
    None
):
    tree = get_shared_tree_fixture("balanced_rooted_ultrametric").path
    ordered_fixture = get_shared_trait_table_fixture("categorical_predictor_match")
    reordered_fixture = get_shared_trait_table_fixture(
        "categorical_predictor_misordered"
    )

    ordered_report = inspect_pgls_inputs(
        tree,
        ordered_fixture.path,
        response="response",
        predictors=["predictor_one", "habitat"],
    )
    reordered_report = inspect_pgls_inputs(
        tree,
        reordered_fixture.path,
        response="response",
        predictors=["predictor_one", "habitat"],
    )

    assert ordered_report.ready is True
    assert reordered_report.ready is True
    assert ordered_report.encoded_columns == reordered_report.encoded_columns
    assert ordered_report.analysis_taxa == ["A", "B", "C", "D"]
    assert reordered_report.analysis_taxa == ["A", "B", "C", "D"]


def test_shared_trait_table_fixture_catalog_supports_phylogenetic_residual_cases() -> (
    None
):
    clean_fixture = get_shared_trait_table_fixture(
        "phylogenetic_residual_allometry_match"
    )
    missing_fixture = get_shared_trait_table_fixture(
        "phylogenetic_residual_allometry_missing"
    )

    assert clean_fixture.tree_fixture_id == "balanced_rooted_six_taxon"
    assert missing_fixture.tree_fixture_id == "balanced_rooted_six_taxon"
    assert clean_fixture.row_count == 6
    assert missing_fixture.row_count == 7
    assert {"body_mass", "brain_mass"} <= set(clean_fixture.primary_trait_columns)
    assert "allometry-case" in clean_fixture.feature_tags
    assert "residual-diagnostics" in clean_fixture.feature_tags
    assert "missing-trait-values" in missing_fixture.feature_tags
    assert "extra-taxa-not-in-tree" in missing_fixture.feature_tags
    assert _read_taxa("phylogenetic_residual_allometry_match") == [
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
    ]
    assert _read_taxa("phylogenetic_residual_allometry_missing")[-1] == "G"


def test_shared_trait_table_fixture_catalog_supports_phylogenetic_anova_cases() -> None:
    clean_fixture = get_shared_trait_table_fixture("phylogenetic_anova_group_effect")
    missing_fixture = get_shared_trait_table_fixture(
        "phylogenetic_anova_group_effect_missing"
    )
    tree = get_shared_tree_fixture(clean_fixture.tree_fixture_id).path

    clean_report = summarize_phylogenetic_anova(
        tree,
        clean_fixture.path,
        response="trait_value",
        group="habitat",
        simulations=16,
        seed=7,
    )
    missing_report = summarize_phylogenetic_anova(
        tree,
        missing_fixture.path,
        response="trait_value",
        group="habitat",
        simulations=16,
        seed=7,
    )

    assert clean_fixture.tree_fixture_id == "balanced_rooted_six_taxon"
    assert missing_fixture.tree_fixture_id == "balanced_rooted_six_taxon"
    assert clean_fixture.row_count == 6
    assert missing_fixture.row_count == 7
    assert {"trait_value", "habitat"} <= set(clean_fixture.primary_trait_columns)
    assert "phylogenetic-anova" in clean_fixture.feature_tags
    assert "unequal-group-sizes" in clean_fixture.feature_tags
    assert "missing-trait-values" in missing_fixture.feature_tags
    assert clean_report.group_count == 2
    assert {row.group: row.taxon_count for row in clean_report.group_rows} == {
        "desert": 2,
        "forest": 4,
    }
    assert {row.taxon: row.reason for row in missing_report.excluded_taxa} == {
        "F": "missing_value",
        "G": "absent_from_tree",
    }


def test_shared_trait_table_fixture_catalog_supports_governed_pic_cases() -> None:
    balanced_fixture = get_shared_trait_table_fixture("pic_continuous_balanced")
    pectinate_fixture = get_shared_trait_table_fixture("pic_continuous_pectinate")
    six_taxon_fixture = get_shared_trait_table_fixture("pic_continuous_six_taxon")
    missing_fixture = get_shared_trait_table_fixture("pic_continuous_missing_values")

    balanced_tree = get_shared_tree_fixture(balanced_fixture.tree_fixture_id).path
    pectinate_tree = get_shared_tree_fixture(pectinate_fixture.tree_fixture_id).path
    six_taxon_tree = get_shared_tree_fixture(six_taxon_fixture.tree_fixture_id).path

    balanced_report = compute_phylogenetic_independent_contrasts(
        balanced_tree,
        balanced_fixture.path,
        trait="response",
    )
    pectinate_report = compute_phylogenetic_independent_contrasts(
        pectinate_tree,
        pectinate_fixture.path,
        trait="response",
    )
    six_taxon_report = compute_phylogenetic_independent_contrasts(
        six_taxon_tree,
        six_taxon_fixture.path,
        trait="response_growth",
    )
    missing_report = compute_phylogenetic_independent_contrasts(
        six_taxon_tree,
        missing_fixture.path,
        trait="response_growth",
    )

    assert balanced_report.taxon_count == 4
    assert [row.node_id for row in balanced_report.contrasts] == [6, 7, 5]
    assert pectinate_report.input_audit.tree_is_ultrametric is False
    assert [row.node_id for row in pectinate_report.contrasts] == [7, 6, 5]
    assert six_taxon_report.taxon_count == 6
    assert len(six_taxon_report.contrasts) == 5
    assert missing_report.taxon_count == 4
    assert missing_report.input_audit.pruned_missing_value_taxa == ["B"]
    assert "nonnumeric-trait-values" in missing_fixture.feature_tags


def test_shared_trait_table_fixture_catalog_supports_governed_continuous_ace_cases() -> (
    None
):
    balanced_fixture = get_shared_trait_table_fixture("ace_continuous_balanced")
    pectinate_fixture = get_shared_trait_table_fixture("ace_continuous_pectinate")
    six_taxon_fixture = get_shared_trait_table_fixture("ace_continuous_six_taxon")
    missing_fixture = get_shared_trait_table_fixture("ace_continuous_missing_values")

    balanced_tree = get_shared_tree_fixture(balanced_fixture.tree_fixture_id).path
    pectinate_tree = get_shared_tree_fixture(pectinate_fixture.tree_fixture_id).path
    six_taxon_tree = get_shared_tree_fixture(six_taxon_fixture.tree_fixture_id).path

    balanced_summary = summarize_continuous_ancestral_report(
        reconstruct_continuous_ancestral_states(
            balanced_tree,
            balanced_fixture.path,
            trait="response",
        )
    )
    pectinate_summary = summarize_continuous_ancestral_report(
        reconstruct_continuous_ancestral_states(
            pectinate_tree,
            pectinate_fixture.path,
            trait="response",
        )
    )
    six_taxon_summary = summarize_continuous_ancestral_report(
        reconstruct_continuous_ancestral_states(
            six_taxon_tree,
            six_taxon_fixture.path,
            trait="response_growth",
        )
    )
    missing_report = reconstruct_continuous_ancestral_states(
        six_taxon_tree,
        missing_fixture.path,
        trait="response_growth",
    )

    assert balanced_summary.analyzed_taxon_count == 4
    assert balanced_summary.tree_is_ultrametric is True
    assert pectinate_summary.analyzed_taxon_count == 4
    assert pectinate_summary.tree_is_ultrametric is False
    assert six_taxon_summary.analyzed_taxon_count == 6
    assert six_taxon_summary.internal_node_count == 5
    assert missing_report.taxon_count == 4
    assert missing_report.dropped_missing_taxa == ["B"]
    assert missing_report.dropped_non_numeric_taxa == ["C"]
    assert missing_report.brownian_fit_diagnostics is not None
    assert missing_report.brownian_fit_diagnostics.tree_is_ultrametric is True


def test_shared_trait_table_fixture_catalog_supports_governed_discrete_ace_cases() -> (
    None
):
    binary_fixture = get_shared_trait_table_fixture("binary_discrete_match")
    multistate_fixture = get_shared_trait_table_fixture("multistate_discrete_match")
    missing_fixture = get_shared_trait_table_fixture("missing_trait_values")

    balanced_tree = get_shared_tree_fixture(binary_fixture.tree_fixture_id).path
    pectinate_tree = get_shared_tree_fixture("pectinate_rooted_non_ultrametric").path

    binary_report = reconstruct_discrete_ancestral_states(
        balanced_tree,
        binary_fixture.path,
        trait="presence",
        model="equal-rates",
    )
    multistate_report = reconstruct_discrete_ancestral_states(
        balanced_tree,
        multistate_fixture.path,
        trait="region",
        model="equal-rates",
    )
    pectinate_report = reconstruct_discrete_ancestral_states(
        pectinate_tree,
        multistate_fixture.path,
        trait="region",
        model="equal-rates",
    )
    missing_report = reconstruct_discrete_ancestral_states(
        balanced_tree,
        missing_fixture.path,
        trait="habitat",
        model="equal-rates",
    )

    assert binary_report.log_likelihood is not None
    assert binary_report.parameter_count == 1
    assert len(binary_report.transition_rate_rows) == 2
    assert multistate_report.log_likelihood is not None
    assert multistate_report.parameter_count == 1
    assert len(multistate_report.transition_rate_rows) == 6
    assert pectinate_report.log_likelihood is not None
    assert pectinate_report.parameter_count == 1
    assert missing_report.taxon_count == 3
    assert missing_report.dropped_missing_taxa == ["D"]


def test_shared_trait_table_fixture_catalog_supports_governed_discrete_sym_ace_cases() -> (
    None
):
    balanced_fixture = get_shared_trait_table_fixture("ace_discrete_sym_balanced")
    pectinate_fixture = get_shared_trait_table_fixture("ace_discrete_sym_pectinate")
    six_taxon_fixture = get_shared_trait_table_fixture("ace_discrete_sym_six_taxon")
    missing_fixture = get_shared_trait_table_fixture("ace_discrete_sym_missing_values")

    balanced_tree = get_shared_tree_fixture(balanced_fixture.tree_fixture_id).path
    pectinate_tree = get_shared_tree_fixture(pectinate_fixture.tree_fixture_id).path
    six_taxon_tree = get_shared_tree_fixture(six_taxon_fixture.tree_fixture_id).path

    balanced_report = reconstruct_discrete_ancestral_states(
        balanced_tree,
        balanced_fixture.path,
        trait="region",
        model="symmetric",
    )
    pectinate_report = reconstruct_discrete_ancestral_states(
        pectinate_tree,
        pectinate_fixture.path,
        trait="region",
        model="symmetric",
    )
    six_taxon_report = reconstruct_discrete_ancestral_states(
        six_taxon_tree,
        six_taxon_fixture.path,
        trait="region",
        model="symmetric",
    )
    missing_report = reconstruct_discrete_ancestral_states(
        six_taxon_tree,
        missing_fixture.path,
        trait="region",
        model="symmetric",
    )

    assert balanced_report.model == "symmetric"
    assert balanced_report.parameter_count == 3
    assert balanced_report.optimizer_diagnostics is not None
    assert pectinate_report.model == "symmetric"
    assert len(pectinate_report.transition_rate_rows) == 6
    assert six_taxon_report.taxon_count == 6
    assert six_taxon_report.baseline_comparison is not None
    assert six_taxon_report.baseline_comparison.delta_aic > 0.0
    assert six_taxon_report.baseline_comparison.preferred_model_by_aic == "equal-rates"
    assert missing_report.taxon_count == 5
    assert missing_report.dropped_missing_taxa == ["F"]


@pytest.mark.slow
def test_shared_trait_table_fixture_catalog_supports_governed_discrete_ard_ace_cases() -> (
    None
):
    binary_fixture = get_shared_trait_table_fixture("ace_discrete_ard_binary_balanced")
    pectinate_fixture = get_shared_trait_table_fixture("ace_discrete_ard_pectinate")
    six_taxon_fixture = get_shared_trait_table_fixture("ace_discrete_ard_six_taxon")
    missing_fixture = get_shared_trait_table_fixture("ace_discrete_ard_missing_values")

    binary_tree = get_shared_tree_fixture(binary_fixture.tree_fixture_id).path
    pectinate_tree = get_shared_tree_fixture(pectinate_fixture.tree_fixture_id).path
    six_taxon_tree = get_shared_tree_fixture(six_taxon_fixture.tree_fixture_id).path

    binary_report = reconstruct_discrete_ancestral_states(
        binary_tree,
        binary_fixture.path,
        trait="habitat",
        model="all-rates-different",
    )
    pectinate_report = reconstruct_discrete_ancestral_states(
        pectinate_tree,
        pectinate_fixture.path,
        trait="region",
        model="all-rates-different",
    )
    six_taxon_report = reconstruct_discrete_ancestral_states(
        six_taxon_tree,
        six_taxon_fixture.path,
        trait="region",
        model="all-rates-different",
    )
    missing_report = reconstruct_discrete_ancestral_states(
        six_taxon_tree,
        missing_fixture.path,
        trait="region",
        model="all-rates-different",
    )

    assert binary_report.model == "all-rates-different"
    assert binary_report.parameter_count == 2
    assert len(binary_report.transition_rate_rows) == 2
    assert pectinate_report.model == "all-rates-different"
    assert pectinate_report.parameter_count == 6
    assert pectinate_report.overparameterized is True
    assert pectinate_report.optimizer_diagnostics is not None
    assert six_taxon_report.taxon_count == 6
    assert six_taxon_report.baseline_comparison is not None
    assert six_taxon_report.baseline_comparison.preferred_model_by_aic == "equal-rates"
    assert missing_report.taxon_count == 5
    assert missing_report.dropped_missing_taxa == ["F"]
