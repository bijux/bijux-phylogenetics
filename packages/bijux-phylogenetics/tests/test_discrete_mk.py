from __future__ import annotations

import math
from pathlib import Path
import random

import pytest

from bijux_phylogenetics.ancestral.common import load_discrete_dataset
from bijux_phylogenetics.comparative.discrete_mk import (
    compare_discrete_mk_model_ranking,
    compare_discrete_mk_model_ranking_from_dataset,
    fit_discrete_mk_model,
    fit_discrete_mk_model_from_dataset,
    write_discrete_mk_summary_table,
)
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_geiger_discrete_fixture,
    get_shared_phytools_comparative_fixture,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.simulation import (
    simulate_discrete_traits,
    write_discrete_trait_table,
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


def _single_allowed_rate(report) -> float:
    allowed_rates = {
        row.rate for row in report.transition_rate_rows if row.transition_allowed
    }
    assert len(allowed_rates) == 1
    return next(iter(allowed_rates))


def _allowed_rate_lookup(report) -> dict[tuple[str, str], float]:
    return {
        (row.source_state, row.target_state): row.rate
        for row in report.transition_rate_rows
        if row.transition_allowed
    }


def _simulate_directional_tip_states(
    tree_path: Path,
    *,
    rates_by_state: dict[str, dict[str, float]],
    root_state: str,
    seed: int,
) -> dict[str, str]:
    tree = load_tree(tree_path)
    rng = random.Random(seed)  # nosec B311
    node_states: dict[int, str] = {}

    def visit(node, state: str) -> None:
        node_states[id(node)] = state
        for child in node.children:
            current_state = state
            branch_length = max(child.branch_length or 0.0, 0.0)
            elapsed = 0.0
            while True:
                outgoing_rates = rates_by_state.get(current_state, {})
                total_rate = sum(outgoing_rates.values())
                if total_rate <= 0.0:
                    break
                wait_time = rng.expovariate(total_rate)
                if elapsed + wait_time >= branch_length:
                    break
                elapsed += wait_time
                threshold = rng.random() * total_rate
                running_total = 0.0
                for target_state, rate in outgoing_rates.items():
                    running_total += rate
                    if threshold <= running_total:
                        current_state = target_state
                        break
            visit(child, current_state)

    visit(tree.root, root_state)
    return {
        node.name: node_states[id(node)] for node in tree.iter_nodes() if node.is_leaf()
    }


def _write_simulated_state_table(
    path: Path,
    *,
    tip_states: dict[str, str],
    trait: str = "state",
) -> Path:
    path.write_text(
        "taxon\t"
        + trait
        + "\n"
        + "".join(f"{taxon}\t{state}\n" for taxon, state in sorted(tip_states.items())),
        encoding="utf-8",
    )
    return path


def test_fit_discrete_mk_model_reports_binary_equal_rates_surface() -> None:
    report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="binary_state",
        taxon_column="taxon",
        model="equal-rates",
    )

    assert report.model == "equal-rates"
    assert report.taxon_count == 24
    assert report.parameter_count == 1
    assert report.input_audit.observed_states == ["0", "1"]
    assert report.input_audit.pruned_missing_value_taxa == []
    assert report.input_audit.missing_value_policy == "prune-overlapping-missing-values"
    assert math.isfinite(report.log_likelihood)
    assert math.isfinite(report.aic)
    assert math.isfinite(report.aicc)
    assert report.optimizer_diagnostics.converged is True
    assert len(report.transition_rate_rows) == 2
    assert _single_allowed_rate(report) > 0.0


def test_fit_discrete_mk_model_reports_multistate_equal_rates_surface() -> None:
    report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="region_state",
        taxon_column="taxon",
        model="equal-rates",
    )

    assert report.model == "equal-rates"
    assert report.taxon_count == 24
    assert report.parameter_count == 1
    assert report.input_audit.observed_states == ["north", "south", "west"]
    assert len(report.transition_rate_rows) == 6
    assert all(row.transition_allowed for row in report.transition_rate_rows)
    assert {row.step_distance for row in report.transition_rate_rows} == {1}
    assert _single_allowed_rate(report) > 0.0


@pytest.mark.slow
def test_fit_discrete_mk_model_reports_multistate_symmetric_surface() -> None:
    report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="region_state",
        taxon_column="taxon",
        model="symmetric",
    )

    rate_lookup = _allowed_rate_lookup(report)

    assert report.model == "symmetric"
    assert report.taxon_count == 24
    assert report.parameter_count == 3
    assert report.input_audit.observed_states == ["north", "south", "west"]
    assert report.baseline_comparison is not None
    assert report.baseline_comparison.baseline_model == "equal-rates"
    assert report.baseline_comparison.preferred_model_by_aic == "equal-rates"
    assert len(report.transition_rate_rows) == 6
    assert math.isclose(
        rate_lookup[("north", "south")],
        rate_lookup[("south", "north")],
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        rate_lookup[("north", "west")],
        rate_lookup[("west", "north")],
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        rate_lookup[("south", "west")],
        rate_lookup[("west", "south")],
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_fit_discrete_mk_model_reports_binary_ard_surface() -> None:
    fixture_entry = get_shared_phytools_comparative_fixture(
        "phytools_discrete_ard_binary_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="all-rates-different",
    )

    rate_lookup = _allowed_rate_lookup(report)

    assert report.model == "all-rates-different"
    assert report.parameter_count == 2
    assert report.baseline_comparison is not None
    assert report.optimizer_diagnostics.converged is True
    assert report.optimizer_diagnostics.hit_lower_parameter_bound is False
    assert report.optimizer_diagnostics.hit_upper_parameter_bound is False
    assert rate_lookup[("0", "1")] > 0.0
    assert rate_lookup[("1", "0")] > 0.0
    assert not math.isclose(
        rate_lookup[("0", "1")],
        rate_lookup[("1", "0")],
        rel_tol=0.0,
        abs_tol=1e-9,
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_reports_multistate_ard_surface() -> None:
    fixture_entry = get_shared_phytools_comparative_fixture(
        "phytools_discrete_ard_multistate_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="all-rates-different",
    )

    rate_lookup = _allowed_rate_lookup(report)

    assert report.model == "all-rates-different"
    assert report.parameter_count == 12
    assert report.baseline_comparison is not None
    assert report.optimizer_diagnostics.converged is False
    assert report.optimizer_diagnostics.hit_lower_parameter_bound is True
    assert report.optimizer_diagnostics.hit_upper_parameter_bound is False
    assert any(
        "weakly identified" in warning for warning in report.input_audit.warnings
    )
    assert rate_lookup[("north", "south")] < 1e-3
    assert rate_lookup[("north", "west")] < 1e-3
    assert rate_lookup[("west", "south")] > 1.0
    assert rate_lookup[("east", "west")] > 1.0
    assert not math.isclose(
        rate_lookup[("east", "west")],
        rate_lookup[("west", "east")],
        rel_tol=0.0,
        abs_tol=1e-9,
    )


def test_fit_discrete_mk_model_from_dataset_matches_path_surface() -> None:
    tree = fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")
    traits = fixture("example_traits_phytools_signal_twenty_four_taxa.tsv")
    dataset = load_discrete_dataset(
        tree,
        traits,
        trait="binary_state",
        taxon_column="taxon",
    )

    path_report = fit_discrete_mk_model(
        tree,
        traits,
        trait="binary_state",
        taxon_column="taxon",
        model="equal-rates",
    )
    dataset_report = fit_discrete_mk_model_from_dataset(
        dataset,
        model="equal-rates",
    )

    assert dataset_report.log_likelihood == path_report.log_likelihood
    assert dataset_report.aic == path_report.aic
    assert dataset_report.aicc == path_report.aicc
    assert dataset_report.transition_rate_rows == path_report.transition_rate_rows


def test_fit_discrete_mk_model_reports_pruned_missing_values() -> None:
    report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_discrete_missing_twenty_four_taxa.tsv"),
        trait="region_state",
        taxon_column="taxon",
        model="equal-rates",
    )

    assert report.taxon_count == 23
    assert report.input_audit.pruned_missing_value_taxa == ["Phy14"]
    assert (
        "one or more taxa were excluded because the discrete trait state was missing"
        in report.input_audit.warnings
    )
    assert len(report.transition_rate_rows) == 6


@pytest.mark.slow
def test_fit_discrete_mk_model_reports_symmetric_pruned_missing_values() -> None:
    report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_discrete_missing_twenty_four_taxa.tsv"),
        trait="region_state",
        taxon_column="taxon",
        model="symmetric",
    )

    assert report.model == "symmetric"
    assert report.taxon_count == 23
    assert report.parameter_count == 3
    assert report.input_audit.pruned_missing_value_taxa == ["Phy14"]
    assert report.baseline_comparison is not None


def test_fit_discrete_mk_model_rejects_meristic_parity_claim() -> None:
    with pytest.raises(ValueError) as excinfo:
        fit_discrete_mk_model(
            fixture("example_tree.nwk"),
            fixture("example_traits_geography.tsv"),
            trait="region",
            taxon_column="taxon",
            model="meristic",
            state_ordering="ordered",
            ordered_states=["north", "south", "island"],
        )

    assert "explicitly excluded this round" in str(excinfo.value)
    assert "ordered-state Mk support is not claimed as meristic parity" in str(
        excinfo.value
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_reports_ard_pruned_missing_values() -> None:
    binary_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_ard_binary_missing_twenty_four_taxa"
    )
    multistate_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_ard_multistate_missing_twenty_four_taxa"
    )

    binary_report = fit_discrete_mk_model(
        binary_fixture.tree_path,
        binary_fixture.traits_path,
        trait=binary_fixture.trait_name,
        taxon_column=binary_fixture.taxon_column,
        model="all-rates-different",
    )
    multistate_report = fit_discrete_mk_model(
        multistate_fixture.tree_path,
        multistate_fixture.traits_path,
        trait=multistate_fixture.trait_name,
        taxon_column=multistate_fixture.taxon_column,
        model="all-rates-different",
    )

    assert binary_report.taxon_count == 23
    assert binary_report.input_audit.pruned_missing_value_taxa == ["Phy10"]
    assert binary_report.parameter_count == 2
    assert multistate_report.taxon_count == 23
    assert multistate_report.input_audit.pruned_missing_value_taxa == ["Phy14"]
    assert multistate_report.parameter_count == 12
    assert multistate_report.optimizer_diagnostics.hit_lower_parameter_bound is True


def test_fit_discrete_mk_model_recovers_binary_er_known_truth(tmp_path: Path) -> None:
    simulation = simulate_discrete_traits(
        fixture("example_tree_phytools_ultrametric_one_hundred_twenty_eight_taxa.nwk"),
        states=["0", "1"],
        transition_rate=0.35,
        root_state="0",
        seed=7128,
    )
    traits_path = write_discrete_trait_table(tmp_path / "binary-er.tsv", simulation)

    report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_one_hundred_twenty_eight_taxa.nwk"),
        traits_path,
        trait="state",
        model="equal-rates",
    )

    assert math.isclose(_single_allowed_rate(report), 0.35, rel_tol=0.0, abs_tol=0.18)


def test_fit_discrete_mk_model_matches_governed_geiger_er_binary_surface() -> None:
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_er_binary_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="equal-rates",
    )

    assert report.optimizer_diagnostics.optimizer_name == "golden-section-search"
    assert math.isclose(
        report.log_likelihood,
        -9.078105640476831,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        _single_allowed_rate(report),
        0.39352316673030907,
        rel_tol=0.0,
        abs_tol=1e-6,
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_matches_governed_geiger_er_lambda_surface() -> None:
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_er_binary_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="equal-rates",
        transform="lambda",
    )

    assert report.transform_fit is not None
    assert report.parameter_count == 2
    assert math.isclose(
        report.log_likelihood,
        -9.078105725221098,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        report.transform_fit.parameter_value,
        1.0,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert report.transform_fit.starting_parameter_policy == (
        "lower-bound-first-evaluation"
    )
    assert math.isclose(
        report.transform_fit.starting_parameter_value,
        0.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert report.transform_fit.refinement_start_count == 1


@pytest.mark.slow
def test_fit_discrete_mk_model_marks_governed_geiger_lambda_weak_signal_surface() -> (
    None
):
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_transform_weak_signal_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="equal-rates",
        transform="lambda",
    )

    assert report.transform_fit is not None
    assert report.parameter_count == 2
    assert report.transform_fit.parameter_value <= 1e-6
    assert math.isclose(
        report.log_likelihood,
        -16.635532333438686,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert any(
        warning.kind == "weak_phylogenetic_signal"
        for warning in report.transform_fit.warnings
    )


@pytest.mark.slow
def test_write_discrete_mk_summary_table_reports_shared_transform_search_audit(
    tmp_path: Path,
) -> None:
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_er_binary_twenty_four_taxa"
    )
    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="equal-rates",
        transform="lambda",
    )

    output_path = tmp_path / "discrete-summary.tsv"
    write_discrete_mk_summary_table(output_path, report)
    lines = output_path.read_text(encoding="utf-8").strip().splitlines()

    assert "transform_starting_parameter_policy" in lines[0]
    assert "transform_refinement_start_count" in lines[0]
    assert "likelihood_constant_policy" in lines[0]
    assert "lower-bound-first-evaluation" in lines[1]
    assert "\t1\t" in lines[1]
    assert (
        "continuous-time-markov-pruning-loglikelihood-has-no-extra-normalizing-constant"
        in lines[1]
    )


@pytest.mark.slow
def test_compare_discrete_mk_model_ranking_reports_ranked_model_surface() -> None:
    report = compare_discrete_mk_model_ranking(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="region_state",
        taxon_column="taxon",
    )

    assert report.better_model == "equal-rates"
    assert report.likelihood_constant_policy == (
        "continuous-time-markov-pruning-loglikelihood-has-no-extra-normalizing-constant"
    )
    assert report.likelihood_comparison_policy == (
        "relative-aic-and-aicc-ranking-is-permitted-only-when-all-candidate-discrete-mk-models-share-one-pruning-likelihood-policy"
    )
    assert report.model_confidence_weight_basis == "AICc"
    assert report.selected_model_akaike_weight is not None
    assert report.selected_model_akaike_weight > 0.0
    assert report.models_within_delta_aicc_threshold == ["equal-rates"]
    assert [row.model for row in report.rows] == [
        "equal-rates",
        "symmetric",
        "all-rates-different",
    ]
    assert [row.rank for row in report.rows] == [1, 2, 3]
    assert report.rows[0].selected is True
    assert report.rows[0].akaike_weight is not None
    assert "equal-rates" in report.uncertainty_language


@pytest.mark.slow
def test_compare_discrete_mk_model_ranking_from_dataset_matches_path_surface() -> None:
    dataset = load_discrete_dataset(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="binary_state",
        taxon_column="taxon",
    )

    path_report = compare_discrete_mk_model_ranking(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.trait,
        taxon_column=dataset.taxon_column,
    )
    dataset_report = compare_discrete_mk_model_ranking_from_dataset(dataset)

    assert dataset_report.better_model == path_report.better_model
    assert dataset_report.likelihood_constant_policy == (
        path_report.likelihood_constant_policy
    )
    assert [row.model for row in dataset_report.rows] == [
        row.model for row in path_report.rows
    ]
    assert [row.rank for row in dataset_report.rows] == [
        row.rank for row in path_report.rows
    ]


@pytest.mark.slow
def test_compare_discrete_mk_model_ranking_marks_infinite_aicc_rows_noncomparable(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "small-tree.nwk"
    traits_path = tmp_path / "small-traits.tsv"
    tree_path.write_text("((a:1,b:1):1,(c:1,d:1):1);\n", encoding="utf-8")
    traits_path.write_text(
        "taxon\tstate\na\tnorth\nb\tsouth\nc\twest\nd\tnorth\n",
        encoding="utf-8",
    )

    report = compare_discrete_mk_model_ranking(
        tree_path,
        traits_path,
        trait="state",
        taxon_column="taxon",
    )

    row_by_model = {row.model: row for row in report.rows}

    assert report.better_model == "equal-rates"
    assert row_by_model["equal-rates"].selected is True
    assert row_by_model["equal-rates"].rank == 1
    assert row_by_model["symmetric"].comparable is False
    assert row_by_model["all-rates-different"].comparable is False
    assert row_by_model["symmetric"].akaike_weight is None
    assert row_by_model["all-rates-different"].akaike_weight is None
    assert (
        row_by_model["symmetric"].comparability_note
        == "sample size is too small to compute finite AICc for this parameter count"
    )
    assert (
        row_by_model["all-rates-different"].comparability_note
        == "sample size is too small to compute finite AICc for this parameter count"
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_matches_governed_geiger_er_kappa_surface() -> None:
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_kappa_branch_sensitive_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="equal-rates",
        transform="kappa",
    )

    assert report.transform_fit is not None
    assert report.parameter_count == 2
    assert math.isclose(
        report.log_likelihood,
        -9.066235312319336,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        report.transform_fit.parameter_value,
        0.9011763252454394,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        _single_allowed_rate(report),
        0.34869516275650086,
        rel_tol=0.0,
        abs_tol=1e-6,
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_marks_governed_geiger_kappa_weak_signal_surface() -> (
    None
):
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_kappa_weak_signal_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="equal-rates",
        transform="kappa",
    )

    assert report.transform_fit is not None
    assert report.parameter_count == 2
    assert report.transform_fit.parameter_value <= 1e-6
    assert math.isclose(
        report.log_likelihood,
        -16.635532333438686,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert any(
        warning.kind == "branch_length_flattening_limit"
        for warning in report.transform_fit.warnings
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_matches_governed_geiger_sym_kappa_missing_surface() -> (
    None
):
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_kappa_missing_three_state_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="symmetric",
        transform="kappa",
    )
    rate_lookup = _allowed_rate_lookup(report)

    assert report.transform_fit is not None
    assert report.parameter_count == 4
    assert report.input_audit.pruned_missing_value_taxa == ["Phy14"]
    assert math.isclose(
        report.log_likelihood,
        -15.153000877556716,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        report.transform_fit.parameter_value,
        1.0,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        rate_lookup[("central", "north")],
        0.733706223323755,
        rel_tol=0.0,
        abs_tol=1e-6,
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_matches_governed_geiger_er_delta_boundary_surface() -> (
    None
):
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_delta_late_change_binary_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="equal-rates",
        transform="delta",
    )

    assert report.transform_fit is not None
    assert report.parameter_count == 2
    assert report.transform_baseline_comparison is not None
    assert report.transform_baseline_comparison.baseline_transform == "untransformed"
    assert math.isclose(
        report.log_likelihood,
        -8.52405,
        rel_tol=0.0,
        abs_tol=1e-5,
    )
    assert math.isclose(
        report.transform_fit.parameter_value,
        3.0,
        rel_tol=0.0,
        abs_tol=1e-6,
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_reports_delta_time_sensitive_review_surface() -> None:
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_delta_time_sensitive_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="equal-rates",
        transform="delta",
    )

    assert report.transform_fit is not None
    assert report.parameter_count == 2
    assert report.transform_baseline_comparison is not None
    assert report.transform_fit.parameter_value >= math.exp(-5.0)
    assert report.transform_fit.parameter_value <= 3.0
    assert report.transform_baseline_comparison.preferred_transform_by_aic == (
        "untransformed"
    )
    assert math.isfinite(report.log_likelihood)


@pytest.mark.slow
def test_fit_discrete_mk_model_matches_governed_geiger_sym_delta_boundary_surface() -> (
    None
):
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_sym_three_state_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="symmetric",
        transform="delta",
    )

    assert report.transform_fit is not None
    assert report.parameter_count == 4
    assert report.transform_baseline_comparison is not None
    assert math.isclose(
        report.log_likelihood,
        -15.03642,
        rel_tol=0.0,
        abs_tol=1e-5,
    )
    assert math.isclose(
        report.transform_fit.parameter_value,
        math.exp(-5.0),
        rel_tol=0.0,
        abs_tol=1e-9,
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_matches_governed_geiger_sym_delta_missing_surface() -> (
    None
):
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_missing_three_state_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="symmetric",
        transform="delta",
    )

    assert report.transform_fit is not None
    assert report.parameter_count == 4
    assert report.input_audit.pruned_missing_value_taxa == ["Phy14"]
    assert report.transform_baseline_comparison is not None
    assert math.isclose(
        report.log_likelihood,
        -15.0279,
        rel_tol=0.0,
        abs_tol=1e-5,
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_matches_governed_geiger_er_early_burst_surface() -> None:
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_early_burst_early_change_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="equal-rates",
        transform="early-burst",
    )

    assert report.transform_fit is not None
    assert report.parameter_count == 2
    assert report.transform_fit.parameter_name == "a"
    assert math.isclose(
        report.log_likelihood,
        -8.532873913569691,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        report.transform_fit.parameter_value,
        2.2291236000336485,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        _single_allowed_rate(report),
        0.018860745234985883,
        rel_tol=0.0,
        abs_tol=1e-6,
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_marks_governed_geiger_early_burst_weak_signal_surface() -> (
    None
):
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_early_burst_weak_signal_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="equal-rates",
        transform="early-burst",
    )

    assert report.transform_fit is not None
    assert report.parameter_count == 2
    assert report.transform_fit.parameter_value == 0.0
    assert math.isclose(
        report.log_likelihood,
        -16.635532333438686,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert any(
        warning.kind == "brownian_like_rate_change"
        for warning in report.transform_fit.warnings
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_reports_governed_geiger_early_burst_late_change_surface() -> (
    None
):
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_early_burst_late_change_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="equal-rates",
        transform="early-burst",
    )

    assert report.transform_fit is not None
    assert report.parameter_count == 2
    assert report.transform_fit.parameter_value < 0.0
    assert report.transform_baseline_comparison is not None
    assert report.transform_baseline_comparison.preferred_transform_by_aic == (
        "untransformed"
    )
    assert math.isfinite(report.log_likelihood)


@pytest.mark.slow
def test_fit_discrete_mk_model_reports_governed_geiger_early_burst_missing_surface() -> (
    None
):
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_early_burst_missing_binary_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="equal-rates",
        transform="early-burst",
    )

    assert report.transform_fit is not None
    assert report.parameter_count == 2
    assert report.input_audit.pruned_missing_value_taxa == ["Phy10"]
    assert report.transform_fit.parameter_value < 0.0
    assert math.isclose(
        report.log_likelihood,
        -15.219263187720406,
        rel_tol=0.0,
        abs_tol=1e-6,
    )


def test_fit_discrete_mk_model_matches_governed_geiger_sym_three_state_surface() -> (
    None
):
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_sym_three_state_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="symmetric",
    )
    rate_lookup = _allowed_rate_lookup(report)

    assert report.optimizer_diagnostics.optimizer_name == "nelder-mead"
    assert math.isclose(
        report.log_likelihood,
        -15.1632793478894,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        rate_lookup[("central", "north")],
        0.7274457297910687,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        rate_lookup[("central", "south")],
        0.9701467744758534,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        rate_lookup[("north", "south")],
        0.3642863446871228,
        rel_tol=0.0,
        abs_tol=1e-6,
    )


def test_fit_discrete_mk_model_matches_governed_geiger_ard_binary_surface() -> None:
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_ard_binary_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="all-rates-different",
    )
    rate_lookup = _allowed_rate_lookup(report)

    assert report.optimizer_diagnostics.optimizer_name == "nelder-mead"
    assert math.isclose(
        report.log_likelihood,
        -10.750446676724945,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        rate_lookup[("0", "1")],
        1.3142787012012467,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        rate_lookup[("1", "0")],
        2.824830939872932,
        rel_tol=0.0,
        abs_tol=1e-6,
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_matches_governed_geiger_ard_missing_surface() -> None:
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_missing_three_state_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="all-rates-different",
    )
    rate_lookup = _allowed_rate_lookup(report)

    assert report.taxon_count == 23
    assert report.input_audit.pruned_missing_value_taxa == ["Phy14"]
    assert report.optimizer_diagnostics.converged is True
    assert report.optimizer_diagnostics.hit_lower_parameter_bound is False
    assert math.isclose(
        report.log_likelihood,
        -14.31297558479869,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        rate_lookup[("central", "north")],
        0.58934415590246,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        rate_lookup[("central", "south")],
        2.3728847457672133,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        rate_lookup[("south", "north")],
        0.8075051753472205,
        rel_tol=0.0,
        abs_tol=1e-6,
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_marks_governed_geiger_ard_four_state_review_surface() -> (
    None
):
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_ard_four_state_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="all-rates-different",
    )
    rate_lookup = _allowed_rate_lookup(report)

    assert report.parameter_count == 12
    assert report.baseline_comparison is not None
    assert report.baseline_comparison.preferred_model_by_aic == "equal-rates"
    assert report.optimizer_diagnostics.converged is False
    assert report.optimizer_diagnostics.hit_lower_parameter_bound is True
    assert report.optimizer_diagnostics.hit_upper_parameter_bound is False
    assert any(
        "weakly identified" in warning for warning in report.input_audit.warnings
    )
    assert any(
        "equal-rates baseline remains preferred" in warning
        for warning in report.input_audit.warnings
    )
    assert math.isclose(
        rate_lookup[("east", "north")],
        0.3807734820867481,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        rate_lookup[("south", "east")],
        1.996892420461163,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(
        rate_lookup[("west", "south")],
        7.648584079856097,
        rel_tol=0.0,
        abs_tol=1e-6,
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_recovers_binary_ard_known_truth(
    tmp_path: Path,
) -> None:
    tree_path = fixture(
        "example_tree_phytools_ultrametric_one_hundred_twenty_eight_taxa.nwk"
    )
    traits_path = _write_simulated_state_table(
        tmp_path / "binary-ard.tsv",
        tip_states=_simulate_directional_tip_states(
            tree_path,
            rates_by_state={
                "0": {"1": 0.8},
                "1": {"0": 0.2},
            },
            root_state="0",
            seed=71000,
        ),
    )

    report = fit_discrete_mk_model(
        tree_path,
        traits_path,
        trait="state",
        model="all-rates-different",
    )

    rate_lookup = _allowed_rate_lookup(report)

    assert report.parameter_count == 2
    assert report.optimizer_diagnostics.converged is True
    assert report.optimizer_diagnostics.hit_lower_parameter_bound is False
    assert rate_lookup[("0", "1")] > rate_lookup[("1", "0")]
    assert math.isclose(rate_lookup[("0", "1")], 0.8, rel_tol=0.0, abs_tol=0.5)
    assert math.isclose(rate_lookup[("1", "0")], 0.2, rel_tol=0.0, abs_tol=0.5)


def test_fit_discrete_mk_model_recovers_multistate_er_known_truth() -> None:
    report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="region_state",
        taxon_column="taxon",
        model="equal-rates",
    )

    assert math.isclose(_single_allowed_rate(report), 0.15, rel_tol=0.0, abs_tol=0.05)


@pytest.mark.slow
def test_fit_discrete_mk_model_reports_multistate_sym_known_truth_limits(
    tmp_path: Path,
) -> None:
    simulation = simulate_discrete_traits(
        fixture("example_tree_phytools_ultrametric_one_hundred_twenty_eight_taxa.nwk"),
        states=["north", "south", "east", "west"],
        transition_rate=0.45,
        root_state="north",
        seed=51024,
    )
    traits_path = write_discrete_trait_table(
        tmp_path / "multistate-sym.tsv", simulation
    )

    report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_one_hundred_twenty_eight_taxa.nwk"),
        traits_path,
        trait="state",
        model="symmetric",
    )

    allowed_rates = _allowed_rate_lookup(report)

    assert report.parameter_count == 6
    assert report.optimizer_diagnostics.converged is False
    assert report.optimizer_diagnostics.hit_lower_parameter_bound is True
    assert report.baseline_comparison is not None
    assert any(
        "weakly identified" in warning for warning in report.input_audit.warnings
    )
    assert any(
        "equal-rates baseline remains preferred" in warning
        for warning in report.input_audit.warnings
    )
    for (source_state, target_state), rate in allowed_rates.items():
        reverse_rate = allowed_rates[(target_state, source_state)]
        assert math.isclose(rate, reverse_rate, rel_tol=0.0, abs_tol=1e-12)
    assert allowed_rates[("east", "south")] > 1.0
    assert allowed_rates[("north", "south")] > 0.1
    assert allowed_rates[("east", "north")] < 1e-3


def test_fit_discrete_mk_model_marks_overparameterized_symmetric_surface(
    tmp_path: Path,
) -> None:
    traits_path = tmp_path / "overparameterized-symmetric.tsv"
    traits_path.write_text(
        "\n".join(
            [
                "taxon\tstate",
                "A\talpha",
                "B\tbeta",
                "C\tgamma",
                "D\tdelta",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = fit_discrete_mk_model(
        fixture("example_tree.nwk"),
        traits_path,
        trait="state",
        taxon_column="taxon",
        model="symmetric",
    )

    assert report.taxon_count == 4
    assert report.parameter_count == 6
    assert report.overparameterized is True
    assert any(
        warning.startswith("the discrete Mk likelihood fit is likely overparameterized")
        for warning in report.input_audit.warnings
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_marks_sparse_shared_symmetric_surface() -> None:
    fixture_entry = get_shared_geiger_discrete_fixture(
        "geiger_discrete_sparse_six_state_twenty_four_taxa"
    )

    report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="symmetric",
    )

    assert report.parameter_count == 15
    assert report.input_audit.sparse_states == ["f"]
    assert report.optimizer_diagnostics.converged is False
    assert report.optimizer_diagnostics.hit_lower_parameter_bound is True
    assert report.baseline_comparison is not None
    assert report.baseline_comparison.preferred_model_by_aic == "equal-rates"
    assert any(
        "represented by fewer than two taxa" in warning
        for warning in report.input_audit.warnings
    )
    assert any(
        "equal-rates baseline remains preferred" in warning
        for warning in report.input_audit.warnings
    )
    assert any(
        "weakly identified" in warning for warning in report.input_audit.warnings
    )


@pytest.mark.slow
def test_fit_discrete_mk_model_marks_overparameterized_ard_surface(
    tmp_path: Path,
) -> None:
    traits_path = tmp_path / "overparameterized-ard.tsv"
    traits_path.write_text(
        "\n".join(
            [
                "taxon\tstate",
                "A\talpha",
                "B\tbeta",
                "C\tgamma",
                "D\tdelta",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = fit_discrete_mk_model(
        fixture("example_tree.nwk"),
        traits_path,
        trait="state",
        taxon_column="taxon",
        model="all-rates-different",
    )

    assert report.taxon_count == 4
    assert report.parameter_count == 12
    assert report.overparameterized is True
    assert any(
        warning.startswith("the discrete Mk likelihood fit is likely overparameterized")
        for warning in report.input_audit.warnings
    )
