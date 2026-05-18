from __future__ import annotations

import math
from pathlib import Path
import random

import pytest

from bijux_phylogenetics.ancestral.common import load_discrete_dataset
from bijux_phylogenetics.comparative.discrete_mk import (
    fit_discrete_mk_model,
    fit_discrete_mk_model_from_dataset,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.fixtures import (
    get_shared_geiger_discrete_fixture,
    get_shared_phytools_comparative_fixture,
)
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
