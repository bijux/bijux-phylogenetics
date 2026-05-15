from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import load_discrete_dataset
from bijux_phylogenetics.comparative.discrete_mk import (
    fit_discrete_mk_model,
    fit_discrete_mk_model_from_dataset,
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


def test_fit_discrete_mk_model_recovers_multistate_er_known_truth() -> None:
    report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="region_state",
        taxon_column="taxon",
        model="equal-rates",
    )

    assert math.isclose(_single_allowed_rate(report), 0.15, rel_tol=0.0, abs_tol=0.05)


def test_fit_discrete_mk_model_recovers_multistate_sym_known_truth(
    tmp_path: Path,
) -> None:
    simulation = simulate_discrete_traits(
        fixture("example_tree_phytools_ultrametric_one_hundred_twenty_eight_taxa.nwk"),
        states=["north", "south", "east", "west"],
        transition_rate=0.45,
        root_state="north",
        seed=51024,
    )
    traits_path = write_discrete_trait_table(tmp_path / "multistate-sym.tsv", simulation)

    report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_one_hundred_twenty_eight_taxa.nwk"),
        traits_path,
        trait="state",
        model="symmetric",
    )

    allowed_rates = _allowed_rate_lookup(report)

    assert report.parameter_count == 6
    assert math.isclose(
        allowed_rates[("north", "south")],
        allowed_rates[("south", "north")],
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert all(
        math.isclose(rate, 0.45, rel_tol=0.0, abs_tol=0.25)
        for rate in allowed_rates.values()
    )


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
        warning.startswith(
            "the discrete Mk likelihood fit is likely overparameterized"
        )
        for warning in report.input_audit.warnings
    )
