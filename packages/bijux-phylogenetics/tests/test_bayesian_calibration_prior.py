from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.calibration_priors import (
    evaluate_calibration_tree_log_prior,
    load_calibration_prior_definitions,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import (
    NonUltrametricTreeError,
    PhylogeneticsError,
)

FIXTURES = Path(__file__).parent / "fixtures"
_NORMAL_QUANTILE_97_5 = 1.959963984540054


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def _normal_cdf(value: float, mean: float, standard_deviation: float) -> float:
    return 0.5 * (
        1.0 + math.erf((value - mean) / (standard_deviation * math.sqrt(2.0)))
    )


def _truncated_normal_log_density(
    *,
    value: float,
    lower_bound: float,
    upper_bound: float,
) -> float:
    mean = (lower_bound + upper_bound) / 2.0
    standard_deviation = (upper_bound - lower_bound) / (2.0 * _NORMAL_QUANTILE_97_5)
    normalization_mass = _normal_cdf(
        upper_bound, mean, standard_deviation
    ) - _normal_cdf(
        lower_bound,
        mean,
        standard_deviation,
    )
    return (
        -math.log(standard_deviation)
        - 0.5 * math.log(2.0 * math.pi)
        - math.log(normalization_mass)
    )


def _truncated_lognormal_log_density(
    *,
    value: float,
    lower_bound: float,
    upper_bound: float,
) -> float:
    log_mean = (math.log(lower_bound) + math.log(upper_bound)) / 2.0
    log_standard_deviation = (math.log(upper_bound) - math.log(lower_bound)) / (
        2.0 * _NORMAL_QUANTILE_97_5
    )
    normalization_mass = _normal_cdf(
        math.log(upper_bound),
        log_mean,
        log_standard_deviation,
    ) - _normal_cdf(
        math.log(lower_bound),
        log_mean,
        log_standard_deviation,
    )
    return (
        -math.log(value)
        - math.log(log_standard_deviation)
        - 0.5 * math.log(2.0 * math.pi)
        - math.log(normalization_mass)
    )


def test_load_calibration_prior_definitions_resolves_mrca_nodes_and_prior_families() -> (
    None
):
    prior_definitions = load_calibration_prior_definitions(
        fixture("trees", "calibration_prior_time_tree_6_taxa.nwk"),
        fixture("metadata", "calibration_priors_6_taxa.tsv"),
    )

    assert [row.calibration_id for row in prior_definitions] == [
        "cal-ab",
        "cal-abc",
        "cal-abcde",
        "cal-de",
        "cal-root",
    ]
    definitions_by_id = {row.calibration_id: row for row in prior_definitions}
    assert definitions_by_id["cal-root"].family == "fixed"
    assert definitions_by_id["cal-root"].translated is True
    assert definitions_by_id["cal-abcde"].family == "uniform"
    assert definitions_by_id["cal-abc"].family == "normal"
    assert definitions_by_id["cal-de"].family == "lognormal"
    assert definitions_by_id["cal-ab"].family == "offset-exponential"
    assert definitions_by_id["cal-ab"].translated is True
    assert definitions_by_id["cal-ab"].parameter_values() == {
        "offset_age": 0.75,
        "exponential_mean": 1.0,
    }


def test_evaluate_calibration_tree_log_prior_matches_bounded_and_offset_fixture() -> (
    None
):
    tree = load_tree(fixture("trees", "calibration_prior_time_tree_6_taxa.nwk"))
    prior_definitions = load_calibration_prior_definitions(
        fixture("trees", "calibration_prior_time_tree_6_taxa.nwk"),
        fixture("metadata", "calibration_priors_6_taxa.tsv"),
    )

    report = evaluate_calibration_tree_log_prior(tree, prior_definitions)

    expected_total_log_prior = (
        0.0
        - math.log(4.5 - 3.5)
        + _truncated_normal_log_density(
            value=2.0,
            lower_bound=1.0,
            upper_bound=3.0,
        )
        + _truncated_lognormal_log_density(
            value=1.5,
            lower_bound=1.0,
            upper_bound=2.25,
        )
        + (-math.log(1.0) - ((1.0 - 0.75) / 1.0))
    )

    assert report.tip_count == 6
    assert report.internal_node_count == 5
    assert report.calibration_count == 5
    assert report.translated_calibration_count == 2
    assert math.isclose(
        report.total_log_prior,
        expected_total_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    rows_by_id = {row.calibration_id: row for row in report.calibration_rows}
    assert rows_by_id["cal-root"].node_age == pytest.approx(5.0, abs=1e-12)
    assert rows_by_id["cal-abcde"].log_prior_contribution == pytest.approx(
        0.0,
        abs=1e-12,
    )
    assert rows_by_id["cal-abc"].node_age == pytest.approx(2.0, abs=1e-12)
    assert rows_by_id["cal-de"].node_age == pytest.approx(1.5, abs=1e-12)
    assert rows_by_id["cal-ab"].node_age == pytest.approx(1.0, abs=1e-12)
    assert rows_by_id["cal-ab"].translated is True


def test_load_calibration_prior_definitions_rejects_contradictory_fixture_before_scoring() -> (
    None
):
    with pytest.raises(
        PhylogeneticsError,
        match="calibration prior constraints are infeasible",
    ):
        load_calibration_prior_definitions(
            fixture("trees", "calibration_prior_time_tree_6_taxa.nwk"),
            fixture("metadata", "calibration_priors_contradictory_6_taxa.tsv"),
        )


def test_evaluate_calibration_tree_log_prior_rejects_non_ultrametric_tree() -> None:
    tree = load_tree(fixture("trees", "calibration_prior_time_tree_6_taxa.nwk"))
    tree.root.children[0].children[1].branch_length = 2.1
    prior_definitions = load_calibration_prior_definitions(
        fixture("trees", "calibration_prior_time_tree_6_taxa.nwk"),
        fixture("metadata", "calibration_priors_6_taxa.tsv"),
    )

    with pytest.raises(
        NonUltrametricTreeError,
        match="requires an ultrametric tree",
    ):
        evaluate_calibration_tree_log_prior(tree, prior_definitions)
