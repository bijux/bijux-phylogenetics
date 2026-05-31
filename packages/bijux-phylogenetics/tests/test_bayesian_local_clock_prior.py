from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.clock_models import (
    build_local_clock_rate_model,
    evaluate_local_clock_tree_log_prior,
    load_local_clock_regime_definitions,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    PhylogeneticsError,
    UnrootedTreeError,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def load_rooted_tree_fixture(name: str):
    tree = load_tree(fixture("trees", name))
    tree.rooted = True
    return tree


def build_local_clock_substitution_tree():
    dated_tree = load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk")
    substitution_tree = dated_tree.copy()
    rate_by_descendant_taxa = {
        ("A", "B", "C"): 0.5,
        ("A", "B"): 0.2,
        ("A",): 0.2,
        ("B",): 0.2,
        ("C",): 0.1,
        ("D",): 0.1,
    }
    for _parent, child in substitution_tree.iter_edges():
        child.branch_length = (
            float(child.branch_length or 0.0)
            * rate_by_descendant_taxa[tuple(child.descendant_taxa)]
        )
    return dated_tree, substitution_tree


def test_local_clock_prior_evaluates_background_and_foreground_regimes() -> None:
    dated_tree, substitution_tree = build_local_clock_substitution_tree()
    regime_definitions = load_local_clock_regime_definitions(
        dated_tree,
        fixture("metadata", "local_clock_regimes_4_taxa.tsv"),
    )
    report = evaluate_local_clock_tree_log_prior(
        substitution_tree,
        dated_tree,
        regime_definitions,
        build_local_clock_rate_model(
            background_clock_rate=0.1,
            regime_clock_rates={"ab_clade": 0.2, "abc_stem": 0.5},
            log_standard_deviation=0.5,
        ),
    )

    assert report.family == "local-clock"
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.regime_count == 3
    assert report.branch_count == 6
    assert math.isfinite(report.total_log_prior)

    regime_rows_by_id = {row.regime_id: row for row in report.regime_rows}
    assert set(regime_rows_by_id) == {"background", "ab_clade", "abc_stem"}
    assert regime_rows_by_id["background"].branch_count == 2
    assert regime_rows_by_id["ab_clade"].branch_count == 3
    assert regime_rows_by_id["abc_stem"].branch_count == 1
    assert regime_rows_by_id["background"].target_kind == "background"
    assert regime_rows_by_id["ab_clade"].target_kind == "clade"
    assert regime_rows_by_id["abc_stem"].target_kind == "branch"
    assert math.isclose(
        regime_rows_by_id["background"].class_clock_rate,
        0.1,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        regime_rows_by_id["ab_clade"].class_clock_rate,
        0.2,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        regime_rows_by_id["abc_stem"].class_clock_rate,
        0.5,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    rows_by_descendant_taxa = {
        tuple(row.descendant_taxa): row for row in report.branch_rows
    }
    assert rows_by_descendant_taxa[("A", "B", "C")].regime_id == "abc_stem"
    assert rows_by_descendant_taxa[("A", "B")].regime_id == "ab_clade"
    assert rows_by_descendant_taxa[("A",)].regime_id == "ab_clade"
    assert rows_by_descendant_taxa[("B",)].regime_id == "ab_clade"
    assert rows_by_descendant_taxa[("C",)].regime_id == "background"
    assert rows_by_descendant_taxa[("D",)].regime_id == "background"


def test_local_clock_prior_changes_when_foreground_rate_changes() -> None:
    dated_tree, substitution_tree = build_local_clock_substitution_tree()
    regime_definitions = load_local_clock_regime_definitions(
        dated_tree,
        fixture("metadata", "local_clock_regimes_4_taxa.tsv"),
    )
    baseline_report = evaluate_local_clock_tree_log_prior(
        substitution_tree,
        dated_tree,
        regime_definitions,
        build_local_clock_rate_model(
            background_clock_rate=0.1,
            regime_clock_rates={"ab_clade": 0.2, "abc_stem": 0.5},
            log_standard_deviation=0.5,
        ),
    )
    shifted_report = evaluate_local_clock_tree_log_prior(
        substitution_tree,
        dated_tree,
        regime_definitions,
        build_local_clock_rate_model(
            background_clock_rate=0.1,
            regime_clock_rates={"ab_clade": 0.3, "abc_stem": 0.5},
            log_standard_deviation=0.5,
        ),
    )

    assert not math.isclose(
        baseline_report.total_log_prior,
        shifted_report.total_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    baseline_rows = {row.regime_id: row for row in baseline_report.regime_rows}
    shifted_rows = {row.regime_id: row for row in shifted_report.regime_rows}
    assert not math.isclose(
        baseline_rows["ab_clade"].log_prior_contribution,
        shifted_rows["ab_clade"].log_prior_contribution,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        baseline_rows["background"].log_prior_contribution,
        shifted_rows["background"].log_prior_contribution,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_local_clock_prior_rejects_ambiguous_overlapping_regimes() -> None:
    dated_tree, substitution_tree = build_local_clock_substitution_tree()
    conflicting_regime_definitions = load_local_clock_regime_definitions(
        dated_tree,
        fixture("metadata", "local_clock_conflicting_regimes_4_taxa.tsv"),
    )

    with pytest.raises(PhylogeneticsError, match="ambiguous"):
        evaluate_local_clock_tree_log_prior(
            substitution_tree,
            dated_tree,
            conflicting_regime_definitions,
            build_local_clock_rate_model(
                background_clock_rate=0.1,
                regime_clock_rates={"ab_alpha": 0.2, "ab_beta": 0.25},
                log_standard_deviation=0.5,
            ),
        )


@pytest.mark.parametrize(
    ("builder_kwargs", "message"),
    [
        (
            {
                "background_clock_rate": 0.0,
                "regime_clock_rates": {"ab_clade": 0.2},
                "log_standard_deviation": 0.5,
            },
            "requires a strictly positive finite background_clock_rate",
        ),
        (
            {
                "background_clock_rate": 0.1,
                "regime_clock_rates": {"background": 0.2},
                "log_standard_deviation": 0.5,
            },
            "reserves regime_id 'background'",
        ),
        (
            {
                "background_clock_rate": 0.1,
                "regime_clock_rates": {"ab_clade": 0.2},
                "log_standard_deviation": 0.0,
            },
            "requires a strictly positive finite log standard deviation",
        ),
    ],
)
def test_local_clock_prior_rejects_invalid_model_parameters(
    builder_kwargs: dict[str, object],
    message: str,
) -> None:
    with pytest.raises(PhylogeneticsError, match=message):
        build_local_clock_rate_model(**builder_kwargs)


def test_local_clock_prior_rejects_unrooted_missing_rates_and_zero_duration() -> None:
    dated_tree, substitution_tree = build_local_clock_substitution_tree()
    regime_definitions = load_local_clock_regime_definitions(
        dated_tree,
        fixture("metadata", "local_clock_regimes_4_taxa.tsv"),
    )
    unrooted_substitution_tree = substitution_tree.copy()
    unrooted_substitution_tree.rooted = False
    zero_duration_dated_tree = dated_tree.copy()
    zero_duration_dated_tree.root.children[1].branch_length = 0.0

    with pytest.raises(
        UnrootedTreeError, match="requires one rooted substitution tree"
    ):
        evaluate_local_clock_tree_log_prior(
            unrooted_substitution_tree,
            dated_tree,
            regime_definitions,
            build_local_clock_rate_model(
                background_clock_rate=0.1,
                regime_clock_rates={"ab_clade": 0.2, "abc_stem": 0.5},
                log_standard_deviation=0.5,
            ),
        )

    with pytest.raises(
        PhylogeneticsError,
        match="requires one clock rate for every named regime",
    ):
        evaluate_local_clock_tree_log_prior(
            substitution_tree,
            dated_tree,
            regime_definitions,
            build_local_clock_rate_model(
                background_clock_rate=0.1,
                regime_clock_rates={"ab_clade": 0.2},
                log_standard_deviation=0.5,
            ),
        )

    with pytest.raises(
        InvalidBranchLengthError,
        match="requires strictly positive dated branch durations",
    ):
        evaluate_local_clock_tree_log_prior(
            substitution_tree,
            zero_duration_dated_tree,
            regime_definitions,
            build_local_clock_rate_model(
                background_clock_rate=0.1,
                regime_clock_rates={"ab_clade": 0.2, "abc_stem": 0.5},
                log_standard_deviation=0.5,
            ),
        )
