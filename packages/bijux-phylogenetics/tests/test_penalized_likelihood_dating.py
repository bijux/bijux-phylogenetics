from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.phylo.dating as dating_api
from bijux_phylogenetics.phylo.dating import (
    fit_least_squares_dating_from_metadata,
    fit_penalized_likelihood_dating_from_metadata,
)
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    PhylogeneticsError,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "metadata")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_package_dating_gateway_exports_penalized_likelihood_surface() -> None:
    assert (
        dating_api.fit_penalized_likelihood_dating_from_metadata
        is fit_penalized_likelihood_dating_from_metadata
    )


def test_penalized_likelihood_dating_separates_data_and_smoothing_scores() -> None:
    low_smoothing_report = fit_penalized_likelihood_dating_from_metadata(
        fixture("penalized_likelihood_dating_substitution_tree_4_taxa.nwk"),
        fixture("penalized_likelihood_dating_tip_dates_4_taxa.tsv"),
        smoothing_parameter=0.01,
    )
    high_smoothing_report = fit_penalized_likelihood_dating_from_metadata(
        fixture("penalized_likelihood_dating_substitution_tree_4_taxa.nwk"),
        fixture("penalized_likelihood_dating_tip_dates_4_taxa.tsv"),
        smoothing_parameter=10.0,
    )
    least_squares_report = fit_least_squares_dating_from_metadata(
        fixture("penalized_likelihood_dating_substitution_tree_4_taxa.nwk"),
        fixture("penalized_likelihood_dating_tip_dates_4_taxa.tsv"),
    )

    assert low_smoothing_report.tree_newick == (
        "(((A:0.339411254969543,B:0.329848450049413):0.554256258422041,"
        "C:1.03691851174526):0.831384387633061,D:1.79097738679192);"
    )
    assert low_smoothing_report.tip_count == 4
    assert low_smoothing_report.internal_node_count == 3
    assert low_smoothing_report.branch_count == 6
    assert low_smoothing_report.parameter_count == 10
    assert low_smoothing_report.tree_path == str(
        fixture("penalized_likelihood_dating_substitution_tree_4_taxa.nwk")
    )
    assert low_smoothing_report.metadata_path == str(
        fixture("penalized_likelihood_dating_tip_dates_4_taxa.tsv")
    )
    assert low_smoothing_report.smoothing_parameter == 0.01
    assert low_smoothing_report.total_score == pytest.approx(
        low_smoothing_report.data_score + low_smoothing_report.penalty_score,
        abs=1e-15,
    )
    assert high_smoothing_report.total_score == pytest.approx(
        high_smoothing_report.data_score + high_smoothing_report.penalty_score,
        abs=1e-15,
    )
    assert low_smoothing_report.data_score < high_smoothing_report.data_score
    low_rate_spread = max(
        row.estimated_branch_rate for row in low_smoothing_report.branch_rows
    ) - min(row.estimated_branch_rate for row in low_smoothing_report.branch_rows)
    high_rate_spread = max(
        row.estimated_branch_rate for row in high_smoothing_report.branch_rows
    ) - min(row.estimated_branch_rate for row in high_smoothing_report.branch_rows)
    assert low_rate_spread > high_rate_spread
    assert low_rate_spread > 0.0
    assert low_smoothing_report.root_date > least_squares_report.root_date
    assert (
        low_smoothing_report.dated_tree_newick != least_squares_report.dated_tree_newick
    )
    assert low_smoothing_report.optimizer_name == (
        "bounded-coordinate-search with closed-form penalized log-rate solve"
    )
    assert low_smoothing_report.optimization_pass_count >= 1
    assert low_smoothing_report.function_evaluation_count >= 1
    assert low_smoothing_report.converged is True
    assert len(low_smoothing_report.node_rows) == 7
    assert len(low_smoothing_report.branch_rows) == 6
    assert all(row.estimated_rate > 0.0 for row in low_smoothing_report.node_rows)
    assert all(
        row.estimated_branch_rate > 0.0 for row in low_smoothing_report.branch_rows
    )


def test_penalized_likelihood_dating_requires_positive_smoothing_parameter() -> None:
    with pytest.raises(
        PhylogeneticsError,
        match="requires a strictly positive smoothing parameter",
    ):
        fit_penalized_likelihood_dating_from_metadata(
            fixture("penalized_likelihood_dating_substitution_tree_4_taxa.nwk"),
            fixture("penalized_likelihood_dating_tip_dates_4_taxa.tsv"),
            smoothing_parameter=0.0,
        )


def test_penalized_likelihood_dating_rejects_zero_branch_lengths(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "zero-branch-tree.nwk"
    tree_path.write_text(
        "(((A:0,B:0.329848450049413):0.554256258422041,C:1.03691851174526):0.831384387633061,D:1.79097738679192);",
        encoding="utf-8",
    )

    with pytest.raises(
        InvalidBranchLengthError,
        match="requires strictly positive branch lengths",
    ):
        fit_penalized_likelihood_dating_from_metadata(
            tree_path,
            fixture("penalized_likelihood_dating_tip_dates_4_taxa.tsv"),
        )
