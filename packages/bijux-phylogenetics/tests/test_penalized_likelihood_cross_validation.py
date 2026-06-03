from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.phylo.dating as dating_api
from bijux_phylogenetics.phylo.dating import (
    cross_validate_penalized_likelihood_smoothing_from_metadata,
)
from bijux_phylogenetics.phylo.dating import (
    cross_validation as cross_validation_module,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_package_dating_gateway_exports_penalized_cross_validation_surface() -> None:
    assert (
        dating_api.cross_validate_penalized_likelihood_smoothing_from_metadata
        is cross_validate_penalized_likelihood_smoothing_from_metadata
    )


def test_penalized_cross_validation_selects_candidate_with_lowest_rmse() -> None:
    report = cross_validate_penalized_likelihood_smoothing_from_metadata(
        fixture(
            "trees",
            "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
        ),
        fixture(
            "metadata",
            "penalized_likelihood_cross_validation_tip_dates_5_taxa.tsv",
        ),
        fixture(
            "metadata",
            "penalized_likelihood_cross_validation_calibrations_5_taxa.tsv",
        ),
        smoothing_parameters=[0.01, 0.1, 1.0, 10.0, 100.0],
    )

    assert report.tip_count == 5
    assert report.internal_node_count == 4
    assert report.branch_count == 8
    assert report.usable_calibration_count == 4
    assert report.candidate_count == 5
    assert report.selected_smoothing_parameter == pytest.approx(0.01, abs=1e-12)
    assert len(report.calibration_rows) == 4
    assert len(report.candidate_rows) == 5
    assert len(report.prediction_rows) == 20
    assert min(
        row.root_mean_squared_error for row in report.candidate_rows
    ) == pytest.approx(report.selected_root_mean_squared_error, abs=1e-15)
    selected_rows = [row for row in report.candidate_rows if row.selected]
    assert [row.smoothing_parameter for row in selected_rows] == [0.01]
    rows_by_taxa = {
        tuple(row.descendant_taxa): row for row in report.selected_fit.node_rows
    }
    assert rows_by_taxa[("A", "B")].estimated_date == pytest.approx(1996.0, abs=1e-9)
    assert rows_by_taxa[("A", "B", "C")].estimated_date == pytest.approx(
        1990.0,
        abs=1e-9,
    )
    assert rows_by_taxa[("D", "E")].estimated_date == pytest.approx(1994.0, abs=1e-9)
    assert rows_by_taxa[("A", "B", "C", "D", "E")].estimated_date == pytest.approx(
        1980.0,
        abs=1e-9,
    )


def test_penalized_cross_validation_rejects_infeasible_calibrations_before_optimization(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_if_called(*args, **kwargs):  # pragma: no cover - defensive assertion
        raise AssertionError(
            "penalized optimizer should not run for infeasible calibrations"
        )

    monkeypatch.setattr(
        cross_validation_module,
        "fit_penalized_likelihood_dating",
        fail_if_called,
    )

    with pytest.raises(
        PhylogeneticsError,
        match="dating calibrations are infeasible",
    ):
        cross_validate_penalized_likelihood_smoothing_from_metadata(
            fixture(
                "trees",
                "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
            ),
            fixture(
                "metadata",
                "penalized_likelihood_cross_validation_tip_dates_5_taxa.tsv",
            ),
            fixture(
                "metadata",
                "penalized_likelihood_cross_validation_contradictory_calibrations_5_taxa.tsv",
            ),
            smoothing_parameters=[0.01, 0.1, 1.0],
        )
