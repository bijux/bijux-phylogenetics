from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.dating import load_fixed_dating_calibrations
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_load_fixed_dating_calibrations_resolves_internal_nodes() -> None:
    calibration_rows = load_fixed_dating_calibrations(
        fixture(
            "trees",
            "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
        ),
        fixture(
            "metadata",
            "penalized_likelihood_cross_validation_calibrations_5_taxa.tsv",
        ),
    )

    assert [row.calibration_id for row in calibration_rows] == [
        "cal-ab",
        "cal-abc",
        "cal-de",
        "cal-root",
    ]
    assert [row.fixed_date for row in calibration_rows] == [
        1996.0,
        1990.0,
        1994.0,
        1980.0,
    ]
    rows_by_taxa = {tuple(row.descendant_taxa): row for row in calibration_rows}
    assert rows_by_taxa[("A", "B")].node_kind == "internal"
    assert rows_by_taxa[("A", "B", "C", "D", "E")].node_kind == "root"


def test_load_fixed_dating_calibrations_rejects_interval_constraints(
    tmp_path: Path,
) -> None:
    calibration_path = tmp_path / "interval-calibrations.tsv"
    calibration_path.write_text(
        "\n".join(
            [
                "calibration_id\tclade_name\ttaxa\tminimum_age\tmaximum_age\tdistribution",
                "cal-ab\t\tA|B\t1995\t1996\tuniform",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(
        PhylogeneticsError,
        match="requires fixed calibrations with identical minimum_age and maximum_age",
    ):
        load_fixed_dating_calibrations(
            fixture(
                "trees",
                "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
            ),
            calibration_path,
        )


def test_load_fixed_dating_calibrations_rejects_contradictory_fixed_constraints() -> (
    None
):
    with pytest.raises(
        PhylogeneticsError,
        match="dating calibrations are infeasible",
    ):
        load_fixed_dating_calibrations(
            fixture(
                "trees",
                "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
            ),
            fixture(
                "metadata",
                "penalized_likelihood_cross_validation_contradictory_calibrations_5_taxa.tsv",
            ),
        )
