from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.compare as compare_api
from bijux_phylogenetics.compare import (
    CladeAgeComparisonRow,
    DateAwareTreeComparisonReport,
    compare_clade_ages,
)
from bijux_phylogenetics.runtime.errors import NonUltrametricTreeError


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_package_compare_gateway_exports_date_aware_tree_comparison_surface() -> None:
    assert compare_api.DateAwareTreeComparisonReport is DateAwareTreeComparisonReport
    assert compare_api.CladeAgeComparisonRow is CladeAgeComparisonRow
    assert compare_api.compare_clade_ages is compare_clade_ages


def test_compare_clade_ages_reports_age_differences_when_rf_distance_is_zero() -> None:
    report = compare_clade_ages(
        fixture("strict_clock_time_tree_4_taxa.nwk"),
        fixture("relaxed_rate_summary_dated_tree_4_taxa.nwk"),
    )

    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.comparison_scope == "full-taxa"
    assert report.topology.topology_equal is True
    assert report.topology.robinson_foulds_distance == 0
    assert report.left_root_age == pytest.approx(3.0, abs=1e-12)
    assert report.right_root_age == pytest.approx(8.0, abs=1e-12)
    assert report.age_rmse == pytest.approx(3.16227766016838, abs=1e-12)
    assert report.mean_absolute_age_difference == pytest.approx(
        2.666666666666667,
        abs=1e-12,
    )
    assert report.max_absolute_age_difference == pytest.approx(5.0, abs=1e-12)
    assert report.unstable_clade_count == 1
    assert [
        (
            row.clade_id,
            row.left_age,
            row.right_age,
            row.age_difference,
            row.unstable_age,
        )
        for row in report.clade_rows
    ] == [
        ("A|B", 1.0, 2.0, 1.0, False),
        ("A|B|C", 2.0, 4.0, 2.0, False),
        ("A|B|C|D", 3.0, 8.0, 5.0, True),
    ]


def test_compare_clade_ages_prunes_to_shared_taxa_before_age_matching(
    tmp_path: Path,
) -> None:
    left_path = tmp_path / "left-extra-taxon.nwk"
    right_path = tmp_path / "right-shared-taxa.nwk"
    left_path.write_text("((((A:1,B:1):1,C:2):1,D:3):1,E:4);\n", encoding="utf-8")
    right_path.write_text("(((A:1,B:1):1,C:2):1,D:3);\n", encoding="utf-8")

    report = compare_clade_ages(left_path, right_path)

    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.left_only_taxa == ["E"]
    assert report.right_only_taxa == []
    assert report.comparison_scope == "pruned-to-shared-taxa"
    assert report.topology.topology_equal is True
    assert report.age_rmse == pytest.approx(0.0, abs=1e-12)
    assert report.matched_clade_count == 3
    assert [row.clade_id for row in report.clade_rows] == [
        "A|B",
        "A|B|C",
        "A|B|C|D",
    ]


def test_compare_clade_ages_rejects_non_ultrametric_trees() -> None:
    with pytest.raises(NonUltrametricTreeError) as error:
        compare_clade_ages(
            fixture("strict_clock_nonclock_tree_4_taxa.nwk"),
            fixture("strict_clock_nonclock_tree_4_taxa.nwk"),
        )

    assert error.value.code == "date_aware_tree_comparison_requires_ultrametric_tree"
