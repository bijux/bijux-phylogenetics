from __future__ import annotations

import math
from pathlib import Path
import random

from bijux_phylogenetics.comparative.common import (
    build_brownian_covariance_matrix,
    stable_covariance,
)
from bijux_phylogenetics.comparative.evolutionary_modes import (
    transform_tree_for_evolutionary_mode,
)
from bijux_phylogenetics.comparative.traits.rate_through_time import (
    summarize_trait_rate_through_time,
    write_trait_rate_through_time_exclusion_table,
    write_trait_rate_through_time_interval_table,
    write_trait_rate_through_time_summary_table,
)
from bijux_phylogenetics.io.trees import load_tree

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


def test_summarize_trait_rate_through_time_reports_interval_rates() -> None:
    report = summarize_trait_rate_through_time(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        interval_count=4,
    )

    assert report.tree_taxon_count == 4
    assert report.analyzed_taxa == ["A", "B", "C", "D"]
    assert report.analyzed_taxon_count == 4
    assert report.interval_count == 4
    assert len(report.interval_rows) == 4
    assert report.nonempty_interval_count >= 2
    assert report.tree_depth > 0.0
    assert report.ancestral_model == "brownian"
    assert report.trend_direction in {"slowdown", "acceleration", "stable"}
    assert any(row.estimated_rate is not None for row in report.interval_rows)


def test_summarize_trait_rate_through_time_tracks_excluded_taxa() -> None:
    report = summarize_trait_rate_through_time(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_continuous_evolution_missing.tsv"),
        trait="response_growth",
        interval_count=4,
    )

    assert report.analyzed_taxa == ["A", "D", "E", "F"]
    assert {row.taxon: row.reason for row in report.excluded_taxa} == {
        "B": "missing_trait_value",
        "C": "non_numeric_trait_value",
        "G": "absent_from_tree",
    }


def test_trait_rate_through_time_writers_emit_review_ledgers(tmp_path: Path) -> None:
    report = summarize_trait_rate_through_time(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_continuous_evolution_missing.tsv"),
        trait="response_growth",
        interval_count=4,
    )
    summary_out = tmp_path / "trait-rate-through-time-summary.tsv"
    intervals_out = tmp_path / "trait-rate-through-time-intervals.tsv"
    excluded_out = tmp_path / "trait-rate-through-time-excluded.tsv"

    write_trait_rate_through_time_summary_table(summary_out, report)
    write_trait_rate_through_time_interval_table(intervals_out, report)
    write_trait_rate_through_time_exclusion_table(excluded_out, report)

    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    interval_rows = intervals_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("trait\ttaxon_column\ttree_taxon_count")
    assert interval_rows[0].startswith("interval_index\tstart_depth\tend_depth")
    assert excluded_rows == [
        "taxon\treason",
        "B\tmissing_trait_value",
        "C\tnon_numeric_trait_value",
        "G\tabsent_from_tree",
    ]


def test_summarize_trait_rate_through_time_detects_simulated_slowdown(
    tmp_path: Path,
) -> None:
    tree_path, traits_path = _write_simulated_trait_dataset(
        tmp_path,
        mode="early-burst",
        seed=1,
    )

    report = summarize_trait_rate_through_time(
        tree_path,
        traits_path,
        trait="trait",
        interval_count=5,
    )

    assert report.trend_direction == "slowdown"
    assert report.earliest_interval_rate is not None
    assert report.latest_interval_rate is not None
    assert report.latest_interval_rate < report.earliest_interval_rate
    assert report.normalized_rate_slope is not None
    assert report.normalized_rate_slope < 0.0


def test_summarize_trait_rate_through_time_reports_stable_on_simulated_brownian(
    tmp_path: Path,
) -> None:
    tree_path, traits_path = _write_simulated_trait_dataset(
        tmp_path,
        mode="brownian",
        seed=19,
    )

    report = summarize_trait_rate_through_time(
        tree_path,
        traits_path,
        trait="trait",
        interval_count=5,
    )

    assert report.trend_direction in {"stable", "insufficient_data"}
    assert report.nonempty_interval_count >= 2
    assert all(
        row.estimated_rate is None or row.estimated_rate >= 0.0
        for row in report.interval_rows
    )


def _write_simulated_trait_dataset(
    tmp_path: Path,
    *,
    mode: str,
    seed: int,
) -> tuple[Path, Path]:
    tree_path = tmp_path / f"simulated-{mode}-tree.nwk"
    traits_path = tmp_path / f"simulated-{mode}-traits.tsv"
    tree_path.write_text(
        (
            "((((A:0.8,B:0.4):0.7,(C:0.6,D:0.3):0.5):1.2,"
            "((E:0.7,F:0.35):0.8,(G:0.5,H:0.25):0.9):1.0):0.6,"
            "((I:0.9,J:0.45):0.6,(K:0.55,L:0.2):1.1):0.85);"
        ),
        encoding="utf-8",
    )
    base_tree = load_tree(tree_path)
    working_tree = (
        base_tree
        if mode == "brownian"
        else transform_tree_for_evolutionary_mode(
            base_tree,
            mode="early-burst",
            parameter_value=4.0,
        )
    )
    covariance = stable_covariance(
        build_brownian_covariance_matrix(
            working_tree,
            base_tree.tip_names,
        )
    )
    cholesky = _cholesky(covariance)
    rng = random.Random(seed)
    standard_normal = [rng.gauss(0.0, 1.0) for _ in base_tree.tip_names]
    values = [
        3.0
        + sum(
            cholesky[row][column] * standard_normal[column] for column in range(row + 1)
        )
        for row in range(len(base_tree.tip_names))
    ]
    traits_path.write_text(
        "taxon\ttrait\n"
        + "".join(
            f"{taxon}\t{value}\n"
            for taxon, value in zip(base_tree.tip_names, values, strict=True)
        ),
        encoding="utf-8",
    )
    return tree_path, traits_path


def _cholesky(matrix: list[list[float]]) -> list[list[float]]:
    size = len(matrix)
    lower = [[0.0] * size for _ in range(size)]
    for row in range(size):
        for column in range(row + 1):
            partial = sum(
                lower[row][index] * lower[column][index] for index in range(column)
            )
            if row == column:
                lower[row][column] = math.sqrt(max(matrix[row][row] - partial, 1e-12))
            else:
                lower[row][column] = (matrix[row][column] - partial) / lower[column][
                    column
                ]
    return lower
