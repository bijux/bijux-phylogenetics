from __future__ import annotations

import math
from pathlib import Path
import random

from bijux_phylogenetics.comparative.common import (
    build_brownian_covariance_matrix,
    stable_covariance,
)
from bijux_phylogenetics.comparative.continuous import (
    summarize_early_burst_trait_evolution,
    write_early_burst_rate_change_profile_table,
    write_early_burst_trait_evolution_comparison_table,
    write_early_burst_trait_evolution_exclusion_table,
    write_early_burst_trait_evolution_summary_table,
)
from bijux_phylogenetics.comparative.evolutionary_modes import (
    transform_tree_for_evolutionary_mode,
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


def test_summarize_early_burst_trait_evolution_reports_comparison_context() -> None:
    report = summarize_early_burst_trait_evolution(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )

    assert report.tree_taxon_count == 4
    assert report.analyzed_taxon_count == 4
    assert report.excluded_taxa == []
    assert report.rate_change == 0.0
    assert report.better_model == "brownian"
    assert len(report.profile_rows) == 161
    assert [row.model for row in report.comparison_rows] == [
        "brownian",
        "early-burst",
        "ornstein-uhlenbeck",
    ]
    assert [warning.kind for warning in report.identifiability_warnings] == [
        "boundary_rate_change",
        "brownian_like_rate_change",
        "comparison_not_preferred",
    ]
    assert report.confidence_intervals[0].name == "rate_change"
    assert report.confidence_intervals[0].lower_95 == 0.0
    assert report.confidence_intervals[0].upper_95 == 13.4375


def test_summarize_early_burst_trait_evolution_records_missing_and_invalid_taxa() -> (
    None
):
    report = summarize_early_burst_trait_evolution(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_continuous_evolution_missing.tsv"),
        trait="response_growth",
    )

    assert report.analyzed_taxa == ["A", "D", "E", "F"]
    assert {row.taxon: row.reason for row in report.excluded_taxa} == {
        "B": "missing_trait_value",
        "C": "non_numeric_trait_value",
        "G": "absent_from_tree",
    }
    assert (
        "one or more overlapping taxa have missing trait values and will be pruned"
        in report.warnings
    )
    assert (
        "one or more overlapping taxa have non-numeric trait values and will be pruned"
        in report.warnings
    )
    assert "trait table contains taxa absent from the tree" in report.warnings


def test_early_burst_trait_evolution_writers_emit_summary_review_ledgers(
    tmp_path: Path,
) -> None:
    report = summarize_early_burst_trait_evolution(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_continuous_evolution_missing.tsv"),
        trait="response_growth",
    )
    summary_out = tmp_path / "early-burst-summary.tsv"
    excluded_out = tmp_path / "early-burst-excluded.tsv"
    comparison_out = tmp_path / "early-burst-comparison.tsv"
    profile_out = tmp_path / "early-burst-profile.tsv"

    write_early_burst_trait_evolution_summary_table(summary_out, report)
    write_early_burst_trait_evolution_exclusion_table(excluded_out, report)
    write_early_burst_trait_evolution_comparison_table(comparison_out, report)
    write_early_burst_rate_change_profile_table(profile_out, report)

    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    comparison_rows = comparison_out.read_text(encoding="utf-8").splitlines()
    profile_rows = profile_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("trait\ttaxon_column\ttree_taxon_count")
    assert excluded_rows == [
        "taxon\treason",
        "B\tmissing_trait_value",
        "C\tnon_numeric_trait_value",
        "G\tabsent_from_tree",
    ]
    assert comparison_rows[0].startswith("row_kind\tmodel\tcomparison_id")
    assert comparison_rows[1].startswith("model_fit\tbrownian\t")
    assert any(
        row.startswith("likelihood_ratio_test\t\tbrownian-vs-early-burst\t")
        for row in comparison_rows[1:]
    )
    assert profile_rows[0].startswith("trait\trate_change\tlog_likelihood")
    assert len(profile_rows) == 162


def test_summarize_early_burst_trait_evolution_prefers_early_burst_on_simulated_data(
    tmp_path: Path,
) -> None:
    tree_path, traits_path = _write_simulated_early_burst_dataset(tmp_path)

    report = summarize_early_burst_trait_evolution(
        tree_path,
        traits_path,
        trait="trait",
        rate_change_bounds=(0.0, 10.0),
    )

    assert report.better_model == "early-burst"
    assert report.rate_change > 0.0
    assert (
        3.0
        <= report.confidence_intervals[0].lower_95
        <= report.confidence_intervals[0].upper_95
    )
    assert report.identifiability_warnings == []
    assert next(
        row for row in report.comparison_rows if row.model == "early-burst"
    ).selected


def _write_simulated_early_burst_dataset(tmp_path: Path) -> tuple[Path, Path]:
    tree_path = tmp_path / "simulated-early-burst-tree.nwk"
    traits_path = tmp_path / "simulated-early-burst-traits.tsv"
    tree_path.write_text(
        (
            "((((A:0.8,B:0.4):0.7,(C:0.6,D:0.3):0.5):1.2,"
            "((E:0.7,F:0.35):0.8,(G:0.5,H:0.25):0.9):1.0):0.6,"
            "((I:0.9,J:0.45):0.6,(K:0.55,L:0.2):1.1):0.85);"
        ),
        encoding="utf-8",
    )
    base_tree = load_tree(tree_path)
    transformed_tree = transform_tree_for_evolutionary_mode(
        base_tree,
        mode="early-burst",
        parameter_value=6.0,
    )
    covariance = stable_covariance(
        build_brownian_covariance_matrix(
            transformed_tree,
            base_tree.tip_names,
        )
    )
    cholesky = _cholesky(covariance)
    rng = random.Random(1)
    standard_normal = [rng.gauss(0.0, 1.0) for _ in base_tree.tip_names]
    values = [
        3.0 + sum(cholesky[row][col] * standard_normal[col] for col in range(row + 1))
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
