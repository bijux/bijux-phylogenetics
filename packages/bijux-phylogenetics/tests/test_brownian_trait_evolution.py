from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.comparative.continuous import (
    summarize_brownian_trait_evolution,
    write_brownian_trait_evolution_exclusion_table,
    write_brownian_trait_evolution_summary_table,
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


def test_summarize_brownian_trait_evolution_matches_reference_example() -> None:
    report = summarize_brownian_trait_evolution(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )

    reference = _brownian_reference_observation()
    assert report.tree_taxon_count == 4
    assert report.analyzed_taxon_count == 4
    assert report.excluded_taxa == []
    assert math.isclose(
        report.root_state,
        reference["root_state"],
        rel_tol=1e-6,
        abs_tol=1e-6,
    )
    assert math.isclose(
        report.sigma_squared,
        reference["rate"],
        rel_tol=1e-6,
        abs_tol=1e-6,
    )
    assert math.isclose(
        report.log_likelihood,
        reference["log_likelihood"],
        rel_tol=1e-6,
        abs_tol=1e-6,
    )
    expected_aic = (2.0 * 2.0) - (2.0 * report.log_likelihood)
    expected_aicc = expected_aic + ((2.0 * 2.0 * 3.0) / (4 - 2 - 1))
    assert math.isclose(report.aic, expected_aic, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(report.aicc, expected_aicc, rel_tol=0.0, abs_tol=1e-12)


def test_summarize_brownian_trait_evolution_records_missing_and_invalid_taxa() -> None:
    report = summarize_brownian_trait_evolution(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_brownian_missing.tsv"),
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


def test_brownian_trait_evolution_writers_emit_summary_and_exclusion_ledgers(
    tmp_path: Path,
) -> None:
    report = summarize_brownian_trait_evolution(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_brownian_missing.tsv"),
        trait="response_growth",
    )
    summary_out = tmp_path / "brownian-trait-summary.tsv"
    excluded_out = tmp_path / "brownian-trait-excluded.tsv"

    write_brownian_trait_evolution_summary_table(summary_out, report)
    write_brownian_trait_evolution_exclusion_table(excluded_out, report)

    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("trait\ttaxon_column\ttree_taxon_count")
    assert summary_rows[1].split("\t")[3] == "4"
    assert excluded_rows == [
        "taxon\treason",
        "B\tmissing_trait_value",
        "C\tnon_numeric_trait_value",
        "G\tabsent_from_tree",
    ]


def _brownian_reference_observation() -> dict[str, float]:
    payload = json.loads(
        fixture("comparative_reference_validation.json").read_text(encoding="utf-8")
    )
    observation = next(
        row for row in payload["observations"] if row["case"] == "brownian-example-tree"
    )
    return {
        key: float(value) for key, value in observation["expected_parameters"].items()
    }
