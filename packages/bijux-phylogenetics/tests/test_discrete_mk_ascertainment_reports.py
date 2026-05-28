from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.discrete_mk import (
    fit_discrete_mk_model,
    write_discrete_mk_summary_table,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_discrete_mk_summary_table_reports_ascertainment_terms(
    tmp_path: Path,
) -> None:
    report = fit_discrete_mk_model(
        fixture("trees", "example_tree.nwk"),
        fixture("metadata", "example_traits_discrete_mk_variable_only_four_taxa.tsv"),
        trait="state",
        taxon_column="taxon",
        model="equal-rates",
        ascertainment_policy="lewis-variable-only",
    )

    output_path = tmp_path / "discrete-mk-summary.tsv"
    write_discrete_mk_summary_table(output_path, report)
    lines = output_path.read_text(encoding="utf-8").strip().splitlines()

    assert "ascertainment_policy" in lines[0]
    assert "ascertainment_conditioning_log_probability" in lines[0]
    assert "invariant_pattern_log_probability" in lines[0]
    assert "lewis-variable-only" in lines[1]
