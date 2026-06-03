from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.pgls import inspect_pgls_inputs, run_pgls

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


def test_pgls_prunes_missing_numeric_cells_instead_of_rejecting_the_column() -> None:
    input_report = inspect_pgls_inputs(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogenetic_residuals_missing.tsv"),
        response="brain_mass",
        predictors=["body_mass"],
    )
    model = run_pgls(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogenetic_residuals_missing.tsv"),
        response="brain_mass",
        predictors=["body_mass"],
        lambda_value=1.0,
    )

    assert input_report.ready is True
    assert input_report.analysis_taxa == ["A", "B", "C", "D", "F"]
    assert [
        (row.taxon, row.reason) for row in input_report.formula_audit.excluded_taxa
    ] == [
        ("E", "missing_value"),
    ]
    assert model.taxa == ["A", "B", "C", "D", "F"]
    assert [row.taxon for row in model.diagnostics.fitted_observed_rows] == model.taxa
