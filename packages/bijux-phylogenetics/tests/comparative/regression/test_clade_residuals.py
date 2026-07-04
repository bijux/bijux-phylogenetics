from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.clades.residuals import (
    analyze_comparative_residual_clades,
    write_comparative_residual_clade_table,
    write_comparative_residual_taxon_table,
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


def test_analyze_comparative_residual_clades_identifies_residual_heavy_clade() -> None:
    report = analyze_comparative_residual_clades(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        formula="response_growth ~ predictor_two",
        lambda_value=0.0,
    )

    assert report.model_family == "pgls"
    assert report.standardized_residual_method == "leveraged-gls-residual"
    assert len(report.taxon_rows) == 6
    assert len(report.clade_rows) == 4
    assert report.residual_heavy_clades == ["E|F"]

    ranked = sorted(report.clade_rows, key=lambda row: row.rank)
    assert ranked[0].clade_id == "E|F"
    assert ranked[0].residual_heavy is True
    assert ranked[0].influence_score > ranked[1].influence_score


def test_analyze_comparative_residual_clades_auto_detects_logistic_family() -> None:
    report = analyze_comparative_residual_clades(
        fixture("example_tree_eight_taxa.nwk"),
        fixture("example_traits_phylogenetic_logistic_model_selection.tsv"),
        formula="presence ~ habitat",
        lambda_value=1.0,
    )

    assert report.model_family == "logistic"
    assert report.standardized_residual_method == "pearson-binomial-residual"
    assert len(report.taxon_rows) == 8
    assert len(report.clade_rows) == 6


def test_write_comparative_residual_tables_write_expected_rows(tmp_path: Path) -> None:
    report = analyze_comparative_residual_clades(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        formula="response_growth ~ predictor_two",
        lambda_value=0.0,
    )
    taxon_out = tmp_path / "comparative-residual-taxa.tsv"
    clade_out = tmp_path / "comparative-residual-clades.tsv"

    write_comparative_residual_taxon_table(taxon_out, report)
    write_comparative_residual_clade_table(clade_out, report)

    taxon_rows = taxon_out.read_text(encoding="utf-8").splitlines()
    clade_rows = clade_out.read_text(encoding="utf-8").splitlines()
    assert taxon_rows[0].startswith("taxon\tobserved_value\tfitted_value\tresidual")
    assert clade_rows[0].startswith("clade_id\tnode_label\ttaxon_count\ttaxa")
    assert len(taxon_rows) == 7
    assert len(clade_rows) == 5
