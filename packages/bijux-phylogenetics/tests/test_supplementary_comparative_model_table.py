from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.reports import write_supplementary_comparative_model_table

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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_write_supplementary_comparative_model_table_writes_pgls_coefficient_rows(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "supplementary-comparative-model.tsv"

    result = write_supplementary_comparative_model_table(
        output_path,
        tree_path=fixture("example_tree_six_taxa.nwk"),
        traits_path=fixture("example_traits_comparative_multiple.tsv"),
        formulas=[
            "response_growth ~ predictor_one",
            "response_growth ~ predictor_two",
            "response_growth ~ predictor_one + predictor_two",
        ],
        lambda_value=0.0,
    )

    assert result.output_path == output_path
    assert result.model_count == 3
    assert result.selected_formula == "response_growth ~ predictor_one"
    assert result.selected_criterion == "AICc"
    assert result.excluded_taxon_count == 0
    assert result.row_count == 7

    selected_rows = [row for row in result.rows if row.selected]
    assert {row.formula for row in selected_rows} == {"response_growth ~ predictor_one"}
    assert {row.coefficient_name for row in selected_rows} == {
        "intercept",
        "predictor_one",
    }
    predictor_row = next(
        row
        for row in result.rows
        if row.formula == "response_growth ~ predictor_one"
        and row.coefficient_name == "predictor_one"
    )
    assert predictor_row.model_family == "pgls"
    assert predictor_row.phylogenetic_parameter_name == "lambda"
    assert predictor_row.phylogenetic_parameter_value == 0.0
    assert predictor_row.outlier_taxon_count == 0
    assert predictor_row.max_leverage is not None
    assert predictor_row.converged is None

    written = read_tsv(output_path)
    assert written[0]["response"] == "response_growth"
    assert written[0]["model_family"] == "pgls"


def test_write_supplementary_comparative_model_table_carries_shared_exclusions(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "supplementary-comparative-model.tsv"

    result = write_supplementary_comparative_model_table(
        output_path,
        tree_path=fixture("example_tree_six_taxa.nwk"),
        traits_path=fixture("example_traits_comparative_missing_predictor.tsv"),
        formulas=[
            "response ~ predictor_one",
            "response ~ predictor_two",
        ],
        lambda_value=0.0,
    )

    assert result.excluded_taxon_count == 3
    assert all(row.excluded_taxon_count == 3 for row in result.rows)
    assert any(
        "B:missing_required_values:predictor_one" in row.excluded_taxa
        for row in result.rows
    )
    assert any(
        "E:missing_from_trait_table:predictor_one,predictor_two,response"
        in row.excluded_taxa
        for row in result.rows
    )

    written = read_tsv(output_path)
    assert any(
        "B:missing_required_values:predictor_one" in row["excluded_taxa"]
        for row in written
    )
