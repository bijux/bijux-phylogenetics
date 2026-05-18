from __future__ import annotations

import csv
from pathlib import Path

import pytest

from bijux_phylogenetics.reports import write_supplementary_model_selection_table

FIXTURES = Path(__file__).parent / "fixtures"


def _write_iqtree_model_fixture(report_path: Path, sidecar_path: Path) -> None:
    report_path.write_text(
        "\n".join(
            [
                " No. Model         -LnL         df  AIC          AICc         BIC",
                "  1  GTR+F         20.415       13  66.831       430.831      67.863",
                "  2  GTR+F+I       19.209       14  66.417       486.417      67.530",
                "  3  JC            23.429       5   56.857       86.857       57.255",
                "  4  K2P+I         21.104       7   56.208       168.208      56.764",
                "Akaike Information Criterion:           K2P+I",
                "Corrected Akaike Information Criterion: JC",
                "Bayesian Information Criterion:         K2P+I",
                "Best-fit model according to BIC: K2P+I",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    sidecar_path.write_text(
        "\n".join(
            [
                "--- # IQ-TREE Checkpoint ver >= 1.6",
                "best_model_AIC: K2P+I",
                "best_model_AICc: JC",
                "best_model_BIC: K2P+I",
                "best_score_AIC: 56.20822137",
                "best_score_AICc: 86.85740671",
                "best_score_BIC: 56.76431216",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_write_supplementary_model_selection_table_writes_candidate_rows(
    tmp_path: Path,
) -> None:
    report_path = tmp_path / "example.iqtree"
    sidecar_path = tmp_path / "example.model"
    output_path = tmp_path / "supplementary-model-selection.tsv"
    _write_iqtree_model_fixture(report_path, sidecar_path)

    result = write_supplementary_model_selection_table(
        output_path,
        iqtree_report_path=report_path,
        model_sidecar_path=sidecar_path,
    )

    assert result.output_path == output_path
    assert result.row_count == 4
    assert result.candidate_count == 4
    assert result.selected_model == "K2P+I"
    assert result.selected_criterion == "BIC"
    assert result.rows[-1].selected_model is True
    assert result.rows[-1].best_aic is True
    assert result.rows[-1].best_bic is True
    assert result.rows[2].best_aicc is True

    written = _read_tsv(output_path)
    assert written[0]["model"] == "GTR+F"
    assert written[-1]["selected_model"] == "True"
    assert written[-1]["selected_model_name"] == "K2P+I"
    assert written[2]["best_aicc"] == "True"


def test_write_supplementary_model_selection_table_auto_resolves_sidecar(
    tmp_path: Path,
) -> None:
    report_path = tmp_path / "example.iqtree"
    sidecar_path = tmp_path / "example.model"
    output_path = tmp_path / "supplementary-model-selection.tsv"
    _write_iqtree_model_fixture(report_path, sidecar_path)

    result = write_supplementary_model_selection_table(
        output_path,
        iqtree_report_path=report_path,
    )

    assert result.row_count == 4
    assert result.rows[0].model_sidecar_source == str(sidecar_path)


def test_write_supplementary_model_selection_table_requires_candidate_rows(
    tmp_path: Path,
) -> None:
    report_path = tmp_path / "version-variant.iqtree"
    output_path = tmp_path / "supplementary-model-selection.tsv"
    fixture_root = FIXTURES / "engine_outputs" / "iqtree"
    report_path.write_text(
        (fixture_root / "model-selection-version-variant.iqtree").read_text(
            encoding="utf-8"
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="do not expose candidate model rows",
    ):
        write_supplementary_model_selection_table(
            output_path,
            iqtree_report_path=report_path,
            model_sidecar_path=fixture_root / "model-selection-version-variant.model",
        )
