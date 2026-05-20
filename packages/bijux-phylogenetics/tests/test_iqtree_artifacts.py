from __future__ import annotations

import gzip
from pathlib import Path

from bijux_phylogenetics.engines.artifacts.iqtree import (
    parse_iqtree_model_selection_summary,
    parse_log_likelihood_text,
    parse_selected_model_decision_text,
    resolve_iqtree_model_sidecar,
    write_iqtree_model_candidates_table,
)


def test_parse_selected_model_decision_text_supports_iqtree_variants() -> None:
    assert parse_selected_model_decision_text(
        "Best-fit model according to BIC: GTR+G\n"
    ) == ("GTR+G", "BIC")
    assert parse_selected_model_decision_text(
        "Best-fit model: K2P+I chosen according to BIC\n"
    ) == ("K2P+I", "BIC")
    assert parse_selected_model_decision_text(
        "Best-fit model according to BIC score = K2P+I\n"
    ) == ("K2P+I", "BIC")


def test_parse_log_likelihood_text_supports_report_and_log_variants() -> None:
    assert (
        parse_log_likelihood_text("Log-likelihood of the tree: -123.456\n") == -123.456
    )
    assert parse_log_likelihood_text("BEST SCORE FOUND : -21.104\n") == -21.104
    assert parse_log_likelihood_text("Log likelihood of tree = -19.209\n") == -19.209


def test_parse_iqtree_model_selection_summary_reads_candidates_and_criteria(
    tmp_path: Path,
) -> None:
    iqtree_report = tmp_path / "example.iqtree"
    iqtree_report.write_text(
        "\n".join(
            [
                "ModelFinder will test up to 4 DNA models ...",
                " No. Model         -LnL         df  AIC          AICc         BIC",
                "  1  GTR+F         20.415       13  66.831       430.831      67.863",
                "  2  GTR+F+I       19.209       14  66.417       486.417      67.530",
                "  3  JC            23.429       5   56.857       86.857       57.255",
                "  4  K2P+I         21.104       7   56.208       168.208      56.764",
                "Akaike Information Criterion:           K2P+I",
                "Corrected Akaike Information Criterion: JC",
                "Bayesian Information Criterion:         K2P+I",
                "Best-fit model: K2P+I chosen according to BIC",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    model_sidecar = tmp_path / "example.model.gz"
    model_sidecar.write_bytes(
        gzip.compress(
            "\n".join(
                [
                    "--- # IQ-TREE Checkpoint ver >= 1.6",
                    "best_model_AIC: K2P+I",
                    "best_model_AICc: JC",
                    "best_model_BIC: K2P+I",
                    "best_model_list_BIC: K2P+I JC GTR+F+I",
                    "best_score_AIC: 56.20822137",
                    "best_score_AICc: 86.85740671",
                    "best_score_BIC: 56.76431216",
                ]
            ).encode("utf-8")
        )
    )

    summary = parse_iqtree_model_selection_summary(
        iqtree_report_path=iqtree_report,
        model_sidecar_path=model_sidecar,
    )

    assert summary is not None
    assert summary.selected_model == "K2P+I"
    assert summary.selected_criterion == "BIC"
    assert summary.best_model_aic == "K2P+I"
    assert summary.best_model_aicc == "JC"
    assert summary.best_model_bic == "K2P+I"
    assert summary.best_score_aic == 56.20822137
    assert summary.best_score_aicc == 86.85740671
    assert summary.best_score_bic == 56.76431216
    assert summary.candidate_count == 4
    assert summary.bic_near_best_models == ["K2P+I", "JC", "GTR+F+I"]
    assert summary.candidates[0].model == "GTR+F"
    assert summary.candidates[0].log_likelihood == -20.415
    assert summary.candidates[-1].model == "K2P+I"
    assert summary.candidates[-1].bic == 56.764


def test_write_iqtree_model_candidates_table_marks_selected_and_criteria(
    tmp_path: Path,
) -> None:
    iqtree_report = tmp_path / "example.iqtree"
    iqtree_report.write_text(
        "\n".join(
            [
                " No. Model         -LnL         df  AIC          AICc         BIC",
                "  1  JC            23.429       5   56.857       86.857       57.255",
                "  2  K2P+I         21.104       7   56.208       168.208      56.764",
                "Akaike Information Criterion:           K2P+I",
                "Corrected Akaike Information Criterion: JC",
                "Bayesian Information Criterion:         K2P+I",
                "Best-fit model according to BIC: K2P+I",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    summary = parse_iqtree_model_selection_summary(iqtree_report_path=iqtree_report)
    assert summary is not None

    table_path = write_iqtree_model_candidates_table(
        tmp_path / "model-candidates.tsv",
        summary,
    )
    lines = table_path.read_text(encoding="utf-8").splitlines()

    assert lines[0].startswith("rank\tmodel\tlog_likelihood")
    assert lines[1].endswith("\tfalse\ttrue\tfalse\tfalse")
    assert lines[2].endswith("\ttrue\tfalse\ttrue\ttrue")


def test_resolve_iqtree_model_sidecar_prefers_plain_text_before_gzip(
    tmp_path: Path,
) -> None:
    prefix = tmp_path / "example"
    plain = prefix.with_suffix(".model")
    plain.write_text("Best-fit model: GTR+G\n", encoding="utf-8")
    gzip_path = prefix.with_suffix(".model.gz")
    gzip_path.write_bytes(gzip.compress(b"best_model_BIC: GTR+G\n"))

    assert resolve_iqtree_model_sidecar(prefix) == plain


def test_parse_iqtree_model_selection_summary_accepts_version_variant_fixture() -> None:
    fixture_root = Path(__file__).parent / "fixtures" / "engine_outputs" / "iqtree"
    summary = parse_iqtree_model_selection_summary(
        iqtree_report_path=fixture_root / "model-selection-version-variant.iqtree",
        model_sidecar_path=fixture_root / "model-selection-version-variant.model",
    )

    assert summary is not None
    assert summary.selected_model == "K2P+I"
    assert summary.selected_criterion == "BIC"
    assert summary.best_model_aic == "K2P+I"
    assert summary.best_model_aicc == "JC"
    assert summary.best_model_bic == "K2P+I"
    assert summary.best_score_aic == 56.20822137
    assert summary.best_score_aicc == 86.85740671
    assert summary.best_score_bic == 56.76431216
    assert summary.bic_near_best_models == ["K2P+I", "JC", "GTR+F+I"]
    assert summary.candidate_count == 0
