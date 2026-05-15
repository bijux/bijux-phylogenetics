from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.io.artifact_schema import (
    validate_comparative_report_manifest_schema,
    validate_comparative_summary_table_schema,
    validate_fasta_to_tree_model_table_schema,
)


def test_model_schema_rejects_column_order_drift(tmp_path: Path) -> None:
    path = tmp_path / "example.model.tsv"
    path.write_text(
        "\t".join(
            [
                "engine_name",
                "workflow",
                "sequence_type",
                "selected_model",
                "report_selected_model",
                "artifact_selected_model",
                "model_consistent",
                "alignment_path",
                "trimmed_alignment_path",
                "manifest_path",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "iqtree",
                "model-selection",
                "dna",
                "GTR+G",
                "GTR+G",
                "GTR+G",
                "true",
                "example.aln",
                "example.trimmed.aln",
                "example.manifest.json",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = validate_fasta_to_tree_model_table_schema(path)

    assert not report.valid
    assert not report.order_matches
    assert report.missing_fields == ()
    assert report.unexpected_fields == ()


def test_comparative_summary_schema_rejects_missing_field(tmp_path: Path) -> None:
    path = tmp_path / "comparative-summary.tsv"
    path.write_text(
        "\t".join(
            [
                "response",
                "formula",
                "predictor_count",
                "analysis_taxa",
                "excluded_taxa",
                "pgls_lambda",
                "pgls_log_likelihood",
                "pgls_r_squared",
                "phylogenetic_signal_k",
                "phylogenetic_signal_lambda",
                "independent_contrast_count",
                "better_model_aicc_delta",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "region_longitude",
                "region_longitude ~ host_group",
                "1",
                "9",
                "0",
                "1.0",
                "-12.3",
                "0.8",
                "2.4",
                "1.0",
                "8",
                "4.9",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = validate_comparative_summary_table_schema(path)

    assert not report.valid
    assert report.missing_fields == ("selected_model",)


def test_comparative_manifest_schema_rejects_unexpected_key(tmp_path: Path) -> None:
    path = tmp_path / "comparative-report.manifest.json"
    path.write_text(
        json.dumps(
            {
                "input_checksums": {},
                "input_paths": [],
                "limitations": [],
                "metrics": {},
                "outputs": {},
                "report_kind": "comparative_package",
                "summary": {},
                "transitional_key": True,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    report = validate_comparative_report_manifest_schema(path)

    assert not report.valid
    assert report.unexpected_fields == ("transitional_key",)
