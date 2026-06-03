from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.io.artifact_schema import (
    artifact_schema_names,
    assert_artifact_schema_valid,
    validate_artifact_schema,
    validate_clade_table_schema,
    validate_run_manifest_schema,
)


def test_artifact_schema_names_list_supported_profiles() -> None:
    assert artifact_schema_names() == (
        "clade_table_tsv",
        "comparative_report_manifest_json",
        "comparative_summary_tsv",
        "comparative_traits_tsv",
        "fasta_to_tree_manifest_json",
        "fasta_to_tree_model_tsv",
        "fasta_to_tree_support_tsv",
        "geographic_event_tsv",
        "host_switch_branch_tsv",
        "run_manifest_json",
    )


def test_validate_artifact_schema_rejects_unknown_profile(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.tsv"
    artifact.write_text("column\nvalue\n", encoding="utf-8")

    with pytest.raises(ValueError, match="unknown artifact schema"):
        validate_artifact_schema(artifact, "unknown_schema")


def test_run_manifest_schema_reports_missing_required_key(tmp_path: Path) -> None:
    manifest_path = tmp_path / "run.manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "arguments": [],
                "command": "demo",
                "dependency_versions": {},
                "host_platform": "test",
                "input_checksums": {},
                "input_paths": [],
                "output_checksums": {},
                "output_paths": [],
                "package_version": "0.0.0",
                "python_version": "3.11.0",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    report = validate_run_manifest_schema(manifest_path)

    assert not report.valid
    assert report.missing_fields == ("timestamp_utc",)
    with pytest.raises(ValueError, match="timestamp_utc"):
        assert_artifact_schema_valid(report)


def test_clade_schema_rejects_invalid_metadata_triplets(tmp_path: Path) -> None:
    clade_path = tmp_path / "clades.tsv"
    clade_path.write_text(
        "\t".join(
            [
                "source_path",
                "tree_index",
                "node_kind",
                "clade_id",
                "node_label",
                "taxon_count",
                "taxa",
                "support",
                "support_fraction",
                "branch_length",
                "root_depth",
                "descendant_tip_depth_min",
                "descendant_tip_depth_max",
                "node_age",
                "host_group_values",
                "host_group_missing_taxa",
                "host_group_distinct_values",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "tree.nwk",
                "",
                "root",
                "A|B",
                "",
                "2",
                "A|B",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "A=bat|B=dog",
                "",
                "bat|dog",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = validate_clade_table_schema(clade_path)

    assert not report.valid
    assert report.unexpected_fields == (
        "host_group_values",
        "host_group_missing_taxa",
        "host_group_distinct_values",
    )
