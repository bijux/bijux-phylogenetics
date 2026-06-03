from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.reports import (
    write_publication_package_revalidation_report,
)

from .support.rabies_cross_host_geography_package import (
    build_stub_rabies_cross_host_geography_package,
)


def test_publication_package_revalidation_passes_on_preserved_rabies_package(
    tmp_path: Path,
    monkeypatch,
) -> None:
    package_result = build_stub_rabies_cross_host_geography_package(
        tmp_path / "package",
        monkeypatch,
    )

    result = write_publication_package_revalidation_report(
        tmp_path / "revalidation",
        package_result.package_manifest_path,
    )

    assert result.report_kind == "rabies_cross_host_geography_package"
    assert result.all_original_artifacts_match is True
    assert result.overall_revalidation_status == "pass"
    assert result.unexpected_file_count == 0
    assert result.missing_artifact_count == 0
    assert result.checksum_mismatch_count == 0
    assert result.size_mismatch_count == 0
    assert result.artifact_table_path.is_file()
    assert result.check_table_path.is_file()
    assert result.summary_path.is_file()
    assert result.report_path.is_file()
    artifact_rows = list(
        csv.DictReader(
            result.artifact_table_path.open("r", encoding="utf-8", newline=""),
            delimiter="\t",
        )
    )
    assert any(
        row["relative_path"] == "workflow/rabies-cross-host-geography-report.html"
        and row["status"] == "pass"
        for row in artifact_rows
    )
    assert any(
        row["relative_path"] == "rabies-cross-host-geography-artifacts.tsv"
        and row["artifact_scope"] == "package_control"
        for row in artifact_rows
    )


def test_publication_package_revalidation_blocks_on_checksum_drift(
    tmp_path: Path,
    monkeypatch,
) -> None:
    package_result = build_stub_rabies_cross_host_geography_package(
        tmp_path / "package",
        monkeypatch,
    )
    drifted_path = package_result.dataset_export.sequences_path
    drifted_path.write_text(
        drifted_path.read_text(encoding="utf-8") + "\n>drift\nACGT\n",
        encoding="utf-8",
    )

    result = write_publication_package_revalidation_report(
        tmp_path / "revalidation",
        package_result.package_manifest_path,
    )

    assert result.all_original_artifacts_match is False
    assert result.overall_revalidation_status == "blocked"
    assert result.checksum_mismatch_count > 0
    assert any(
        row.relative_path == "dataset/sequences.fasta" and row.status == "blocked"
        for row in result.artifact_rows
    )
    assert any(
        row.check_id == "inventory-listed-artifacts-match" and row.status == "blocked"
        for row in result.check_rows
    )


def test_publication_package_revalidation_marks_unexpected_files_as_risk(
    tmp_path: Path,
    monkeypatch,
) -> None:
    package_result = build_stub_rabies_cross_host_geography_package(
        tmp_path / "package",
        monkeypatch,
    )
    (package_result.output_root / "reviewer-notes.txt").write_text(
        "manual note\n",
        encoding="utf-8",
    )

    result = write_publication_package_revalidation_report(
        tmp_path / "revalidation",
        package_result.package_manifest_path,
    )

    assert result.all_original_artifacts_match is True
    assert result.overall_revalidation_status == "risk"
    assert result.unexpected_file_count == 1
    assert any(
        row.check_id == "unexpected-package-files" and row.status == "risk"
        for row in result.check_rows
    )
