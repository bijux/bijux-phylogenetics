from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.reports import (
    write_publication_package_comparison_report,
)

from .support.rabies_cross_host_geography_package import (
    build_stub_rabies_cross_host_geography_package,
    refresh_stub_rabies_cross_host_geography_package,
)


def test_publication_package_comparison_passes_on_identical_rabies_packages(
    tmp_path: Path,
    monkeypatch,
) -> None:
    left = build_stub_rabies_cross_host_geography_package(
        tmp_path / "left", monkeypatch
    )
    right = build_stub_rabies_cross_host_geography_package(
        tmp_path / "right",
        monkeypatch,
    )

    result = write_publication_package_comparison_report(
        tmp_path / "comparison",
        left.package_manifest_path,
        right.package_manifest_path,
    )

    assert result.report_kind == "rabies_cross_host_geography_package"
    assert result.dataset_id == "rabies_cross_host_geography_panel"
    assert result.overall_comparison_status == "pass"
    assert result.changed_artifact_count == 0
    assert result.left_only_artifact_count == 0
    assert result.right_only_artifact_count == 0
    assert result.config_difference_count == 0
    assert result.alignment_difference_count == 0
    assert result.figure_or_report_difference_count == 0
    assert result.scientific_finding_difference_count == 0
    assert result.artifact_table_path.is_file()
    assert result.check_table_path.is_file()
    assert result.summary_path.is_file()
    assert result.report_path.is_file()


def test_publication_package_comparison_flags_input_tree_model_figure_and_conclusion_drift(
    tmp_path: Path,
    monkeypatch,
) -> None:
    left = build_stub_rabies_cross_host_geography_package(
        tmp_path / "left", monkeypatch
    )
    right = build_stub_rabies_cross_host_geography_package(
        tmp_path / "right",
        monkeypatch,
    )
    with right.dataset_export.sequences_path.open("a", encoding="utf-8") as handle:
        handle.write(">wolf_extra_rv999\nACGT\n")
    with right.dataset_export.accession_table_path.open(
        "a", encoding="utf-8"
    ) as handle:
        handle.write("RV999\thttps://example.org/RV999\n")
    right.workflow_bundle.selected_model = "GTR+F+I"
    right.workflow_bundle.root_host = "canid"
    right.workflow_bundle.root_region = "south_america"
    right.workflow_bundle.tree_path.write_text(
        "(bat_chile_rv108:0.1,fox_canada_rv241:0.2,wolf_extra_rv999:0.3)root;\n",
        encoding="utf-8",
    )
    (
        right.workflow_bundle.output_root / "rabies-cross-host-geography-panel.aln"
    ).write_text(
        ">taxon_a\nACGT\n>taxon_b\nACGT\n>taxon_c\nACGT\n",
        encoding="utf-8",
    )
    (
        right.workflow_bundle.output_root / "biogeography" / "biogeography-map.svg"
    ).write_text("<svg><text>updated-map</text></svg>\n", encoding="utf-8")
    right.workflow_bundle.scientific_findings_path.write_text(
        "finding_id\tquestion\tclaim\tevidence\tcaution\tsource_artifact\n"
        "comparative_longitude\tquestion\tchanged claim\tevidence\tcaution\tscientific-findings.tsv\n",
        encoding="utf-8",
    )
    right = refresh_stub_rabies_cross_host_geography_package(right)

    result = write_publication_package_comparison_report(
        tmp_path / "comparison",
        left.package_manifest_path,
        right.package_manifest_path,
    )

    assert result.overall_comparison_status == "risk"
    assert result.changed_artifact_count > 0
    assert result.sequence_right_only_count == 1
    assert result.accession_right_only_count == 1
    assert result.config_difference_count == 0
    assert result.alignment_difference_count > 0
    assert result.figure_or_report_difference_count > 0
    assert result.scientific_finding_difference_count == 1
    assert any(
        row.check_id == "taxa-and-accessions" and row.status == "risk"
        for row in result.check_rows
    )
    assert any(
        row.check_id == "rooted-tree" and row.status == "risk"
        for row in result.check_rows
    )
    assert any(
        row.relative_path == "workflow/biogeography/biogeography-map.svg"
        and row.status == "changed"
        for row in result.artifact_rows
    )
    artifact_rows = list(
        csv.DictReader(
            result.artifact_table_path.open("r", encoding="utf-8", newline=""),
            delimiter="\t",
        )
    )
    assert any(
        row["relative_path"] == "dataset/sequences.fasta" and row["status"] == "changed"
        for row in artifact_rows
    )
