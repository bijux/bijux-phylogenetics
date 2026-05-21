from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.reports.publication.alignment import (
    build_alignment_figure_package,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "alignments" / name


def test_build_alignment_figure_package_writes_publication_bundle(
    tmp_path: Path,
) -> None:
    result = build_alignment_figure_package(
        fixture("example_alignment.fasta"),
        out_dir=tmp_path / "alignment-figure-package",
    )

    assert result.heatmap_figure_path.exists()
    assert result.site_summary_figure_path.exists()
    assert result.sequence_panel_figure_path.exists()
    assert result.heatmap_table_path.exists()
    assert result.window_table_path.exists()
    assert result.ranking_table_path.exists()
    assert result.legend_path.exists()
    assert result.caption_path.exists()
    assert result.review_path.exists()
    assert result.manifest_path.exists()
    assert result.reproducibility_manifest_path.exists()
    assert result.reviewer_audit_checklist_path.exists()
    assert result.audit.publication_ready is True
    assert result.audit.heatmap_visible is True
    assert result.audit.site_summary_visible is True
    assert result.audit.sequence_panel_visible is True
    assert result.audit.heatmap_row_count == 4
    assert result.audit.heatmap_bin_count == 8
    assert result.audit.plotted_window_count == 1
    assert result.audit.plotted_sequence_count == 4

    html = result.review_path.read_text(encoding="utf-8")
    assert "Bijux Alignment Quality Review" in html
    assert "Missingness Heatmap" in html
    assert "Site-Quality Summary" in html
    assert "Sequence-Quality Panel" in html
    assert "Reviewer Audit Checklist" in html

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    reproducibility = json.loads(
        result.reproducibility_manifest_path.read_text(encoding="utf-8")
    )
    checklist_rows = result.reviewer_audit_checklist_path.read_text(
        encoding="utf-8"
    ).splitlines()
    assert manifest["report_kind"] == "alignment_quality_figure_package"
    assert manifest["metrics"]["publication_ready"] is True
    assert manifest["reproducibility_manifest_path"] == str(
        result.reproducibility_manifest_path
    )
    assert manifest["reviewer_audit_checklist_path"].endswith(
        "reviewer-audit-checklist.tsv"
    )
    assert len(manifest["reviewer_audit_checklist"]["items"]) == 5
    assert checklist_rows[0] == "section\tstatus\tsummary\tevidence\tartifact_paths"
    assert any(
        line.startswith("publication_readiness\t") for line in checklist_rows[1:]
    )
    assert reproducibility["report_kind"] == "alignment_quality_figure_package"
    assert reproducibility["settings"]["maximum_site_bins"] == 120


def test_build_alignment_figure_package_blocks_suspicious_alignment_publication(
    tmp_path: Path,
) -> None:
    result = build_alignment_figure_package(
        fixture("example_alignment_missingness.fasta"),
        out_dir=tmp_path / "alignment-figure-package",
    )

    assert result.audit.publication_ready is False
    assert result.audit.suspicious_alignment is True
    assert result.audit.heatmap_visible is True
    assert result.audit.site_summary_visible is True
    assert result.audit.sequence_panel_visible is True
    assert result.audit.heatmap_row_count == 3
    assert result.audit.heatmap_bin_count == 6
    assert result.audit.plotted_window_count == 1
    assert result.audit.plotted_sequence_count == 3
    assert any("missing" in item for item in result.audit.limitations)
    assert len({cell.uncertainty_fraction for cell in result.heatmap_cells}) > 1
