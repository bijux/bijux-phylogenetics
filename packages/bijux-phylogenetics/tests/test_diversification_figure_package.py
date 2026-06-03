from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.comparative import build_diversification_figure_package

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


def test_build_diversification_figure_package_writes_publication_bundle(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "diversification-figure-package"
    result = build_diversification_figure_package(
        fixture("example_tree.nwk"),
        metadata_path=fixture("example_sampling_fractions.tsv"),
        out_dir=output_dir,
    )

    assert result.audit.publication_ready is True
    assert result.audit.lineage_curve_visible is True
    assert result.audit.clade_outlier_surface_visible is True
    assert result.audit.model_comparison_visible is True
    assert result.audit.sampling_metadata_complete is True
    assert result.audit.plotted_ltt_point_count == 4
    assert result.audit.plotted_clade_count == 3
    assert result.audit.highlighted_outlier_count == 2
    assert result.audit.plotted_model_count == 2
    assert result.audit.better_model in {"yule", "birth-death"}
    assert result.lineage_figure_path.exists()
    assert result.clade_figure_path.exists()
    assert result.model_figure_path.exists()
    assert result.legend_path.exists()
    assert result.caption_path.exists()
    assert result.methods_summary_path.exists()
    assert result.review_path.exists()
    assert result.manifest_path.exists()
    assert result.reproducibility_manifest_path.exists()
    assert "Lineage-through-time curve" in result.lineage_figure_path.read_text(
        encoding="utf-8"
    )
    assert "Clade diversification outliers" in result.clade_figure_path.read_text(
        encoding="utf-8"
    )
    assert "Diversification model comparison" in result.model_figure_path.read_text(
        encoding="utf-8"
    )
    assert (
        "Diversification Analysis Methods Summary"
        in result.methods_summary_path.read_text(encoding="utf-8")
    )
    assert "publication_ready" in result.review_path.read_text(encoding="utf-8")
    assert "Methods Summary" in result.review_path.read_text(encoding="utf-8")
    assert (
        result.machine_manifest["metrics"]["publication_ready"]
        == result.audit.publication_ready
    )
    assert result.machine_manifest["outputs"]["methods_summary_path"].endswith(
        "diversification-methods-summary.md"
    )
    assert result.machine_manifest["metrics"]["methods_summary_warning_count"] >= 0
    reproducibility = json.loads(
        result.reproducibility_manifest_path.read_text(encoding="utf-8")
    )
    assert result.machine_manifest["reproducibility_manifest_path"] == str(
        result.reproducibility_manifest_path
    )
    assert reproducibility["report_kind"] == "diversification_figure_package"
    assert reproducibility["model"]["kind"] == "diversification"


def test_build_diversification_figure_package_blocks_incomplete_sampling_metadata(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "diversification-figure-package-incomplete-sampling"
    result = build_diversification_figure_package(
        fixture("example_tree.nwk"),
        metadata_path=fixture("example_sampling_fractions_incomplete.tsv"),
        out_dir=output_dir,
    )

    assert result.audit.publication_ready is False
    assert result.audit.sampling_metadata_complete is False
    assert result.audit.lineage_curve_visible is True
    assert result.audit.clade_outlier_surface_visible is True
    assert result.audit.model_comparison_visible is True
    assert any(
        "sampling metadata is incomplete or invalid" in limitation
        for limitation in result.audit.limitations
    )
    assert result.review_path.exists()
    assert result.manifest_path.exists()
    assert result.sampling_report is not None
    assert result.sampling_report.complete is False
