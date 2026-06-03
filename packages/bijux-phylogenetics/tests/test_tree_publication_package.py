from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.render.tree_figure_package import build_tree_figure_package
from bijux_phylogenetics.render.tree_svg import AnnotationStrip


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_build_tree_figure_package_writes_publication_legend_and_caption_artifacts(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "tree-publication-package"
    result = build_tree_figure_package(
        tree_fixture("example_tree_support_left.nwk"),
        out_dir=output_dir,
        title="Rabies host phylogeny",
        show_support_values=True,
        categorical_traits={
            "A": "bat",
            "B": "canid",
            "C": "canid",
            "D": "livestock",
        },
        metadata_strips=[
            AnnotationStrip(
                "region",
                {"A": "north", "B": "north", "C": "south", "D": "south"},
            )
        ],
        heatmap_columns=[
            AnnotationStrip(
                "height_cm",
                {"A": "1.2", "B": "1.4", "C": "1.8", "D": "2.0"},
            )
        ],
    )

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    reproducibility = json.loads(
        result.reproducibility_manifest_path.read_text(encoding="utf-8")
    )
    legend_lines = result.legend_path.read_text(encoding="utf-8").splitlines()
    caption = result.caption_path.read_text(encoding="utf-8")

    assert result.figure_path.exists()
    assert result.caption_path.exists()
    assert result.legend_path.exists()
    assert result.reproducibility_manifest_path.exists()
    assert result.caption_draft.caption_ready is True
    assert result.legibility_audit.legible is True
    assert len(result.legend_entries) >= 5
    assert legend_lines[0] == "surface\tlabel\tswatch\tdetail"
    assert any(
        line.startswith("branch-length\tscale bar\t") for line in legend_lines[1:]
    )
    assert any("validated support labels" in line for line in legend_lines[1:])
    assert "## Draft Caption" in caption
    assert "## Figure Specifications" in caption
    assert "## Legibility" in caption
    assert manifest["legend_path"] == str(result.legend_path)
    assert manifest["reproducibility_manifest_path"] == str(
        result.reproducibility_manifest_path
    )
    assert manifest["caption_draft"]["caption_ready"] is True
    assert manifest["legibility_audit"]["legible"] is True
    assert reproducibility["report_kind"] == "tree_figure_package"
    assert reproducibility["input_files"][0]["label"] == "tree"
    assert {row["label"] for row in reproducibility["generated_tables"]} == {
        "figure_legend",
        "tip_annotations",
    }


def test_build_tree_figure_package_flags_long_labels_against_publication_lane(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "tree-publication-package"
    result = build_tree_figure_package(
        tree_fixture("example_tree.nwk"),
        out_dir=output_dir,
        labels={
            "A": "Artibeus_jamaicensis_long_publication_label",
            "B": "Beta",
            "C": "Gamma",
            "D": "Delta",
        },
        metadata_strips=[
            AnnotationStrip(
                "region",
                {"A": "north", "B": "north", "C": "south", "D": "south"},
            )
        ],
    )

    assert result.legibility_audit.legible is False
    assert any(
        "reserved label lane" in warning for warning in result.legibility_audit.warnings
    )


def test_build_tree_figure_package_treats_taxon_names_as_default_labels(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "tree-publication-package"
    result = build_tree_figure_package(
        tree_fixture("example_tree.nwk"),
        out_dir=output_dir,
    )

    label_coverage = next(
        row for row in result.audit.annotation_coverage if row.surface == "labels"
    )

    assert label_coverage.aligned is True
    assert label_coverage.covered_taxa == 4
    assert label_coverage.missing_taxa == []
