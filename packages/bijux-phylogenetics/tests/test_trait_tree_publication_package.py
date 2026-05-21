from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.render.annotated_trait_tree_package import (
    build_annotated_trait_tree_package,
)


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def metadata_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "metadata" / name


def traits_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "metadata" / name


def test_build_annotated_trait_tree_package_writes_publication_review_bundle(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "annotated-trait-tree-package"
    result = build_annotated_trait_tree_package(
        tree_fixture("example_tree_support_left.nwk"),
        out_dir=output_dir,
        metadata_path=metadata_fixture("example_metadata.tsv"),
        traits_path=traits_fixture("example_traits_validate.tsv"),
        label_column="species",
        categorical_column="habitat",
        continuous_column="height_cm",
        metadata_strip_columns=["location"],
        heatmap_columns=["height_cm"],
        title="Rabies host trait tree",
    )

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    reproducibility = json.loads(
        result.reproducibility_manifest_path.read_text(encoding="utf-8")
    )
    coverage_lines = result.coverage_path.read_text(encoding="utf-8").splitlines()
    summary_lines = result.summary_path.read_text(encoding="utf-8").splitlines()
    review_html = result.review_path.read_text(encoding="utf-8")

    assert result.audit.publication_ready is True
    assert result.audit.required_surface_count == 5
    assert result.audit.complete_surface_count == 5
    assert result.figure_package.caption_draft.caption_ready is True
    assert result.figure_package.legibility_audit.legible is True
    assert result.coverage_path.exists()
    assert result.summary_path.exists()
    assert result.review_path.exists()
    assert result.reproducibility_manifest_path.exists()
    assert coverage_lines[0].startswith("surface\tsource_kind\trequired")
    assert summary_lines[0].startswith("surface\tsource_kind\tvalue_kind")
    assert "Rabies host trait tree" in review_html
    assert manifest["report_kind"] == "annotated_trait_tree_package"
    assert manifest["audit"]["publication_ready"] is True
    assert manifest["reproducibility_manifest_path"] == str(
        result.reproducibility_manifest_path
    )
    assert reproducibility["report_kind"] == "annotated_trait_tree_package"
    assert {row["label"] for row in reproducibility["generated_tables"]} == {
        "tree_legend",
        "tree_annotations",
        "annotation_coverage",
        "annotation_surface_summary",
    }


def test_build_annotated_trait_tree_package_flags_incomplete_required_surfaces(
    tmp_path: Path,
) -> None:
    metadata_path = tmp_path / "metadata.tsv"
    metadata_path.write_text(
        "\n".join(
            [
                "taxon\tspecies\tlocation",
                "A\tAlpha species\tSweden",
                "B\tBeta species\tNorway",
                "C\tGamma species\t",
                "D\tDelta species\tFinland",
                "",
            ]
        ),
        encoding="utf-8",
    )

    result = build_annotated_trait_tree_package(
        tree_fixture("example_tree_support_left.nwk"),
        out_dir=tmp_path / "annotated-trait-tree-package",
        metadata_path=metadata_path,
        traits_path=traits_fixture("example_traits_validate.tsv"),
        label_column="species",
        categorical_column="habitat",
        metadata_strip_columns=["location"],
    )

    location_row = next(
        row for row in result.coverage_rows if row.surface == "location"
    )

    assert result.audit.publication_ready is False
    assert result.audit.missing_surface_count == 1
    assert location_row.complete is False
    assert location_row.missing_taxa == ["C"]
    assert any(
        "location is incomplete for publication review" in limitation
        for limitation in result.audit.limitations
    )
