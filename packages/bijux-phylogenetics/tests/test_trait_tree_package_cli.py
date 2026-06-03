from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def metadata_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "metadata" / name


def traits_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "metadata" / name


def test_cli_report_trait_tree_package_writes_publication_review_bundle(
    tmp_path: Path, capsys
) -> None:
    output_dir = tmp_path / "annotated-trait-tree-package"
    exit_code = main(
        [
            "report",
            "trait-tree-package",
            str(tree_fixture("example_tree_support_left.nwk")),
            "--metadata",
            str(metadata_fixture("example_metadata.tsv")),
            "--traits",
            str(traits_fixture("example_traits_validate.tsv")),
            "--label-column",
            "species",
            "--categorical-column",
            "habitat",
            "--continuous-column",
            "height_cm",
            "--metadata-strip-columns",
            "location",
            "--heatmap-columns",
            "height_cm",
            "--support-labels",
            "--out-dir",
            str(output_dir),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["publication_ready"] is True
    assert payload["metrics"]["required_surface_count"] == 5
    assert payload["metrics"]["complete_surface_count"] == 5
    assert payload["metrics"]["caption_ready"] is True
    assert payload["metrics"]["legible"] is True
    assert payload["data"]["audit"]["publication_ready"] is True
    assert payload["data"]["figure_package"]["caption_draft"]["caption_ready"] is True
    assert (output_dir / "annotated-trait-tree-review.html").exists()
    assert (output_dir / "annotation-coverage.tsv").exists()
    assert (output_dir / "annotation-surface-summary.tsv").exists()
    assert (output_dir / "annotated-trait-tree-reproducibility.manifest.json").exists()
