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


def test_cli_render_reports_publication_package_metrics(tmp_path: Path, capsys) -> None:
    output = tmp_path / "annotated.svg"
    package_dir = tmp_path / "tree-publication-package"
    exit_code = main(
        [
            "render",
            str(tree_fixture("example_tree_support_left.nwk")),
            "--metadata",
            str(metadata_fixture("example_metadata.tsv")),
            "--label-column",
            "species",
            "--metadata-strip-columns",
            "location",
            "--traits",
            str(traits_fixture("example_traits_validate.tsv")),
            "--layout",
            "phylogram",
            "--support-labels",
            "--categorical-column",
            "habitat",
            "--continuous-column",
            "height_cm",
            "--heatmap-columns",
            "height_cm,status",
            "--package-dir",
            str(package_dir),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["figure_package_legible"] is True
    assert payload["metrics"]["figure_package_caption_ready"] is True
    assert payload["metrics"]["figure_package_legend_entry_count"] >= 5
    assert payload["data"]["figure_package_caption_draft"]["caption_ready"] is True
    assert payload["data"]["figure_package_legibility_audit"]["legible"] is True
    assert len(payload["data"]["figure_package_legend_entries"]) >= 5
    assert (package_dir / "figure.svg").exists()
    assert (package_dir / "figure-caption.md").exists()
    assert (package_dir / "figure-legend.tsv").exists()
