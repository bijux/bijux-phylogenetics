from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.render.reproducibility import (
    FigureReproducibilityFilter,
    write_figure_reproducibility_manifest,
)


def test_write_figure_reproducibility_manifest_tracks_inputs_outputs_and_settings(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "input.nwk"
    figure_path = tmp_path / "figure.svg"
    table_path = tmp_path / "table.tsv"
    note_path = tmp_path / "caption.md"
    manifest_path = tmp_path / "figure-reproducibility.manifest.json"

    input_path.write_text("(A:1,B:1);\n", encoding="utf-8")
    figure_path.write_text("<svg />\n", encoding="utf-8")
    table_path.write_text("column\nvalue\n", encoding="utf-8")
    note_path.write_text("# Caption\n", encoding="utf-8")

    manifest = write_figure_reproducibility_manifest(
        manifest_path,
        report_kind="tree_figure_package",
        input_files=[("tree", input_path)],
        generated_figures=[("figure", figure_path)],
        generated_tables=[("table", table_path)],
        filters=[
            FigureReproducibilityFilter(
                name="collapsed_clades",
                value="bat-clade",
                detail="collapse one named clade before rendering the figure",
            )
        ],
        model={
            "kind": "none",
            "name": None,
        },
        settings={
            "layout": "phylogram",
            "show_support_values": True,
        },
        linked_artifacts=[("caption", note_path)],
    )

    payload = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert manifest_path.exists()
    assert payload == manifest
    assert payload["report_kind"] == "tree_figure_package"
    assert payload["format_version"] == 1
    assert payload["input_files"][0]["label"] == "tree"
    assert payload["generated_figures"][0]["label"] == "figure"
    assert payload["generated_tables"][0]["label"] == "table"
    assert payload["filters"][0]["name"] == "collapsed_clades"
    assert payload["settings"]["layout"] == "phylogram"
    assert payload["linked_artifacts"][0]["label"] == "caption"
