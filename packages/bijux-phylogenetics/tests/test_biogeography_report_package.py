from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.biogeography.report_package import (
    build_biogeography_report_package,
)

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


def test_build_biogeography_report_package_writes_review_bundle(
    tmp_path: Path,
) -> None:
    result = build_biogeography_report_package(
        tree_path=fixture("example_tree.nwk"),
        traits_path=fixture("example_traits_geography.tsv"),
        centroids_path=fixture("example_geographic_region_centroids.tsv"),
        trait="region",
        out_dir=tmp_path / "biogeography-report",
        model="ard",
    )

    assert result.report_path.exists()
    assert result.tree_figure_path.exists()
    assert result.map_path.exists()
    assert result.summary_table_path.exists()
    assert result.region_count_table_path.exists()
    assert result.node_table_path.exists()
    assert result.transition_matrix_path.exists()
    assert result.event_table_path.exists()
    assert result.map_marker_table_path.exists()
    assert result.map_line_table_path.exists()
    assert result.exclusion_table_path.exists()
    assert result.manifest_path.exists()

    assert "Bijux Biogeography Report" in result.report_path.read_text(encoding="utf-8")
    assert "Regional Transition Map Review" in result.map_path.read_text(
        encoding="utf-8"
    )
    assert result.tree_figure_path.read_text(encoding="utf-8").startswith("<svg")

    assert result.summary_table_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\tinternal_model"
    )
    assert result.region_count_table_path.read_text(encoding="utf-8").startswith(
        "region\ttip_taxon_count\tanalyzed_taxon_fraction"
    )
    assert result.node_table_path.read_text(encoding="utf-8").startswith(
        "node\tnode_name\tdescendant_taxa\tmost_likely_region"
    )
    assert result.transition_matrix_path.read_text(encoding="utf-8").startswith(
        "source_region\ttarget_region\trate\tlower_95_interval"
    )
    assert result.event_table_path.read_text(encoding="utf-8").startswith(
        "branch_id\tparent_node\tchild_node\tchild_descendant_taxa"
    )
    assert result.map_marker_table_path.read_text(encoding="utf-8").startswith(
        "marker_id\tlabel\tmarker_kind\tlatitude"
    )
    assert result.map_line_table_path.read_text(encoding="utf-8").startswith(
        "line_id\tline_kind\tsource_label\ttarget_label"
    )
    assert result.exclusion_table_path.read_text(encoding="utf-8").startswith(
        "surface\tsubject_id\tsubject_kind\traw_left"
    )

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest["report_kind"] == "biogeography_report_package"
    assert manifest["metrics"]["observed_region_count"] == 3
    assert manifest["metrics"]["event_count"] == 2
    assert manifest["metrics"]["visible_map_line_count"] >= 0
