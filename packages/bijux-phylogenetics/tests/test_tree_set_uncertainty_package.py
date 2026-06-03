from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.trees.uncertainty import (
    build_tree_set_uncertainty_figure_package,
)

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_build_tree_set_uncertainty_figure_package_writes_publication_bundle(
    tmp_path: Path,
) -> None:
    result = build_tree_set_uncertainty_figure_package(
        fixture("example_tree_set_left.nwk"),
        out_dir=tmp_path / "tree-set-uncertainty-package",
    )

    assert result.consensus_tree_path.exists()
    assert result.consensus_figure_path.exists()
    assert result.clade_support_plot_path.exists()
    assert result.unstable_taxa_plot_path.exists()
    assert result.topology_clusters_plot_path.exists()
    assert result.unstable_taxa_table_path.exists()
    assert result.topology_clusters_table_path.exists()
    assert result.uncertainty_conclusions_table_path.exists()
    assert result.conclusion_summary_path.exists()
    assert result.legend_path.exists()
    assert result.caption_path.exists()
    assert result.methods_summary_path.exists()
    assert result.review_path.exists()
    assert result.manifest_path.exists()
    assert result.reproducibility_manifest_path.exists()
    assert result.audit.publication_ready is True
    assert result.audit.support_labels_validated is True
    assert result.audit.plotted_unstable_taxon_count > 0
    assert result.audit.plotted_topology_cluster_count > 0
    assert result.consensus_render.rendered_support_count > 0

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    reproducibility = json.loads(
        result.reproducibility_manifest_path.read_text(encoding="utf-8")
    )
    assert manifest["report_kind"] == "tree_set_uncertainty_figure_package"
    assert manifest["audit"]["publication_ready"] is True
    assert manifest["outputs"]["methods_summary_path"].endswith(
        "tree-set-uncertainty-methods-summary.md"
    )
    assert manifest["metrics"]["methods_summary_warning_count"] >= 0
    assert manifest["reproducibility_manifest_path"] == str(
        result.reproducibility_manifest_path
    )
    assert reproducibility["report_kind"] == "tree_set_uncertainty_figure_package"
    assert reproducibility["generated_figures"][0]["label"] == "consensus_tree"
    assert (
        "Tree-Set Uncertainty Methods Summary"
        in result.methods_summary_path.read_text(encoding="utf-8")
    )
    assert "Methods Summary" in result.review_path.read_text(encoding="utf-8")


def test_build_tree_set_uncertainty_figure_package_keeps_empty_instability_panel_explicit(
    tmp_path: Path,
) -> None:
    tree_set_path = tmp_path / "single-topology-tree-set.nwk"
    source_tree = (
        FIXTURES.joinpath("example_tree_set_left.nwk")
        .read_text(encoding="utf-8")
        .splitlines()[0]
    )
    tree_set_path.write_text(
        "\n".join([source_tree, source_tree, source_tree]) + "\n", encoding="utf-8"
    )

    result = build_tree_set_uncertainty_figure_package(
        tree_set_path,
        out_dir=tmp_path / "tree-set-uncertainty-package",
    )

    assert result.audit.publication_ready is True
    assert result.audit.unstable_taxon_count == 0
    assert result.audit.plotted_unstable_taxon_count == 0
    assert result.audit.unstable_taxa_visible is True
    assert result.audit.plotted_topology_cluster_count == 1
    assert "No unstable taxa were detected" in result.unstable_taxa_plot_path.read_text(
        encoding="utf-8"
    )
