from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.bayesian.presentation.time_tree_figure_bundle import (
    build_time_tree_figure_package,
)


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "metadata" / name


def metadata_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "metadata" / name


def test_build_time_tree_figure_package_writes_publication_review_bundle(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "time-tree-package"
    result = build_time_tree_figure_package(
        tree_fixture("beast2_strict_yule_posterior.trees"),
        out_dir=output_dir,
        source_format="beast",
        burnin_fraction=0.25,
        metadata_path=metadata_fixture("example_metadata.tsv"),
        label_column="species",
        title="Rabies time tree",
    )

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    reproducibility = json.loads(
        result.reproducibility_manifest_path.read_text(encoding="utf-8")
    )
    review_html = result.review_path.read_text(encoding="utf-8")
    interval_rows = result.interval_table_path.read_text(encoding="utf-8").splitlines()

    assert result.audit.publication_ready is True
    assert result.retained_tree_count > 0
    assert result.figure_path.exists()
    assert result.retained_tree_set_path.exists()
    assert result.mcc_tree_path.exists()
    assert result.interval_table_path.exists()
    assert result.legend_path.exists()
    assert result.caption_path.exists()
    assert result.review_path.exists()
    assert result.reproducibility_manifest_path.exists()
    assert result.render.rendered_interval_count == result.render.internal_node_count
    assert manifest["report_kind"] == "time_tree_package"
    assert manifest["audit"]["publication_ready"] is True
    assert manifest["reproducibility_manifest_path"] == str(
        result.reproducibility_manifest_path
    )
    assert reproducibility["report_kind"] == "time_tree_package"
    assert reproducibility["model"]["name"] == "beast"
    assert "Rabies time tree" in review_html
    assert interval_rows[0].startswith("clade\tnode_kind\tmean_age")


def test_build_time_tree_figure_package_blocks_non_ultrametric_tree_sets(
    tmp_path: Path,
) -> None:
    tree_set_path = tmp_path / "non-ultrametric-posterior.nwk"
    tree_set_path.write_text(
        "\n".join(
            [
                "((A:0.1,B:0.2):0.2,(C:0.1,D:0.1):0.2);",
                "((A:0.1,B:0.2):0.2,(C:0.1,D:0.1):0.2);",
                "",
            ]
        ),
        encoding="utf-8",
    )

    result = build_time_tree_figure_package(
        tree_set_path,
        out_dir=tmp_path / "time-tree-package",
        source_format="generic",
        burnin_fraction=0.0,
        title="Non-ultrametric time tree",
    )

    assert result.audit.publication_ready is False
    assert result.audit.ultrametric is False
    assert any(
        "not ultrametric" in limitation for limitation in result.audit.limitations
    )
