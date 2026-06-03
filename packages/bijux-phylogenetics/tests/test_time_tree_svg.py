from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.render.time_tree_svg import (
    TimeTreeNodeInterval,
    render_time_tree_svg,
)


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_render_time_tree_svg_writes_age_labels_hpd_intervals_and_axis(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "time-tree.svg"
    result = render_time_tree_svg(
        tree_fixture("example_tree_support_left.nwk"),
        out_path=output_path,
        node_intervals=[
            TimeTreeNodeInterval(
                clade="A|B",
                node_kind="internal",
                mean_age=0.1,
                median_age=0.1,
                minimum_age=0.1,
                maximum_age=0.1,
                lower_95_credible_interval=0.08,
                upper_95_credible_interval=0.12,
                tree_count=12,
            ),
            TimeTreeNodeInterval(
                clade="C|D",
                node_kind="internal",
                mean_age=0.1,
                median_age=0.1,
                minimum_age=0.1,
                maximum_age=0.1,
                lower_95_credible_interval=0.08,
                upper_95_credible_interval=0.12,
                tree_count=12,
            ),
            TimeTreeNodeInterval(
                clade="A|B|C|D",
                node_kind="root",
                mean_age=0.3,
                median_age=0.3,
                minimum_age=0.28,
                maximum_age=0.32,
                lower_95_credible_interval=0.27,
                upper_95_credible_interval=0.33,
                tree_count=12,
            ),
        ],
        title="Rabies time tree",
    )

    svg = output_path.read_text(encoding="utf-8")

    assert result.output_path == output_path
    assert result.tip_count == 4
    assert result.internal_node_count == 3
    assert result.rendered_interval_count == 3
    assert result.rendered_age_label_count == 3
    assert result.axis_tick_count >= 2
    assert result.ultrametric is True
    assert "node age before present" in svg
    assert "Rabies time tree" in svg
    assert "hpd-interval" in svg
    assert "node-age-label" in svg
