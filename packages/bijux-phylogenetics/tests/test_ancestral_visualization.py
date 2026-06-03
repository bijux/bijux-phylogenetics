from __future__ import annotations

from pathlib import Path
import shutil

import pytest

from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.ancestral.presentation import (
    render_ancestral_state_visualization,
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


def test_render_ancestral_state_visualization_writes_continuous_svg_with_regime_branch_colors(
    tmp_path: Path,
) -> None:
    report = reconstruct_continuous_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    output = tmp_path / "ancestral-continuous.svg"

    result = render_ancestral_state_visualization(
        fixture("example_tree.nwk"),
        report,
        out_path=output,
        layout="phylogram",
        branch_coloring="regime",
    )

    svg = output.read_text(encoding="utf-8")
    assert result.format == "svg"
    assert result.branch_coloring == "regime"
    assert result.tree_render.rendered_branch_color_count >= 4
    assert result.tree_render.rendered_internal_annotation_count == 3
    assert 'class="branch branch-colored"' in svg
    assert 'class="internal-annotation-label"' in svg


def test_render_ancestral_state_visualization_writes_discrete_html_with_pies(
    tmp_path: Path,
) -> None:
    report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
    )
    output = tmp_path / "ancestral-discrete.html"

    result = render_ancestral_state_visualization(
        fixture("example_tree.nwk"),
        report,
        out_path=output,
        discrete_node_style="pies",
        branch_coloring="state",
    )

    html = output.read_text(encoding="utf-8")
    svg = result.svg_path.read_text(encoding="utf-8")
    assert result.format == "html"
    assert result.discrete_node_style == "pies"
    assert result.tree_render.rendered_internal_pie_count == 3
    assert result.tree_render.rendered_branch_color_count >= 4
    assert "<figure><svg" in html
    assert "Bijux Ancestral Visualization" in html
    assert 'class="internal-pie-slice"' in svg


def test_render_ancestral_state_visualization_writes_png(tmp_path: Path) -> None:
    if shutil.which("rsvg-convert") is None and shutil.which("sips") is None:
        pytest.skip("no SVG-to-PNG converter is available")
    report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
    )
    output = tmp_path / "ancestral-discrete.png"

    result = render_ancestral_state_visualization(
        fixture("example_tree.nwk"),
        report,
        out_path=output,
        discrete_node_style="pies",
        branch_coloring="state",
    )

    assert result.format == "png"
    assert output.exists()
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert result.svg_path.exists()
