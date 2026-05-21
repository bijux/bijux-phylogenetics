from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.ancestral.comparison import (
    compare_continuous_ancestral_models,
    compare_continuous_ancestral_trees,
    compare_discrete_ancestral_models,
    compare_discrete_ancestral_reconstructions,
    compare_discrete_ancestral_trees,
)
from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.ancestral.presentation.report_rendering import (
    render_ancestral_state_report,
    render_ancestral_state_tree,
    write_ancestral_state_table,
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


def test_compare_continuous_ancestral_models_reports_node_deltas() -> None:
    report = compare_continuous_ancestral_models(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        left_model="brownian",
        right_model="ou",
        right_alpha=1.5,
    )
    assert report.left_model == "brownian"
    assert report.right_model == "ou"
    assert any(row.estimate_delta != 0 for row in report.rows if row.node == "A|B|C|D")


def test_write_ancestral_state_table_exports_continuous_and_discrete_rows(
    tmp_path: Path,
) -> None:
    continuous = reconstruct_continuous_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    discrete = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="habitat",
    )
    continuous_path = tmp_path / "continuous.tsv"
    discrete_path = tmp_path / "discrete.tsv"
    write_ancestral_state_table(continuous_path, continuous)
    write_ancestral_state_table(discrete_path, discrete)
    assert "estimate\tstandard_error" in continuous_path.read_text(encoding="utf-8")
    assert "most_likely_state\tstate_set" in discrete_path.read_text(encoding="utf-8")


@pytest.mark.slow
def test_compare_discrete_ancestral_models_selects_supported_model() -> None:
    report = compare_discrete_ancestral_models(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
    )
    assert len(report.rows) == 3
    assert report.selected_model in {"equal-rates", "symmetric", "all-rates-different"}
    assert any(
        difference.comparison_model == "symmetric"
        for difference in report.node_differences
    )


def test_compare_discrete_ancestral_reconstructions_supports_fitch_baseline() -> None:
    report = compare_discrete_ancestral_reconstructions(
        fixture("example_tree.nwk"),
        fixture("example_traits_ancestral_sparse.tsv"),
        trait="habitat",
        left_model="fitch",
        right_model="equal-rates",
    )
    assert report.left_model == "fitch"
    assert report.right_model == "equal-rates"
    assert report.left_minimal_change_count == 1
    assert report.right_minimal_change_count is None
    assert len(report.rows) == 3
    assert any(row.ambiguity_changed for row in report.rows)


def test_compare_ancestral_reconstructions_across_trees_reports_shared_nodes() -> None:
    continuous = compare_continuous_ancestral_trees(
        fixture("example_tree.nwk"),
        fixture("example_tree_topology_diff.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    discrete = compare_discrete_ancestral_trees(
        fixture("example_tree.nwk"),
        fixture("example_tree_topology_diff.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
    )
    assert continuous.rows
    assert discrete.rows


def test_render_ancestral_state_tree_adds_internal_annotations(tmp_path: Path) -> None:
    report = reconstruct_continuous_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    output = tmp_path / "ancestral.svg"
    result = render_ancestral_state_tree(
        fixture("example_tree.nwk"), report, out_path=output, layout="phylogram"
    )
    svg = output.read_text(encoding="utf-8")
    assert result.rendered_internal_annotation_count == 3
    assert 'class="internal-annotation-label"' in svg
    assert "2.25" in svg


def test_render_ancestral_state_report_writes_html_and_svg(tmp_path: Path) -> None:
    output = tmp_path / "ancestral-report.html"
    result = render_ancestral_state_report(
        tree_path=fixture("example_tree.nwk"),
        traits_path=fixture("example_traits_comparative.tsv"),
        trait="response",
        reconstruction_kind="continuous",
        model="brownian",
        compare_model="ou",
        out_path=output,
    )
    html = output.read_text(encoding="utf-8")
    manifest = json.loads(
        html.split('<script id="bijux-report-manifest" type="application/json">', 1)[
            1
        ].split("</script>", 1)[0]
    )
    assert result.output_path == output
    assert output.with_suffix(".svg").exists()
    assert "ancestral-reconstruction" in html
    assert "ancestral-comparison" in html
    assert "ancestral-sensitivity" in html
    assert "ancestral-node-table" in html
    assert "limitations" in html
    assert manifest["report_kind"] == "ancestral-state"
    assert manifest["supplement_sections"] == [
        "ancestral-methods",
        "ancestral-exclusions",
        "ancestral-node-table",
        "ancestral-uncertainty",
        "ancestral-sensitivity",
        "limitations",
    ]
    assert result.supplement_sections == manifest["supplement_sections"]
    assert result.sensitivity is not None
    assert manifest["limitations"]


def test_render_ancestral_state_report_supports_fitch_comparison(
    tmp_path: Path,
) -> None:
    output = tmp_path / "ancestral-discrete-report.html"
    result = render_ancestral_state_report(
        tree_path=fixture("example_tree.nwk"),
        traits_path=fixture("example_traits_ancestral_sparse.tsv"),
        trait="habitat",
        reconstruction_kind="discrete",
        model="fitch",
        compare_model="equal-rates",
        out_path=output,
    )
    html = output.read_text(encoding="utf-8")
    assert result.output_path == output
    assert "ancestral-comparison" in html
    assert "limitations" in html
