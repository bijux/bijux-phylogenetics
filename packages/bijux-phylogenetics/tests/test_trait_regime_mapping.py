from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.traits.regime_mapping import (
    render_trait_regime_map,
    summarize_trait_regime_mapping,
    write_trait_regime_branch_table,
    write_trait_regime_exclusion_table,
    write_trait_regime_node_table,
    write_trait_regime_summary_table,
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


def test_summarize_trait_regime_mapping_reconstructs_branch_regimes_from_tip_states() -> (
    None
):
    report = summarize_trait_regime_mapping(
        fixture("example_tree.nwk"),
        tip_states_path=fixture("example_traits_geography.tsv"),
        trait="region",
    )

    assert report.source_kind == "tip-state-reconstruction"
    assert report.tree_taxon_count == 4
    assert report.analyzed_taxon_count == 4
    assert report.branch_id_column == "branch_id"
    assert report.regime_column == "regime"
    assert report.observed_regimes == ["island", "north", "south"]
    assert len(report.branch_rows) == 6
    assert len(report.node_rows) == 7
    assert report.ambiguous_branch_count == 1
    branch_lookup = {row.branch_id: row for row in report.branch_rows}
    node_lookup = {row.node_id: row for row in report.node_rows}
    assert branch_lookup["A|B"].regime == "north"
    assert branch_lookup["A|B"].ambiguous_assignment is False
    assert branch_lookup["C|D"].candidate_regimes == ["island", "south"]
    assert branch_lookup["C|D"].ambiguous_assignment is True
    assert node_lookup["A|B|C|D"].candidate_regimes == ["island", "north", "south"]
    assert node_lookup["A|B|C|D"].ambiguous_assignment is True
    assert report.analysis_tree_newick is not None
    assert (
        "one or more discrete ancestral nodes remain unstable across candidate states"
        in report.warnings
    )


def test_summarize_trait_regime_mapping_normalizes_user_map() -> None:
    report = summarize_trait_regime_mapping(
        fixture("example_tree.nwk"),
        regime_map_path=fixture("example_branch_regimes.tsv"),
    )

    assert report.source_kind == "user-provided-map"
    assert report.trait is None
    assert report.analyzed_taxon_count == 4
    assert report.node_rows == []
    assert report.ambiguous_branch_count == 0
    assert [
        row.regime for row in report.branch_rows if row.branch_id in {"A|B", "C|D"}
    ] == [
        "fast",
        "slow",
    ]


def test_trait_regime_mapping_tracks_missing_and_extra_tip_states() -> None:
    report = summarize_trait_regime_mapping(
        fixture("example_tree.nwk"),
        tip_states_path=fixture("example_traits_regime_mapping_missing.tsv"),
        trait="region",
    )

    assert report.analyzed_taxa == ["A", "C", "D"]
    assert {row.taxon: row.reason for row in report.excluded_taxa} == {
        "B": "missing_state_value",
        "E": "absent_from_tree",
    }


def test_trait_regime_mapping_writers_and_renderer_emit_review_outputs(
    tmp_path: Path,
) -> None:
    reconstructed = summarize_trait_regime_mapping(
        fixture("example_tree.nwk"),
        tip_states_path=fixture("example_traits_geography.tsv"),
        trait="region",
    )
    normalized = summarize_trait_regime_mapping(
        fixture("example_tree.nwk"),
        regime_map_path=fixture("example_branch_regimes.tsv"),
    )
    summary_out = tmp_path / "trait-regime-summary.tsv"
    branch_out = tmp_path / "trait-regime-branches.tsv"
    node_out = tmp_path / "trait-regime-nodes.tsv"
    exclusion_out = tmp_path / "trait-regime-excluded.tsv"
    svg_out = tmp_path / "trait-regime.svg"

    write_trait_regime_summary_table(summary_out, reconstructed)
    write_trait_regime_branch_table(branch_out, reconstructed)
    write_trait_regime_node_table(node_out, reconstructed)
    write_trait_regime_exclusion_table(exclusion_out, reconstructed)
    render = render_trait_regime_map(normalized, out_path=svg_out, layout="circular")

    assert (
        summary_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("source_kind\ttrait\ttaxon_column\treconstruction_model")
    )
    assert (
        branch_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("branch_id\tchild_node_name\tis_tip_branch")
    )
    assert (
        node_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("node_id\tnode_name\tis_tip")
    )
    assert exclusion_out.read_text(encoding="utf-8").splitlines() == ["taxon\treason"]
    svg = svg_out.read_text(encoding="utf-8")
    assert render.layout == "circular"
    assert render.rendered_internal_annotation_count >= 2
    assert render.rendered_categorical_trait_count == 4
    assert 'class="internal-annotation-label"' in svg
    assert "categorical trait" in svg
