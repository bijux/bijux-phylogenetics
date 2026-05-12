from __future__ import annotations

import json
from pathlib import Path
import shutil

from bijux_phylogenetics.cli import main

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


def test_ancestral_continuous_cli_can_export_table(tmp_path: Path, capsys) -> None:
    table_path = tmp_path / "ancestral.tsv"
    summary_path = tmp_path / "ancestral-summary.tsv"
    uncertainty_path = tmp_path / "ancestral-uncertainty.tsv"
    exclusions_path = tmp_path / "ancestral-excluded.tsv"
    exit_code = main(
        [
            "ancestral",
            "continuous",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--table-out",
            str(table_path),
            "--summary-out",
            str(summary_path),
            "--uncertainty-out",
            str(uncertainty_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["model"] == "brownian"
    assert payload["metrics"]["internal_node_count"] == 3
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["unstable_node_count"] >= 0
    assert "estimate\tstandard_error" in table_path.read_text(encoding="utf-8")
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\talpha"
    )
    assert uncertainty_path.read_text(encoding="utf-8").startswith(
        "node\tnode_name\tdescendant_taxa\testimate\tstandard_error"
    )
    assert exclusions_path.read_text(encoding="utf-8") == "taxon\treason\n"


def test_ancestral_discrete_cli_reports_sparse_state_warning(capsys) -> None:
    exit_code = main(
        [
            "ancestral",
            "discrete",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_ancestral_sparse.tsv")),
            "--trait",
            "habitat",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert (
        "one or more discrete states are represented by fewer than two taxa and should be interpreted cautiously"
        in payload["warnings"]
    )
    assert (
        "one or more discrete ancestral nodes remain unstable across candidate states"
        in payload["warnings"]
    )


def test_ancestral_discrete_cli_can_export_probability_review(
    tmp_path: Path, capsys
) -> None:
    table_path = tmp_path / "ancestral-discrete.tsv"
    summary_path = tmp_path / "ancestral-discrete-summary.tsv"
    probabilities_path = tmp_path / "ancestral-discrete-probabilities.tsv"
    exclusions_path = tmp_path / "ancestral-discrete-excluded.tsv"
    exit_code = main(
        [
            "ancestral",
            "discrete",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "equal-rates",
            "--table-out",
            str(table_path),
            "--summary-out",
            str(summary_path),
            "--probabilities-out",
            str(probabilities_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["internal_node_count"] == 3
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["unstable_node_count"] >= 0
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\tstate_ordering"
    )
    assert probabilities_path.read_text(encoding="utf-8").startswith(
        "node\tnode_name\tdescendant_taxa\tmost_likely_state\tstate_set"
    )
    assert exclusions_path.read_text(encoding="utf-8") == "taxon\treason\n"
    assert "most_likely_state\tstate_set" in table_path.read_text(encoding="utf-8")


def test_ancestral_discrete_cli_can_export_parsimony_comparison(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-discrete-summary.tsv"
    comparison_path = tmp_path / "ancestral-discrete-comparison.tsv"
    exit_code = main(
        [
            "ancestral",
            "discrete",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_ancestral_sparse.tsv")),
            "--trait",
            "habitat",
            "--model",
            "fitch",
            "--compare-model",
            "equal-rates",
            "--summary-out",
            str(summary_path),
            "--comparison-out",
            str(comparison_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["model"] == "fitch"
    assert payload["metrics"]["minimal_change_count"] == 1
    assert payload["metrics"]["ambiguous_internal_node_count"] == 1
    assert payload["metrics"]["comparison_node_count"] == 3
    assert payload["metrics"]["comparison_differing_node_count"] >= 0
    assert "minimal_change_count" in summary_path.read_text(encoding="utf-8")
    assert comparison_path.read_text(encoding="utf-8").startswith(
        "node\tdescendant_taxa\tleft_model\tright_model"
    )


def test_ancestral_tree_set_cli_can_export_continuous_review(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-tree-set-summary.tsv"
    trees_path = tmp_path / "ancestral-tree-set-trees.tsv"
    nodes_path = tmp_path / "ancestral-tree-set-nodes.tsv"
    clades_path = tmp_path / "ancestral-tree-set-clades.tsv"
    exclusions_path = tmp_path / "ancestral-tree-set-excluded.tsv"
    exit_code = main(
        [
            "ancestral",
            "tree-set",
            str(fixture("example_posterior_tree_set_six_taxa.nwk")),
            str(fixture("example_traits_clade_summary.tsv")),
            "--trait",
            "body_mass",
            "--kind",
            "continuous",
            "--summary-out",
            str(summary_path),
            "--trees-out",
            str(trees_path),
            "--nodes-out",
            str(nodes_path),
            "--clades-out",
            str(clades_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["kind"] == "continuous"
    assert payload["metrics"]["model"] == "brownian"
    assert payload["metrics"]["kept_tree_count"] == 5
    assert payload["metrics"]["clade_summary_count"] == 14
    assert payload["metrics"]["unstable_clade_count"] == 14
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\talpha"
    )
    assert trees_path.read_text(encoding="utf-8").startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id"
    )
    assert nodes_path.read_text(encoding="utf-8").startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id\tunrooted_topology_id\tclade_id"
    )
    assert clades_path.read_text(encoding="utf-8").startswith(
        "clade_id\tclade_taxa\ttree_presence_count\ttree_presence_fraction"
    )
    assert exclusions_path.read_text(encoding="utf-8") == "taxon\treason\n"


def test_ancestral_tree_set_cli_reports_discrete_stability_warnings(capsys) -> None:
    exit_code = main(
        [
            "ancestral",
            "tree-set",
            str(fixture("example_posterior_tree_set_six_taxa.nwk")),
            str(fixture("example_traits_clade_summary.tsv")),
            "--trait",
            "habitat",
            "--kind",
            "discrete",
            "--model",
            "equal-rates",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["kind"] == "discrete"
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["kept_tree_count"] == 5
    assert payload["metrics"]["clade_summary_count"] == 14
    assert payload["metrics"]["unstable_clade_count"] == 14
    assert (
        "one or more comparable ancestral clades are absent from some retained trees"
        in payload["warnings"]
    )
    assert (
        "one or more discrete ancestral clades change state or support profile across retained trees"
        in payload["warnings"]
    )


def test_ancestral_transitions_cli_can_export_branch_review(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-transition-summary.tsv"
    branch_path = tmp_path / "ancestral-transition-branches.tsv"
    count_path = tmp_path / "ancestral-transition-counts.tsv"
    exclusion_path = tmp_path / "ancestral-transition-excluded.tsv"
    exit_code = main(
        [
            "ancestral",
            "transitions",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_ancestral_sparse.tsv")),
            "--trait",
            "habitat",
            "--summary-out",
            str(summary_path),
            "--branches-out",
            str(branch_path),
            "--counts-out",
            str(count_path),
            "--exclusions-out",
            str(exclusion_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["tree_set"] is False
    assert payload["metrics"]["model"] == "fitch"
    assert payload["metrics"]["total_branch_count"] == 6
    assert payload["metrics"]["changed_branch_count"] == 2
    assert payload["metrics"]["uncertain_change_count"] == 2
    assert payload["metrics"]["transition_pair_count"] == 2
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\tstate_ordering"
    )
    assert branch_path.read_text(encoding="utf-8").startswith(
        "parent_node\tchild_node\tchild_descendant_taxa\tbranch_length"
    )
    assert count_path.read_text(encoding="utf-8").startswith(
        "transition\tsource_state\ttarget_state\tcertain_change_count"
    )
    assert exclusion_path.read_text(encoding="utf-8") == "taxon\treason\n"


def test_ancestral_transitions_cli_can_export_tree_set_review(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-transition-tree-set-summary.tsv"
    tree_path = tmp_path / "ancestral-transition-tree-set-trees.tsv"
    branch_path = tmp_path / "ancestral-transition-tree-set-branches.tsv"
    count_path = tmp_path / "ancestral-transition-tree-set-counts.tsv"
    exclusion_path = tmp_path / "ancestral-transition-tree-set-excluded.tsv"
    exit_code = main(
        [
            "ancestral",
            "transitions",
            str(fixture("example_posterior_tree_set_six_taxa.nwk")),
            str(fixture("example_traits_clade_summary.tsv")),
            "--trait",
            "habitat",
            "--model",
            "equal-rates",
            "--tree-set",
            "--summary-out",
            str(summary_path),
            "--trees-out",
            str(tree_path),
            "--branches-out",
            str(branch_path),
            "--counts-out",
            str(count_path),
            "--exclusions-out",
            str(exclusion_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["tree_set"] is True
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["kept_tree_count"] == 5
    assert payload["metrics"]["transition_pair_count"] >= 2
    assert payload["metrics"]["topology_sensitive_transition_pair_count"] >= 1
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\tstate_ordering"
    )
    assert tree_path.read_text(encoding="utf-8").startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id"
    )
    assert branch_path.read_text(encoding="utf-8").startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id\tunrooted_topology_id\tparent_node"
    )
    assert count_path.read_text(encoding="utf-8").startswith(
        "transition\tsource_state\ttarget_state\ttree_presence_count"
    )
    assert exclusion_path.read_text(encoding="utf-8") == "taxon\treason\n"


def test_ancestral_render_cli_writes_svg_with_internal_annotations(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "ancestral.svg"
    exit_code = main(
        [
            "ancestral",
            "render",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["format"] == "svg"
    assert payload["metrics"]["rendered_internal_annotation_count"] == 3
    assert 'class="internal-annotation-label"' in output.read_text(encoding="utf-8")


def test_ancestral_render_cli_can_export_discrete_pie_html(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "ancestral-discrete.html"
    exit_code = main(
        [
            "ancestral",
            "render",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--kind",
            "discrete",
            "--model",
            "equal-rates",
            "--discrete-node-style",
            "pies",
            "--branch-coloring",
            "state",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["format"] == "html"
    assert payload["metrics"]["rendered_internal_pie_count"] == 3
    assert payload["metrics"]["rendered_branch_color_count"] >= 4
    assert "<figure><svg" in output.read_text(encoding="utf-8")
    assert output.with_suffix(".svg").exists()


def test_ancestral_render_cli_can_export_png(tmp_path: Path, capsys) -> None:
    if shutil.which("rsvg-convert") is None and shutil.which("sips") is None:
        return
    output = tmp_path / "ancestral-discrete.png"
    exit_code = main(
        [
            "ancestral",
            "render",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--kind",
            "discrete",
            "--model",
            "equal-rates",
            "--discrete-node-style",
            "pies",
            "--branch-coloring",
            "state",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["format"] == "png"
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert output.with_suffix(".svg").exists()


def test_ancestral_report_cli_writes_html_and_svg(tmp_path: Path, capsys) -> None:
    output = tmp_path / "ancestral-report.html"
    exit_code = main(
        [
            "ancestral",
            "report",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--compare-model",
            "ou",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["report_kind"] == "ancestral-state"
    assert output.exists()
    assert output.with_suffix(".svg").exists()


def test_ancestral_sensitivity_cli_reports_available_comparisons(capsys) -> None:
    exit_code = main(
        [
            "ancestral",
            "sensitivity",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--compare-model",
            "ou",
            "--compare-tree",
            str(fixture("example_tree_topology_diff.nwk")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["has_model_sensitivity"] is True
    assert payload["metrics"]["has_tree_sensitivity"] is True


def test_ancestral_package_cli_writes_publication_bundle(
    tmp_path: Path, capsys
) -> None:
    output_dir = tmp_path / "ancestral-package"
    exit_code = main(
        [
            "ancestral",
            "package",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--out-dir",
            str(output_dir),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["artifact_count"] == 9
    assert (output_dir / "figure-manifest.json").exists()
    assert (output_dir / "ancestral-figure.png").exists()
    assert (output_dir / "ancestral-figure.html").exists()


def test_ancestral_discrete_cli_rejects_ordered_fitch_model(capsys) -> None:
    try:
        main(
            [
                "ancestral",
                "report",
                str(fixture("example_tree.nwk")),
                str(fixture("example_traits_geography.tsv")),
                "--trait",
                "region",
                "--kind",
                "discrete",
                "--state-ordering",
                "ordered",
                "--out",
                "ignored.html",
            ]
        )
    except SystemExit as error:
        assert error.code == 2
    captured = capsys.readouterr()
    assert (
        "ordered ancestral discrete reconstruction requires a likelihood model"
        in captured.err
    )
