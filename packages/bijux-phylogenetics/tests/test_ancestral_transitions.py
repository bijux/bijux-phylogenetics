from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.discrete.review import (
    summarize_ancestral_transition_report,
    summarize_ancestral_transition_tree_set,
    summarize_ancestral_transition_tree_set_report,
    summarize_ancestral_transitions,
    write_ancestral_transition_branch_table,
    write_ancestral_transition_count_table,
    write_ancestral_transition_exclusion_table,
    write_ancestral_transition_summary_table,
    write_ancestral_transition_tree_set_branch_table,
    write_ancestral_transition_tree_set_count_table,
    write_ancestral_transition_tree_set_summary_table,
    write_ancestral_transition_tree_set_tree_table,
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


def test_summarize_ancestral_transitions_separates_certain_and_uncertain_changes(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "tree.nwk"
    table_path = tmp_path / "traits.tsv"
    tree_path.write_text(
        "(((A:1,B:1):1,(C:1,D:1):1):5,(E:0.1,F:0.1):0.1);\n",
        encoding="utf-8",
    )
    table_path.write_text(
        "\n".join(
            [
                "taxon\thabitat",
                "A\tnorth",
                "B\tnorth",
                "C\tnorth",
                "D\tnorth",
                "E\tsouth",
                "F\tsouth",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = summarize_ancestral_transitions(
        tree_path,
        table_path,
        trait="habitat",
        model="equal-rates",
    )
    summary = summarize_ancestral_transition_report(report)
    certain_transition = report.transition_rows[0]

    assert summary.total_branch_count == 10
    assert summary.changed_branch_count == 1
    assert summary.certain_change_count == 1
    assert summary.uncertain_change_count == 0
    assert certain_transition.certain_change_count == 1
    assert certain_transition.uncertain_change_count == 0
    assert certain_transition.transition == "south->north"


def test_summarize_ancestral_transitions_tracks_missing_tip_exclusions(
    tmp_path: Path,
) -> None:
    table_path = tmp_path / "traits.tsv"
    table_path.write_text(
        "\n".join(
            [
                "taxon\thabitat",
                "A\tforest",
                "B\tforest",
                "C\tdesert",
                "D\t",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = summarize_ancestral_transitions(
        fixture("example_tree.nwk"),
        table_path,
        trait="habitat",
        model="fitch",
    )

    assert [(row.taxon, row.reason) for row in report.exclusions] == [
        ("D", "missing_discrete_trait_state")
    ]


def test_write_ancestral_transition_tables_write_expected_rows(tmp_path: Path) -> None:
    report = summarize_ancestral_transitions(
        fixture("example_tree.nwk"),
        fixture("example_traits_ancestral_sparse.tsv"),
        trait="habitat",
        model="fitch",
    )
    summary_path = tmp_path / "ancestral-transition-summary.tsv"
    branch_path = tmp_path / "ancestral-transition-branches.tsv"
    count_path = tmp_path / "ancestral-transition-counts.tsv"
    exclusion_path = tmp_path / "ancestral-transition-excluded.tsv"

    write_ancestral_transition_summary_table(summary_path, report)
    write_ancestral_transition_branch_table(branch_path, report)
    write_ancestral_transition_count_table(count_path, report)
    write_ancestral_transition_exclusion_table(exclusion_path, report)

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    branch_rows = branch_path.read_text(encoding="utf-8").splitlines()
    count_rows = count_path.read_text(encoding="utf-8").splitlines()
    exclusion_rows = exclusion_path.read_text(encoding="utf-8").splitlines()

    assert summary_rows[0].startswith("trait\ttaxon_column\tmodel\tstate_ordering")
    assert branch_rows[0].startswith(
        "parent_node\tchild_node\tchild_descendant_taxa\tbranch_length"
    )
    assert count_rows[0].startswith(
        "transition\tsource_state\ttarget_state\tcertain_change_count"
    )
    assert exclusion_rows == ["taxon\treason"]
    assert len(summary_rows) == 2
    assert len(branch_rows) == 7
    assert len(count_rows) == 3


def test_summarize_ancestral_transition_tree_set_reports_pair_stability() -> None:
    report = summarize_ancestral_transition_tree_set(
        fixture("example_posterior_tree_set_six_taxa.nwk"),
        fixture("example_traits_clade_summary.tsv"),
        trait="habitat",
        model="equal-rates",
    )
    summary = summarize_ancestral_transition_tree_set_report(report)

    assert report.total_tree_count == 5
    assert report.burnin_tree_count == 0
    assert report.kept_tree_count == 5
    assert report.shared_tree_taxa == ["A", "B", "C", "D", "E", "F"]
    assert report.analysis_taxa == ["A", "B", "C", "D", "E", "F"]
    assert report.rooted_topology_count == 5
    assert report.unrooted_topology_count == 4
    assert len(report.tree_rows) == 5
    assert len(report.branch_rows) == 50
    assert len(report.transition_rows) >= 2
    assert summary.transition_pair_count == len(report.transition_rows)
    assert (
        summary.topology_sensitive_transition_pair_count
        + summary.uncertainty_sensitive_transition_pair_count
        >= 1
    )
    assert (
        "one or more inferred transition pairs depend on uncertain branchwise ancestral changes"
        in report.warnings
    )


def test_write_ancestral_transition_tree_set_tables_write_expected_rows(
    tmp_path: Path,
) -> None:
    report = summarize_ancestral_transition_tree_set(
        fixture("example_posterior_tree_set_six_taxa.nwk"),
        fixture("example_traits_clade_summary.tsv"),
        trait="habitat",
        model="equal-rates",
    )
    summary_path = tmp_path / "ancestral-transition-tree-set-summary.tsv"
    tree_path = tmp_path / "ancestral-transition-tree-set-trees.tsv"
    branch_path = tmp_path / "ancestral-transition-tree-set-branches.tsv"
    count_path = tmp_path / "ancestral-transition-tree-set-counts.tsv"
    exclusion_path = tmp_path / "ancestral-transition-tree-set-excluded.tsv"

    write_ancestral_transition_tree_set_summary_table(summary_path, report)
    write_ancestral_transition_tree_set_tree_table(tree_path, report)
    write_ancestral_transition_tree_set_branch_table(branch_path, report)
    write_ancestral_transition_tree_set_count_table(count_path, report)
    write_ancestral_transition_exclusion_table(exclusion_path, report)

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    tree_rows = tree_path.read_text(encoding="utf-8").splitlines()
    branch_rows = branch_path.read_text(encoding="utf-8").splitlines()
    count_rows = count_path.read_text(encoding="utf-8").splitlines()
    exclusion_rows = exclusion_path.read_text(encoding="utf-8").splitlines()

    assert summary_rows[0].startswith("trait\ttaxon_column\tmodel\tstate_ordering")
    assert tree_rows[0].startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id"
    )
    assert branch_rows[0].startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id\tunrooted_topology_id\tparent_node"
    )
    assert count_rows[0].startswith(
        "transition\tsource_state\ttarget_state\ttree_presence_count"
    )
    assert exclusion_rows == ["taxon\treason"]
    assert len(summary_rows) == 2
    assert len(tree_rows) == 6
    assert len(branch_rows) == 51
