from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from bijux_phylogenetics.trees import (
    summarize_posterior_agreement_subtree,
    write_posterior_agreement_subtree_artifacts,
    write_posterior_agreement_subtree_removed_taxa_table,
    write_posterior_agreement_subtree_search_table,
    write_posterior_agreement_subtree_summary_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_summarize_posterior_agreement_subtree_finds_stable_retained_topology() -> None:
    tree, report = summarize_posterior_agreement_subtree(
        fixture("posterior_agreement_subtree_tree_set.nwk")
    )

    assert dumps_newick(tree) == "((A:1,B:1):2,(D:1,E:1):2);"
    assert report.tree_count == 3
    assert report.shared_taxa == ["A", "B", "C", "D", "E"]
    assert report.search_strategy == "exact-descending-retained-subsets"
    assert report.possible_retained_subset_count == 26
    assert report.evaluated_candidate_count == 4
    assert report.retained_taxa == ["A", "B", "D", "E"]
    assert report.agreement_removed_taxa == ["C"]
    assert report.stable_rooted_topology_id == "A|B||D|E"
    assert report.agreement_subtree_newick == "((A:1,B:1):2,(D:1,E:1):2);"
    assert [
        (
            row.candidate_index,
            row.retained_taxa,
            row.rooted_topology_count,
            row.dominant_rooted_topology_frequency,
            row.stable_topology_reached,
        )
        for row in report.candidate_rows
    ] == [
        (1, ["A", "B", "C", "D", "E"], 2, 0.666666666666667, False),
        (2, ["A", "B", "C", "D"], 2, 0.666666666666667, False),
        (3, ["A", "B", "C", "E"], 2, 0.666666666666667, False),
        (4, ["A", "B", "D", "E"], 1, 1.0, True),
    ]


def test_write_posterior_agreement_subtree_tables_write_expected_rows(
    tmp_path: Path,
) -> None:
    _tree, report = summarize_posterior_agreement_subtree(
        fixture("posterior_agreement_subtree_tree_set.nwk")
    )

    summary_path = tmp_path / "posterior-agreement-subtree-summary.tsv"
    removed_path = tmp_path / "posterior-agreement-subtree-removed.tsv"
    search_path = tmp_path / "posterior-agreement-subtree-search.tsv"
    write_posterior_agreement_subtree_summary_table(summary_path, report)
    write_posterior_agreement_subtree_removed_taxa_table(removed_path, report)
    write_posterior_agreement_subtree_search_table(search_path, report)

    assert summary_path.read_text(encoding="utf-8").splitlines() == [
        "tree_count\tsearch_strategy\tpossible_retained_subset_count\tevaluated_candidate_count\tretained_taxa\tagreement_removed_taxa\tstable_rooted_topology_id\tagreement_subtree_newick",
        "3\texact-descending-retained-subsets\t26\t4\tA|B|D|E\tC\tA|B||D|E\t((A:1,B:1):2,(D:1,E:1):2);",
    ]
    assert removed_path.read_text(encoding="utf-8").splitlines() == [
        "taxon\tremoved_for_agreement_subtree",
        "A\tfalse",
        "B\tfalse",
        "C\ttrue",
        "D\tfalse",
        "E\tfalse",
    ]
    assert search_path.read_text(encoding="utf-8").splitlines() == [
        "candidate_index\tretained_taxon_count\tretained_taxa\tremoved_taxa\trooted_topology_count\tdominant_rooted_topology_frequency\tstable_topology_reached",
        "1\t5\tA|B|C|D|E\t\t2\t0.666666666666667\tfalse",
        "2\t4\tA|B|C|D\tE\t2\t0.666666666666667\tfalse",
        "3\t4\tA|B|C|E\tD\t2\t0.666666666666667\tfalse",
        "4\t4\tA|B|D|E\tC\t1\t1\ttrue",
    ]


def test_write_posterior_agreement_subtree_artifacts_writes_bundle(
    tmp_path: Path,
) -> None:
    tree, report = summarize_posterior_agreement_subtree(
        fixture("posterior_agreement_subtree_tree_set.nwk")
    )

    paths = write_posterior_agreement_subtree_artifacts(tmp_path, tree, report)

    assert sorted(paths) == [
        "posterior_agreement_subtree_path",
        "posterior_agreement_subtree_removed_taxa_path",
        "posterior_agreement_subtree_search_path",
        "posterior_agreement_subtree_summary_path",
    ]
    assert paths["posterior_agreement_subtree_path"].read_text(encoding="utf-8") == (
        "((A:1,B:1):2,(D:1,E:1):2);\n"
    )


def test_summarize_posterior_agreement_subtree_requires_exact_taxa() -> None:
    with pytest.raises(InvalidAlignmentError, match="exact same taxon set"):
        summarize_posterior_agreement_subtree(
            fixture("example_tree_set_mismatched.nwk")
        )
