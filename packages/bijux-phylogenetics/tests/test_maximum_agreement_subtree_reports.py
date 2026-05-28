from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare import (
    write_maximum_agreement_subtree_pruning_table,
    write_maximum_agreement_subtree_removed_taxa_table,
    write_maximum_agreement_subtree_search_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_write_maximum_agreement_subtree_removed_taxa_table_marks_conflict_taxa(
    tmp_path: Path,
) -> None:
    output = tmp_path / "maximum-agreement-subtree-removed.tsv"

    write_maximum_agreement_subtree_removed_taxa_table(
        output,
        fixture("agreement_subtree_left.nwk"),
        fixture("agreement_subtree_right.nwk"),
        max_evaluated_candidate_count=13,
    )

    assert output.read_text(encoding="utf-8") == (
        "tree_side\ttree_path\ttaxon\treason\tshared_taxon\tremoved_for_maximum_agreement_subtree\n"
        f"left\t{fixture('agreement_subtree_left.nwk')}\tC\tnot_requested\ttrue\ttrue\n"
        f"right\t{fixture('agreement_subtree_right.nwk')}\tC\tnot_requested\ttrue\ttrue\n"
    )


def test_write_maximum_agreement_subtree_search_table_writes_greedy_trace(
    tmp_path: Path,
) -> None:
    output = tmp_path / "maximum-agreement-subtree-search.tsv"

    write_maximum_agreement_subtree_search_table(
        output,
        fixture("agreement_subtree_left.nwk"),
        fixture("agreement_subtree_right.nwk"),
        max_evaluated_candidate_count=13,
    )

    assert output.read_text(encoding="utf-8") == (
        "evaluation_index\tstep_index\tretained_taxon_count\tretained_taxa\tremoved_taxa\trobinson_foulds_distance\tnormalized_robinson_foulds\ttopology_equal\tselected_for_next_step\n"
        "1\t0\t5\tA|B|C|D|E\t\t2\t0.3333333333333333\tfalse\ttrue\n"
        "2\t1\t4\tB|C|D|E\tA\t2\t0.5\tfalse\tfalse\n"
        "3\t1\t4\tA|C|D|E\tB\t2\t0.5\tfalse\tfalse\n"
        "4\t1\t4\tA|B|D|E\tC\t0\t0.0\ttrue\ttrue\n"
        "5\t1\t4\tA|B|C|E\tD\t2\t0.5\tfalse\tfalse\n"
        "6\t1\t4\tA|B|C|D\tE\t2\t0.5\tfalse\tfalse\n"
    )


def test_write_maximum_agreement_subtree_pruning_table_includes_budget_summary(
    tmp_path: Path,
) -> None:
    output = tmp_path / "maximum-agreement-subtree-pruning.tsv"

    write_maximum_agreement_subtree_pruning_table(
        output,
        fixture("agreement_subtree_left.nwk"),
        fixture("agreement_subtree_right.nwk"),
        max_evaluated_candidate_count=13,
    )

    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("tree_side\ttree_path\trf_mode\tsearch_strategy\t")
    assert "greedy-single-taxon-removal" in lines[1]
    assert "minimize-robinson-foulds-then-normalized-distance" in lines[1]
    assert "heuristic-solution-not-guaranteed-optimal" in lines[1]
    assert "\t26\t13\t6\tA|B|D|E\tC\t" in lines[1]
    assert "\t26\t13\t6\tA|B|D|E\tC\t" in lines[2]
