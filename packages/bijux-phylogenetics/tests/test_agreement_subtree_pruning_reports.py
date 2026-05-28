from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare import (
    write_agreement_subtree_pruning_table,
    write_agreement_subtree_removed_taxa_table,
    write_agreement_subtree_search_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_write_agreement_subtree_removed_taxa_table_marks_conflict_taxa(
    tmp_path: Path,
) -> None:
    output = tmp_path / "agreement-subtree-removed.tsv"

    write_agreement_subtree_removed_taxa_table(
        output,
        fixture("agreement_subtree_left.nwk"),
        fixture("agreement_subtree_right.nwk"),
    )

    assert output.read_text(encoding="utf-8") == (
        "tree_side\ttree_path\ttaxon\treason\tshared_taxon\tremoved_for_agreement_subtree\n"
        f"left\t{fixture('agreement_subtree_left.nwk')}\tC\tnot_requested\ttrue\ttrue\n"
        f"right\t{fixture('agreement_subtree_right.nwk')}\tC\tnot_requested\ttrue\ttrue\n"
    )


def test_write_agreement_subtree_search_table_writes_candidate_trace(
    tmp_path: Path,
) -> None:
    output = tmp_path / "agreement-subtree-search.tsv"

    write_agreement_subtree_search_table(
        output,
        fixture("agreement_subtree_left.nwk"),
        fixture("agreement_subtree_right.nwk"),
    )

    assert output.read_text(encoding="utf-8") == (
        "candidate_index\tretained_taxon_count\tretained_taxa\tremoved_taxa\trobinson_foulds_distance\tnormalized_robinson_foulds\ttopology_equal\n"
        "1\t5\tA|B|C|D|E\t\t2\t0.3333333333333333\tfalse\n"
        "2\t4\tA|B|C|D\tE\t2\t0.5\tfalse\n"
        "3\t4\tA|B|C|E\tD\t2\t0.5\tfalse\n"
        "4\t4\tA|B|D|E\tC\t0\t0.0\ttrue\n"
    )


def test_write_agreement_subtree_pruning_table_includes_search_summary(
    tmp_path: Path,
) -> None:
    output = tmp_path / "agreement-subtree-pruning.tsv"

    write_agreement_subtree_pruning_table(
        output,
        fixture("agreement_subtree_left.nwk"),
        fixture("agreement_subtree_right.nwk"),
    )

    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("tree_side\ttree_path\trf_mode\tsearch_strategy\t")
    assert "exact-descending-retained-subsets" in lines[1]
    assert "\t26\t4\tA|B|D|E\tC\t" in lines[1]
    assert "\t26\t4\tA|B|D|E\tC\t" in lines[2]
