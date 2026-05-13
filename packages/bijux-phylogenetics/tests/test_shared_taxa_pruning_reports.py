from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.topology import (
    write_shared_taxa_pruning_table,
    write_shared_taxa_removed_taxa_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_write_shared_taxa_pruning_table_writes_branch_audit_rows(
    tmp_path: Path,
) -> None:
    output = tmp_path / "shared-taxa-pruning.tsv"

    write_shared_taxa_pruning_table(
        output,
        fixture("example_tree.nwk"),
        fixture("example_tree_overlap.nwk"),
    )

    assert output.read_text(encoding="utf-8") == (
        "tree_side\ttree_path\toriginal_tip_count\tretained_tip_count\tremoved_tip_count\trequested_taxa\tkept_taxa\tremoved_taxa\tabsent_requested_taxa\tremoved_taxa_with_reasons\ttransformation\troot_to_tip_complete\tmin_root_to_tip\tmax_root_to_tip\tunary_internal_nodes\toriginal_total_branch_length\tpruned_total_branch_length\tbranch_length_delta\tlost_taxa_count\tlost_taxa_fraction\tlost_clade_count\tlost_clade_fraction\tlost_branch_length\tlost_branch_length_fraction\n"
        f"left\t{fixture('example_tree.nwk')}\t4\t3\t1\tA|B|C\tA|B|C\tD\t\tD:not_requested\tprune-tree-to-requested-taxa\ttrue\t0.30000000000000004\t0.30000000000000004\t\t0.9\t0.7\t-0.2\t1\t0.25\t1\t0.5\t0.2\t0.22222222222222224\n"
        f"right\t{fixture('example_tree_overlap.nwk')}\t4\t3\t1\tA|B|C\tA|B|C\tE\t\tE:not_requested\tprune-tree-to-requested-taxa\ttrue\t0.30000000000000004\t0.30000000000000004\t\t0.9\t0.7\t-0.2\t1\t0.25\t1\t0.5\t0.2\t0.22222222222222224\n"
    )


def test_write_shared_taxa_removed_taxa_table_writes_one_row_per_taxon(
    tmp_path: Path,
) -> None:
    output = tmp_path / "shared-taxa-removed.tsv"

    write_shared_taxa_removed_taxa_table(
        output,
        fixture("example_tree.nwk"),
        fixture("example_tree_overlap.nwk"),
    )

    assert output.read_text(encoding="utf-8") == (
        "tree_side\ttree_path\ttaxon\treason\n"
        f"left\t{fixture('example_tree.nwk')}\tD\tnot_requested\n"
        f"right\t{fixture('example_tree_overlap.nwk')}\tE\tnot_requested\n"
    )
