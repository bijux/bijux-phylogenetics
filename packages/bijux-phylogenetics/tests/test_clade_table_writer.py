from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.trees import (
    extract_tree_clades,
    extract_tree_set_clades,
    write_clade_table,
)


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def metadata_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "metadata" / name


def test_write_clade_table_writes_tree_rows_with_metadata_columns(
    tmp_path: Path,
) -> None:
    output = tmp_path / "clades.tsv"
    report = extract_tree_clades(
        tree_fixture("example_tree_support_conflict_left.nwk"),
        metadata_path=metadata_fixture("example_metadata.tsv"),
        metadata_columns=["species"],
    )

    write_clade_table(output, report)

    assert output.read_text(encoding="utf-8").splitlines()[0] == (
        "source_path\ttree_index\tnode_kind\tclade_id\tnode_label\ttaxon_count\t"
        "taxa\tsupport\tsupport_fraction\tbranch_length\troot_depth\t"
        "descendant_tip_depth_min\tdescendant_tip_depth_max\tnode_age\t"
        "species_values\tspecies_distinct_values\tspecies_missing_taxa"
    )
    assert (
        f"{tree_fixture('example_tree_support_conflict_left.nwk')}\t\tinternal\tA|B\t95\t2\tA|B\t95.0\t0.95\t0.2\t0.2\t0.1\t0.1\t0.1\tA=Alpha species|B=Beta species\tAlpha species|Beta species\t"
        in output.read_text(encoding="utf-8")
    )


def test_write_clade_table_writes_tree_set_rows_with_tree_indices(
    tmp_path: Path,
) -> None:
    output = tmp_path / "tree-set-clades.tsv"
    report = extract_tree_set_clades(tree_fixture("example_tree_set_left.nwk"))

    write_clade_table(output, report)

    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[0] == (
        "source_path\ttree_index\tnode_kind\tclade_id\tnode_label\ttaxon_count\t"
        "taxa\tsupport\tsupport_fraction\tbranch_length\troot_depth\t"
        "descendant_tip_depth_min\tdescendant_tip_depth_max\tnode_age"
    )
    assert any(
        line.startswith(
            f"{tree_fixture('example_tree_set_left.nwk')}\t3\tinternal\tA|C\t\t2\tA|C\t\t\t0.2\t0.2\t0.1\t0.1\t0.1"
        )
        for line in lines[1:]
    )
