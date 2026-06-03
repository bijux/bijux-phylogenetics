from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.trees import (
    compute_clade_compatibility_graph,
    write_clade_compatibility_edge_table,
    write_clade_compatibility_graph_dot,
    write_clade_compatibility_node_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_write_clade_compatibility_node_table_writes_expected_columns(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "clade-compatibility-nodes.tsv"
    report = compute_clade_compatibility_graph(
        fixture("clade_compatibility_tree_set.nwk")
    )

    write_clade_compatibility_node_table(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()

    assert lines == [
        "clade\ttree_count\tfrequency\tcompatible_neighbor_count\tconflict_neighbor_count",
        "A|B\t1\t0.5\t1\t2",
        "A|C\t1\t0.5\t1\t2",
        "B|D\t1\t0.5\t1\t2",
        "C|D\t1\t0.5\t1\t2",
    ]


def test_write_clade_compatibility_edge_table_writes_expected_columns(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "clade-compatibility-edges.tsv"
    report = compute_clade_compatibility_graph(
        fixture("clade_compatibility_tree_set.nwk")
    )

    write_clade_compatibility_edge_table(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()

    assert lines == [
        (
            "left_clade\tright_clade\tcompatibility_relation\tcompatibility_reason\t"
            "left_tree_count\tright_tree_count\tleft_frequency\tright_frequency"
        ),
        "A|B\tA|C\tconflict\toverlap-without-containment\t1\t1\t0.5\t0.5",
        "A|B\tB|D\tconflict\toverlap-without-containment\t1\t1\t0.5\t0.5",
        "A|B\tC|D\tcompatible\tdisjoint\t1\t1\t0.5\t0.5",
        "A|C\tB|D\tcompatible\tdisjoint\t1\t1\t0.5\t0.5",
        "A|C\tC|D\tconflict\toverlap-without-containment\t1\t1\t0.5\t0.5",
        "B|D\tC|D\tconflict\toverlap-without-containment\t1\t1\t0.5\t0.5",
    ]


def test_write_clade_compatibility_graph_dot_writes_expected_relations(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "clade-compatibility.dot"
    report = compute_clade_compatibility_graph(
        fixture("clade_compatibility_tree_set.nwk")
    )

    write_clade_compatibility_graph_dot(output_path, report)

    dot_text = output_path.read_text(encoding="utf-8")

    assert dot_text.startswith("graph clade_compatibility {\n")
    assert 'label="A|B\\nfrequency=0.5\\ncompatible=1\\nconflict=2"' in dot_text
    assert 'label="C|D\\nfrequency=0.5\\ncompatible=1\\nconflict=2"' in dot_text
    assert 'color="darkgreen", style="solid", label="disjoint"' in dot_text
    assert (
        'color="firebrick", style="dashed", label="overlap-without-containment"'
        in dot_text
    )
