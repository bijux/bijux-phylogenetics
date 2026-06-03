from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.trees import load_tree_set, summarize_bootstrap_tree_set

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def _write_large_tree_set(
    path: Path,
    *,
    tree_count: int,
    malformed_record: str | None = None,
) -> Path:
    source_rows = [
        row.strip()
        for row in fixture("example_tree_set_left.nwk")
        .read_text(encoding="utf-8")
        .splitlines()
        if row.strip()
    ]
    rows = [source_rows[index % len(source_rows)] for index in range(tree_count)]
    if malformed_record is not None:
        rows.insert(tree_count // 2, malformed_record)
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def test_load_tree_set_skips_malformed_newick_records_and_reports_processing(
    tmp_path: Path,
) -> None:
    tree_set_path = _write_large_tree_set(
        tmp_path / "mixed-tree-set.nwk",
        tree_count=6,
        malformed_record="((A:0.1,B:0.1),(C:0.1,D:0.1);",
    )

    report = load_tree_set(tree_set_path)

    assert report.tree_count == 6
    assert report.processing.skipped_malformed_tree_count == 1
    assert report.processing.runtime_seconds >= 0.0
    assert report.processing.peak_memory_bytes >= 0
    assert report.rooted_topology_count == 2


def test_load_tree_set_reads_newick_records_from_trees_suffix(tmp_path: Path) -> None:
    tree_set_path = _write_large_tree_set(
        tmp_path / "posterior.trees",
        tree_count=3,
    )

    report = load_tree_set(tree_set_path)

    assert report.tree_count == 3
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.processing.skipped_malformed_tree_count == 0
    assert report.rooted_topology_count == 2


def test_summarize_bootstrap_tree_set_handles_thousand_tree_inputs(
    tmp_path: Path,
) -> None:
    tree_set_path = _write_large_tree_set(
        tmp_path / "large-tree-set.nwk",
        tree_count=1024,
    )

    report = summarize_bootstrap_tree_set(tree_set_path)

    assert report.tree_count == 1024
    assert report.processing.skipped_malformed_tree_count == 0
    assert report.processing.runtime_seconds >= 0.0
    assert report.processing.peak_memory_bytes >= 0
    assert report.diversity.pair_count == (1024 * 1023) // 2
    assert sum(row.pair_count for row in report.diversity.rf_distribution) == (
        (1024 * 1023) // 2
    )
    assert any(
        row.robinson_foulds_distance == 0 for row in report.diversity.rf_distribution
    )
