from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.trees import (
    summarize_bootstrap_tree_set,
    write_bootstrap_tree_set_artifacts,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_summarize_bootstrap_tree_set_reports_consensus_diversity_and_unstable_branches() -> (
    None
):
    report = summarize_bootstrap_tree_set(fixture("example_tree_set_left.nwk"))

    assert report.tree_count == 3
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.consensus.consensus_newick == (
        "((A:0.1,B:0.1)66.6666666666667:0.2,(C:0.1,D:0.1)66.6666666666667:0.2);"
    )
    assert report.diversity.rooted_topology_count == 2
    assert report.diversity.dominant_topology_frequency == 0.666666666666667
    assert [
        (
            row.clade,
            row.bootstrap_tree_count,
            row.bootstrap_frequency,
            row.conflict_count,
            row.support_classification,
        )
        for row in report.unstable_branches
    ] == [
        ("A|B", 2, 0.666666666666667, 2, "intermediate-support"),
        ("C|D", 2, 0.666666666666667, 2, "intermediate-support"),
    ]
    assert report.warnings == [
        "bootstrap replicate trees contain multiple rooted topologies",
        "consensus tree contains branches below the robust bootstrap threshold or with conflicting alternatives",
    ]


def test_summarize_bootstrap_tree_set_respects_robust_support_threshold() -> None:
    report = summarize_bootstrap_tree_set(
        fixture("example_tree_set_left.nwk"),
        robust_support_threshold=0.6,
    )

    assert report.unstable_branch_count == 2
    assert {row.clade for row in report.unstable_branches} == {"A|B", "C|D"}


def test_write_bootstrap_tree_set_artifacts_writes_governed_tables_and_consensus(
    tmp_path: Path,
) -> None:
    report = write_bootstrap_tree_set_artifacts(
        fixture("example_tree_set_left.nwk"),
        out_dir=tmp_path,
        prefix="bootstrap-review",
    )

    assert sorted(report.output_paths) == [
        "clade_frequencies",
        "consensus_tree",
        "distance_matrix",
        "rf_distribution",
        "summary_table",
        "topology_clusters",
        "unstable_branches",
        "unstable_clades",
    ]
    assert (
        report.output_paths["summary_table"]
        .read_text(encoding="utf-8")
        .startswith(
            "tree_count\truntime_seconds\tpeak_memory_bytes\tskipped_malformed_tree_count\tshared_taxon_count\trooted_topology_count\t"
        )
    )
    assert report.output_paths["consensus_tree"].read_text(
        encoding="utf-8"
    ).strip() == (
        "((A:0.1,B:0.1)66.6666666666667:0.2,(C:0.1,D:0.1)66.6666666666667:0.2);"
    )
    assert (
        report.output_paths["rf_distribution"]
        .read_text(encoding="utf-8")
        .startswith(
            "robinson_foulds_distance\tnormalized_robinson_foulds\tpair_count\tfrequency\n"
        )
    )
    assert (
        report.output_paths["unstable_branches"]
        .read_text(encoding="utf-8")
        .startswith(
            "clade\tbootstrap_tree_count\tbootstrap_frequency\tbootstrap_support_percent\t"
        )
    )
