from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.simulation import (
    simulate_coalescent_trees,
    write_coalescent_skyline_table,
)


def test_write_coalescent_skyline_table_emits_interval_ledgers(
    tmp_path: Path,
) -> None:
    _trees, report = simulate_coalescent_trees(
        tree_count=64,
        tip_count=5,
        population_size=2.5,
        waiting_time_tolerance=0.2,
        seed=19,
    )

    skyline_path = write_coalescent_skyline_table(
        tmp_path / "coalescent-skyline.tsv",
        report,
    )

    assert skyline_path.read_text(encoding="utf-8").splitlines() == [
        "interval\tlineage_count\tduration\teffective_population_size_estimate\tobservation_count\trelative_error\tuncertainty_flag",
        "5->4\t5\t0.22986771303003\t2.2986771303003\t64\t0.08052914787988\tlow",
        "4->3\t4\t0.385252631360545\t2.31151578816327\t64\t0.075393684734692\tlow",
        "3->2\t3\t0.876301079716505\t2.62890323914952\t64\t0.051561295659806\tlow",
        "2->1\t2\t3.11097449418924\t3.11097449418924\t64\t0.244389797675697\thigh",
    ]
