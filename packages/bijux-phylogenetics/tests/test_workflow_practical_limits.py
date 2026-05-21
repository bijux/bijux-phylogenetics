from __future__ import annotations

import json

import pytest

import bijux_phylogenetics.benchmark as benchmark_api
from bijux_phylogenetics.benchmark import benchmark_workflow_practical_limits
from bijux_phylogenetics.command_line import main


def _entries_by_workflow(report):
    return {entry.workflow: entry for entry in report.entries}


@pytest.mark.slow
def test_benchmark_workflow_practical_limits_reports_governed_limits() -> None:
    report = benchmark_workflow_practical_limits(
        replicates=1,
        tree_tip_counts=[8, 16],
        alignment_size_classes=[
            ("sequences-4-sites-16", 4, 16),
            ("sequences-6-sites-24", 6, 24),
        ],
        tree_set_size_classes=[
            ("trees-8-taxa-6", 8, 6),
            ("trees-12-taxa-8", 12, 8),
        ],
        stress_tiers=["small"],
    )

    assert report.replicates == 1
    assert report.stress_tiers == ["small"]
    assert report.limitations
    entries = _entries_by_workflow(report)
    assert set(entries) >= {
        "tree-validation",
        "alignment-diagnostics",
        "tree-set-consensus",
        "large-alignment-inference",
        "posterior-tree-set-consensus",
    }

    tree_validation = entries["tree-validation"]
    assert tree_validation.evidence_source == "large-tree-scaling"
    assert tree_validation.tested_taxon_limit == 16
    assert tree_validation.tested_tree_limit == 1

    alignment = entries["alignment-diagnostics"]
    assert alignment.tested_taxon_limit == 6
    assert alignment.tested_site_limit == 24

    tree_set = entries["tree-set-consensus"]
    assert tree_set.tested_taxon_limit == 8
    assert tree_set.tested_tree_limit == 12
    assert tree_set.tested_posterior_size == 12

    stress_tree_set = entries["posterior-tree-set-consensus"]
    assert stress_tree_set.evidence_source == "stress-suite"
    assert stress_tree_set.tested_taxon_limit == 64
    assert stress_tree_set.tested_tree_limit == 256
    assert stress_tree_set.tested_posterior_size == 256
    assert any("small" in note for note in stress_tree_set.notes)


def test_benchmark_workflow_practical_limits_rejects_empty_stress_tiers() -> None:
    with pytest.raises(ValueError, match="at least one governed tier"):
        benchmark_workflow_practical_limits(
            replicates=1,
            tree_tip_counts=[8],
            alignment_size_classes=[("sequences-4-sites-16", 4, 16)],
            tree_set_size_classes=[("trees-8-taxa-6", 8, 6)],
            stress_tiers=[],
        )


@pytest.mark.slow
def test_cli_benchmark_workflow_practical_limits_reports_entry_metrics(
    capsys: pytest.CaptureFixture[str],
) -> None:
    exit_code = main(
        [
            "benchmark",
            "workflow-practical-limits",
            "--replicates",
            "1",
            "--tree-tip-count",
            "8",
            "--tree-tip-count",
            "16",
            "--sequence-count",
            "4",
            "--alignment-length",
            "16",
            "--sequence-count",
            "6",
            "--alignment-length",
            "24",
            "--posterior-tree-count",
            "8",
            "--tree-set-tip-count",
            "6",
            "--posterior-tree-count",
            "12",
            "--tree-set-tip-count",
            "8",
            "--stress-tier",
            "small",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["replicates"] == 1
    assert payload["metrics"]["entry_count"] == 17
    assert payload["metrics"]["workflow_count"] == 17
    assert payload["metrics"]["stress_tier_count"] == 1
    assert payload["metrics"]["max_taxon_limit"] == 256
    assert payload["metrics"]["max_site_limit"] == 704
    assert payload["metrics"]["max_tree_limit"] == 256
    assert payload["metrics"]["max_posterior_size"] == 256


def test_public_runtime_exports_workflow_practical_limit_surface() -> None:
    assert (
        benchmark_api.benchmark_workflow_practical_limits
        is benchmark_workflow_practical_limits
    )
