from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.phylo.topology.clades import (
    informative_rooted_clades,
    informative_unrooted_splits,
)

from .contracts import TreeSetSplitFrequencyReport, TreeSetSplitFrequencyRow
from .inventory import _analyze_tree_set, _require_exact_taxa
from .topology import _format_clade

TREE_SET_SPLIT_FREQUENCY_POLICIES = ("rooted", "unrooted")


def compute_tree_set_split_frequency_table(
    path: Path,
    *,
    split_policy: str,
) -> TreeSetSplitFrequencyReport:
    """Compute split frequencies across one tree set under one explicit policy."""
    analysis = _analyze_tree_set(path)
    exact_taxa = _require_exact_taxa(analysis)
    validated_split_policy = _validate_split_policy(split_policy)
    shared_taxa = set(exact_taxa)
    split_counts: dict[frozenset[str], int] = {}
    for tree in analysis.trees:
        if validated_split_policy == "rooted":
            informative_splits = informative_rooted_clades(tree, shared_taxa)
        else:
            informative_splits = informative_unrooted_splits(tree, shared_taxa)
        for split in informative_splits:
            split_counts[split] = split_counts.get(split, 0) + 1
    tree_count = len(analysis.trees)
    return TreeSetSplitFrequencyReport(
        path=path,
        tree_count=tree_count,
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        split_policy=validated_split_policy,
        split_frequencies=[
            TreeSetSplitFrequencyRow(
                split=_format_clade(split),
                tree_count=count,
                frequency=round(count / tree_count, 15),
            )
            for split, count in sorted(
                split_counts.items(),
                key=lambda item: _format_clade(item[0]),
            )
        ],
    )


def write_tree_set_split_frequency_table(
    path: Path,
    report: TreeSetSplitFrequencyReport,
) -> Path:
    """Write one tree-set split-frequency table as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["split", "tree_count", "frequency", "split_policy"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.split_frequencies:
            writer.writerow(
                {
                    "split": row.split,
                    "tree_count": row.tree_count,
                    "frequency": format(row.frequency, ".15g"),
                    "split_policy": report.split_policy,
                }
            )
    return path


def _validate_split_policy(split_policy: str) -> str:
    if split_policy not in TREE_SET_SPLIT_FREQUENCY_POLICIES:
        raise ValueError(
            "split_policy must be one of "
            f"{TREE_SET_SPLIT_FREQUENCY_POLICIES}, got {split_policy!r}"
        )
    return split_policy
