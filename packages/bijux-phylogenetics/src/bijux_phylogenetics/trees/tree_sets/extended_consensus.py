from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick, loads_newick, write_newick
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .consensus import _build_consensus_tree_from_selected_clades, _mean
from .contracts import (
    MajorityRuleExtendedAcceptedCladeRow,
    MajorityRuleExtendedConsensusReport,
    MajorityRuleExtendedRejectedCladeRow,
)
from .inventory import _analyze_tree_set, _require_exact_taxa
from .topology import _clades_conflict, _format_clade

_MAJORITY_THRESHOLD = 0.5


def _candidate_sort_key(item: tuple[frozenset[str], int]) -> tuple[int, str]:
    clade, count = item
    return (-count, _format_clade(clade))


def _select_extended_consensus_clades(
    counts: dict[frozenset[str], int],
    *,
    tree_count: int,
) -> tuple[
    list[MajorityRuleExtendedAcceptedCladeRow],
    list[MajorityRuleExtendedRejectedCladeRow],
]:
    majority_candidates = sorted(
        (
            (clade, count)
            for clade, count in counts.items()
            if count / tree_count >= _MAJORITY_THRESHOLD
        ),
        key=_candidate_sort_key,
    )
    accepted_rows: list[MajorityRuleExtendedAcceptedCladeRow] = []
    rejected_rows: list[MajorityRuleExtendedRejectedCladeRow] = []
    included_clades: list[frozenset[str]] = []

    for insertion_rank, (clade, count) in enumerate(majority_candidates, start=1):
        included_clades.append(clade)
        accepted_rows.append(
            MajorityRuleExtendedAcceptedCladeRow(
                insertion_rank=insertion_rank,
                clade=_format_clade(clade),
                tree_count=count,
                frequency=round(count / tree_count, 15),
                inclusion_stage="majority",
            )
        )

    insertion_rank = len(accepted_rows)
    remaining_candidates = sorted(
        (
            (clade, count)
            for clade, count in counts.items()
            if count / tree_count < _MAJORITY_THRESHOLD
        ),
        key=_candidate_sort_key,
    )
    for clade, count in remaining_candidates:
        blocking_clades = sorted(
            _format_clade(included_clade)
            for included_clade in included_clades
            if _clades_conflict(clade, included_clade)
        )
        if blocking_clades:
            rejected_rows.append(
                MajorityRuleExtendedRejectedCladeRow(
                    clade=_format_clade(clade),
                    tree_count=count,
                    frequency=round(count / tree_count, 15),
                    blocking_clades=blocking_clades,
                )
            )
            continue
        insertion_rank += 1
        included_clades.append(clade)
        accepted_rows.append(
            MajorityRuleExtendedAcceptedCladeRow(
                insertion_rank=insertion_rank,
                clade=_format_clade(clade),
                tree_count=count,
                frequency=round(count / tree_count, 15),
                inclusion_stage="compatible-extension",
            )
        )
    return accepted_rows, rejected_rows


def _build_majority_rule_extended_consensus_report(
    path: Path,
) -> tuple[PhyloTree, MajorityRuleExtendedConsensusReport]:
    analysis = _analyze_tree_set(path)
    shared_taxa = _require_exact_taxa(analysis)
    counts = analysis.clade_counts or {}
    accepted_rows, rejected_rows = _select_extended_consensus_clades(
        counts,
        tree_count=len(analysis.trees),
    )
    included_clades = {
        frozenset(row.clade.split("|")) for row in accepted_rows if row.clade
    }
    clade_support = {
        clade: round((counts[clade] / len(analysis.trees)) * 100.0, 15)
        for clade in included_clades
    }
    clade_lengths = {
        clade: _mean(lengths)
        for clade, lengths in analysis.clade_branch_lengths.items()
        if clade in included_clades and lengths
    }
    terminal_lengths = {
        taxon: _mean(lengths)
        for taxon, lengths in analysis.terminal_lengths.items()
        if lengths
    }
    tree = _build_consensus_tree_from_selected_clades(
        shared_taxa=shared_taxa,
        source_format=analysis.source_format,
        included_clades=included_clades,
        clade_support=clade_support,
        clade_lengths=clade_lengths,
        terminal_lengths=terminal_lengths,
    )
    majority_included_clade_count = sum(
        1 for row in accepted_rows if row.inclusion_stage == "majority"
    )
    extension_included_clade_count = sum(
        1 for row in accepted_rows if row.inclusion_stage == "compatible-extension"
    )
    report = MajorityRuleExtendedConsensusReport(
        path=analysis.path,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=shared_taxa,
        consensus_method="majority-rule-extended",
        majority_threshold=_MAJORITY_THRESHOLD,
        included_clade_count=len(accepted_rows),
        majority_included_clade_count=majority_included_clade_count,
        extension_included_clade_count=extension_included_clade_count,
        rejected_conflict_count=len(rejected_rows),
        consensus_newick=dumps_newick(tree),
        accepted_clades=accepted_rows,
        rejected_clades=rejected_rows,
    )
    return tree, report


def compute_majority_rule_extended_consensus(
    path: Path,
) -> tuple[PhyloTree, MajorityRuleExtendedConsensusReport]:
    """Compute a majority-rule extended consensus tree from one exact-taxa tree set."""
    return _build_majority_rule_extended_consensus_report(path)


def write_majority_rule_extended_consensus_inclusion_table(
    path: Path,
    report: MajorityRuleExtendedConsensusReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "insertion_rank",
                "clade",
                "tree_count",
                "frequency",
                "inclusion_stage",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.accepted_clades:
            writer.writerow(
                {
                    "insertion_rank": row.insertion_rank,
                    "clade": row.clade,
                    "tree_count": row.tree_count,
                    "frequency": format(row.frequency, ".15g"),
                    "inclusion_stage": row.inclusion_stage,
                }
            )
    return path


def write_majority_rule_extended_consensus_rejected_conflict_table(
    path: Path,
    report: MajorityRuleExtendedConsensusReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "clade",
                "tree_count",
                "frequency",
                "blocking_clades",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.rejected_clades:
            writer.writerow(
                {
                    "clade": row.clade,
                    "tree_count": row.tree_count,
                    "frequency": format(row.frequency, ".15g"),
                    "blocking_clades": "||".join(row.blocking_clades),
                }
            )
    return path


def write_majority_rule_extended_consensus_artifacts(
    out_dir: Path,
    report: MajorityRuleExtendedConsensusReport,
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    consensus_tree_path = write_newick(
        out_dir / "consensus_tree.nwk",
        loads_newick(report.consensus_newick),
    )
    inclusion_table_path = write_majority_rule_extended_consensus_inclusion_table(
        out_dir / "inclusion_order.tsv",
        report,
    )
    rejected_conflicts_path = (
        write_majority_rule_extended_consensus_rejected_conflict_table(
            out_dir / "rejected_conflicts.tsv",
            report,
        )
    )
    return {
        "consensus_tree_path": consensus_tree_path,
        "inclusion_table_path": inclusion_table_path,
        "rejected_conflicts_path": rejected_conflicts_path,
    }
