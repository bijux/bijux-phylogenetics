from __future__ import annotations

import csv
from itertools import combinations
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import prune_tree_object_to_requested_taxa
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_bipartition,
    informative_unrooted_splits,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

from .contracts import TreeSetQuartetSupportReport, TreeSetQuartetSupportRow
from .inventory import _analyze_tree_set, _require_exact_taxa, _TreeSetAnalysis


def _format_split_id(split: frozenset[str], taxon_scope: set[str]) -> str:
    left_taxa = sorted(split)
    right_taxa = sorted(taxon_scope - split)
    return f"{'|'.join(left_taxa)}::{'|'.join(right_taxa)}"


def _resolve_induced_quartet_split(
    tree: PhyloTree,
    quartet_taxa: tuple[str, str, str, str],
) -> frozenset[str] | None:
    quartet_scope = set(quartet_taxa)
    induced_tree = prune_tree_object_to_requested_taxa(tree, list(quartet_taxa))
    splits = informative_unrooted_splits(induced_tree, quartet_scope)
    if len(splits) != 1:
        return None
    return next(iter(splits))


def _build_reference_tree_quartet_support_report(
    *,
    reference_tree_path: Path,
    reference_tree: PhyloTree,
    analysis: _TreeSetAnalysis,
) -> TreeSetQuartetSupportReport:
    exact_taxa = _require_exact_taxa(analysis)
    taxon_scope = set(exact_taxa)
    quartet_split_cache: list[
        dict[tuple[str, str, str, str], frozenset[str] | None]
    ] = [{} for _ in analysis.trees]
    rows: list[TreeSetQuartetSupportRow] = []
    total_concordant = 0
    total_discordant = 0
    total_uninformative = 0
    total_quartets = 0
    for split in sorted(
        informative_unrooted_splits(reference_tree, taxon_scope),
        key=lambda signature: (len(signature), tuple(sorted(signature))),
    ):
        left_taxa = sorted(split)
        right_taxa = sorted(taxon_scope - split)
        concordant = 0
        discordant = 0
        uninformative = 0
        quartet_count_per_tree = 0
        for left_pair in combinations(left_taxa, 2):
            for right_pair in combinations(right_taxa, 2):
                quartet_count_per_tree += 1
                quartet_taxa = tuple(sorted((*left_pair, *right_pair)))
                reference_quartet_split = canonical_bipartition(
                    set(left_pair),
                    set(quartet_taxa),
                )
                for tree_index, tree in enumerate(analysis.trees):
                    cached_split = quartet_split_cache[tree_index].get(quartet_taxa)
                    if quartet_taxa not in quartet_split_cache[tree_index]:
                        cached_split = _resolve_induced_quartet_split(
                            tree, quartet_taxa
                        )
                        quartet_split_cache[tree_index][quartet_taxa] = cached_split
                    if cached_split is None:
                        uninformative += 1
                    elif cached_split == reference_quartet_split:
                        concordant += 1
                    else:
                        discordant += 1
        total_row_quartets = quartet_count_per_tree * len(analysis.trees)
        total_quartets += total_row_quartets
        total_concordant += concordant
        total_discordant += discordant
        total_uninformative += uninformative
        rows.append(
            TreeSetQuartetSupportRow(
                branch_id=_format_split_id(split, taxon_scope),
                left_taxa=left_taxa,
                right_taxa=right_taxa,
                quartet_count_per_tree=quartet_count_per_tree,
                concordant_quartet_count=concordant,
                discordant_quartet_count=discordant,
                uninformative_quartet_count=uninformative,
                concordant_frequency=round(concordant / total_row_quartets, 15),
                discordant_frequency=round(discordant / total_row_quartets, 15),
                uninformative_frequency=round(uninformative / total_row_quartets, 15),
            )
        )
    return TreeSetQuartetSupportReport(
        reference_tree_path=reference_tree_path,
        comparison_tree_set_path=analysis.path,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        branch_count=len(rows),
        total_quartet_count=total_quartets,
        concordant_quartet_count=total_concordant,
        discordant_quartet_count=total_discordant,
        uninformative_quartet_count=total_uninformative,
        rows=rows,
    )


def compute_reference_tree_quartet_support(
    reference_tree_path: Path,
    comparison_tree_set_path: Path,
) -> TreeSetQuartetSupportReport:
    """Score every informative reference-tree branch by induced quartet agreement."""
    reference_tree = load_tree(reference_tree_path)
    analysis = _analyze_tree_set(comparison_tree_set_path)
    if analysis.exact_taxa is None:
        raise InvalidAlignmentError(
            "reference tree quartet support requires all comparison trees to share the exact same taxon set"
        )
    exact_taxa = analysis.exact_taxa
    if sorted(reference_tree.tip_names) != exact_taxa:
        raise InvalidAlignmentError(
            "reference tree and comparison tree set must share the exact same taxon set"
        )
    return _build_reference_tree_quartet_support_report(
        reference_tree_path=reference_tree_path,
        reference_tree=reference_tree,
        analysis=analysis,
    )


def write_reference_tree_quartet_support_table(
    path: Path,
    report: TreeSetQuartetSupportReport,
) -> Path:
    """Write one reference-tree quartet support table as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "branch_id",
                "left_taxa",
                "right_taxa",
                "quartet_count_per_tree",
                "concordant_quartet_count",
                "discordant_quartet_count",
                "uninformative_quartet_count",
                "concordant_frequency",
                "discordant_frequency",
                "uninformative_frequency",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.rows:
            writer.writerow(
                {
                    "branch_id": row.branch_id,
                    "left_taxa": "|".join(row.left_taxa),
                    "right_taxa": "|".join(row.right_taxa),
                    "quartet_count_per_tree": row.quartet_count_per_tree,
                    "concordant_quartet_count": row.concordant_quartet_count,
                    "discordant_quartet_count": row.discordant_quartet_count,
                    "uninformative_quartet_count": row.uninformative_quartet_count,
                    "concordant_frequency": format(row.concordant_frequency, ".15g"),
                    "discordant_frequency": format(row.discordant_frequency, ".15g"),
                    "uninformative_frequency": format(
                        row.uninformative_frequency,
                        ".15g",
                    ),
                }
            )
    return path
