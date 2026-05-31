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

from .contracts import (
    GeneTreeQuartetConcordanceReport,
    GeneTreeQuartetConcordanceRow,
)
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


def _quartet_topology_splits(
    quartet_taxa: tuple[str, str, str, str],
) -> list[frozenset[str]]:
    quartet_scope = set(quartet_taxa)
    first, second, third, fourth = quartet_taxa
    splits = {
        canonical_bipartition({first, second}, quartet_scope),
        canonical_bipartition({first, third}, quartet_scope),
        canonical_bipartition({first, fourth}, quartet_scope),
    }
    return sorted(splits, key=lambda signature: tuple(sorted(signature)))


def _build_gene_tree_quartet_concordance_report(
    *,
    species_tree_path: Path,
    species_tree: PhyloTree,
    analysis: _TreeSetAnalysis,
) -> GeneTreeQuartetConcordanceReport:
    exact_taxa = _require_exact_taxa(analysis)
    taxon_scope = set(exact_taxa)
    quartet_split_cache: list[
        dict[tuple[str, str, str, str], frozenset[str] | None]
    ] = [{} for _ in analysis.trees]
    rows: list[GeneTreeQuartetConcordanceRow] = []
    total_concordant = 0
    total_discordant_first = 0
    total_discordant_second = 0
    total_uninformative = 0
    total_informative = 0
    total_quartets = 0
    for split in sorted(
        informative_unrooted_splits(species_tree, taxon_scope),
        key=lambda signature: (len(signature), tuple(sorted(signature))),
    ):
        left_taxa = sorted(split)
        right_taxa = sorted(taxon_scope - split)
        concordant = 0
        discordant_first = 0
        discordant_second = 0
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
                discordant_splits = [
                    candidate
                    for candidate in _quartet_topology_splits(quartet_taxa)
                    if candidate != reference_quartet_split
                ]
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
                    elif cached_split == discordant_splits[0]:
                        discordant_first += 1
                    else:
                        discordant_second += 1
        total_row_quartets = quartet_count_per_tree * len(analysis.trees)
        informative = concordant + discordant_first + discordant_second
        total_quartets += total_row_quartets
        total_concordant += concordant
        total_discordant_first += discordant_first
        total_discordant_second += discordant_second
        total_uninformative += uninformative
        total_informative += informative
        rows.append(
            GeneTreeQuartetConcordanceRow(
                branch_id=_format_split_id(split, taxon_scope),
                left_taxa=left_taxa,
                right_taxa=right_taxa,
                quartet_count_per_tree=quartet_count_per_tree,
                concordant_quartet_count=concordant,
                discordant_first_quartet_count=discordant_first,
                discordant_second_quartet_count=discordant_second,
                uninformative_quartet_count=uninformative,
                informative_quartet_count=informative,
                concordance_factor=(
                    round(concordant / informative, 15) if informative else None
                ),
                concordant_frequency=round(concordant / total_row_quartets, 15),
                discordant_first_frequency=round(
                    discordant_first / total_row_quartets,
                    15,
                ),
                discordant_second_frequency=round(
                    discordant_second / total_row_quartets,
                    15,
                ),
                uninformative_frequency=round(uninformative / total_row_quartets, 15),
            )
        )
    return GeneTreeQuartetConcordanceReport(
        species_tree_path=species_tree_path,
        gene_tree_set_path=analysis.path,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        branch_count=len(rows),
        total_quartet_count=total_quartets,
        concordant_quartet_count=total_concordant,
        discordant_first_quartet_count=total_discordant_first,
        discordant_second_quartet_count=total_discordant_second,
        uninformative_quartet_count=total_uninformative,
        informative_quartet_count=total_informative,
        rows=rows,
    )


def compute_gene_tree_quartet_concordance_factors(
    species_tree_path: Path,
    gene_tree_set_path: Path,
) -> GeneTreeQuartetConcordanceReport:
    """Compute concordance factors for every informative species-tree branch from gene trees."""
    species_tree = load_tree(species_tree_path)
    analysis = _analyze_tree_set(gene_tree_set_path)
    if analysis.exact_taxa is None:
        raise InvalidAlignmentError(
            "gene tree quartet concordance requires all gene trees to share the exact same taxon set"
        )
    exact_taxa = analysis.exact_taxa
    if sorted(species_tree.tip_names) != exact_taxa:
        raise InvalidAlignmentError(
            "species tree and gene tree set must share the exact same taxon set"
        )
    return _build_gene_tree_quartet_concordance_report(
        species_tree_path=species_tree_path,
        species_tree=species_tree,
        analysis=analysis,
    )


def write_gene_tree_quartet_concordance_table(
    path: Path,
    report: GeneTreeQuartetConcordanceReport,
) -> Path:
    """Write one species-tree branch concordance-factor table as TSV."""
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
                "discordant_first_quartet_count",
                "discordant_second_quartet_count",
                "uninformative_quartet_count",
                "informative_quartet_count",
                "concordance_factor",
                "concordant_frequency",
                "discordant_first_frequency",
                "discordant_second_frequency",
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
                    "discordant_first_quartet_count": (
                        row.discordant_first_quartet_count
                    ),
                    "discordant_second_quartet_count": (
                        row.discordant_second_quartet_count
                    ),
                    "uninformative_quartet_count": row.uninformative_quartet_count,
                    "informative_quartet_count": row.informative_quartet_count,
                    "concordance_factor": (
                        ""
                        if row.concordance_factor is None
                        else format(row.concordance_factor, ".15g")
                    ),
                    "concordant_frequency": format(
                        row.concordant_frequency,
                        ".15g",
                    ),
                    "discordant_first_frequency": format(
                        row.discordant_first_frequency,
                        ".15g",
                    ),
                    "discordant_second_frequency": format(
                        row.discordant_second_frequency,
                        ".15g",
                    ),
                    "uninformative_frequency": format(
                        row.uninformative_frequency,
                        ".15g",
                    ),
                }
            )
    return path
