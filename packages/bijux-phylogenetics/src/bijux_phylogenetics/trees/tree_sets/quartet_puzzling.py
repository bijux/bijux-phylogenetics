from __future__ import annotations

import csv
from itertools import combinations
import json
from pathlib import Path
from random import Random

from bijux_phylogenetics.io.newick import (
    dumps_newick,
    loads_newick,
    write_newick,
    write_newick_tree_set,
)
from bijux_phylogenetics.phylo.pruning import prune_tree_object_to_requested_taxa
from bijux_phylogenetics.phylo.topology.clades import informative_unrooted_splits
from bijux_phylogenetics.phylo.topology.rooting import _root_tree_by_outgroup_node
from bijux_phylogenetics.phylo.topology.stepwise_addition import (
    apply_stepwise_addition_candidate,
    iter_stepwise_addition_edge_candidates,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

from .consensus import _build_consensus_tree_with_threshold_from_trees
from .contracts import (
    QuartetPuzzlingAssemblyRow,
    QuartetPuzzlingReport,
    QuartetTopologyScoreRow,
)
from .inventory import _analyze_tree_set, _require_exact_taxa, _TreeSetAnalysis
from .topology import _unrooted_topology_id

DEFAULT_QUARTET_PUZZLING_MAX_ORDER_COUNT = 32
DEFAULT_QUARTET_PUZZLING_RANDOM_SEED = 0
CANONICAL_QUARTET_PUZZLING_ROOTING_STRATEGY = "lexicographic-tip-outgroup"


def _canonical_quartet(
    quartet_taxa: tuple[str, str, str, str],
) -> tuple[str, str, str, str]:
    return tuple(sorted(quartet_taxa))


def _quartet_split_options(
    quartet_taxa: tuple[str, str, str, str],
) -> tuple[frozenset[str], frozenset[str], frozenset[str]]:
    first, second, third, fourth = quartet_taxa
    return (
        frozenset({first, second}),
        frozenset({first, third}),
        frozenset({first, fourth}),
    )


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


def _format_quartet_score_row(
    quartet_taxa: tuple[str, str, str, str],
    split_counts: dict[frozenset[str], int],
    *,
    uninformative_tree_count: int,
    tree_count: int,
) -> QuartetTopologyScoreRow:
    first_split, second_split, third_split = _quartet_split_options(quartet_taxa)
    scored_splits = [
        (first_split, split_counts.get(first_split, 0)),
        (second_split, split_counts.get(second_split, 0)),
        (third_split, split_counts.get(third_split, 0)),
    ]
    best_count = max(count for _split, count in scored_splits)
    tied_best_splits = [
        sorted(split)
        for split, count in scored_splits
        if count == best_count and count > 0
    ]
    best_split_taxa = None if len(tied_best_splits) != 1 else tied_best_splits[0]
    best_split_support_frequency = (
        None if best_split_taxa is None else round(best_count / tree_count, 15)
    )
    return QuartetTopologyScoreRow(
        quartet_taxa=list(quartet_taxa),
        first_split_taxa=sorted(first_split),
        first_split_tree_count=split_counts.get(first_split, 0),
        second_split_taxa=sorted(second_split),
        second_split_tree_count=split_counts.get(second_split, 0),
        third_split_taxa=sorted(third_split),
        third_split_tree_count=split_counts.get(third_split, 0),
        uninformative_tree_count=uninformative_tree_count,
        best_split_taxa=best_split_taxa,
        best_split_support_frequency=best_split_support_frequency,
        tied_best_split_taxa=sorted(tied_best_splits),
    )


def _score_tree_set_quartets(
    analysis: _TreeSetAnalysis,
) -> tuple[
    list[QuartetTopologyScoreRow],
    dict[tuple[str, str, str, str], dict[frozenset[str], float]],
]:
    exact_taxa = _require_exact_taxa(analysis)
    quartet_rows: list[QuartetTopologyScoreRow] = []
    quartet_frequency_lookup: dict[
        tuple[str, str, str, str], dict[frozenset[str], float]
    ] = {}
    for quartet_taxa in combinations(exact_taxa, 4):
        canonical_quartet = _canonical_quartet(tuple(quartet_taxa))
        split_counts: dict[frozenset[str], int] = {}
        uninformative_tree_count = 0
        for tree in analysis.trees:
            induced_split = _resolve_induced_quartet_split(tree, canonical_quartet)
            if induced_split is None:
                uninformative_tree_count += 1
                continue
            split_counts[induced_split] = split_counts.get(induced_split, 0) + 1
        quartet_rows.append(
            _format_quartet_score_row(
                canonical_quartet,
                split_counts,
                uninformative_tree_count=uninformative_tree_count,
                tree_count=len(analysis.trees),
            )
        )
        quartet_frequency_lookup[canonical_quartet] = {
            split: round(count / len(analysis.trees), 15)
            for split, count in split_counts.items()
        }
    return quartet_rows, quartet_frequency_lookup


def _build_quartet_start_tree(
    quartet_taxa: tuple[str, str, str, str],
    quartet_lookup: dict[tuple[str, str, str, str], dict[frozenset[str], float]],
) -> PhyloTree:
    split_options = _quartet_split_options(quartet_taxa)
    frequencies = quartet_lookup[quartet_taxa]
    scored_options = sorted(
        (
            (-frequencies.get(split, 0.0), tuple(sorted(split)), split)
            for split in split_options
        ),
    )
    selected_split = scored_options[0][2]
    complement = frozenset(quartet_taxa) - selected_split
    left_taxa = sorted(selected_split)
    right_taxa = sorted(complement)
    return PhyloTree(
        root=TreeNode(
            children=[
                TreeNode(
                    children=[TreeNode(name=left_taxa[0]), TreeNode(name=left_taxa[1])]
                ),
                TreeNode(
                    children=[
                        TreeNode(name=right_taxa[0]),
                        TreeNode(name=right_taxa[1]),
                    ]
                ),
            ]
        ),
        rooted=True,
    ).refresh()


def _score_quartet_puzzling_tree(
    tree: PhyloTree,
    quartet_lookup: dict[tuple[str, str, str, str], dict[frozenset[str], float]],
) -> float:
    if tree.tip_count < 4:
        return 0.0
    score = 0.0
    for quartet_taxa in combinations(sorted(tree.tip_names), 4):
        canonical_quartet = _canonical_quartet(tuple(quartet_taxa))
        induced_split = _resolve_induced_quartet_split(tree, canonical_quartet)
        if induced_split is None:
            continue
        score += quartet_lookup.get(canonical_quartet, {}).get(induced_split, 0.0)
    return round(score, 15)


def _assemble_quartet_puzzling_tree(
    taxon_order: list[str],
    quartet_lookup: dict[tuple[str, str, str, str], dict[frozenset[str], float]],
) -> tuple[PhyloTree, float]:
    current_tree = _build_quartet_start_tree(
        _canonical_quartet(tuple(taxon_order[:4])),
        quartet_lookup,
    )
    current_score = _score_quartet_puzzling_tree(current_tree, quartet_lookup)
    for taxon in taxon_order[4:]:
        best_tree: PhyloTree | None = None
        best_score: float | None = None
        best_newick: str | None = None
        for candidate in iter_stepwise_addition_edge_candidates(current_tree):
            candidate_tree = apply_stepwise_addition_candidate(
                current_tree, candidate, taxon
            )
            candidate_score = _score_quartet_puzzling_tree(
                candidate_tree, quartet_lookup
            )
            candidate_newick = candidate_tree.to_newick()
            if (
                (
                    best_score is None
                    or candidate_score > best_score
                    or (
                        candidate_score == best_score
                        and best_newick is not None
                        and candidate_newick < best_newick
                    )
                )
                or best_score is None
                and best_newick is None
            ):
                best_tree = candidate_tree
                best_score = candidate_score
                best_newick = candidate_newick
        if best_tree is None or best_score is None:
            raise AssertionError(
                "quartet puzzling insertion must evaluate at least one edge"
            )
        current_tree = best_tree
        current_score = best_score
    return current_tree, current_score


def _find_named_tip(tree: PhyloTree, tip_name: str) -> TreeNode:
    for node in tree.iter_leaves():
        if node.name == tip_name:
            return node
    raise ValueError(f"tree does not contain tip '{tip_name}'")


def _canonicalize_quartet_puzzling_tree(
    tree: PhyloTree,
    *,
    canonical_root_taxon: str,
) -> PhyloTree:
    """Root one assembled unrooted topology on a deterministic singleton outgroup."""
    return _root_tree_by_outgroup_node(
        tree,
        outgroup_node=_find_named_tip(tree, canonical_root_taxon),
    )


def _generate_taxon_orders(
    taxa: list[str],
    *,
    max_order_count: int,
    random_seed: int,
) -> list[list[str]]:
    if max_order_count <= 0:
        raise ValueError("quartet puzzling requires at least one taxon order")
    ordered_taxa = sorted(taxa)
    candidate_orders: list[list[str]] = []
    seen: set[tuple[str, ...]] = set()

    def add(order: list[str]) -> None:
        signature = tuple(order)
        if signature in seen or len(candidate_orders) >= max_order_count:
            return
        seen.add(signature)
        candidate_orders.append(order)

    add(list(ordered_taxa))
    add(list(reversed(ordered_taxa)))
    for offset in range(len(ordered_taxa)):
        rotated = ordered_taxa[offset:] + ordered_taxa[:offset]
        add(rotated)
        add(list(reversed(rotated)))
    rng = Random(random_seed)  # nosec B311
    while len(candidate_orders) < max_order_count:
        shuffled = list(ordered_taxa)
        rng.shuffle(shuffled)
        add(shuffled)
        if len(seen) >= math_factorial_bound(len(ordered_taxa)):
            break
    return candidate_orders


def math_factorial_bound(value: int) -> int:
    product = 1
    for factor in range(2, value + 1):
        product *= factor
    return product


def build_quartet_puzzling_consensus(
    tree_set_path: Path,
    *,
    max_order_count: int = DEFAULT_QUARTET_PUZZLING_MAX_ORDER_COUNT,
    random_seed: int = DEFAULT_QUARTET_PUZZLING_RANDOM_SEED,
    consensus_threshold: float = 0.5,
) -> tuple[PhyloTree, QuartetPuzzlingReport]:
    """Assemble quartet-scored trees across taxon orders and summarize consensus support."""
    if not 0.0 < consensus_threshold <= 1.0:
        raise ValueError(
            "quartet puzzling consensus threshold must be greater than 0 and at most 1"
        )
    analysis = _analyze_tree_set(tree_set_path)
    exact_taxa = _require_exact_taxa(analysis)
    if len(exact_taxa) < 4:
        raise InvalidAlignmentError(
            "quartet puzzling requires at least four shared taxa"
        )
    quartet_rows, quartet_lookup = _score_tree_set_quartets(analysis)
    taxon_orders = _generate_taxon_orders(
        exact_taxa,
        max_order_count=max_order_count,
        random_seed=random_seed,
    )
    canonical_root_taxon = exact_taxa[0]
    shared_taxa_set = set(exact_taxa)
    assembled_trees: list[PhyloTree] = []
    assembly_rows: list[QuartetPuzzlingAssemblyRow] = []
    for order_index, taxon_order in enumerate(taxon_orders, start=1):
        assembled_tree, quartet_score = _assemble_quartet_puzzling_tree(
            taxon_order,
            quartet_lookup,
        )
        assembled_topology_id = _unrooted_topology_id(assembled_tree, shared_taxa_set)
        canonical_tree = _canonicalize_quartet_puzzling_tree(
            assembled_tree,
            canonical_root_taxon=canonical_root_taxon,
        )
        assembled_trees.append(canonical_tree)
        assembly_rows.append(
            QuartetPuzzlingAssemblyRow(
                order_index=order_index,
                taxon_order=list(taxon_order),
                quartet_score=quartet_score,
                assembled_topology_id=assembled_topology_id,
                assembled_tree_newick=dumps_newick(canonical_tree),
            )
        )
    consensus_tree, included_clade_count = (
        _build_consensus_tree_with_threshold_from_trees(
            assembled_trees,
            threshold=consensus_threshold,
        )
    )
    report = QuartetPuzzlingReport(
        tree_set_path=tree_set_path,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        quartet_count=len(quartet_rows),
        assembly_count=len(assembly_rows),
        unique_assembled_topology_count=len(
            {row.assembled_topology_id for row in assembly_rows}
        ),
        canonical_root_taxon=canonical_root_taxon,
        canonical_rooting_strategy=CANONICAL_QUARTET_PUZZLING_ROOTING_STRATEGY,
        consensus_method="majority-rule"
        if consensus_threshold == 0.5
        else "thresholded",
        consensus_threshold=consensus_threshold,
        included_clade_count=included_clade_count,
        consensus_newick=dumps_newick(consensus_tree),
        quartet_rows=quartet_rows,
        assembly_rows=assembly_rows,
    )
    return consensus_tree, report


def write_quartet_topology_scores_table(
    path: Path,
    report: QuartetPuzzlingReport,
) -> Path:
    """Write one deterministic table of scored quartet topologies."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "quartet_taxa",
                "first_split_taxa",
                "first_split_tree_count",
                "second_split_taxa",
                "second_split_tree_count",
                "third_split_taxa",
                "third_split_tree_count",
                "uninformative_tree_count",
                "best_split_taxa",
                "best_split_support_frequency",
                "tied_best_split_taxa",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.quartet_rows:
            writer.writerow(
                {
                    "quartet_taxa": "|".join(row.quartet_taxa),
                    "first_split_taxa": "|".join(row.first_split_taxa),
                    "first_split_tree_count": row.first_split_tree_count,
                    "second_split_taxa": "|".join(row.second_split_taxa),
                    "second_split_tree_count": row.second_split_tree_count,
                    "third_split_taxa": "|".join(row.third_split_taxa),
                    "third_split_tree_count": row.third_split_tree_count,
                    "uninformative_tree_count": row.uninformative_tree_count,
                    "best_split_taxa": ""
                    if row.best_split_taxa is None
                    else "|".join(row.best_split_taxa),
                    "best_split_support_frequency": ""
                    if row.best_split_support_frequency is None
                    else format(row.best_split_support_frequency, ".15g"),
                    "tied_best_split_taxa": ";".join(
                        "|".join(split) for split in row.tied_best_split_taxa
                    ),
                }
            )
    return path


def write_quartet_puzzling_assembly_table(
    path: Path,
    report: QuartetPuzzlingReport,
) -> Path:
    """Write one deterministic quartet-puzzling assembly ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "order_index",
                "taxon_order",
                "quartet_score",
                "assembled_topology_id",
                "assembled_tree_newick",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.assembly_rows:
            writer.writerow(
                {
                    "order_index": row.order_index,
                    "taxon_order": "|".join(row.taxon_order),
                    "quartet_score": format(row.quartet_score, ".15g"),
                    "assembled_topology_id": row.assembled_topology_id,
                    "assembled_tree_newick": row.assembled_tree_newick,
                }
            )
    return path


def write_quartet_puzzling_run_json(
    path: Path,
    report: QuartetPuzzlingReport,
) -> Path:
    """Write one machine-readable quartet-puzzling payload."""
    payload = {
        "tree_set_path": str(report.tree_set_path),
        "tree_count": report.tree_count,
        "processing": {
            "runtime_seconds": report.processing.runtime_seconds,
            "peak_memory_bytes": report.processing.peak_memory_bytes,
            "skipped_malformed_tree_count": (
                report.processing.skipped_malformed_tree_count
            ),
        },
        "shared_taxa": report.shared_taxa,
        "quartet_count": report.quartet_count,
        "assembly_count": report.assembly_count,
        "unique_assembled_topology_count": report.unique_assembled_topology_count,
        "canonical_root_taxon": report.canonical_root_taxon,
        "canonical_rooting_strategy": report.canonical_rooting_strategy,
        "consensus_method": report.consensus_method,
        "consensus_threshold": report.consensus_threshold,
        "included_clade_count": report.included_clade_count,
        "consensus_newick": report.consensus_newick,
        "quartet_rows": [
            {
                "quartet_taxa": row.quartet_taxa,
                "first_split_taxa": row.first_split_taxa,
                "first_split_tree_count": row.first_split_tree_count,
                "second_split_taxa": row.second_split_taxa,
                "second_split_tree_count": row.second_split_tree_count,
                "third_split_taxa": row.third_split_taxa,
                "third_split_tree_count": row.third_split_tree_count,
                "uninformative_tree_count": row.uninformative_tree_count,
                "best_split_taxa": row.best_split_taxa,
                "best_split_support_frequency": row.best_split_support_frequency,
                "tied_best_split_taxa": row.tied_best_split_taxa,
            }
            for row in report.quartet_rows
        ],
        "assembly_rows": [
            {
                "order_index": row.order_index,
                "taxon_order": row.taxon_order,
                "quartet_score": row.quartet_score,
                "assembled_topology_id": row.assembled_topology_id,
                "assembled_tree_newick": row.assembled_tree_newick,
            }
            for row in report.assembly_rows
        ],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_quartet_puzzling_artifacts(
    out_dir: Path,
    report: QuartetPuzzlingReport,
) -> dict[str, Path]:
    """Write the governed artifact bundle for one quartet-puzzling run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    consensus_tree_path = write_newick(
        out_dir / "consensus_tree.nwk",
        loads_newick(report.consensus_newick),
    )
    assembled_trees_path = write_newick_tree_set(
        out_dir / "assembled_trees.nwk",
        [loads_newick(row.assembled_tree_newick) for row in report.assembly_rows],
    )
    quartet_scores_path = write_quartet_topology_scores_table(
        out_dir / "quartet_scores.tsv",
        report,
    )
    assembly_scores_path = write_quartet_puzzling_assembly_table(
        out_dir / "assembly_scores.tsv",
        report,
    )
    run_json_path = write_quartet_puzzling_run_json(out_dir / "run.json", report)
    return {
        "consensus_tree_path": consensus_tree_path,
        "assembled_trees_path": assembled_trees_path,
        "quartet_scores_path": quartet_scores_path,
        "assembly_scores_path": assembly_scores_path,
        "run_json_path": run_json_path,
    }
