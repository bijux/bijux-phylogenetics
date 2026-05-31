from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.io.newick import loads_newick, write_newick
from bijux_phylogenetics.phylo.topology.models import RandomBifurcatingTreeReport
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.simulation.trees import simulate_random_tree

_SUPPORTED_RANDOM_BIFURCATING_BRANCH_LENGTH_POLICIES = frozenset({"none", "uniform"})
_RANDOM_BIFURCATING_TAXON_PREFIX = "RandomBifurcatingTaxon"


def generate_random_bifurcating_tree(
    taxa: list[str],
    *,
    seed: int = 1,
    branch_length_policy: str = "none",
) -> tuple[PhyloTree, RandomBifurcatingTreeReport]:
    """Generate one rooted seeded bifurcating tree over an explicit taxon set."""
    ordered_taxa = validate_random_bifurcating_taxa(taxa)
    validated_branch_length_policy = validate_random_bifurcating_branch_length_policy(
        branch_length_policy
    )
    random_tree, _simulation_report = simulate_random_tree(
        tip_count=len(ordered_taxa),
        seed=seed,
        taxon_prefix=_RANDOM_BIFURCATING_TAXON_PREFIX,
    )
    label_map = {
        f"{_RANDOM_BIFURCATING_TAXON_PREFIX}{index}": taxon
        for index, taxon in enumerate(ordered_taxa, start=1)
    }
    for leaf in random_tree.iter_leaves():
        if leaf.name is None:
            raise ValueError(
                "random bifurcating tree generation produced an unnamed tip"
            )
        leaf.name = label_map[leaf.name]
        if validated_branch_length_policy == "none":
            leaf.branch_length = None
    if validated_branch_length_policy == "none":
        for node in random_tree.iter_internal_nodes(order="preorder"):
            if node is not random_tree.root:
                node.branch_length = None
    random_tree.refresh()
    report = summarize_random_bifurcating_tree(
        random_tree,
        requested_taxa=ordered_taxa,
        seed=seed,
        branch_length_policy=validated_branch_length_policy,
    )
    return random_tree, report


def validate_random_bifurcating_taxa(taxa: list[str]) -> list[str]:
    """Require at least two distinct non-empty taxa and canonicalize their order."""
    if len(taxa) < 2:
        raise ValueError(
            "random bifurcating tree generation requires at least two taxa"
        )
    blank_taxa = sorted({taxon for taxon in taxa if not taxon.strip()})
    if blank_taxa:
        raise ValueError(
            "random bifurcating tree generation does not allow blank taxon labels"
        )
    duplicate_taxa = sorted({taxon for taxon in taxa if taxa.count(taxon) > 1})
    if duplicate_taxa:
        raise ValueError(
            "random bifurcating tree generation requires distinct taxa; duplicates: "
            + ", ".join(duplicate_taxa)
        )
    return sorted(taxa)


def validate_random_bifurcating_branch_length_policy(branch_length_policy: str) -> str:
    """Validate the supported branch-length treatment for random bifurcating trees."""
    normalized_policy = branch_length_policy.strip().lower()
    if normalized_policy not in _SUPPORTED_RANDOM_BIFURCATING_BRANCH_LENGTH_POLICIES:
        raise ValueError(
            "branch_length_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_RANDOM_BIFURCATING_BRANCH_LENGTH_POLICIES))
        )
    return normalized_policy


def summarize_random_bifurcating_tree(
    tree: PhyloTree,
    *,
    requested_taxa: list[str],
    seed: int,
    branch_length_policy: str,
) -> RandomBifurcatingTreeReport:
    """Summarize one generated tree against the requested taxon set."""
    tip_order = tree.tip_names
    tip_name_counts: dict[str, int] = {}
    for taxon in tip_order:
        tip_name_counts[taxon] = tip_name_counts.get(taxon, 0) + 1
    duplicate_generated_taxa = sorted(
        taxon for taxon, count in tip_name_counts.items() if count > 1
    )
    generated_taxa = set(tip_order)
    requested_taxon_set = set(requested_taxa)
    missing_requested_taxa = sorted(requested_taxon_set - generated_taxa)
    unexpected_generated_taxa = sorted(generated_taxa - requested_taxon_set)
    validation_errors = tree.validation_errors()
    strictly_bifurcating = _tree_is_strictly_bifurcating(tree)
    return RandomBifurcatingTreeReport(
        algorithm="random-bifurcating-tree-generation",
        seed=seed,
        branch_length_policy=branch_length_policy,
        requested_taxa=list(requested_taxa),
        tip_order=tip_order,
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        rooted=tree.rooted,
        strictly_bifurcating=strictly_bifurcating,
        all_requested_taxa_present_once=not (
            missing_requested_taxa
            or duplicate_generated_taxa
            or unexpected_generated_taxa
        ),
        missing_requested_taxa=missing_requested_taxa,
        duplicate_generated_taxa=duplicate_generated_taxa,
        unexpected_generated_taxa=unexpected_generated_taxa,
        validation_errors=validation_errors,
        tree_newick=tree.to_newick(),
    )


def write_random_bifurcating_tree_report(
    path: Path,
    report: RandomBifurcatingTreeReport,
) -> Path:
    """Write one durable TSV summary for a seeded random bifurcating tree."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "algorithm",
                "seed",
                "branch_length_policy",
                "requested_taxa",
                "tip_order",
                "tip_count",
                "internal_node_count",
                "rooted",
                "strictly_bifurcating",
                "all_requested_taxa_present_once",
                "missing_requested_taxa",
                "duplicate_generated_taxa",
                "unexpected_generated_taxa",
                "validation_errors",
                "tree_newick",
            ]
        ),
        "\t".join(
            [
                report.algorithm,
                str(report.seed),
                report.branch_length_policy,
                ",".join(report.requested_taxa),
                ",".join(report.tip_order),
                str(report.tip_count),
                str(report.internal_node_count),
                "" if report.rooted is None else ("true" if report.rooted else "false"),
                "true" if report.strictly_bifurcating else "false",
                "true" if report.all_requested_taxa_present_once else "false",
                ",".join(report.missing_requested_taxa),
                ",".join(report.duplicate_generated_taxa),
                ",".join(report.unexpected_generated_taxa),
                " | ".join(report.validation_errors),
                report.tree_newick,
            ]
        ),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_random_bifurcating_tree_run_json(
    path: Path,
    report: RandomBifurcatingTreeReport,
) -> Path:
    """Write one machine-readable random bifurcating tree payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "algorithm": report.algorithm,
        "seed": report.seed,
        "branch_length_policy": report.branch_length_policy,
        "requested_taxa": report.requested_taxa,
        "tip_order": report.tip_order,
        "tip_count": report.tip_count,
        "internal_node_count": report.internal_node_count,
        "rooted": report.rooted,
        "strictly_bifurcating": report.strictly_bifurcating,
        "all_requested_taxa_present_once": report.all_requested_taxa_present_once,
        "missing_requested_taxa": report.missing_requested_taxa,
        "duplicate_generated_taxa": report.duplicate_generated_taxa,
        "unexpected_generated_taxa": report.unexpected_generated_taxa,
        "validation_errors": report.validation_errors,
        "tree_newick": report.tree_newick,
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def write_random_bifurcating_tree_artifacts(
    out_dir: Path,
    report: RandomBifurcatingTreeReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one random bifurcating tree run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    tree_path = write_newick(out_dir / "tree.nwk", loads_newick(report.tree_newick))
    summary_path = write_random_bifurcating_tree_report(out_dir / "summary.tsv", report)
    run_json_path = write_random_bifurcating_tree_run_json(out_dir / "run.json", report)
    return {
        "tree_path": tree_path,
        "summary_path": summary_path,
        "run_json_path": run_json_path,
    }


def _tree_is_strictly_bifurcating(tree: PhyloTree) -> bool:
    return all(node.is_leaf() or len(node.children) == 2 for node in tree.iter_nodes())
