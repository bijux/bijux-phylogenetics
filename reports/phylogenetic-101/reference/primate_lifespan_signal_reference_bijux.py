from __future__ import annotations

import csv
from dataclasses import asdict
import json
from math import erfc, sqrt
from pathlib import Path
import sys

from bijux_phylogenetics import (
    estimate_pagels_lambda,
    inspect_tree_path,
    reconstruct_continuous_ancestral_states,
    summarize_numeric_trait,
    summarize_numeric_trait_readiness,
    validate_tree_path,
)
from bijux_phylogenetics.ancestral.common import node_descendant_taxa, node_signature
from bijux_phylogenetics.comparative.common import (
    build_brownian_covariance_matrix,
    lambda_transform_covariance,
    load_comparative_dataset,
)
from bijux_phylogenetics.core.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.core.topology import (
    rotate_all_internal_nodes,
    rotate_named_node,
    unroot_tree,
)
from bijux_phylogenetics.io.trees import load_tree


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def _tip_node_id_lookup(tree_path: Path) -> dict[str, int]:
    tree = load_tree(tree_path)
    tip_lookup = {name: index for index, name in enumerate(tree.tip_names, start=1)}
    for node in tree.iter_nodes():
        if node.name and node.name.startswith("Node"):
            tip_lookup[node.name] = len(tree.tip_names) + int(node.name.removeprefix("Node"))
    return tip_lookup


def _find_smallest_covering_node(tree, taxa: set[str]):
    best_node = None
    best_descendants: list[str] | None = None
    for node in tree.iter_nodes():
        descendants = node_descendant_taxa(node)
        if taxa <= set(descendants):
            if best_descendants is None or len(descendants) < len(best_descendants):
                best_node = node
                best_descendants = descendants
    if best_node is None or best_descendants is None:
        raise ValueError(f"no covering node found for taxa: {sorted(taxa)}")
    return best_node, best_descendants


def _count_branch_increases(tree, estimate_lookup: dict[str, float]) -> dict[str, int]:
    increase_count = 0
    increase_gt12_count = 0

    def visit(node, parent_signature: str | None = None) -> None:
        nonlocal increase_count, increase_gt12_count
        signature = node_signature(node)
        estimate = estimate_lookup[signature]
        if parent_signature is not None:
            diff = estimate - estimate_lookup[parent_signature]
            if diff > 0.0:
                increase_count += 1
            if diff > 12.0:
                increase_gt12_count += 1
        for child in node.children:
            visit(child, signature)

    visit(tree.root)
    return {
        "increase_count": increase_count,
        "increase_gt12_count": increase_gt12_count,
    }


def _matrix_top3(matrix: list[list[float]]) -> list[list[float]]:
    return [row[:3] for row in matrix[:3]]


def _internal_clade_catalog(tree) -> list[dict[str, object]]:
    catalog: list[dict[str, object]] = []
    for node in tree.iter_nodes():
        if node.is_leaf():
            continue
        descendants = sorted(node_descendant_taxa(node))
        catalog.append(
            {
                "node_name": node.name,
                "signature": "|".join(descendants),
                "tip_count": len(descendants),
                "taxa": descendants,
            }
        )
    return catalog


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: primate_lifespan_signal_reference_bijux.py <r_repo_root> <out_dir>"
        )

    r_repo_root = Path(sys.argv[1]).resolve()
    out_dir = Path(sys.argv[2]).resolve()
    data_dir = r_repo_root / "PCM1_plots_signal" / "Lecture" / "R" / "Data"

    original_tree_path = data_dir / "primatetree.nex"
    trimmed_tree_path = data_dir / "trimmed_primatetree.nex"
    traits_path = data_dir / "primate.csv"

    processed_rows = _read_csv_rows(traits_path)
    tree = load_tree(trimmed_tree_path)
    original_tree = load_tree(original_tree_path)

    inspection = inspect_tree_path(trimmed_tree_path)
    validation = validate_tree_path(trimmed_tree_path)
    readiness = summarize_numeric_trait_readiness(
        trimmed_tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
    )
    summary = summarize_numeric_trait(
        trimmed_tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
    )

    pruned_tree, pruning_report = prune_tree_to_requested_taxa(
        original_tree_path,
        requested_taxa=summary.taxa,
    )
    unrooted_tree, unroot_report = unroot_tree(trimmed_tree_path)
    rotated_tree, rotated_report = rotate_named_node(
        trimmed_tree_path,
        clade_name="Node56",
    )
    rotated_all_tree, rotated_all_report = rotate_all_internal_nodes(trimmed_tree_path)
    internal_clades = _internal_clade_catalog(tree)

    tip_lookup = _tip_node_id_lookup(trimmed_tree_path)
    aligned_species = [row["species"] for row in processed_rows if row["species"] in tree.tip_names]
    aligned_species = sorted(aligned_species, key=lambda species: tree.tip_names.index(species))

    random_tree_path = out_dir / "random_tree_seed1.nwk"
    random_examples: list[dict[str, object]] = []
    for name in [
        "random_data",
        "random_data2",
        "random_data3",
        "random_data4",
        "random_data5",
    ]:
        report = estimate_pagels_lambda(
            random_tree_path,
            out_dir / f"{name}.csv",
            trait="value",
            taxon_column="species",
            fine_step=0.001,
        )
        random_examples.append(
            {
                "name": name,
                "lambda_value": report.lambda_value,
                "log_likelihood": report.log_likelihood,
                "tip_count": report.taxon_count,
            }
        )

    lambda_report = estimate_pagels_lambda(
        trimmed_tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
        fine_step=0.001,
    )
    ll_diff0 = -2.0 * (lambda_report.null_log_likelihood - lambda_report.log_likelihood)
    lambda_p_value = erfc(sqrt(ll_diff0 / 2.0))

    dataset = load_comparative_dataset(
        trimmed_tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
    )
    covariance = build_brownian_covariance_matrix(dataset.tree, dataset.taxa)
    lambda0_covariance = lambda_transform_covariance(covariance, 0.0)

    ancestral = reconstruct_continuous_ancestral_states(
        trimmed_tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
        model="brownian",
    )
    internal_nodes = [
        {
            "signature": estimate.node,
            "estimate": estimate.estimate,
            "lower_95": estimate.lower_95_interval,
            "upper_95": estimate.upper_95_interval,
        }
        for estimate in ancestral.estimates
        if not estimate.is_tip
    ]
    estimate_lookup = {estimate.node: estimate.estimate for estimate in ancestral.estimates}
    mrca_node, mrca_descendants = _find_smallest_covering_node(
        tree, {"Pan_paniscus", "Hylobates_lar"}
    )
    mrca_signature = node_signature(mrca_node)

    result = {
        "data_processing": {
            "processed_row_count": len(processed_rows),
            "processed_species_count": len({row["species"] for row in processed_rows}),
            "analysis_taxa_count": len(summary.taxa),
        },
        "tree_processing": {
            "original_tree": {
                "inspect": asdict(inspect_tree_path(original_tree_path)),
                "validate": asdict(validate_tree_path(original_tree_path)),
            },
            "trimmed_tree": {
                "inspect": asdict(inspection),
                "validate": asdict(validation),
            },
            "pruning": asdict(pruning_report),
        },
        "tree_examples": {
            "extract_clade": {
                "internal_clades": internal_clades,
            },
            "unroot_tree": {
                "rooted_flag": unrooted_tree.rooted,
                "tip_count": unrooted_tree.tip_count,
                "root_child_count": len(unrooted_tree.root.children),
                "report": asdict(unroot_report),
            },
            "rotate_node": {
                "source_node_label": "Node56",
                "tip_count": rotated_tree.tip_count,
                "tip_order": rotated_tree.tip_names,
                "same_tip_set": sorted(rotated_tree.tip_names) == sorted(tree.tip_names),
                "report": asdict(rotated_report),
            },
            "rotate_all": {
                "tip_count": rotated_all_tree.tip_count,
                "tip_order": rotated_all_tree.tip_names,
                "same_tip_set": sorted(rotated_all_tree.tip_names) == sorted(tree.tip_names),
                "report": asdict(rotated_all_report),
            },
        },
        "data_tree_alignment": {
            "aligned_species_equals_tip_order": aligned_species == tree.tip_names,
            "aligned_species_first_6": aligned_species[:6],
            "tip_order_first_6": tree.tip_names[:6],
            "nodeid_examples": {
                "Pan_paniscus": tip_lookup["Pan_paniscus"],
                "Hylobates_lar": tip_lookup["Hylobates_lar"],
                "Node32": tip_lookup["Node32"],
            },
            "readiness": asdict(readiness),
            "summary": asdict(summary),
        },
        "random_signal": {
            "seed": 1,
            "random_tree_tip_count": load_tree(random_tree_path).tip_count,
            "examples": random_examples,
        },
        "primate_lambda": {
            "lambda_value": lambda_report.lambda_value,
            "log_likelihood": lambda_report.log_likelihood,
        },
        "primate_lambda_zero": {
            "lambda0_log_likelihood": lambda_report.null_log_likelihood,
            "likelihood_ratio": ll_diff0,
            "p_value": lambda_p_value,
            "lambda0_vcv_top3": _matrix_top3(lambda0_covariance),
            "real_vcv_top3": _matrix_top3(covariance),
        },
        "ancestral": {
            "nodewise": internal_nodes,
            "mrca_signature": mrca_signature,
            "mrca_descendant_taxa": mrca_descendants,
            "mrca_estimate": next(
                estimate.estimate
                for estimate in ancestral.estimates
                if estimate.node == mrca_signature
            ),
            **_count_branch_increases(tree, estimate_lookup),
        },
    }

    (out_dir / "bijux_reference_results.json").write_text(
        json.dumps(result, indent=2, default=str) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
