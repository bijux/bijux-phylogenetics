from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, replace
from importlib import metadata
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess  # nosec B404 - parity helpers invoke repository-owned reference commands
import tempfile

from bijux_phylogenetics.ancestral.common import (
    load_continuous_dataset,
    load_discrete_dataset,
    node_signature,
)
from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete import (
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.trees import extract_tree_clades, extract_tree_set_clades
from bijux_phylogenetics.comparative.brownian_covariance import (
    summarize_brownian_covariance,
)
from bijux_phylogenetics.comparative.signal import (
    compute_phylogenetic_independent_contrasts,
)
from bijux_phylogenetics.compare.structural_parity import (
    compare_tree_sets_structurally,
    compare_tree_structurally,
)
from bijux_phylogenetics.compare.topology_distance import (
    compare_topology_distance_trees,
)
from bijux_phylogenetics.core._node_identity import build_ape_internal_node_map
from bijux_phylogenetics.core.branching_times import compute_tree_branching_times
from bijux_phylogenetics.core.clade_sets import informative_unrooted_splits
from bijux_phylogenetics.core.node_depth import compute_tree_node_depths
from bijux_phylogenetics.core.pruning import (
    drop_tree_taxa,
    prune_tree_to_requested_taxa,
)
from bijux_phylogenetics.core.topology import (
    assess_tree_monophyly,
    extract_tree_clade_by_node_id,
    find_tree_mrca,
    root_tree_on_outgroup,
    unroot_tree,
)
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.core.tree_distance import compute_tree_tip_distance_matrix
from bijux_phylogenetics.core.ultrametric import assess_tree_ultrametricity
from bijux_phylogenetics.diagnostics.validation import inspect_tree_path
from bijux_phylogenetics.distance import (
    build_tree_from_imported_distance_matrix,
    compute_pairwise_genetic_distance_matrix,
)
from bijux_phylogenetics.fixtures import get_shared_tree_simulation_fixture
from bijux_phylogenetics.comparative import (
    compute_diversification_gamma_statistic,
)
from bijux_phylogenetics.io.fasta.coding import translate_coding_alignment
from bijux_phylogenetics.io.fasta.matrix import (
    compute_alignment_base_frequency_report,
    compute_alignment_segregating_site_report,
    load_dna_bin_alignment,
)
from bijux_phylogenetics.io.newick import (
    dumps_newick,
    load_newick_tree_set,
    write_newick,
    write_newick_tree_set,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.simulation import (
    simulate_coalescent_trees,
    simulate_random_trees,
)
from bijux_phylogenetics.trees import (
    compute_clade_frequency_table,
    compute_consensus_tree,
    compute_reference_tree_clade_support,
    compute_strict_consensus_tree,
)
from .registry import ApeParityCase, _selected_cases, _write_case_file


@dataclass(frozen=True, slots=True)
class ApeParityObservation:
    """One live parity comparison between Bijux and `ape`."""

    case_id: str
    fixture_kind: str
    fixture_id: str
    function_name: str
    python_function_name: str
    input_fixture: Path
    tolerance: float
    r_version: str | None
    ape_version: str | None
    bijux_version: str
    bijux_commit: str | None
    status: str
    passed: bool
    mismatch_reason: str | None
    reproducible_artifact_root: Path | None
    reference_summary: dict[str, object] | None
    bijux_summary: dict[str, object] | None
    reference_error: dict[str, object] | None
    bijux_error: dict[str, object] | None


@dataclass(frozen=True, slots=True)
class ApeParitySummaryRow:
    """One function-level summary across governed `ape` parity cases."""

    function_name: str
    case_count: int
    passed_case_count: int
    failed_case_count: int
    skipped_case_count: int


@dataclass(slots=True)
class ApeParityReport:
    """Aggregate report for governed live `ape` parity cases."""

    observations: list[ApeParityObservation]
    summary_rows: list[ApeParitySummaryRow]
    case_count: int
    passed_case_count: int
    failed_case_count: int
    skipped_case_count: int
    all_passed: bool
    limitations: list[str]


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[6]


def _ape_runner_path() -> Path:
    return (
        Path(__file__).resolve().parents[2]
        / "resources"
        / "reference"
        / "ape_parity_runner.R"
    )


def _failure_root() -> Path:
    return _repository_root() / "artifacts" / "ape-parity-failures"


def _reference_environment() -> dict[str, str]:
    environment = dict(os.environ)
    r_library = _repository_root() / "artifacts" / "r-lib"
    if "R_LIBS_USER" not in environment and r_library.is_dir():
        environment["R_LIBS_USER"] = str(r_library)
    return environment


def _bijux_version() -> str:
    try:
        return metadata.version("bijux-phylogenetics")
    except metadata.PackageNotFoundError:
        return "0.1.0"


def _bijux_commit() -> str | None:
    # Fixed repository git metadata probe.
    result = subprocess.run(  # nosec
        ["git", "rev-parse", "--short", "HEAD"],
        capture_output=True,
        check=False,
        cwd=_repository_root(),
        text=True,
    )
    if result.returncode != 0:
        return None
    commit = result.stdout.strip()
    return commit or None


def _node_kind_order(node_kind: str) -> int:
    return {"root": 0, "internal": 1, "tip": 2}.get(node_kind, 9)


def _clade_rows_to_parity_rows(rows) -> list[dict[str, object]]:
    parity_rows = [
        {
            "tree_index": "" if row.tree_index is None else row.tree_index,
            "node_kind": row.node_kind,
            "clade_id": row.clade_id,
            "node_label": "" if row.node_label is None else row.node_label,
            "taxon_count": row.taxon_count,
            "taxa": "|".join(row.taxa),
            "support": "" if row.support is None else row.support,
            "branch_length": "" if row.branch_length is None else row.branch_length,
        }
        for row in rows
    ]
    return _sort_parity_rows(parity_rows)


def _build_bijux_tree_structure(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    tree = load_tree(input_fixture)
    inspection = inspect_tree_path(input_fixture)
    clades = extract_tree_clades(input_fixture)
    return _tree_structure_payload(tree, inspection.rooted, clades.rows)


def _tree_structure_payload(
    tree: TreeNode | PhyloTree,
    rooted: bool | None,
    clade_rows,
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    phylo_tree = (
        tree if isinstance(tree, PhyloTree) else PhyloTree(root=tree, rooted=rooted)
    )
    summary = {
        "tree_count": 1,
        "tip_count": phylo_tree.tip_count,
        "internal_node_count": phylo_tree.internal_node_count,
        "edge_count": phylo_tree.tip_count + phylo_tree.internal_node_count - 1,
        "rooted": rooted,
        "tip_labels": phylo_tree.tip_names,
        "branch_length_count": sum(
            1
            for branch_length in phylo_tree.branch_lengths()
            if branch_length is not None
        ),
    }
    return summary, _clade_rows_to_parity_rows(clade_rows), dumps_newick(phylo_tree)


def _build_bijux_tree_set_structure(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]], None]:
    clades = extract_tree_set_clades(input_fixture)
    parity_rows = _clade_rows_to_parity_rows(clades.rows)
    tree_indices = sorted(
        {
            tree_index
            for row in parity_rows
            for tree_index in [row.get("tree_index")]
            if isinstance(tree_index, int)
        }
    )
    first_tree_tip_labels = [
        node_label
        for row in parity_rows
        for node_label in [row.get("node_label")]
        if row.get("tree_index") == 1
        and row.get("node_kind") == "tip"
        and isinstance(node_label, str)
    ]
    summary = {
        "tree_count": clades.tree_count,
        "source_format": clades.source_format,
        "tree_indices": tree_indices,
        "shared_tip_labels": sorted(first_tree_tip_labels),
        "unique_tip_label_count": len(first_tree_tip_labels),
    }
    return summary, parity_rows, None


def _build_bijux_root_outgroup_structure(
    input_fixture: Path,
    *,
    outgroup_taxa: tuple[str, ...],
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    rooted_tree, _report = root_tree_on_outgroup(
        input_fixture,
        outgroup_taxa=list(outgroup_taxa),
    )
    with tempfile.TemporaryDirectory(prefix="bijux-ape-root-") as tmpdir:
        rooted_path = Path(tmpdir) / "rooted.nwk"
        write_newick(rooted_path, rooted_tree)
        clades = extract_tree_clades(rooted_path)
    return _tree_structure_payload(rooted_tree, rooted_tree.rooted, clades.rows)


def _build_bijux_unroot_structure(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    unrooted_tree, _report = unroot_tree(input_fixture)
    with tempfile.TemporaryDirectory(prefix="bijux-ape-unroot-") as tmpdir:
        unrooted_path = Path(tmpdir) / "unrooted.nwk"
        write_newick(unrooted_path, unrooted_tree)
        clades = extract_tree_clades(unrooted_path)
    return _tree_structure_payload(unrooted_tree, unrooted_tree.rooted, clades.rows)


def _build_bijux_drop_tip_structure(
    input_fixture: Path,
    *,
    excluded_taxa: tuple[str, ...],
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    pruned_tree, report = drop_tree_taxa(input_fixture, list(excluded_taxa))
    with tempfile.TemporaryDirectory(prefix="bijux-ape-drop-tip-") as tmpdir:
        pruned_path = Path(tmpdir) / "drop-tip.nwk"
        write_newick(pruned_path, pruned_tree)
        clades = extract_tree_clades(pruned_path)
    summary, rows, normalized_text = _tree_structure_payload(
        pruned_tree,
        pruned_tree.rooted,
        clades.rows,
    )
    summary["dropped_taxa"] = report.removed_taxa
    summary["absent_requested_taxa"] = report.absent_requested_taxa
    return summary, rows, normalized_text


def _build_bijux_keep_tip_structure(
    input_fixture: Path,
    *,
    requested_taxa: tuple[str, ...],
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    pruned_tree, report = prune_tree_to_requested_taxa(
        input_fixture, list(requested_taxa)
    )
    with tempfile.TemporaryDirectory(prefix="bijux-ape-keep-tip-") as tmpdir:
        pruned_path = Path(tmpdir) / "keep-tip.nwk"
        write_newick(pruned_path, pruned_tree)
        clades = extract_tree_clades(pruned_path)
    summary, rows, normalized_text = _tree_structure_payload(
        pruned_tree,
        pruned_tree.rooted,
        clades.rows,
    )
    summary["requested_taxa"] = report.requested_taxa
    summary["dropped_taxa"] = report.removed_taxa
    return summary, rows, normalized_text


def _build_bijux_extract_clade_structure(
    input_fixture: Path,
    *,
    node_id: int,
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    subtree, report = extract_tree_clade_by_node_id(input_fixture, node_id=node_id)
    with tempfile.TemporaryDirectory(prefix="bijux-ape-extract-clade-") as tmpdir:
        subtree_path = Path(tmpdir) / "extract-clade.nwk"
        write_newick(subtree_path, subtree)
        clades = extract_tree_clades(subtree_path)
    summary, rows, normalized_text = _tree_structure_payload(
        subtree,
        subtree.rooted,
        clades.rows,
    )
    summary["requested_node_id"] = report.requested_node_id
    summary["matched_node_id"] = report.matched_node_id
    summary["matched_node_name"] = report.matched_node_name
    return summary, rows, normalized_text


def _build_bijux_mrca_summary(
    input_fixture: Path,
    *,
    mrca_taxa: tuple[str, ...],
) -> dict[str, object]:
    report = find_tree_mrca(input_fixture, taxa=list(mrca_taxa))
    return {
        "requested_taxa": report.requested_taxa,
        "unique_requested_taxa": report.unique_requested_taxa,
        "duplicate_requested_taxa": report.duplicate_requested_taxa,
        "matched_node_id": report.matched_node_id,
        "matched_node_name": report.matched_node_name or "",
        "matched_taxa": report.matched_taxa,
        "matched_extra_taxa": report.matched_extra_taxa,
        "matched_tip_count": report.matched_tip_count,
        "is_root": report.is_root,
    }


def _build_bijux_monophyly_summary(
    input_fixture: Path,
    *,
    requested_taxa: tuple[str, ...],
    reroot: bool,
) -> dict[str, object]:
    report = assess_tree_monophyly(
        input_fixture,
        taxa=list(requested_taxa),
        reroot=reroot,
    )
    return {
        "requested_taxa": report.requested_taxa,
        "unique_requested_taxa": report.unique_requested_taxa,
        "duplicate_requested_taxa": report.duplicate_requested_taxa,
        "missing_requested_taxa": report.missing_requested_taxa,
        "present_requested_taxa": report.present_requested_taxa,
        "reroot": report.reroot,
        "rooted": report.rooted,
        "monophyletic": report.monophyletic,
        "complementary_clade_used": report.complementary_clade_used,
        "matched_node_id": report.matched_node_id,
        "matched_node_name": report.matched_node_name or "",
        "matched_taxa": report.matched_taxa,
        "matched_extra_taxa": report.matched_extra_taxa,
        "matched_tip_count": report.matched_tip_count,
        "is_root": report.is_root,
    }


def _materialize_reference_input(case: ApeParityCase, working_root: Path) -> Path:
    reference_input_path = working_root / "bijux-reference-input.nwk"
    if case.operation == "write-tree-structure":
        tree = load_tree(case.input_fixture)
        write_newick(reference_input_path, tree)
        return reference_input_path
    if case.operation == "write-tree-set-structure":
        trees = load_newick_tree_set(case.input_fixture)
        write_newick_tree_set(reference_input_path, trees)
        return reference_input_path
    return case.input_fixture


def _build_bijux_dnabin_rows(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    matrix = load_dna_bin_alignment(input_fixture)
    return {
        "sequence_count": matrix.sequence_count,
        "alignment_length": matrix.alignment_length,
        "state_count": len(matrix.rows),
    }, [
        {
            "identifier": row.identifier,
            "position": row.position,
            "state": row.state,
        }
        for row in matrix.rows
    ]


def _ape_base_frequency_rows(input_fixture: Path) -> list[dict[str, object]]:
    report = compute_alignment_base_frequency_report(input_fixture)
    return [
        {"state": row.state, "frequency": row.frequency}
        for row in report.alignment_rows
    ]


def _build_bijux_base_frequency_summary(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_alignment_base_frequency_report(input_fixture)
    rows = [
        {"state": row.state, "frequency": row.frequency}
        for row in report.alignment_rows
    ]
    return {
        "sequence_count": report.sequence_count,
        "alignment_length": report.alignment_length,
        "state_count": len(rows),
    }, rows


def _build_bijux_segregating_site_rows(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_alignment_segregating_site_report(input_fixture)
    return {
        "sequence_count": report.sequence_count,
        "alignment_length": report.alignment_length,
        "segregating_site_count": len(report.segregating_site_positions),
    }, [{"position": row.position} for row in report.rows]


def _build_bijux_distance_rows(
    input_fixture: Path,
    *,
    pairwise_deletion: bool,
    distance_model: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_pairwise_genetic_distance_matrix(
        input_fixture,
        model=distance_model,
        gap_handling="pairwise-deletion" if pairwise_deletion else "complete-deletion",
        ambiguity_policy="ignore",
    )
    pair_lookup = {
        (row.left_identifier, row.right_identifier): row for row in report.pairs
    }
    rows: list[dict[str, object]] = []
    finite_distance_count = 0
    undefined_distance_count = 0
    infinite_distance_count = 0
    for left_identifier in report.identifiers:
        for right_identifier in report.identifiers:
            pair = pair_lookup.get((left_identifier, right_identifier))
            if pair is None:
                pair = pair_lookup.get((right_identifier, left_identifier))
            if pair is None:
                raise ValueError(
                    "distance parity rows require a complete symmetric pair lookup"
                )
            if pair.distance is None:
                if (
                    pair.saturation_reason
                    and "tends to infinity" in pair.saturation_reason
                ):
                    distance = ""
                    distance_status = "infinite"
                    infinite_distance_count += 1
                else:
                    distance = ""
                    distance_status = "undefined"
                    undefined_distance_count += 1
            else:
                distance = pair.distance
                distance_status = "finite"
                finite_distance_count += 1
            rows.append(
                {
                    "left_identifier": left_identifier,
                    "right_identifier": right_identifier,
                    "distance": distance,
                    "distance_status": distance_status,
                }
            )
    return {
        "sequence_count": len(report.identifiers),
        "alignment_length": report.alignment_length,
        "pairwise_deletion": pairwise_deletion,
        "distance_model": report.model,
        "finite_distance_count": finite_distance_count,
        "undefined_distance_count": undefined_distance_count,
        "infinite_distance_count": infinite_distance_count,
    }, rows


def _build_bijux_tree_tip_distance_rows(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_tree_tip_distance_matrix(input_fixture)
    return {
        "tip_count": len(report.identifiers),
        "rooted": report.rooted,
        "tip_labels": report.identifiers,
        "pair_count": report.pair_count,
        "diagonal_zero": report.diagonal_zero,
        "symmetric": report.symmetric,
        "complete_branch_lengths": report.complete_branch_lengths,
        "missing_branch_length_policy": report.missing_branch_length_policy,
    }, [
        {
            "left_identifier": row.left_identifier,
            "right_identifier": row.right_identifier,
            "distance": row.distance,
        }
        for row in report.pairs
    ]


def _build_bijux_neighbor_joining_structure(
    input_fixture: Path,
) -> tuple[dict[str, object], None, str]:
    tree, _report = build_tree_from_imported_distance_matrix(
        input_fixture,
        method="neighbor-joining",
    )
    summary, _rows, normalized_text = _tree_structure_payload(tree, False, [])
    return summary, None, normalized_text


def _build_bijux_consensus_rows(
    input_fixture: Path,
    *,
    consensus_method: str,
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    if consensus_method == "strict":
        tree, report = compute_strict_consensus_tree(input_fixture)
    elif consensus_method == "majority-rule":
        tree, report = compute_consensus_tree(input_fixture)
    else:
        raise ValueError(f"unsupported consensus method '{consensus_method}'")
    tree.rooted = False
    frequency_report = compute_clade_frequency_table(input_fixture)
    included_split_count = len(informative_unrooted_splits(tree, set(tree.tip_names)))
    return (
        {
            "tree_count": report.tree_count,
            "shared_taxa": report.shared_taxa,
            "shared_taxon_count": len(report.shared_taxa),
            "tip_count": len(tree.tip_names),
            "rooted": False,
            "consensus_method": report.consensus_method,
            "consensus_threshold": report.consensus_threshold,
            "included_clade_count": included_split_count,
            "clade_frequency_count": len(frequency_report.clade_frequencies),
        },
        [
            {
                "clade": row.clade,
                "tree_count": row.tree_count,
                "frequency": row.frequency,
            }
            for row in frequency_report.clade_frequencies
        ],
        dumps_newick(tree),
    )


def _build_bijux_prop_clades_rows(
    reference_tree_path: Path,
    comparison_tree_set_path: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_reference_tree_clade_support(
        reference_tree_path,
        comparison_tree_set_path,
    )
    return {
        "tree_count": report.tree_count,
        "shared_taxa": report.shared_taxa,
        "shared_taxon_count": len(report.shared_taxa),
        "internal_node_count": len(report.rows),
        "supported_clade_count": report.supported_clade_count,
        "absent_clade_count": report.absent_clade_count,
        "unscored_clade_count": report.unscored_clade_count,
    }, [
        {
            "node_id": row.node_id,
            "node_kind": row.node_kind,
            "node_label": "" if row.node_label is None else row.node_label,
            "descendant_taxa": "|".join(row.descendant_taxa),
            "supporting_tree_count": (
                "" if row.supporting_tree_count is None else row.supporting_tree_count
            ),
            "clade_frequency": ""
            if row.clade_frequency is None
            else row.clade_frequency,
            "support_percent": ""
            if row.support_percent is None
            else row.support_percent,
            "support_status": row.support_status,
            "explanation": row.explanation,
            "reference_branch_length": (
                ""
                if row.reference_branch_length is None
                else row.reference_branch_length
            ),
            "reference_root_depth": (
                "" if row.reference_root_depth is None else row.reference_root_depth
            ),
        }
        for row in report.rows
    ]


def _inspect_tree_set_rooted_flags(input_fixture: Path) -> tuple[bool, bool]:
    lines = [
        line.strip()
        for line in input_fixture.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(lines) != 2:
        raise ValueError(
            "ape topology-distance parity fixtures must contain exactly two trees"
        )
    rooted_flags: list[bool] = []
    with tempfile.TemporaryDirectory(
        prefix="bijux-topology-distance-rootedness-"
    ) as tmpdir:
        temporary_root = Path(tmpdir)
        for index, line in enumerate(lines, start=1):
            temporary_path = temporary_root / f"tree-{index}.nwk"
            temporary_path.write_text(f"{line}\n", encoding="utf-8")
            rooted_flags.append(inspect_tree_path(temporary_path).rooted)
    return rooted_flags[0], rooted_flags[1]


def _build_bijux_topology_distance_rows(
    input_fixture: Path,
    *,
    rf_mode: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    tree_set = load_newick_tree_set(input_fixture)
    if len(tree_set) != 2:
        raise ValueError(
            "ape topology-distance parity fixtures must contain exactly two trees"
        )
    report = compare_topology_distance_trees(
        tree_set[0],
        tree_set[1],
        left_path=input_fixture,
        right_path=input_fixture,
        rf_mode=rf_mode,
        taxon_overlap_policy="require-identical",
    )
    rooted_left, rooted_right = _inspect_tree_set_rooted_flags(input_fixture)
    return {
        "tip_count": len(report.shared_taxa),
        "shared_taxa": report.shared_taxa,
        "left_only_taxa": report.left_only_taxa,
        "right_only_taxa": report.right_only_taxa,
        "taxon_overlap_policy": report.taxon_overlap_policy,
        "rf_mode": report.rf_mode,
        "rooted_left": rooted_left,
        "rooted_right": rooted_right,
        "polytomy_present_left": report.polytomy_present_left,
        "polytomy_present_right": report.polytomy_present_right,
        "left_split_count": report.left_split_count,
        "right_split_count": report.right_split_count,
        "shared_split_count": report.shared_split_count,
        "left_only_split_count": report.left_only_split_count,
        "right_only_split_count": report.right_only_split_count,
        "robinson_foulds_distance": report.robinson_foulds_distance,
        "normalized_robinson_foulds": report.normalized_robinson_foulds,
        "topology_equal": report.topology_equal,
    }, [
        {
            "split_id": row.split_id,
            "split_kind": row.split_kind,
            "comparison_status": row.comparison_status,
            "taxon_count": row.taxon_count,
            "descendant_taxa": "|".join(row.descendant_taxa),
            "left_present": row.left_present,
            "right_present": row.right_present,
        }
        for row in report.split_rows
    ]


def _build_bijux_brownian_covariance_rows(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = summarize_brownian_covariance(input_fixture)
    return {
        "tip_count": len(report.taxa),
        "rooted": report.tree_is_rooted,
        "tip_labels": report.taxa,
        "pair_count": len(report.rows),
        "tree_is_ultrametric": report.tree_is_ultrametric,
        "minimum_root_to_tip_depth": report.minimum_root_to_tip_depth,
        "maximum_root_to_tip_depth": report.maximum_root_to_tip_depth,
        "minimum_branch_length": report.minimum_branch_length,
        "maximum_branch_length": report.maximum_branch_length,
        "matrix_dimension": report.matrix_dimension,
        "matrix_rank": report.matrix_rank,
        "singular": report.singular,
        "near_singular": report.near_singular,
        "positive_definite": report.positive_definite,
        "condition_number": (
            None if math.isinf(report.condition_number) else report.condition_number
        ),
        "raw_log_determinant": report.raw_log_determinant,
    }, [
        {
            "left_taxon": row.left_taxon,
            "right_taxon": row.right_taxon,
            "shared_ancestry_covariance": row.shared_ancestry_covariance,
        }
        for row in report.rows
    ]


def _build_bijux_independent_contrast_rows(
    input_fixture: Path,
    *,
    trait_table_path: Path,
    trait_name: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_phylogenetic_independent_contrasts(
        input_fixture,
        trait_table_path,
        trait=trait_name,
    )
    rows = sorted(
        [
            {
                "node_id": row.node_id,
                "node": row.node,
                "left_taxa": "|".join(row.left_taxa),
                "right_taxa": "|".join(row.right_taxa),
                "contrast": row.contrast,
                "expected_variance": row.expected_variance,
            }
            for row in report.contrasts
        ],
        key=lambda row: int(row["node_id"]),
    )
    return {
        "trait": report.trait,
        "taxon_count": report.taxon_count,
        "contrast_count": len(report.contrasts),
        "tree_is_ultrametric": report.input_audit.tree_is_ultrametric,
        "minimum_root_to_tip_depth": report.input_audit.minimum_root_to_tip_depth,
        "maximum_root_to_tip_depth": report.input_audit.maximum_root_to_tip_depth,
    }, rows


def _build_bijux_continuous_ancestral_rows(
    input_fixture: Path,
    *,
    trait_table_path: Path,
    trait_name: str,
    trait_taxon_column: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    dataset = load_continuous_dataset(
        input_fixture,
        trait_table_path,
        trait=trait_name,
        taxon_column=trait_taxon_column,
    )
    report = reconstruct_continuous_ancestral_states(
        input_fixture,
        trait_table_path,
        trait=trait_name,
        taxon_column=trait_taxon_column,
        model="brownian",
    )
    internal_node_map = {
        node_signature(node): node_id
        for node_id, node in build_ape_internal_node_map(dataset.tree).items()
    }
    rows = sorted(
        [
            {
                "node_id": internal_node_map[estimate.node],
                "node": estimate.node,
                "estimate": estimate.estimate,
                "standard_error": estimate.standard_error,
                "lower_95_interval": estimate.lower_95_interval,
                "upper_95_interval": estimate.upper_95_interval,
            }
            for estimate in report.estimates
            if not estimate.is_tip
        ],
        key=lambda row: int(row["node_id"]),
    )
    diagnostics = report.brownian_fit_diagnostics
    return {
        "trait": report.trait,
        "taxon_count": report.taxon_count,
        "excluded_taxon_count": len(report.dropped_missing_taxa)
        + len(report.dropped_non_numeric_taxa),
        "dropped_missing_taxa": report.dropped_missing_taxa,
        "dropped_non_numeric_taxa": report.dropped_non_numeric_taxa,
        "internal_node_count": len(rows),
        "method": "pic",
        "tree_is_ultrametric": (
            None if diagnostics is None else diagnostics.tree_is_ultrametric
        ),
        "minimum_root_to_tip_depth": (
            None if diagnostics is None else diagnostics.minimum_root_to_tip_depth
        ),
        "maximum_root_to_tip_depth": (
            None if diagnostics is None else diagnostics.maximum_root_to_tip_depth
        ),
    }, rows


def _build_bijux_discrete_ancestral_rows(
    input_fixture: Path,
    *,
    trait_table_path: Path,
    trait_name: str,
    trait_taxon_column: str,
    ancestral_model: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    dataset = load_discrete_dataset(
        input_fixture,
        trait_table_path,
        trait=trait_name,
        taxon_column=trait_taxon_column,
    )
    report = reconstruct_discrete_ancestral_states(
        input_fixture,
        trait_table_path,
        trait=trait_name,
        taxon_column=trait_taxon_column,
        model=ancestral_model,
    )
    internal_node_map = {
        node_signature(node): node_id
        for node_id, node in build_ape_internal_node_map(dataset.tree).items()
    }
    rows = sorted(
        [
            {
                "node_id": internal_node_map[estimate.node],
                "node": estimate.node,
                "state": _coerce_table_cell(state),
                "posterior_probability": probability,
                "most_likely_state": _coerce_table_cell(estimate.most_likely_state),
                "max_posterior_probability": estimate.confidence,
            }
            for estimate in report.estimates
            if not estimate.is_tip
            for state, probability in sorted(estimate.state_probabilities.items())
        ],
        key=lambda row: (int(row["node_id"]), str(row["state"])),
    )
    transition_rows = [
        {
            "source_state": row.source_state,
            "target_state": row.target_state,
            "transition_allowed": row.transition_allowed,
            "step_distance": row.step_distance,
            "rate": row.rate,
        }
        for row in report.transition_rate_rows
    ]
    return {
        "trait": report.trait,
        "taxon_count": report.taxon_count,
        "excluded_taxon_count": len(report.dropped_missing_taxa),
        "dropped_missing_taxa": report.dropped_missing_taxa,
        "internal_node_count": len(
            [estimate for estimate in report.estimates if not estimate.is_tip]
        ),
        "model": report.model,
        "state_count": len(report.observed_states),
        "state_labels": report.observed_states,
        "log_likelihood": report.log_likelihood,
        "parameter_count": report.parameter_count,
        "aic": report.aic,
        "overparameterized": report.overparameterized,
        "baseline_model": (
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.baseline_model
        ),
        "baseline_delta_aic": (
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.delta_aic
        ),
        "preferred_model_by_aic": (
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.preferred_model_by_aic
        ),
        "transition_rate_rows": transition_rows,
    }, rows


def _build_bijux_tree_node_depth_rows(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_tree_node_depths(input_fixture)
    return {
        "node_count": report.node_count,
        "tip_count": report.tip_count,
        "internal_node_count": report.internal_node_count,
        "rooted": report.rooted,
        "tip_labels": report.tip_labels,
        "tree_is_ultrametric": report.tree_is_ultrametric,
        "zero_branch_length_count": report.zero_branch_length_count,
        "minimum_tip_depth": report.minimum_tip_depth,
        "maximum_tip_depth": report.maximum_tip_depth,
        "minimum_internal_depth": report.minimum_internal_depth,
        "maximum_internal_depth": report.maximum_internal_depth,
    }, [
        {
            "node_id": row.node_id,
            "node_kind": row.node_kind,
            "node_label": row.node_label or "",
            "descendant_taxa": "|".join(row.descendant_taxa),
            "branch_length_depth": row.branch_length_depth,
            "branch_length": "" if row.branch_length is None else row.branch_length,
        }
        for row in report.rows
    ]


def _build_bijux_branching_time_rows(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_tree_branching_times(input_fixture)
    return {
        "internal_node_count": report.internal_node_count,
        "rooted": report.rooted,
        "tip_labels": report.tip_labels,
        "tree_is_ultrametric": report.tree_is_ultrametric,
        "root_age": report.root_age,
        "zero_branch_length_count": report.zero_branch_length_count,
        "minimum_tip_depth": report.minimum_tip_depth,
        "maximum_tip_depth": report.maximum_tip_depth,
        "max_tip_depth_deviation": report.max_tip_depth_deviation,
        "tolerance": report.tolerance,
    }, [
        {
            "node_id": row.node_id,
            "node_kind": row.node_kind,
            "node_label": row.node_label or "",
            "descendant_taxa": "|".join(row.descendant_taxa),
            "node_depth": row.node_depth,
            "branching_time": row.branching_time,
        }
        for row in report.rows
    ]


def _build_bijux_diversification_gamma_rows(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_diversification_gamma_statistic(input_fixture)
    return {
        "tip_count": report.tip_count,
        "rooted": report.rooted,
        "ultrametric": report.ultrametric,
        "bifurcating": report.bifurcating,
        "root_age": report.root_age,
        "branching_time_count": report.branching_time_count,
        "interval_count": report.interval_count,
        "minimum_branching_time": report.minimum_branching_time,
        "maximum_branching_time": report.maximum_branching_time,
        "gamma_statistic": report.gamma_statistic,
    }, [
        {
            "tip_count": report.tip_count,
            "rooted": report.rooted,
            "ultrametric": report.ultrametric,
            "bifurcating": report.bifurcating,
            "root_age": report.root_age,
            "branching_time_count": report.branching_time_count,
            "interval_count": report.interval_count,
            "minimum_branching_time": report.minimum_branching_time,
            "maximum_branching_time": report.maximum_branching_time,
            "gamma_statistic": report.gamma_statistic,
        }
    ]


def _build_bijux_tree_simulation_envelope_rows(
    fixture_id: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    fixture = get_shared_tree_simulation_fixture(fixture_id)
    if fixture.simulation_model == "random-tree":
        _trees, report = simulate_random_trees(
            tree_count=fixture.replicate_count,
            tip_count=fixture.tip_count,
            seed=fixture.seed,
            branch_length_model=fixture.branch_length_model or "uniform",
        )
    elif fixture.simulation_model == "coalescent":
        _trees, report = simulate_coalescent_trees(
            tree_count=fixture.replicate_count,
            tip_count=fixture.tip_count,
            population_size=fixture.population_size or 1.0,
            seed=fixture.seed,
        )
    else:
        raise ValueError(
            f"unsupported governed simulation model {fixture.simulation_model!r}"
        )
    return {
        "simulation_model": report.model,
        "reference_function": fixture.reference_function,
        "tree_count": report.tree_count,
        "tip_count": report.tip_count,
        "seed": report.seed,
        "branch_length_model": report.branch_length_model,
        "population_size": report.population_size,
        "rooted": report.rooted,
        "binary": report.binary,
        "pooled_branch_count": report.pooled_branch_count,
        "envelope_metric_count": len(report.envelope_metrics),
    }, [
        {
            "metric": row.metric,
            "sample_scope": row.sample_scope,
            "observation_count": row.observation_count,
            "mean": row.mean,
            "standard_deviation": row.standard_deviation,
            "minimum": row.minimum,
            "median": row.median,
            "maximum": row.maximum,
        }
        for row in report.envelope_metrics
    ]


def _build_bijux_tree_ultrametric_rows(
    input_fixture: Path,
    *,
    tolerance: float,
    option: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = assess_tree_ultrametricity(
        input_fixture,
        tolerance=tolerance,
        option=option,
    )
    return {
        "tip_count": len(report.tip_labels),
        "rooted": report.rooted,
        "tip_labels": report.tip_labels,
        "ultrametric": report.ultrametric,
        "criterion_name": report.criterion_name,
        "criterion_value": report.criterion_value,
        "tolerance": report.tolerance,
        "option": report.option,
        "minimum_tip_depth": report.minimum_tip_depth,
        "maximum_tip_depth": report.maximum_tip_depth,
        "mean_tip_depth": report.mean_tip_depth,
        "max_tip_depth_deviation": report.max_tip_depth_deviation,
        "root_age": report.root_age,
        "offending_taxa": report.offending_taxa,
    }, [
        {
            "node_id": row.node_id,
            "tip_label": row.tip_label,
            "root_to_tip_depth": row.root_to_tip_depth,
            "deviation_from_mean_depth": row.deviation_from_mean_depth,
            "deviation_from_min_depth": row.deviation_from_min_depth,
            "deviation_from_max_depth": row.deviation_from_max_depth,
            "is_offending_taxon": row.is_offending_taxon,
        }
        for row in report.rows
    ]


def _build_bijux_translation_rows(
    input_fixture: Path,
    *,
    genetic_code_id: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    translated, report = translate_coding_alignment(
        input_fixture, genetic_code=genetic_code_id
    )
    return {
        "sequence_count": report.translated_sequence_count,
        "translated_length": report.translated_alignment_length,
        "stop_codon_count": report.stop_codon_count,
        "dropped_trailing_nucleotide_count": report.dropped_trailing_nucleotide_count,
        "warning_count": len(report.warnings),
        "warnings": report.warnings,
    }, [
        {
            "identifier": row.identifier,
            "amino_acid_sequence": row.sequence,
        }
        for row in translated
    ]


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _normalize_newick_label(label: str) -> str:
    if len(label) >= 2 and label.startswith("'") and label.endswith("'"):
        return label[1:-1].replace("''", "'")
    return label


def _normalize_expected_label(
    label: str,
    *,
    expected_tip_labels: set[str] | None,
) -> str:
    normalized = _normalize_newick_label(label)
    if (
        expected_tip_labels
        and normalized not in expected_tip_labels
        and normalized.replace("_", " ") in expected_tip_labels
    ):
        return normalized.replace("_", " ")
    return normalized


def _normalize_joined_labels(
    value: str,
    *,
    expected_tip_labels: set[str] | None,
) -> str:
    if value == "":
        return value
    labels = [
        _normalize_expected_label(label, expected_tip_labels=expected_tip_labels)
        for label in value.split("|")
    ]
    return "|".join(sorted(labels))


def _normalize_reference_summary(summary: dict[str, object]) -> dict[str, object]:
    normalized = dict(summary)
    tip_labels = normalized.get("tip_labels")
    if isinstance(tip_labels, list):
        expected_tip_labels = {
            _normalize_newick_label(str(label)) for label in tip_labels
        }
        normalized["tip_labels"] = [
            _normalize_expected_label(
                str(label),
                expected_tip_labels=expected_tip_labels,
            )
            for label in tip_labels
        ]
    return normalized


def _summary_rooted_flag(summary: dict[str, object]) -> bool:
    rooted = summary.get("rooted")
    if isinstance(rooted, bool):
        return rooted
    raise ValueError("reference summary must include a boolean rooted flag")


def _optional_payload_string(payload: dict[str, object], key: str) -> str | None:
    value = payload.get(key)
    return value if isinstance(value, str) else None


def _coerce_table_cell(value: str) -> object:
    if value == "":
        return ""
    if value.lower() in {"true", "false"}:
        return value.lower() == "true"
    if re.fullmatch(r"-?\d+", value):
        return int(value)
    if re.fullmatch(r"-?(?:\d+\.\d*|\d*\.\d+)(?:[eE][+-]?\d+)?", value):
        return float(value)
    return value


def _sort_parity_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return sorted(
        rows,
        key=lambda row: (
            0 if row["tree_index"] == "" else int(row["tree_index"]),
            _node_kind_order(str(row["node_kind"])),
            str(row["clade_id"]),
            str(row["node_label"]),
        ),
    )


def _load_rows_table(
    path: Path,
    *,
    expected_tip_labels: set[str] | None = None,
    sort_rows: bool = False,
) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    normalized_rows = [
        {
            key: _normalize_expected_label(
                value, expected_tip_labels=expected_tip_labels
            )
            if key.endswith("label")
            else _normalize_joined_labels(
                value,
                expected_tip_labels=expected_tip_labels,
            )
            if key in {"clade_id", "taxa", "descendant_taxa", "shared_taxa"}
            else _coerce_table_cell(value)
            for key, value in row.items()
        }
        for row in rows
    ]
    if sort_rows:
        return _sort_parity_rows(normalized_rows)
    return normalized_rows


def _compare_scalar(expected: object, observed: object, *, tolerance: float) -> bool:
    if isinstance(expected, (int, float)) and isinstance(observed, (int, float)):
        return abs(float(expected) - float(observed)) <= tolerance
    return expected == observed


def _compare_json(expected: object, observed: object, *, tolerance: float) -> bool:
    if isinstance(expected, dict) and isinstance(observed, dict):
        if set(expected) != set(observed):
            return False
        return all(
            _compare_json(expected[key], observed[key], tolerance=tolerance)
            for key in expected
        )
    if isinstance(expected, list) and isinstance(observed, list):
        if len(expected) != len(observed):
            return False
        return all(
            _compare_json(left, right, tolerance=tolerance)
            for left, right in zip(expected, observed, strict=True)
        )
    return _compare_scalar(expected, observed, tolerance=tolerance)


def _normalize_tree_labels(
    node: TreeNode,
    *,
    expected_tip_labels: set[str] | None,
) -> None:
    if node.name is not None:
        normalized = _normalize_newick_label(node.name)
        if (
            expected_tip_labels
            and normalized not in expected_tip_labels
            and normalized.replace("_", " ") in expected_tip_labels
        ):
            normalized = normalized.replace("_", " ")
        node.name = normalized
    for child in node.children:
        _normalize_tree_labels(child, expected_tip_labels=expected_tip_labels)


def _clear_branch_lengths(node: TreeNode) -> None:
    node.branch_length = None
    for child in node.children:
        _clear_branch_lengths(child)


def _canonical_newick(
    path: Path,
    *,
    expected_tip_labels: set[str] | None = None,
) -> str:
    tree = load_tree(path)
    _normalize_tree_labels(tree.root, expected_tip_labels=expected_tip_labels)
    return dumps_newick(tree)


def _copy_if_exists(source: Path, destination: Path) -> None:
    if source.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )


def _write_rows_table(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def _persist_failure_bundle(
    *,
    failure_root: Path,
    case: ApeParityCase,
    case_file: Path,
    execution_root: Path,
    execution_payload: dict[str, object] | None,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
    reference_error: dict[str, object] | None,
    bijux_error: dict[str, object] | None,
    reference_rows: list[dict[str, object]] | None,
    bijux_rows: list[dict[str, object]] | None,
    bijux_normalized_text: str | None,
    mismatch_reason: str,
) -> Path:
    artifact_root = failure_root / case.case_id
    artifact_root.mkdir(parents=True, exist_ok=True)
    _copy_if_exists(case_file, artifact_root / "case.json")
    _copy_if_exists(
        execution_root.parent / "bijux-reference-input.nwk",
        artifact_root / "bijux-reference-input.nwk",
    )
    _copy_if_exists(
        execution_root / "reference-execution.json",
        artifact_root / "reference-execution.json",
    )
    if execution_payload is not None:
        outputs = execution_payload.get("outputs")
        if isinstance(outputs, dict):
            for path_string in outputs.values():
                if isinstance(path_string, str):
                    source = Path(path_string)
                    _copy_if_exists(source, artifact_root / f"reference-{source.name}")
    if execution_payload is not None:
        _write_json(
            artifact_root / "reference-execution.observed.json", execution_payload
        )
    if reference_summary is not None:
        _write_json(
            artifact_root / "reference-summary.observed.json", reference_summary
        )
    if reference_rows is not None:
        _write_json(artifact_root / "reference-rows.observed.json", reference_rows)
        _write_rows_table(artifact_root / "reference-rows.observed.tsv", reference_rows)
    if bijux_summary is not None:
        _write_json(artifact_root / "bijux-summary.json", bijux_summary)
    if reference_error is not None:
        _write_json(artifact_root / "reference-error.observed.json", reference_error)
    if bijux_error is not None:
        _write_json(artifact_root / "bijux-error.json", bijux_error)
    if bijux_rows is not None:
        _write_json(artifact_root / "bijux-rows.json", bijux_rows)
        _write_rows_table(artifact_root / "bijux-rows.tsv", bijux_rows)
    if bijux_normalized_text is not None:
        (artifact_root / "bijux-normalized.txt").write_text(
            f"{bijux_normalized_text}\n",
            encoding="utf-8",
        )
    _write_json(
        artifact_root / "comparison.json",
        {
            "case_id": case.case_id,
            "function_name": case.function_name,
            "mismatch_reason": mismatch_reason,
        },
    )
    return artifact_root


def _build_bijux_case_payload(
    case: ApeParityCase,
) -> tuple[dict[str, object], list[dict[str, object]] | None, str | None]:
    if case.operation in {"read-tree-structure", "write-tree-structure"}:
        summary, rows, normalized_text = _build_bijux_tree_structure(case.input_fixture)
        return summary, rows, normalized_text
    if case.operation == "root-tree-outgroup":
        summary, rows, normalized_text = _build_bijux_root_outgroup_structure(
            case.input_fixture,
            outgroup_taxa=case.outgroup_taxa,
        )
        return summary, rows, normalized_text
    if case.operation == "unroot-tree":
        summary, rows, normalized_text = _build_bijux_unroot_structure(
            case.input_fixture,
        )
        return summary, rows, normalized_text
    if case.operation == "drop-tree-taxa":
        summary, rows, normalized_text = _build_bijux_drop_tip_structure(
            case.input_fixture,
            excluded_taxa=case.excluded_taxa,
        )
        return summary, rows, normalized_text
    if case.operation == "keep-tree-taxa":
        summary, rows, normalized_text = _build_bijux_keep_tip_structure(
            case.input_fixture,
            requested_taxa=case.requested_taxa,
        )
        return summary, rows, normalized_text
    if case.operation == "extract-tree-clade":
        if case.node_id is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing an extraction node id"
            )
        summary, rows, normalized_text = _build_bijux_extract_clade_structure(
            case.input_fixture,
            node_id=case.node_id,
        )
        return summary, rows, normalized_text
    if case.operation == "get-tree-mrca":
        if not case.mrca_taxa:
            raise ValueError(f"ape parity case '{case.case_id}' is missing MRCA taxa")
        summary = _build_bijux_mrca_summary(
            case.input_fixture,
            mrca_taxa=case.mrca_taxa,
        )
        return summary, None, None
    if case.operation == "assess-tree-monophyly":
        if case.monophyly_reroot is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a monophyly reroot policy"
            )
        summary = _build_bijux_monophyly_summary(
            case.input_fixture,
            requested_taxa=case.requested_taxa,
            reroot=case.monophyly_reroot,
        )
        return summary, None, None
    if case.operation in {"read-tree-set-structure", "write-tree-set-structure"}:
        summary, rows, normalized_text = _build_bijux_tree_set_structure(
            case.input_fixture
        )
        return summary, rows, normalized_text
    if case.operation == "tree-consensus":
        if case.consensus_method is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a consensus method"
            )
        summary, rows, normalized_text = _build_bijux_consensus_rows(
            case.input_fixture,
            consensus_method=case.consensus_method,
        )
        return summary, rows, normalized_text
    if case.operation == "tree-clade-support":
        if case.reference_tree_path is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a reference tree path"
            )
        summary, rows = _build_bijux_prop_clades_rows(
            case.reference_tree_path,
            case.input_fixture,
        )
        return summary, rows, None
    if case.operation == "tree-tip-distance":
        summary, rows = _build_bijux_tree_tip_distance_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "distance-matrix-neighbor-joining":
        return _build_bijux_neighbor_joining_structure(case.input_fixture)
    if case.operation == "tree-topology-distance":
        if case.rf_mode is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a topology rf mode"
            )
        summary, rows = _build_bijux_topology_distance_rows(
            case.input_fixture,
            rf_mode=case.rf_mode,
        )
        return summary, rows, None
    if case.operation == "tree-brownian-covariance":
        summary, rows = _build_bijux_brownian_covariance_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "tree-continuous-ancestral-states":
        if (
            case.trait_table_path is None
            or case.trait_name is None
            or case.trait_taxon_column is None
        ):
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a trait table path, trait name, or taxon column"
            )
        summary, rows = _build_bijux_continuous_ancestral_rows(
            case.input_fixture,
            trait_table_path=case.trait_table_path,
            trait_name=case.trait_name,
            trait_taxon_column=case.trait_taxon_column,
        )
        return summary, rows, None
    if case.operation == "tree-discrete-ancestral-states":
        if (
            case.trait_table_path is None
            or case.trait_name is None
            or case.trait_taxon_column is None
        ):
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a trait table path, trait name, or taxon column"
            )
        summary, rows = _build_bijux_discrete_ancestral_rows(
            case.input_fixture,
            trait_table_path=case.trait_table_path,
            trait_name=case.trait_name,
            trait_taxon_column=case.trait_taxon_column,
            ancestral_model=case.ancestral_model or "equal-rates",
        )
        return summary, rows, None
    if case.operation == "tree-independent-contrasts":
        if case.trait_table_path is None or case.trait_name is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a trait table path or trait name"
            )
        summary, rows = _build_bijux_independent_contrast_rows(
            case.input_fixture,
            trait_table_path=case.trait_table_path,
            trait_name=case.trait_name,
        )
        return summary, rows, None
    if case.operation == "tree-node-depth":
        summary, rows = _build_bijux_tree_node_depth_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "tree-branching-times":
        summary, rows = _build_bijux_branching_time_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "tree-diversification-gamma-statistic":
        summary, rows = _build_bijux_diversification_gamma_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "tree-simulation-envelope":
        summary, rows = _build_bijux_tree_simulation_envelope_rows(case.fixture_id)
        return summary, rows, None
    if case.operation == "tree-ultrametricity":
        if case.ultrametric_option is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing an ultrametric option"
            )
        summary, rows = _build_bijux_tree_ultrametric_rows(
            case.input_fixture,
            tolerance=case.tolerance,
            option=case.ultrametric_option,
        )
        return summary, rows, None
    if case.operation == "dna-dnabin-structure":
        summary, rows = _build_bijux_dnabin_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "dna-base-frequency":
        summary, rows = _build_bijux_base_frequency_summary(case.input_fixture)
        return summary, rows, None
    if case.operation == "dna-segregating-sites":
        summary, rows = _build_bijux_segregating_site_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "dna-distance":
        if case.pairwise_deletion is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing pairwise deletion policy"
            )
        if case.distance_model is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a distance model"
            )
        summary, rows = _build_bijux_distance_rows(
            case.input_fixture,
            pairwise_deletion=case.pairwise_deletion,
            distance_model=case.distance_model,
        )
        return summary, rows, None
    if case.operation == "dna-translation":
        if case.genetic_code_id is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a genetic code id"
            )
        summary, rows = _build_bijux_translation_rows(
            case.input_fixture,
            genetic_code_id=case.genetic_code_id,
        )
        return summary, rows, None
    raise ValueError(f"unsupported ape parity operation '{case.operation}'")


def _load_reference_case_payload(
    case: ApeParityCase,
    execution_root: Path,
) -> tuple[dict[str, object], list[dict[str, object]] | None, str | None]:
    if case.operation in {
        "read-tree-structure",
        "write-tree-structure",
        "root-tree-outgroup",
        "unroot-tree",
        "drop-tree-taxa",
        "keep-tree-taxa",
        "extract-tree-clade",
    }:
        summary = _normalize_reference_summary(
            _load_json(execution_root / "summary.json")
        )
        expected_tip_labels = {str(label) for label in summary.get("tip_labels", [])}
        rows = _load_rows_table(
            execution_root / "clades.tsv",
            expected_tip_labels=expected_tip_labels,
            sort_rows=True,
        )
        normalized_text = _canonical_newick(
            execution_root / "normalized-tree.nwk",
            expected_tip_labels=expected_tip_labels,
        )
        return summary, rows, normalized_text
    if case.operation == "get-tree-mrca":
        summary = _normalize_reference_summary(
            _load_json(execution_root / "summary.json")
        )
        return summary, None, None
    if case.operation == "assess-tree-monophyly":
        summary = _normalize_reference_summary(
            _load_json(execution_root / "summary.json")
        )
        return summary, None, None
    if case.operation in {"read-tree-set-structure", "write-tree-set-structure"}:
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "clades.tsv", sort_rows=True)
        return summary, rows, None
    if case.operation == "tree-consensus":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "clade-frequencies.tsv")
        normalized_text = _canonical_newick(execution_root / "normalized-tree.nwk")
        return summary, rows, normalized_text
    if case.operation == "tree-clade-support":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "support-table.tsv")
        return summary, rows, None
    if case.operation == "tree-tip-distance":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "tip-distance-long.tsv")
        return summary, rows, None
    if case.operation == "distance-matrix-neighbor-joining":
        summary = _normalize_reference_summary(
            _load_json(execution_root / "summary.json")
        )
        normalized_text = _canonical_newick(execution_root / "normalized-tree.nwk")
        return summary, None, normalized_text
    if case.operation == "tree-topology-distance":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "split-table.tsv")
        return summary, rows, None
    if case.operation == "tree-brownian-covariance":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "covariance-long.tsv")
        return summary, rows, None
    if case.operation == "tree-continuous-ancestral-states":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "continuous-ancestral.tsv")
        return summary, rows, None
    if case.operation == "tree-discrete-ancestral-states":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "discrete-ancestral.tsv")
        return summary, rows, None
    if case.operation == "tree-independent-contrasts":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "independent-contrasts.tsv")
        return summary, rows, None
    if case.operation == "tree-node-depth":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "node-depths.tsv")
        return summary, rows, None
    if case.operation == "tree-branching-times":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "branching-times.tsv")
        return summary, rows, None
    if case.operation == "tree-diversification-gamma-statistic":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "gamma-statistic.tsv")
        return summary, rows, None
    if case.operation == "tree-simulation-envelope":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "simulation-envelope.tsv")
        return summary, rows, None
    if case.operation == "tree-ultrametricity":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "ultrametric-diagnostics.tsv")
        return summary, rows, None
    if case.operation == "dna-dnabin-structure":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "dnabin.tsv")
        return summary, rows, None
    if case.operation == "dna-base-frequency":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "base-frequency.tsv")
        return summary, rows, None
    if case.operation == "dna-segregating-sites":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "segregating-sites.tsv")
        return summary, rows, None
    if case.operation == "dna-distance":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "distance-matrix.tsv")
        return summary, rows, None
    if case.operation == "dna-translation":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "translation.tsv")
        return summary, rows, None
    raise ValueError(f"unsupported ape parity operation '{case.operation}'")


def _tree_structure_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    expected_tree = load_tree(case.input_fixture)
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _tree_set_structure_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree_set = load_newick_tree_set(
        execution_root / "normalized-tree-set.nwk"
    )
    expected_tree_set = load_newick_tree_set(case.input_fixture)
    expected_tip_labels = {
        tip_name for tree in expected_tree_set for tip_name in tree.tip_names
    }
    for tree in reference_tree_set:
        _normalize_tree_labels(tree.root, expected_tip_labels=expected_tip_labels)
    report = compare_tree_sets_structurally(
        expected_tree_set,
        reference_tree_set,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _root_tree_outgroup_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    # Canonical Newick does not preserve rootedness metadata, but ape::root
    # produced this record explicitly as a rooted output for these governed cases.
    reference_tree.rooted = True
    expected_tree, _report = root_tree_on_outgroup(
        case.input_fixture,
        outgroup_taxa=list(case.outgroup_taxa),
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _unroot_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = False
    expected_tree, _report = unroot_tree(case.input_fixture)
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _drop_tip_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_summary = _normalize_reference_summary(
        _load_json(execution_root / "summary.json")
    )
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = _summary_rooted_flag(reference_summary)
    expected_tree, _report = drop_tree_taxa(
        case.input_fixture, list(case.excluded_taxa)
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _keep_tip_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_summary = _normalize_reference_summary(
        _load_json(execution_root / "summary.json")
    )
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = _summary_rooted_flag(reference_summary)
    expected_tree, _report = prune_tree_to_requested_taxa(
        case.input_fixture,
        list(case.requested_taxa),
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _extract_clade_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_summary = _normalize_reference_summary(
        _load_json(execution_root / "summary.json")
    )
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = _summary_rooted_flag(reference_summary)
    if case.node_id is None:
        raise ValueError(
            f"ape parity case '{case.case_id}' is missing an extraction node id"
        )
    expected_tree, _report = extract_tree_clade_by_node_id(
        case.input_fixture,
        node_id=case.node_id,
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _consensus_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    if case.consensus_method == "strict":
        expected_tree, _report = compute_strict_consensus_tree(case.input_fixture)
    elif case.consensus_method == "majority-rule":
        expected_tree, _report = compute_consensus_tree(case.input_fixture)
    else:
        raise ValueError(
            f"ape parity case '{case.case_id}' has unsupported consensus method "
            f"{case.consensus_method!r}"
        )
    expected_tree.rooted = False
    _clear_branch_lengths(expected_tree.root)
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = False
    _clear_branch_lengths(reference_tree.root)
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=False,
    )
    return report.mismatch_reason


def _neighbor_joining_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = False
    expected_tree, _report = build_tree_from_imported_distance_matrix(
        case.input_fixture,
        method="neighbor-joining",
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=False,
    )
    return report.mismatch_reason


def _supports_ard_rate_multiset_equivalence(
    *,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
) -> bool:
    if reference_summary is None or bijux_summary is None:
        return False
    return (
        reference_summary.get("model") == "all-rates-different"
        and bijux_summary.get("model") == "all-rates-different"
        and reference_summary.get("overparameterized") is True
        and bijux_summary.get("overparameterized") is True
    )


def _group_transition_rate_rows(
    rows: list[dict[str, object]],
) -> dict[tuple[bool, int], list[float]]:
    grouped: dict[tuple[bool, int], list[float]] = {}
    for row in rows:
        grouped.setdefault(
            (
                bool(row["transition_allowed"]),
                int(row["step_distance"]),
            ),
            [],
        ).append(float(row["rate"]))
    for values in grouped.values():
        values.sort()
    return grouped


def _transition_rate_rows_match(
    *,
    reference_rows: list[dict[str, object]],
    bijux_rows: list[dict[str, object]],
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
    tolerance: float,
) -> bool:
    if _compare_json(reference_rows, bijux_rows, tolerance=tolerance):
        return True
    if not _supports_ard_rate_multiset_equivalence(
        reference_summary=reference_summary,
        bijux_summary=bijux_summary,
    ):
        return False
    return _compare_json(
        _group_transition_rate_rows(reference_rows),
        _group_transition_rate_rows(bijux_rows),
        tolerance=tolerance,
    )


def _summary_rows(
    observations: list[ApeParityObservation],
) -> list[ApeParitySummaryRow]:
    rows: list[ApeParitySummaryRow] = []
    for function_name in sorted({item.function_name for item in observations}):
        selected = [
            item for item in observations if item.function_name == function_name
        ]
        rows.append(
            ApeParitySummaryRow(
                function_name=function_name,
                case_count=len(selected),
                passed_case_count=sum(
                    1 for item in selected if item.status == "passed"
                ),
                failed_case_count=sum(
                    1 for item in selected if item.status == "failed"
                ),
                skipped_case_count=sum(
                    1 for item in selected if item.status == "skipped"
                ),
            )
        )
    return rows


def run_ape_parity_cases(
    *,
    case_ids: list[str] | None = None,
    rscript_executable: str = "Rscript",
    failure_root: Path | None = None,
    fixtures_root: Path | None = None,
) -> ApeParityReport:
    """Run governed live `ape` parity cases through the checked-in R runner."""
    selected = _selected_cases(case_ids=case_ids, fixtures_root=fixtures_root)
    observations: list[ApeParityObservation] = []
    active_failure_root = _failure_root() if failure_root is None else failure_root
    bijux_version = _bijux_version()
    bijux_commit = _bijux_commit()
    for case in selected:
        with tempfile.TemporaryDirectory(
            prefix=f"bijux-ape-parity-{case.case_id}-"
        ) as tmpdir:
            working_root = Path(tmpdir)
            reference_input_path = _materialize_reference_input(case, working_root)
            reference_case = replace(case, input_fixture=reference_input_path)
            case_file = _write_case_file(working_root / "case.json", reference_case)
            execution_root = working_root / "reference"
            execution_root.mkdir(parents=True, exist_ok=True)
            bijux_summary: dict[str, object] | None = None
            bijux_rows: list[dict[str, object]] | None = None
            bijux_normalized_text: str | None = None
            bijux_error: dict[str, object] | None = None
            try:
                (
                    bijux_summary,
                    bijux_rows,
                    bijux_normalized_text,
                ) = _build_bijux_case_payload(case)
            except Exception as error:
                bijux_error = {
                    "error_type": type(error).__name__,
                    "message": str(error),
                }
            execution_payload: dict[str, object] | None = None
            reference_summary: dict[str, object] | None = None
            reference_error: dict[str, object] | None = None
            reference_rows: list[dict[str, object]] | None = None
            reference_normalized_text: str | None = None
            status = "failed"
            mismatch_reason: str | None = None
            artifact_root: Path | None = None
            r_version: str | None = None
            ape_version: str | None = None
            process_stdout = ""
            process_stderr = ""
            try:
                # Repository-owned R parity runner.
                process = subprocess.run(  # nosec
                    [
                        rscript_executable,
                        str(_ape_runner_path()),
                        str(case_file),
                        str(execution_root),
                    ],
                    capture_output=True,
                    check=False,
                    cwd=_repository_root(),
                    env=_reference_environment(),
                    text=True,
                )
                process_stdout = process.stdout
                process_stderr = process.stderr
            except FileNotFoundError:
                process = None
                status = "skipped"
                mismatch_reason = "rscript_unavailable"
            if process is None:
                pass
            elif process.returncode != 0:
                mismatch_reason = "reference_execution_failed"
            else:
                execution_path = execution_root / "reference-execution.json"
                if not execution_path.exists():
                    mismatch_reason = "reference_execution_failed"
                else:
                    execution_payload = _load_json(execution_path)
                    r_version = _optional_payload_string(execution_payload, "r_version")
                    ape_version = _optional_payload_string(
                        execution_payload, "ape_version"
                    )
                    execution_status = execution_payload.get("status")
                    if execution_status == "unavailable":
                        status = "skipped"
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason", "ape_package_unavailable"
                            )
                        )
                    elif execution_status != "ok":
                        reference_error = {
                            "error_type": str(
                                execution_payload.get(
                                    "error_type",
                                    execution_payload.get(
                                        "mismatch_reason", "reference_execution_failed"
                                    ),
                                )
                            ),
                            "message": str(execution_payload.get("message", "")),
                        }
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason", "reference_execution_failed"
                            )
                        )
                    else:
                        (
                            reference_summary,
                            reference_rows,
                            reference_normalized_text,
                        ) = _load_reference_case_payload(case, execution_root)
                        if case.operation in {
                            "read-tree-structure",
                            "write-tree-structure",
                        }:
                            mismatch_reason = _tree_structure_mismatch_reason(
                                case,
                                execution_root,
                            )
                        elif case.operation == "root-tree-outgroup":
                            mismatch_reason = _root_tree_outgroup_mismatch_reason(
                                case,
                                execution_root,
                            )
                        elif case.operation == "unroot-tree":
                            mismatch_reason = _unroot_tree_mismatch_reason(
                                case,
                                execution_root,
                            )
                        elif case.operation == "drop-tree-taxa":
                            mismatch_reason = _drop_tip_tree_mismatch_reason(
                                case,
                                execution_root,
                            )
                            if mismatch_reason is None and not _compare_json(
                                reference_summary,
                                bijux_summary,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "summary_mismatch"
                        elif case.operation == "keep-tree-taxa":
                            mismatch_reason = _keep_tip_tree_mismatch_reason(
                                case,
                                execution_root,
                            )
                            if mismatch_reason is None and not _compare_json(
                                reference_summary,
                                bijux_summary,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "summary_mismatch"
                        elif case.operation == "extract-tree-clade":
                            mismatch_reason = _extract_clade_tree_mismatch_reason(
                                case,
                                execution_root,
                            )
                            if mismatch_reason is None and not _compare_json(
                                reference_summary,
                                bijux_summary,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "summary_mismatch"
                        elif (
                            case.operation == "get-tree-mrca"
                            or case.operation == "assess-tree-monophyly"
                        ):
                            if not _compare_json(
                                reference_summary,
                                bijux_summary,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "summary_mismatch"
                            else:
                                status = "passed"
                        elif case.operation in {
                            "read-tree-set-structure",
                            "write-tree-set-structure",
                        }:
                            mismatch_reason = _tree_set_structure_mismatch_reason(
                                case,
                                execution_root,
                            )
                        elif case.operation == "tree-consensus":
                            mismatch_reason = _consensus_tree_mismatch_reason(
                                case,
                                execution_root,
                            )
                            if mismatch_reason is None and not _compare_json(
                                reference_summary,
                                bijux_summary,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "summary_mismatch"
                            elif mismatch_reason is None and not _compare_json(
                                reference_rows,
                                bijux_rows,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "rows_mismatch"
                        elif case.operation == "distance-matrix-neighbor-joining":
                            mismatch_reason = _neighbor_joining_tree_mismatch_reason(
                                case,
                                execution_root,
                            )
                            if mismatch_reason is None and not _compare_json(
                                reference_summary,
                                bijux_summary,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "summary_mismatch"
                        elif case.operation == "tree-discrete-ancestral-states":
                            reference_transition_rows = (
                                []
                                if reference_summary is None
                                else reference_summary.get("transition_rate_rows", [])
                            )
                            bijux_transition_rows = (
                                []
                                if bijux_summary is None
                                else bijux_summary.get("transition_rate_rows", [])
                            )
                            reference_summary_without_transition_rows = (
                                {}
                                if reference_summary is None
                                else {
                                    key: value
                                    for key, value in reference_summary.items()
                                    if key != "transition_rate_rows"
                                }
                            )
                            bijux_summary_without_transition_rows = (
                                {}
                                if bijux_summary is None
                                else {
                                    key: value
                                    for key, value in bijux_summary.items()
                                    if key != "transition_rate_rows"
                                }
                            )
                            if not _compare_json(
                                reference_summary_without_transition_rows,
                                bijux_summary_without_transition_rows,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "summary_mismatch"
                            elif not _transition_rate_rows_match(
                                reference_rows=reference_transition_rows,
                                bijux_rows=bijux_transition_rows,
                                reference_summary=reference_summary,
                                bijux_summary=bijux_summary,
                                tolerance=(
                                    case.transition_rate_tolerance
                                    if case.transition_rate_tolerance is not None
                                    else case.tolerance
                                ),
                            ):
                                mismatch_reason = "transition_rate_rows_mismatch"
                        elif not _compare_json(
                            reference_summary, bijux_summary, tolerance=case.tolerance
                        ):
                            mismatch_reason = "summary_mismatch"
                        elif not _compare_json(
                            reference_rows,
                            bijux_rows,
                            tolerance=case.tolerance,
                        ):
                            mismatch_reason = "rows_mismatch"
                        elif reference_normalized_text != bijux_normalized_text:
                            mismatch_reason = "normalized_text_mismatch"
                        else:
                            status = "passed"
                        if mismatch_reason is None:
                            status = "passed"
            if case.expected_status == "parse-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_parse_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_parse_error_missing"
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "parse_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "rooting-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_rooting_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_rooting_error_missing"
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "rooting_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "clade-extraction-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_clade_extraction_error_missing"
                elif reference_error is None:
                    mismatch_reason = (
                        "reference_expected_clade_extraction_error_missing"
                    )
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "clade_extraction_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "mrca-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_mrca_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_mrca_error_missing"
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "mrca_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "monophyly-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_monophyly_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_monophyly_error_missing"
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "monophyly_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "consensus-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_consensus_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_consensus_error_missing"
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "consensus_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "prop-clades-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_prop_clades_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_prop_clades_error_missing"
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "prop_clades_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "dna-distance-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_dna_distance_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_dna_distance_error_missing"
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "dna_distance_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if status != "passed":
                artifact_root = _persist_failure_bundle(
                    failure_root=active_failure_root,
                    case=case,
                    case_file=case_file,
                    execution_root=execution_root,
                    execution_payload=execution_payload,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    reference_error=reference_error,
                    bijux_error=bijux_error,
                    reference_rows=reference_rows,
                    bijux_rows=bijux_rows,
                    bijux_normalized_text=bijux_normalized_text,
                    mismatch_reason=mismatch_reason or "reference_execution_failed",
                )
                if process_stdout:
                    (artifact_root / "reference-stdout.txt").write_text(
                        process_stdout, encoding="utf-8"
                    )
                if process_stderr:
                    (artifact_root / "reference-stderr.txt").write_text(
                        process_stderr, encoding="utf-8"
                    )
            observations.append(
                ApeParityObservation(
                    case_id=case.case_id,
                    fixture_kind=case.fixture_kind,
                    fixture_id=case.fixture_id,
                    function_name=case.function_name,
                    python_function_name=case.python_function_name,
                    input_fixture=case.input_fixture,
                    tolerance=case.tolerance,
                    r_version=r_version,
                    ape_version=ape_version,
                    bijux_version=bijux_version,
                    bijux_commit=bijux_commit,
                    status=status,
                    passed=status == "passed",
                    mismatch_reason=mismatch_reason,
                    reproducible_artifact_root=artifact_root,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    reference_error=reference_error,
                    bijux_error=bijux_error,
                )
            )
    case_count = len(observations)
    passed_case_count = sum(1 for item in observations if item.status == "passed")
    failed_case_count = sum(1 for item in observations if item.status == "failed")
    skipped_case_count = sum(1 for item in observations if item.status == "skipped")
    return ApeParityReport(
        observations=observations,
        summary_rows=_summary_rows(observations),
        case_count=case_count,
        passed_case_count=passed_case_count,
        failed_case_count=failed_case_count,
        skipped_case_count=skipped_case_count,
        all_passed=case_count > 0
        and passed_case_count == case_count
        and failed_case_count == 0
        and skipped_case_count == 0,
        limitations=[
            "The governed live `ape` parity registry is intentionally narrow until later rounds expand the shared fixture surface.",
            "This harness requires Rscript plus the `ape` and `jsonlite` R packages for live reference execution.",
        ],
    )


def write_ape_parity_summary_table(path: Path, report: ApeParityReport) -> Path:
    """Write one row per governed `ape` function summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "function_name",
                "case_count",
                "passed_case_count",
                "failed_case_count",
                "skipped_case_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.summary_rows:
            writer.writerow(asdict(row))
    return path


def write_ape_parity_observation_table(path: Path, report: ApeParityReport) -> Path:
    """Write one row per governed `ape` parity observation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "fixture_kind",
                "fixture_id",
                "function_name",
                "python_function_name",
                "input_fixture",
                "tolerance",
                "r_version",
                "ape_version",
                "bijux_version",
                "bijux_commit",
                "status",
                "passed",
                "mismatch_reason",
                "reproducible_artifact_root",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for observation in report.observations:
            writer.writerow(
                {
                    "case_id": observation.case_id,
                    "fixture_kind": observation.fixture_kind,
                    "fixture_id": observation.fixture_id,
                    "function_name": observation.function_name,
                    "python_function_name": observation.python_function_name,
                    "input_fixture": str(observation.input_fixture),
                    "tolerance": format(observation.tolerance, ".12g"),
                    "r_version": observation.r_version or "",
                    "ape_version": observation.ape_version or "",
                    "bijux_version": observation.bijux_version,
                    "bijux_commit": observation.bijux_commit or "",
                    "status": observation.status,
                    "passed": str(observation.passed).lower(),
                    "mismatch_reason": observation.mismatch_reason or "",
                    "reproducible_artifact_root": ""
                    if observation.reproducible_artifact_root is None
                    else str(observation.reproducible_artifact_root),
                }
            )
    return path
