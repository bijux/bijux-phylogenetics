from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_support_values,
    compare_tree_paths,
)
from bijux_phylogenetics.engines.common import load_engine_manifest

from .contracts import InferenceTreeComparisonReport


def compare_ml_trees_across_models(
    left_manifest_path: Path,
    right_manifest_path: Path,
) -> InferenceTreeComparisonReport:
    """Compare maximum-likelihood trees produced under different model choices."""
    return _compare_inference_trees(
        left_manifest_path,
        right_manifest_path,
        comparison_kind="model",
        left_label=_manifest_comparison_label(
            left_manifest_path, fallback="left-model"
        ),
        right_label=_manifest_comparison_label(
            right_manifest_path, fallback="right-model"
        ),
    )


def compare_inferred_trees_across_engines(
    left_manifest_path: Path,
    right_manifest_path: Path,
) -> InferenceTreeComparisonReport:
    """Compare inferred trees across two engine workflows."""
    left_manifest = load_engine_manifest(left_manifest_path)
    right_manifest = load_engine_manifest(right_manifest_path)
    return _compare_inference_trees(
        left_manifest_path,
        right_manifest_path,
        comparison_kind="engine",
        left_label=_display_engine_name(str(left_manifest["engine_name"])),
        right_label=_display_engine_name(str(right_manifest["engine_name"])),
    )


def _compare_inference_trees(
    left_manifest_path: Path,
    right_manifest_path: Path,
    *,
    comparison_kind: str,
    left_label: str,
    right_label: str,
) -> InferenceTreeComparisonReport:
    left_manifest = load_engine_manifest(left_manifest_path)
    right_manifest = load_engine_manifest(right_manifest_path)
    left_tree_path = _manifest_tree_output_path(left_manifest)
    right_tree_path = _manifest_tree_output_path(right_manifest)
    topology = compare_tree_paths(left_tree_path, right_tree_path)
    support = compare_support_values(left_tree_path, right_tree_path)
    branch_lengths = compare_branch_lengths(left_tree_path, right_tree_path)
    warnings: list[str] = []
    if not topology.topology_equal:
        warnings.append("inferred topologies differ across compared workflows")
    if topology.same_unrooted_topology and not topology.topology_equal:
        warnings.append(
            "compared workflows agree on unrooted splits but differ in rooting"
        )
    if topology.same_topology_different_branch_lengths:
        warnings.append(
            "compared workflows preserve topology but change branch-length interpretation"
        )
    return InferenceTreeComparisonReport(
        comparison_kind=comparison_kind,
        left_manifest_path=left_manifest_path,
        right_manifest_path=right_manifest_path,
        left_label=left_label,
        right_label=right_label,
        left_tree_path=left_tree_path,
        right_tree_path=right_tree_path,
        left_engine_name=str(left_manifest["engine_name"]),
        right_engine_name=str(right_manifest["engine_name"]),
        left_selected_model=_manifest_selected_model(left_manifest),
        right_selected_model=_manifest_selected_model(right_manifest),
        topology=topology,
        support=support,
        branch_lengths=branch_lengths,
        warnings=warnings,
    )


def _manifest_tree_output_path(manifest: dict[str, object]) -> Path:
    output_paths = {
        key: Path(value) for key, value in dict(manifest["output_paths"]).items()
    }
    tree_path = (
        output_paths.get("tree")
        or output_paths.get("support_tree")
        or output_paths.get("consensus_tree")
    )
    if tree_path is None:
        raise ValueError("manifest does not expose a tree output")
    return tree_path


def _manifest_selected_model(manifest: dict[str, object]) -> str | None:
    selected_model = manifest.get("selected_model")
    return None if selected_model is None else str(selected_model)


def _manifest_comparison_label(manifest_path: Path, *, fallback: str) -> str:
    manifest = load_engine_manifest(manifest_path)
    selected_model = _manifest_selected_model(manifest)
    if selected_model is not None:
        return selected_model
    return fallback


def _display_engine_name(raw: str) -> str:
    mapping = {
        "iqtree": "IQ-TREE",
        "iqtree2": "IQ-TREE",
        "fasttree": "FastTree",
        "mafft": "MAFFT",
        "trimal": "trimAl",
    }
    return mapping.get(raw.lower(), raw)
