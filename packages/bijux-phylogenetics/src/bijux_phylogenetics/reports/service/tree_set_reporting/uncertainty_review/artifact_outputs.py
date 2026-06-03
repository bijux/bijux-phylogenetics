from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.trees import (
    write_clade_frequency_table,
    write_consensus_tree,
    write_topology_cluster_table,
    write_tree_distance_distribution_table,
    write_unstable_clade_table,
)

from ...artifacts import write_json_artifact, write_tabular_artifact


def write_tree_uncertainty_artifacts(
    *,
    artifact_root: Path,
    summary,
    methods_summary_result,
    consensus_tree,
    clade_frequencies,
    diversity,
    clusters,
    unstable_taxa,
    unstable_clades,
    clade_conflicts,
    conclusion_summary,
    thinning_sensitivity,
    consensus_sensitivity,
    benchmark,
    multimodality,
    storage_risk,
    maturity,
    scaled_report_note: dict[str, object],
) -> dict[str, Path]:
    return {
        "tree_set_summary": write_json_artifact(
            artifact_root / "tree-set-summary.json", asdict(summary)
        ),
        "methods_summary": methods_summary_result.output_path,
        "consensus_tree": write_consensus_tree(
            artifact_root / "consensus-tree.nwk", consensus_tree
        ),
        "clade_frequencies": write_clade_frequency_table(
            artifact_root / "clade-frequencies.tsv", clade_frequencies
        ),
        "rf_distance_distribution": write_tree_distance_distribution_table(
            artifact_root / "rf-distance-distribution.tsv", diversity
        ),
        "topology_clusters": write_topology_cluster_table(
            artifact_root / "topology-clusters.tsv", clusters
        ),
        "unstable_taxa": write_tabular_artifact(
            artifact_root / "unstable-taxa.tsv",
            [asdict(row) for row in unstable_taxa.taxa],
        ),
        "unstable_clades": write_unstable_clade_table(
            artifact_root / "unstable-clades.tsv", unstable_clades
        ),
        "clade_credibility_conflicts": write_json_artifact(
            artifact_root / "clade-credibility-conflicts.json",
            scaled_report_note if clade_conflicts is None else asdict(clade_conflicts),
        ),
        "uncertainty_aware_conclusions": write_json_artifact(
            artifact_root / "uncertainty-aware-conclusions.json",
            scaled_report_note
            if conclusion_summary is None
            else asdict(conclusion_summary),
        ),
        "thinning_sensitivity": write_tabular_artifact(
            artifact_root / "thinning-sensitivity.tsv",
            []
            if thinning_sensitivity is None
            else [asdict(row) for row in thinning_sensitivity.rows],
        ),
        "consensus_threshold_sensitivity": write_tabular_artifact(
            artifact_root / "consensus-threshold-sensitivity.tsv",
            []
            if consensus_sensitivity is None
            else [asdict(row) for row in consensus_sensitivity.rows],
        ),
        "tree_set_benchmark": write_tabular_artifact(
            artifact_root / "tree-set-benchmark.tsv",
            [] if benchmark is None else [asdict(row) for row in benchmark.rows],
        ),
        "topological_diversity": write_json_artifact(
            artifact_root / "topological-diversity.json", asdict(diversity)
        ),
        "topology_multimodality": write_json_artifact(
            artifact_root / "topology-multimodality.json",
            scaled_report_note if multimodality is None else asdict(multimodality),
        ),
        "storage_risk": write_json_artifact(
            artifact_root / "storage-risk.json", asdict(storage_risk)
        ),
        "maturity_gate": write_json_artifact(
            artifact_root / "maturity-gate.json",
            scaled_report_note if maturity is None else asdict(maturity),
        ),
    }


def finalize_tree_uncertainty_outputs(
    *,
    title: str,
    out_path: Path,
    sections,
    summary_metrics,
    artifact_links: list[tuple[str, str, str | None]],
    artifact_manifest_path: Path,
    artifact_paths: dict[str, Path],
    machine_manifest: dict[str, object],
) -> tuple[dict[str, object], int, int]:
    write_json_artifact(artifact_manifest_path, machine_manifest)
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
        summary_metrics=summary_metrics,
        artifact_links=[
            *artifact_links,
            (
                "tree-uncertainty-manifest",
                artifact_manifest_path.relative_to(out_path.parent).as_posix(),
                None,
            ),
        ],
    )
    html_size_bytes = out_path.stat().st_size
    manifest_size_bytes = artifact_manifest_path.stat().st_size
    linked_artifact_bytes = (
        sum(path.stat().st_size for path in artifact_paths.values())
        + manifest_size_bytes
    )
    total_output_bytes = html_size_bytes + linked_artifact_bytes
    machine_manifest["linked_artifacts"]["tree_uncertainty_manifest"] = {
        "path": artifact_manifest_path.relative_to(out_path.parent).as_posix(),
        "byte_count": manifest_size_bytes,
    }
    machine_manifest["html_size_bytes"] = html_size_bytes
    machine_manifest["linked_artifact_bytes"] = linked_artifact_bytes
    machine_manifest["total_output_bytes"] = total_output_bytes
    write_json_artifact(artifact_manifest_path, machine_manifest)
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
        summary_metrics=summary_metrics,
        artifact_links=[
            *artifact_links,
            (
                "tree-uncertainty-manifest",
                artifact_manifest_path.relative_to(out_path.parent).as_posix(),
                f"{artifact_manifest_path.stat().st_size} bytes",
            ),
        ],
    )
    html_size_bytes = out_path.stat().st_size
    total_output_bytes = html_size_bytes + linked_artifact_bytes
    machine_manifest["html_size_bytes"] = html_size_bytes
    machine_manifest["total_output_bytes"] = total_output_bytes
    write_json_artifact(artifact_manifest_path, machine_manifest)
    final_manifest_size_bytes = artifact_manifest_path.stat().st_size
    if final_manifest_size_bytes != manifest_size_bytes:
        linked_artifact_bytes = (
            sum(path.stat().st_size for path in artifact_paths.values())
            + final_manifest_size_bytes
        )
        total_output_bytes = html_size_bytes + linked_artifact_bytes
        machine_manifest["linked_artifacts"]["tree_uncertainty_manifest"] = {
            "path": artifact_manifest_path.relative_to(out_path.parent).as_posix(),
            "byte_count": final_manifest_size_bytes,
        }
        machine_manifest["linked_artifact_bytes"] = linked_artifact_bytes
        machine_manifest["total_output_bytes"] = total_output_bytes
        write_json_artifact(artifact_manifest_path, machine_manifest)
    return machine_manifest, html_size_bytes, linked_artifact_bytes
