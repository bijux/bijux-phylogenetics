from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.core.manifest import build_run_manifest, write_run_manifest


def _finalize_outputs(
    args: Any,
    *,
    command: str,
    inputs: list[Path | str],
    outputs: list[Path | str] | None = None,
) -> list[Path | str]:
    finalized_outputs = list(outputs or [])
    manifest_path = getattr(args, "manifest", None)
    if manifest_path is None:
        return finalized_outputs
    manifest = build_run_manifest(
        command=command,
        arguments=list(getattr(args, "_argv", [])),
        input_paths=inputs,
        output_paths=finalized_outputs,
    )
    write_run_manifest(manifest_path, manifest)
    finalized_outputs.append(manifest_path)
    return finalized_outputs


def _command_inputs(args: Any) -> list[Path | str]:
    if args.command == "commands":
        return []
    if args.command == "env":
        return []
    if args.command == "metadata":
        return [args.table]
    if args.command == "traits":
        if args.traits_command == "validate":
            return [args.table]
        if args.traits_command == "missing":
            return [args.table]
        if args.traits_command == "prune":
            return [args.tree, args.table, args.out]
        return [args.tree, args.table]
    if args.command == "prune":
        inputs = [args.tree, args.out]
        if getattr(args, "keep_from", None) is not None:
            inputs.append(args.keep_from)
        if args.pruned_taxa_out is not None:
            inputs.append(args.pruned_taxa_out)
        return inputs
    if args.command == "alignment":
        if args.alignment_command == "alphabet":
            return [args.alignment]
        if args.alignment_command == "gc":
            return [args.alignment]
        if args.alignment_command == "inspect":
            return [args.alignment]
        if args.alignment_command == "composition":
            return [args.alignment]
        if args.alignment_command == "outliers":
            return [args.alignment]
        if args.alignment_command == "duplicates":
            return [args.alignment]
        if args.alignment_command == "invalid":
            return [args.alignment]
        if args.alignment_command == "quality":
            return [args.alignment]
        if args.alignment_command == "validate-input":
            return [args.alignment]
        if args.alignment_command == "repair-input":
            return [args.alignment, args.out]
        if args.alignment_command == "concatenate":
            inputs: list[Path | str] = [
                *args.alignments,
                args.out,
                args.partitions_out,
                args.matrix_out,
            ]
            for optional_path in (args.taxa_out, args.loci_out):
                if optional_path is not None:
                    inputs.append(optional_path)
            return inputs
        if args.alignment_command == "occupancy":
            inputs: list[Path | str] = [args.alignment, args.partitions]
            for optional_path in (
                args.taxa_out,
                args.loci_out,
                args.matrix_out,
                args.filtered_alignment_out,
                args.filtered_partitions_out,
            ):
                if optional_path is not None:
                    inputs.append(optional_path)
            return inputs
        if args.alignment_command == "trim":
            return [args.alignment, args.out]
        if args.alignment_command == "identity-matrix":
            inputs = [args.alignment]
            if args.out is not None:
                inputs.append(args.out)
            return inputs
        if args.alignment_command == "distance-matrix":
            inputs = [args.alignment]
            if args.out is not None:
                inputs.append(args.out)
            return inputs
        if args.alignment_command == "distance-quality":
            return [args.alignment]
        if args.alignment_command == "distance-saturation":
            return [args.alignment]
        if args.alignment_command == "distance-additivity":
            return [args.alignment, args.out_dir]
        if args.alignment_command == "distance-suitability":
            return [args.alignment]
        if args.alignment_command == "distance-assumptions":
            return [args.alignment]
        if args.alignment_command == "bootstrap-tree":
            inputs = [args.alignment]
            for optional_path in (
                args.support_out,
                args.tree_set_out,
                args.draws_out,
            ):
                if optional_path is not None:
                    inputs.append(optional_path)
            return inputs
        if args.alignment_command == "distance-support-summary":
            return [args.alignment]
        if args.alignment_command == "build-tree":
            return [args.alignment, args.out]
        if args.alignment_command == "compare-distance-trees":
            return [args.alignment]
        if args.alignment_command == "coding":
            return [args.alignment]
        if args.alignment_command == "translate":
            return [args.alignment, args.out]
        return [args.tree, args.alignment]
    if args.command == "distance":
        if args.distance_command == "validate":
            return [args.matrix]
        if args.distance_command == "additivity":
            return [args.matrix, args.out_dir]
        if args.distance_command == "build-tree":
            return [args.matrix, args.out]
        if args.distance_command == "minimum-evolution":
            return [args.matrix, args.tree, args.out]
        if args.distance_command == "fitch-margoliash":
            return [args.matrix, args.tree, args.out]
        if args.distance_command == "ordinary-least-squares":
            return [args.matrix, args.tree, args.out]
        if args.distance_command == "nonnegative-least-squares":
            return [args.matrix, args.tree, args.out]
        if args.distance_command == "patristic-residuals":
            return [args.matrix, args.tree, args.out_dir]
        if args.distance_command == "taxon-influence":
            return [args.matrix, args.reference_tree, args.out_dir]
        if args.distance_command == "taxon-jackknife":
            return [args.matrix, args.out_dir]
        if args.distance_command == "method-comparison":
            return [args.matrix, args.out_dir]
        if args.distance_command == "bme-nni-search":
            return [args.matrix, args.out_dir]
        if args.distance_command == "report":
            return [args.matrix, args.out]
        return [args.matrix]
    if args.command == "comparative":
        return [args.tree, args.table]
    if args.command == "ancestral":
        inputs = [args.tree, args.table]
        if getattr(args, "out", None) is not None:
            inputs.append(args.out)
        if getattr(args, "table_out", None) is not None:
            inputs.append(args.table_out)
        return inputs
    if args.command == "tree-set":
        if args.tree_set_command == "compare":
            return [args.left, args.right]
        if args.tree_set_command == "bootstrap-summary":
            return [args.tree_set, args.out_dir]
        inputs = [args.tree_set]
        if getattr(args, "out", None) is not None:
            inputs.append(args.out)
        return inputs
    if args.command == "simulate":
        if getattr(args, "tree", None) is not None:
            return [args.tree, args.out]
        return [args.out]
    if args.command == "benchmark":
        return []
    if args.command == "parity":
        return []
    if args.command in {"validate", "inspect"}:
        return [args.tree]
    if args.command == "diagnose":
        if getattr(args, "tree", None) is not None:
            return [args.tree]
        return [Path(args.target)]
    if args.command == "normalize":
        return [args.tree, args.out]
    if args.command == "normalize-taxa":
        inputs = [args.tree, args.out]
        if args.mapping_out is not None:
            inputs.append(args.mapping_out)
        return inputs
    if args.command == "topology":
        return [args.tree, args.out]
    if args.command == "compare":
        if args.left == "clades":
            inputs = [Path(args.right)]
            if getattr(args, "third", None) is not None:
                inputs.append(Path(args.third))
            for path in getattr(args, "extra_trees", []) or []:
                inputs.append(path)
            if getattr(args, "out", None) is not None:
                inputs.append(args.out)
            return inputs
        if args.left == "support":
            inputs = [Path(args.right)]
            if getattr(args, "third", None) is not None:
                inputs.append(Path(args.third))
            if getattr(args, "out", None) is not None:
                inputs.append(args.out)
            return inputs
        if args.left == "clade-ages":
            inputs = [Path(args.right)]
            if getattr(args, "third", None) is not None:
                inputs.append(Path(args.third))
            if getattr(args, "out", None) is not None:
                inputs.append(args.out)
            return inputs
        if getattr(args, "third", None) is not None:
            return [Path(args.right), Path(args.third)]
        return [Path(args.left), Path(args.right)]
    if args.command == "annotate":
        inputs = [args.tree, args.metadata]
        if args.out is not None:
            inputs.append(args.out)
        if getattr(args, "joined_out", None) is not None:
            inputs.append(args.joined_out)
        return inputs
    if args.command == "render":
        inputs = [args.tree, args.out]
        if args.metadata is not None:
            inputs.append(args.metadata)
        return inputs
    if args.command == "evidence":
        if args.evidence_command == "bundle":
            return [*args.inputs, *args.outputs]
        bundle_root = getattr(args, "bundle_root", None)
        return [] if bundle_root is None else [bundle_root]
    if args.command == "demo":
        return []
    if args.command == "report":
        if args.report_command == "tree":
            return [args.tree, args.out]
        if args.report_command == "alignment":
            return [args.alignment, args.out]
        if args.report_command == "dataset":
            inputs = [args.tree, args.metadata, args.out]
            if args.traits is not None:
                inputs.append(args.traits)
            return inputs
        if args.report_command == "phylo-inputs":
            return [args.tree, args.alignment, args.out]
        if args.report_command == "taxonomy":
            inputs = [args.tree, args.out]
            if args.synonym_table is not None:
                inputs.append(args.synonym_table)
            if args.metadata is not None:
                inputs.append(args.metadata)
            if args.traits is not None:
                inputs.append(args.traits)
            if args.alignment is not None:
                inputs.append(args.alignment)
            if args.filtered_alignment is not None:
                inputs.append(args.filtered_alignment)
            if args.inference_tree is not None:
                inputs.append(args.inference_tree)
            if args.reported_taxa is not None:
                inputs.append(args.reported_taxa)
            return inputs
        inputs = [args.tree, args.out]
        if getattr(args, "alignment", None) is not None:
            inputs.append(args.alignment)
        if getattr(args, "traits", None) is not None:
            inputs.append(args.traits)
        if getattr(args, "metadata", None) is not None:
            inputs.append(args.metadata)
        return inputs
    if args.command == "adapter":
        if args.adapter_command == "inspect":
            return [args.engine_name]
        if args.adapter_command == "report":
            return [args.manifest_path, args.out]
        if args.adapter_command == "compare":
            inputs = [args.fast_tree, args.ml_tree]
            if args.out is not None:
                inputs.append(args.out)
            return inputs
        if getattr(args, "out", None) is not None:
            return [args.input_path, args.out]
        if getattr(args, "out_dir", None) is not None:
            return [args.input_path, args.out_dir]
        return [args.input_path]
    if args.command == "phylo":
        if getattr(args, "phylo_command", None) == "run":
            return [args.config_path]
        if getattr(args, "phylo_command", None) == "replay":
            inputs: list[Path | str] = [args.manifest_path]
            if getattr(args, "out_dir", None) is not None:
                inputs.append(args.out_dir)
            return inputs
        if getattr(args, "phylo_command", None) == "bundle":
            return [args.manifest_path, args.out_dir]
        if getattr(args, "phylo_command", None) == "validate-bundle":
            return [args.bundle_root]
        return [] if getattr(args, "workflow", None) is None else [args.workflow]
    return []
