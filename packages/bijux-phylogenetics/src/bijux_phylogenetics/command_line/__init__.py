from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any

from bijux_phylogenetics import __version__
from bijux_phylogenetics.bayesian import (
    assess_beast_convergence,
    assess_mrbayes_convergence,
    build_bayesian_evidence_package,
    compute_mrbayes_effective_sample_sizes,
    parse_beast_log,
    parse_mrbayes_parameter_traces,
    prepare_beast_time_tree_analysis,
    prepare_mrbayes_analysis,
    render_calibration_audit_report,
    render_bayesian_posterior_report,
    run_mrbayes_posterior_inference,
    summarize_mrbayes_posterior_trees,
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
)
from bijux_phylogenetics.benchmark import (
    benchmark_alignment_diagnostics,
    benchmark_tree_comparison,
    benchmark_tree_validation,
)
from bijux_phylogenetics.ancestral.continuous import reconstruct_continuous_ancestral_states
from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.ancestral.package import build_ancestral_figure_package
from bijux_phylogenetics.ancestral.sensitivity import build_ancestral_sensitivity_report
from bijux_phylogenetics.ancestral.service import (
    compare_continuous_ancestral_models,
    render_ancestral_state_report,
    render_ancestral_state_tree,
    write_ancestral_state_table,
)
from bijux_phylogenetics.comparative.common import summarize_numeric_trait, summarize_numeric_trait_readiness
from bijux_phylogenetics.comparative.models import (
    compare_brownian_and_ou_models,
    fit_brownian_motion_model,
    fit_ornstein_uhlenbeck_model,
    run_comparative_sensitivity_analysis,
    validate_comparative_reference_examples,
)
from bijux_phylogenetics.comparative.pgls import inspect_pgls_inputs, run_pgls, run_pgls_multiple_testing
from bijux_phylogenetics.comparative.reporting import (
    build_comparative_method_report,
    build_trait_influence_report,
    compare_comparative_results_across_pruning,
    compare_comparative_results_across_trees,
    write_comparative_method_report,
)
from bijux_phylogenetics.comparative.signal import (
    compute_blombergs_k,
    compute_phylogenetic_independent_contrasts,
    compute_phylogenetic_signal_test,
    estimate_pagels_lambda,
)
from bijux_phylogenetics.command_line.registry import COMMAND_SPECS, get_command_spec
from bijux_phylogenetics.core.environment import inspect_environment
from bijux_phylogenetics.core.manifest import build_run_manifest, write_run_manifest
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.diagnostics.root_to_tip import (
    compute_root_to_tip_distances,
    diagnose_ultrametricity,
    write_root_to_tip_tsv,
)
from bijux_phylogenetics.io.fasta import (
    build_alignment_forensic_report,
    build_alignment_quality_report,
    clean_alignment_with_profile,
    classify_alignment_sequences,
    compare_alignment_versions,
    infer_alignment_alphabet,
    link_alignment_to_tree,
    list_alignment_filter_profiles,
    load_fasta_alignment,
    summarise_fasta,
    summarize_alignment_readiness,
    summarize_alignment_windows,
)
from bijux_phylogenetics.io.fasta import (
    compute_pairwise_sequence_identity_matrix,
    detect_identical_duplicate_sequences,
    detect_invalid_alignment_characters,
    detect_near_duplicate_sequences,
    detect_composition_outlier_sequences,
    detect_over_aligned_regions,
    detect_sequence_length_outliers,
    detect_under_aligned_regions,
    inspect_coding_alignment,
    translate_coding_alignment,
    trim_alignment,
    write_sequence_identity_matrix,
    write_fasta_alignment,
)
from bijux_phylogenetics.core.metadata import inspect_metadata_table
from bijux_phylogenetics.core.pruning import (
    drop_tree_taxa,
    prune_tree_to_requested_taxa,
    prune_tree_to_taxa,
    write_pruned_taxa,
)
from bijux_phylogenetics.core.traits import link_tree_to_traits, prune_traits_to_tree, validate_traits_table
from bijux_phylogenetics.core.traits import detect_missing_trait_values
from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_clade_sets,
    compare_support_values,
    compare_tree_paths,
    detect_clade_changes,
    prune_trees_to_shared_taxa,
    write_tree_comparison_table,
)
from bijux_phylogenetics.compare.reports import build_tree_comparison_report
from bijux_phylogenetics.core.demo import run_capability_demo
from bijux_phylogenetics.diagnostics.validation import diagnose_tree_path, inspect_tree_path, validate_tree_path
from bijux_phylogenetics.diagnostics.assumptions import assess_tree_assumptions
from bijux_phylogenetics.distance import (
    assess_distance_method_assumptions,
    assess_imported_distance_method_assumptions,
    bootstrap_distance_trees,
    build_distance_tree,
    build_tree_from_imported_distance_matrix,
    compare_distance_tree_topologies,
    compute_pairwise_genetic_distance_matrix,
    inspect_distance_matrix_quality,
    validate_distance_reference_examples,
    validate_imported_distance_matrix,
    write_distance_bootstrap_support,
    write_distance_reproducibility_bundle,
    write_genetic_distance_matrix,
)
from bijux_phylogenetics.discrete_evolution import (
    compare_discrete_state_models,
    detect_state_imbalance_problems,
    estimate_ancestral_geographic_states,
    render_discrete_state_evolution_report,
    render_tree_with_geographic_states,
    validate_discrete_transition_reference_examples,
    validate_discrete_state_coding,
    write_discrete_model_comparison_table,
    write_node_state_probability_table,
    write_transition_summary_table,
)
from bijux_phylogenetics.diversification import (
    compare_diversification_models,
    compute_lineage_through_time_curve,
    detect_diversification_outlier_clades,
    detect_incomplete_taxon_sampling_metadata,
    estimate_diversification_rate,
    render_diversification_report,
    run_trait_dependent_diversification_analysis,
    write_clade_diversification_table,
    write_lineage_through_time_table,
    write_trait_dependent_diversification_table,
)
from bijux_phylogenetics.engines import (
    compare_fast_and_ml_trees,
    render_inference_workflow_report,
    read_engine_version,
    run_alignment_trimming,
    run_bootstrap_consensus_tree,
    run_bootstrap_support_estimation,
    run_fast_tree_inference,
    run_maximum_likelihood_tree_inference,
    run_model_selection,
    run_multiple_sequence_alignment,
)
from bijux_phylogenetics.evidence.bundles import bundle_directory, validate_bundle
from bijux_phylogenetics.errors import EngineUnavailableError, EvidenceContractError, MetadataJoinError, PhylogeneticsError
from bijux_phylogenetics.core.taxonomy import (
    audit_tree_taxon_synonyms,
    inspect_tree_taxon_namespaces,
    normalize_tree_taxa,
    resolve_tree_taxon_synonyms,
    write_synonym_resolution_mapping,
    write_taxon_mapping,
)
from bijux_phylogenetics.core.taxon_workflows import (
    build_taxon_stability_report,
    build_taxon_workflow_loss_report,
    load_taxon_run_source,
)
from bijux_phylogenetics.core.topology import reroot_tree_by_midpoint, root_tree_on_outgroup, unroot_tree
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.render.package import build_tree_figure_package
from bijux_phylogenetics.render.svg import AnnotationStrip, render_tree_svg
from bijux_phylogenetics.reports.service import (
    annotate_tree_against_table,
    distance_method_limitations,
    render_dataset_report,
    render_distance_report,
    render_phylo_inputs_report,
    render_tree_uncertainty_report,
    render_tree_report,
    write_annotation_report,
)
from bijux_phylogenetics.simulation import (
    simulate_birth_death_trees,
    simulate_brownian_traits,
    simulate_coalescent_trees,
    simulate_discrete_traits,
    simulate_dna_alignment,
    simulate_ou_traits,
    simulate_protein_alignment,
    write_continuous_trait_table,
    write_discrete_trait_table,
    write_simulated_alignment,
    write_tree_set,
)
from bijux_phylogenetics.tree_set import (
    cluster_trees_by_topology,
    compare_posterior_tree_sets,
    compute_clade_frequency_table,
    compute_consensus_tree,
    compute_tree_distance_matrix,
    detect_unstable_clades,
    detect_unstable_taxa,
    load_tree_set,
    write_clade_frequency_table,
    write_consensus_tree,
    write_tree_distance_matrix,
)
from bijux_phylogenetics.results import build_command_result, build_error_result


def _json_ready(value: Any) -> Any:
    if is_dataclass(value):
        return {key: _json_ready(item) for key, item in asdict(value).items()}
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _json_ready(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    return value


def _print_result(result: Any, *, json_output: bool) -> None:
    if json_output:
        print(json.dumps(_json_ready(result), indent=2, sort_keys=True))
        return
    if isinstance(result, str):
        print(result)
        return
    print(json.dumps(_json_ready(result), indent=2, sort_keys=True))


def _split_csv_values(raw: str | None) -> list[str]:
    if raw is None:
        return []
    return [item.strip() for item in raw.split(",") if item.strip()]


def _parse_assignment_map(raw: str | None) -> dict[str, str]:
    assignments: dict[str, str] = {}
    for item in _split_csv_values(raw):
        if "=" not in item:
            raise ValueError(f"mapping item must be KEY=VALUE, got '{item}'")
        key, value = item.split("=", 1)
        if not key.strip() or not value.strip():
            raise ValueError(f"mapping item must include both KEY and VALUE, got '{item}'")
        assignments[key.strip()] = value.strip()
    return assignments


def _parse_labelled_run(raw: str) -> tuple[str, Path]:
    if "=" not in raw:
        raise ValueError(f"run source must be in LABEL=PATH form, got '{raw}'")
    label, raw_path = raw.split("=", 1)
    if not label.strip() or not raw_path.strip():
        raise ValueError(f"run source must include both LABEL and PATH, got '{raw}'")
    return label.strip(), Path(raw_path.strip())


def _validate_ancestral_discrete_model_arguments(args: Any, parser: argparse.ArgumentParser) -> None:
    if getattr(args, "kind", None) == "discrete" and getattr(args, "state_ordering", "unordered") == "ordered":
        resolved_model = getattr(args, "model", None) or "fitch"
        if resolved_model == "fitch":
            parser.error("ordered ancestral discrete reconstruction requires a likelihood model")


def _build_annotation_strips(table, columns: list[str]) -> list[AnnotationStrip]:
    missing_columns = [column for column in columns if column not in table.columns]
    if missing_columns:
        raise MetadataJoinError(f"table does not contain columns: {', '.join(missing_columns)}")
    return [
        AnnotationStrip(
            name=column,
            values={row[table.taxon_column]: row[column] for row in table.rows if row[column]},
        )
        for column in columns
    ]


def _build_numeric_trait_map(table, column: str) -> dict[str, float]:
    if column not in table.columns:
        raise MetadataJoinError(f"table does not contain column '{column}'")
    values: dict[str, float] = {}
    for row in table.rows:
        raw = row[column]
        if not raw:
            continue
        try:
            values[row[table.taxon_column]] = float(raw)
        except ValueError as error:
            raise MetadataJoinError(f"column '{column}' contains a non-numeric value for taxon '{row[table.taxon_column]}'") from error
    return values


def _build_string_trait_map(table, column: str) -> dict[str, str]:
    if column not in table.columns:
        raise MetadataJoinError(f"table does not contain column '{column}'")
    return {row[table.taxon_column]: row[column] for row in table.rows if row[column]}


def _print_commands(*, output_format: str) -> None:
    payload = build_command_result(
        command="commands",
        inputs=[],
        outputs=[],
        metrics={"command_count": len(COMMAND_SPECS)},
        data={"commands": list(COMMAND_SPECS)},
    )
    if output_format == "json":
        print(json.dumps(_json_ready(payload), indent=2, sort_keys=True))
        return
    for command in _json_ready(payload.data)["commands"]:
        print(f"{command['name']}: {command['domain']} - {command['summary']}")


def _json_requested(args: Any) -> bool:
    return bool(getattr(args, "json", False) or getattr(args, "format", "") == "json")


def _add_manifest_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--manifest", type=Path, help="Write a reproducibility manifest to this JSON path.")


def _adapter_version_args(engine_name: str) -> tuple[str, ...]:
    normalized = engine_name.lower()
    if normalized == "fasttree":
        return ("-help",)
    return ("--version",)


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
        if args.distance_command == "build-tree":
            return [args.matrix, args.out]
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
        return [args.bundle_root]
    if args.command == "demo":
        return []
    if args.command == "report":
        if args.report_command == "tree":
            return [args.tree, args.out]
        if args.report_command == "dataset":
            inputs = [args.tree, args.metadata, args.out]
            if args.traits is not None:
                inputs.append(args.traits)
            return inputs
        if args.report_command == "phylo-inputs":
            return [args.tree, args.alignment, args.out]
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
    return []


def build_parser() -> argparse.ArgumentParser:
    """Build the repository CLI parser."""
    parser = argparse.ArgumentParser(prog="bijux-phylogenetics")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    commands = subparsers.add_parser("commands", help="List the registered command taxonomy.")
    commands.add_argument("--format", choices=("text", "json"), default="text")

    env = subparsers.add_parser(get_command_spec("env").name, help=get_command_spec("env").summary)
    env_subparsers = env.add_subparsers(dest="env_command", required=True)
    env_inspect = env_subparsers.add_parser("inspect", help="Inspect runtime dependency availability.")
    env_inspect.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(env_inspect)

    metadata = subparsers.add_parser(get_command_spec("metadata").name, help=get_command_spec("metadata").summary)
    metadata_subparsers = metadata.add_subparsers(dest="metadata_command", required=True)
    metadata_inspect = metadata_subparsers.add_parser("inspect", help="Inspect a metadata table keyed by taxon.")
    metadata_inspect.add_argument("table", type=Path)
    metadata_inspect.add_argument("--taxon-column")
    metadata_inspect.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(metadata_inspect)

    traits = subparsers.add_parser(get_command_spec("traits").name, help=get_command_spec("traits").summary)
    traits_subparsers = traits.add_subparsers(dest="traits_command", required=True)
    traits_validate = traits_subparsers.add_parser("validate", help="Validate a traits table keyed by taxon.")
    traits_validate.add_argument("table", type=Path)
    traits_validate.add_argument("--taxon-column")
    traits_validate.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(traits_validate)
    traits_missing = traits_subparsers.add_parser("missing", help="List missing trait values by taxon and column.")
    traits_missing.add_argument("table", type=Path)
    traits_missing.add_argument("--taxon-column")
    traits_missing.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(traits_missing)
    traits_link = traits_subparsers.add_parser("link", help="Link tree tips to a traits table.")
    traits_link.add_argument("tree", type=Path)
    traits_link.add_argument("table", type=Path)
    traits_link.add_argument("--taxon-column")
    traits_link.add_argument("--strict", action="store_true")
    traits_link.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(traits_link)
    traits_prune = traits_subparsers.add_parser("prune", help="Prune a traits table to tree taxa.")
    traits_prune.add_argument("tree", type=Path)
    traits_prune.add_argument("table", type=Path)
    traits_prune.add_argument("--taxon-column")
    traits_prune.add_argument("--out", required=True, type=Path)
    traits_prune.add_argument("--json", action="store_true", help="Emit the pruning report as JSON.")
    _add_manifest_argument(traits_prune)

    prune = subparsers.add_parser(get_command_spec("prune").name, help=get_command_spec("prune").summary)
    prune.add_argument("tree", type=Path)
    prune_targets = prune.add_mutually_exclusive_group(required=True)
    prune_targets.add_argument("--keep-from", type=Path)
    prune_targets.add_argument("--taxa", nargs="+")
    prune_targets.add_argument("--exclude-taxa", nargs="+")
    prune.add_argument("--taxon-column")
    prune.add_argument("--out", required=True, type=Path)
    prune.add_argument("--pruned-taxa-out", type=Path)
    prune.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(prune)

    alignment = subparsers.add_parser(get_command_spec("alignment").name, help=get_command_spec("alignment").summary)
    alignment_subparsers = alignment.add_subparsers(dest="alignment_command", required=True)
    alignment_alphabet = alignment_subparsers.add_parser("alphabet", help="Infer the alignment sequence alphabet.")
    alignment_alphabet.add_argument("alignment", type=Path)
    alignment_alphabet.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_alphabet)
    alignment_profiles = alignment_subparsers.add_parser(
        "profiles",
        help="List the supported named alignment-cleaning profiles.",
    )
    alignment_profiles.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_profiles)
    alignment_gc = alignment_subparsers.add_parser("gc", help="Report per-sequence and whole-alignment GC content.")
    alignment_gc.add_argument("alignment", type=Path)
    alignment_gc.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_gc)
    alignment_inspect = alignment_subparsers.add_parser("inspect", help="Inspect an aligned FASTA file.")
    alignment_inspect.add_argument("alignment", type=Path)
    alignment_inspect.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_inspect)
    alignment_classify = alignment_subparsers.add_parser(
        "classify",
        help="Classify whether a FASTA input is aligned, raw, or shape-ambiguous.",
    )
    alignment_classify.add_argument("alignment", type=Path)
    alignment_classify.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_classify)
    alignment_quality = alignment_subparsers.add_parser("quality", help="Generate a higher-level alignment quality report.")
    alignment_quality.add_argument("alignment", type=Path)
    alignment_quality.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_quality)
    alignment_windows = alignment_subparsers.add_parser(
        "windows",
        help="Summarize sliding-window alignment quality and suspicious regions.",
    )
    alignment_windows.add_argument("alignment", type=Path)
    alignment_windows.add_argument("--window-size", type=int, default=30)
    alignment_windows.add_argument("--step-size", type=int, default=10)
    alignment_windows.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_windows)
    alignment_readiness = alignment_subparsers.add_parser(
        "readiness",
        help="Classify whether an alignment is ready for distance, ML, Bayesian, coding, or protein workflows.",
    )
    alignment_readiness.add_argument("alignment", type=Path)
    alignment_readiness.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_readiness)
    alignment_forensic = alignment_subparsers.add_parser(
        "forensic",
        help="Build a reviewer-facing alignment forensic report.",
    )
    alignment_forensic.add_argument("alignment", type=Path)
    alignment_forensic.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_forensic)
    alignment_composition = alignment_subparsers.add_parser(
        "composition",
        help="Inspect inferred alphabet, composition, and GC content.",
    )
    alignment_composition.add_argument("alignment", type=Path)
    alignment_composition.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_composition)
    alignment_invalid = alignment_subparsers.add_parser(
        "invalid",
        help="List alignment characters invalid for a declared alphabet.",
    )
    alignment_invalid.add_argument("alignment", type=Path)
    alignment_invalid.add_argument("--alphabet", choices=("dna", "rna", "protein"), required=True)
    alignment_invalid.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_invalid)
    alignment_duplicates = alignment_subparsers.add_parser(
        "duplicates",
        help="Report identical and near-duplicate aligned sequences.",
    )
    alignment_duplicates.add_argument("alignment", type=Path)
    alignment_duplicates.add_argument("--identity-threshold", type=float, default=0.95)
    alignment_duplicates.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_duplicates)
    alignment_outliers = alignment_subparsers.add_parser(
        "outliers",
        help="Report composition outlier sequences from an alignment.",
    )
    alignment_outliers.add_argument("alignment", type=Path)
    alignment_outliers.add_argument("--deviation-threshold", type=float, default=0.25)
    alignment_outliers.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_outliers)
    alignment_length_outliers = alignment_subparsers.add_parser(
        "length-outliers",
        help="Report raw sequence length outliers before alignment assumptions are imposed.",
    )
    alignment_length_outliers.add_argument("alignment", type=Path)
    alignment_length_outliers.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_length_outliers)
    alignment_filter = alignment_subparsers.add_parser(
        "filter",
        help="Clean an alignment through one named profile and report what changed.",
    )
    alignment_filter.add_argument("alignment", type=Path)
    alignment_filter.add_argument("--profile", required=True, choices=tuple(profile.name for profile in list_alignment_filter_profiles()))
    alignment_filter.add_argument("--out", required=True, type=Path)
    alignment_filter.add_argument("--group-table", type=Path)
    alignment_filter.add_argument("--group-columns")
    alignment_filter.add_argument("--json", action="store_true", help="Emit the cleaning report as JSON.")
    _add_manifest_argument(alignment_filter)
    alignment_compare = alignment_subparsers.add_parser(
        "compare",
        help="Compare two alignment versions for taxa, sites, missingness, gaps, signal, and composition.",
    )
    alignment_compare.add_argument("left_alignment", type=Path)
    alignment_compare.add_argument("right_alignment", type=Path)
    alignment_compare.add_argument("--json", action="store_true", help="Emit the comparison report as JSON.")
    _add_manifest_argument(alignment_compare)
    alignment_trim = alignment_subparsers.add_parser(
        "trim",
        help="Trim all-gap or all-missing sites and optionally drop high-missingness sequences.",
    )
    alignment_trim.add_argument("alignment", type=Path)
    alignment_trim.add_argument("--out", required=True, type=Path)
    alignment_trim.add_argument("--keep-all-gap-sites", action="store_true")
    alignment_trim.add_argument("--keep-all-missing-sites", action="store_true")
    alignment_trim.add_argument("--site-missingness-threshold", type=float)
    alignment_trim.add_argument("--sequence-missingness-threshold", type=float)
    alignment_trim.add_argument("--json", action="store_true", help="Emit the trimming report as JSON.")
    _add_manifest_argument(alignment_trim)
    alignment_identity = alignment_subparsers.add_parser(
        "identity-matrix",
        help="Compute a pairwise sequence identity matrix.",
    )
    alignment_identity.add_argument("alignment", type=Path)
    alignment_identity.add_argument("--out", type=Path, help="Write the matrix as TSV.")
    alignment_identity.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_identity)
    alignment_distance = alignment_subparsers.add_parser(
        "distance-matrix",
        help="Compute a pairwise DNA genetic distance matrix.",
    )
    alignment_distance.add_argument("alignment", type=Path)
    alignment_distance.add_argument(
        "--model",
        choices=("p-distance", "jukes-cantor", "kimura-2-parameter", "amino-acid-p-distance"),
        default="p-distance",
    )
    alignment_distance.add_argument(
        "--gap-handling",
        choices=("pairwise-deletion", "complete-deletion"),
        default="pairwise-deletion",
    )
    alignment_distance.add_argument(
        "--ambiguity-policy",
        choices=("ignore", "partial-match", "strict-mismatch"),
        default="ignore",
    )
    alignment_distance.add_argument("--out", type=Path, help="Write the matrix as TSV.")
    alignment_distance.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_distance)
    alignment_distance_quality = alignment_subparsers.add_parser(
        "distance-quality",
        help="Inspect saturation, divergence, and low-information risks in a computed distance matrix.",
    )
    alignment_distance_quality.add_argument("alignment", type=Path)
    alignment_distance_quality.add_argument(
        "--model",
        choices=("p-distance", "jukes-cantor", "kimura-2-parameter", "amino-acid-p-distance"),
        default="p-distance",
    )
    alignment_distance_quality.add_argument(
        "--gap-handling",
        choices=("pairwise-deletion", "complete-deletion"),
        default="pairwise-deletion",
    )
    alignment_distance_quality.add_argument(
        "--ambiguity-policy",
        choices=("ignore", "partial-match", "strict-mismatch"),
        default="ignore",
    )
    alignment_distance_quality.add_argument("--json", action="store_true", help="Emit the diagnostics as JSON.")
    _add_manifest_argument(alignment_distance_quality)
    alignment_distance_assumptions = alignment_subparsers.add_parser(
        "distance-assumptions",
        help="Audit NJ and UPGMA assumptions, including UPGMA ultrametric compatibility.",
    )
    alignment_distance_assumptions.add_argument("alignment", type=Path)
    alignment_distance_assumptions.add_argument(
        "--model",
        choices=("p-distance", "jukes-cantor", "kimura-2-parameter", "amino-acid-p-distance"),
        default="p-distance",
    )
    alignment_distance_assumptions.add_argument(
        "--gap-handling",
        choices=("pairwise-deletion", "complete-deletion"),
        default="pairwise-deletion",
    )
    alignment_distance_assumptions.add_argument(
        "--ambiguity-policy",
        choices=("ignore", "partial-match", "strict-mismatch"),
        default="ignore",
    )
    alignment_distance_assumptions.add_argument("--json", action="store_true", help="Emit the assumption audit as JSON.")
    _add_manifest_argument(alignment_distance_assumptions)
    alignment_build_tree = alignment_subparsers.add_parser(
        "build-tree",
        help="Build a neighbor-joining or UPGMA tree from a DNA distance matrix.",
    )
    alignment_build_tree.add_argument("alignment", type=Path)
    alignment_build_tree.add_argument("--method", choices=("neighbor-joining", "upgma"), required=True)
    alignment_build_tree.add_argument(
        "--model",
        choices=("p-distance", "jukes-cantor", "kimura-2-parameter", "amino-acid-p-distance"),
        default="p-distance",
    )
    alignment_build_tree.add_argument(
        "--gap-handling",
        choices=("pairwise-deletion", "complete-deletion"),
        default="pairwise-deletion",
    )
    alignment_build_tree.add_argument(
        "--ambiguity-policy",
        choices=("ignore", "partial-match", "strict-mismatch"),
        default="ignore",
    )
    alignment_build_tree.add_argument("--out", required=True, type=Path)
    alignment_build_tree.add_argument("--json", action="store_true", help="Emit the build report as JSON.")
    _add_manifest_argument(alignment_build_tree)
    alignment_compare_distance_trees = alignment_subparsers.add_parser(
        "compare-distance-trees",
        help="Compare NJ and UPGMA topologies built from the same DNA alignment.",
    )
    alignment_compare_distance_trees.add_argument("alignment", type=Path)
    alignment_compare_distance_trees.add_argument(
        "--model",
        choices=("p-distance", "jukes-cantor", "kimura-2-parameter", "amino-acid-p-distance"),
        default="p-distance",
    )
    alignment_compare_distance_trees.add_argument(
        "--gap-handling",
        choices=("pairwise-deletion", "complete-deletion"),
        default="pairwise-deletion",
    )
    alignment_compare_distance_trees.add_argument(
        "--ambiguity-policy",
        choices=("ignore", "partial-match", "strict-mismatch"),
        default="ignore",
    )
    alignment_compare_distance_trees.add_argument("--json", action="store_true", help="Emit the comparison as JSON.")
    _add_manifest_argument(alignment_compare_distance_trees)
    alignment_bootstrap_tree = alignment_subparsers.add_parser(
        "bootstrap-tree",
        help="Bootstrap a distance tree by resampling alignment sites with replacement.",
    )
    alignment_bootstrap_tree.add_argument("alignment", type=Path)
    alignment_bootstrap_tree.add_argument("--method", choices=("neighbor-joining", "upgma"), required=True)
    alignment_bootstrap_tree.add_argument(
        "--model",
        choices=("p-distance", "jukes-cantor", "kimura-2-parameter", "amino-acid-p-distance"),
        default="p-distance",
    )
    alignment_bootstrap_tree.add_argument(
        "--gap-handling",
        choices=("pairwise-deletion", "complete-deletion"),
        default="pairwise-deletion",
    )
    alignment_bootstrap_tree.add_argument(
        "--ambiguity-policy",
        choices=("ignore", "partial-match", "strict-mismatch"),
        default="ignore",
    )
    alignment_bootstrap_tree.add_argument("--replicates", type=int, default=100)
    alignment_bootstrap_tree.add_argument("--seed", type=int, default=1)
    alignment_bootstrap_tree.add_argument("--support-out", type=Path, help="Write bootstrap clade support as TSV.")
    alignment_bootstrap_tree.add_argument("--tree-set-out", type=Path, help="Write bootstrap replicate trees as Newick.")
    alignment_bootstrap_tree.add_argument("--json", action="store_true", help="Emit the bootstrap report as JSON.")
    _add_manifest_argument(alignment_bootstrap_tree)
    alignment_distance_bundle = alignment_subparsers.add_parser(
        "distance-bundle",
        help="Write a reproducibility bundle for one distance-analysis workflow.",
    )
    alignment_distance_bundle.add_argument("alignment", type=Path)
    alignment_distance_bundle.add_argument("--method", choices=("neighbor-joining", "upgma"), required=True)
    alignment_distance_bundle.add_argument(
        "--model",
        choices=("p-distance", "jukes-cantor", "kimura-2-parameter", "amino-acid-p-distance"),
        default="p-distance",
    )
    alignment_distance_bundle.add_argument(
        "--gap-handling",
        choices=("pairwise-deletion", "complete-deletion"),
        default="pairwise-deletion",
    )
    alignment_distance_bundle.add_argument(
        "--ambiguity-policy",
        choices=("ignore", "partial-match", "strict-mismatch"),
        default="ignore",
    )
    alignment_distance_bundle.add_argument("--replicates", type=int, default=100)
    alignment_distance_bundle.add_argument("--seed", type=int, default=1)
    alignment_distance_bundle.add_argument("--out-dir", required=True, type=Path)
    alignment_distance_bundle.add_argument("--json", action="store_true", help="Emit the bundle report as JSON.")
    _add_manifest_argument(alignment_distance_bundle)
    alignment_coding = alignment_subparsers.add_parser(
        "coding",
        help="Inspect a nucleotide coding alignment for frameshift-like lengths and stop codons.",
    )
    alignment_coding.add_argument("alignment", type=Path)
    alignment_coding.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_coding)
    alignment_translate = alignment_subparsers.add_parser(
        "translate",
        help="Translate a nucleotide coding alignment to amino acids.",
    )
    alignment_translate.add_argument("alignment", type=Path)
    alignment_translate.add_argument("--out", required=True, type=Path)
    alignment_translate.add_argument("--json", action="store_true", help="Emit the translation report as JSON.")
    _add_manifest_argument(alignment_translate)
    alignment_link = alignment_subparsers.add_parser("link", help="Link tree tips to an aligned FASTA file.")
    alignment_link.add_argument("tree", type=Path)
    alignment_link.add_argument("alignment", type=Path)
    alignment_link.add_argument("--strict", action="store_true")
    alignment_link.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_link)

    comparative = subparsers.add_parser(
        get_command_spec("comparative").name,
        help=get_command_spec("comparative").summary,
    )
    comparative_subparsers = comparative.add_subparsers(dest="comparative_command", required=True)
    comparative_readiness = comparative_subparsers.add_parser(
        "readiness",
        help="Check whether a rooted tree and numeric trait are ready for comparative analysis.",
    )
    comparative_readiness.add_argument("tree", type=Path)
    comparative_readiness.add_argument("table", type=Path)
    comparative_readiness.add_argument("--trait", required=True)
    comparative_readiness.add_argument("--taxon-column")
    comparative_readiness.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(comparative_readiness)
    comparative_summarize = comparative_subparsers.add_parser(
        "summarize",
        help="Summarize a numeric trait after pruning to overlapping phylogenetic taxa.",
    )
    comparative_summarize.add_argument("tree", type=Path)
    comparative_summarize.add_argument("table", type=Path)
    comparative_summarize.add_argument("--trait", required=True)
    comparative_summarize.add_argument("--taxon-column")
    comparative_summarize.add_argument("--json", action="store_true", help="Emit the summary as JSON.")
    _add_manifest_argument(comparative_summarize)
    comparative_contrasts = comparative_subparsers.add_parser(
        "contrasts",
        help="Compute phylogenetic independent contrasts for one numeric trait.",
    )
    comparative_contrasts.add_argument("tree", type=Path)
    comparative_contrasts.add_argument("table", type=Path)
    comparative_contrasts.add_argument("--trait", required=True)
    comparative_contrasts.add_argument("--taxon-column")
    comparative_contrasts.add_argument("--json", action="store_true", help="Emit the contrast report as JSON.")
    _add_manifest_argument(comparative_contrasts)
    comparative_signal = comparative_subparsers.add_parser(
        "signal",
        help="Estimate phylogenetic signal metrics for one numeric trait.",
    )
    comparative_signal.add_argument("tree", type=Path)
    comparative_signal.add_argument("table", type=Path)
    comparative_signal.add_argument("--trait", required=True)
    comparative_signal.add_argument("--taxon-column")
    comparative_signal.add_argument("--permutations", type=int, default=199)
    comparative_signal.add_argument("--seed", type=int, default=1)
    comparative_signal.add_argument("--json", action="store_true", help="Emit the signal report as JSON.")
    _add_manifest_argument(comparative_signal)
    comparative_brownian = comparative_subparsers.add_parser(
        "brownian",
        help="Fit a standalone Brownian-motion continuous-trait model.",
    )
    comparative_brownian.add_argument("tree", type=Path)
    comparative_brownian.add_argument("table", type=Path)
    comparative_brownian.add_argument("--trait", required=True)
    comparative_brownian.add_argument("--taxon-column")
    comparative_brownian.add_argument("--json", action="store_true", help="Emit the Brownian model fit as JSON.")
    _add_manifest_argument(comparative_brownian)
    comparative_ou = comparative_subparsers.add_parser(
        "ou",
        help="Fit a standalone Ornstein-Uhlenbeck continuous-trait model.",
    )
    comparative_ou.add_argument("tree", type=Path)
    comparative_ou.add_argument("table", type=Path)
    comparative_ou.add_argument("--trait", required=True)
    comparative_ou.add_argument("--taxon-column")
    comparative_ou.add_argument("--json", action="store_true", help="Emit the OU model fit as JSON.")
    _add_manifest_argument(comparative_ou)
    comparative_compare_models = comparative_subparsers.add_parser(
        "compare-models",
        help="Compare standalone Brownian-motion and OU models for one continuous trait.",
    )
    comparative_compare_models.add_argument("tree", type=Path)
    comparative_compare_models.add_argument("table", type=Path)
    comparative_compare_models.add_argument("--trait", required=True)
    comparative_compare_models.add_argument("--taxon-column")
    comparative_compare_models.add_argument("--json", action="store_true", help="Emit the model comparison as JSON.")
    _add_manifest_argument(comparative_compare_models)
    comparative_validate_reference = comparative_subparsers.add_parser(
        "validate-reference",
        help="Validate built-in Brownian-motion and OU reference examples.",
    )
    comparative_validate_reference.add_argument("--json", action="store_true", help="Emit the validation report as JSON.")
    _add_manifest_argument(comparative_validate_reference)
    comparative_sensitivity = comparative_subparsers.add_parser(
        "sensitivity",
        help="Run leave-one-taxon-out sensitivity for a standalone BM or OU model.",
    )
    comparative_sensitivity.add_argument("tree", type=Path)
    comparative_sensitivity.add_argument("table", type=Path)
    comparative_sensitivity.add_argument("--trait", required=True)
    comparative_sensitivity.add_argument("--model", choices=("brownian", "ou"), required=True)
    comparative_sensitivity.add_argument("--taxon-column")
    comparative_sensitivity.add_argument("--json", action="store_true", help="Emit the sensitivity report as JSON.")
    _add_manifest_argument(comparative_sensitivity)
    comparative_pgls = comparative_subparsers.add_parser(
        "pgls",
        help="Fit a phylogenetic generalized least-squares model.",
    )
    comparative_pgls.add_argument("tree", type=Path)
    comparative_pgls.add_argument("table", type=Path)
    comparative_pgls.add_argument("--response")
    comparative_pgls.add_argument("--predictors", nargs="+")
    comparative_pgls.add_argument("--formula", help="Formula-style specification such as 'response ~ body_mass * habitat'.")
    comparative_pgls.add_argument("--taxon-column")
    comparative_pgls.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_pgls.add_argument("--json", action="store_true", help="Emit the model result as JSON.")
    _add_manifest_argument(comparative_pgls)
    comparative_multiple_testing = comparative_subparsers.add_parser(
        "multiple-testing",
        help="Adjust PGLS coefficient p-values across many response traits.",
    )
    comparative_multiple_testing.add_argument("tree", type=Path)
    comparative_multiple_testing.add_argument("table", type=Path)
    comparative_multiple_testing.add_argument("--responses", nargs="+", required=True)
    comparative_multiple_testing.add_argument("--predictors", nargs="+", required=True)
    comparative_multiple_testing.add_argument("--taxon-column")
    comparative_multiple_testing.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_multiple_testing.add_argument("--json", action="store_true", help="Emit the correction report as JSON.")
    _add_manifest_argument(comparative_multiple_testing)
    comparative_report = comparative_subparsers.add_parser(
        "report",
        help="Build an integrated comparative-method report.",
    )
    comparative_report.add_argument("tree", type=Path)
    comparative_report.add_argument("table", type=Path)
    comparative_report.add_argument("--response")
    comparative_report.add_argument("--predictors", nargs="+")
    comparative_report.add_argument("--formula")
    comparative_report.add_argument("--taxon-column")
    comparative_report.add_argument("--out", type=Path)
    comparative_report.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_report.add_argument("--json", action="store_true", help="Emit the comparative report as JSON.")
    _add_manifest_argument(comparative_report)
    comparative_influence = comparative_subparsers.add_parser(
        "influence",
        help="Identify predictor terms and taxa driving one comparative result.",
    )
    comparative_influence.add_argument("tree", type=Path)
    comparative_influence.add_argument("table", type=Path)
    comparative_influence.add_argument("--response")
    comparative_influence.add_argument("--predictors", nargs="+")
    comparative_influence.add_argument("--formula")
    comparative_influence.add_argument("--taxon-column")
    comparative_influence.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_influence.add_argument("--json", action="store_true", help="Emit the influence report as JSON.")
    _add_manifest_argument(comparative_influence)
    comparative_compare_trees = comparative_subparsers.add_parser(
        "compare-trees",
        help="Compare comparative results across two alternative trees.",
    )
    comparative_compare_trees.add_argument("left_tree", type=Path)
    comparative_compare_trees.add_argument("right_tree", type=Path)
    comparative_compare_trees.add_argument("table", type=Path)
    comparative_compare_trees.add_argument("--response")
    comparative_compare_trees.add_argument("--predictors", nargs="+")
    comparative_compare_trees.add_argument("--formula")
    comparative_compare_trees.add_argument("--taxon-column")
    comparative_compare_trees.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_compare_trees.add_argument("--json", action="store_true", help="Emit the comparison as JSON.")
    _add_manifest_argument(comparative_compare_trees)
    comparative_compare_pruning = comparative_subparsers.add_parser(
        "compare-pruning",
        help="Compare comparative results before and after explicit pruning.",
    )
    comparative_compare_pruning.add_argument("tree", type=Path)
    comparative_compare_pruning.add_argument("table", type=Path)
    comparative_compare_pruning.add_argument("--response")
    comparative_compare_pruning.add_argument("--predictors", nargs="+")
    comparative_compare_pruning.add_argument("--formula")
    comparative_compare_pruning.add_argument("--drop-taxa", nargs="+")
    comparative_compare_pruning.add_argument("--keep-taxa", nargs="+")
    comparative_compare_pruning.add_argument("--taxon-column")
    comparative_compare_pruning.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_compare_pruning.add_argument("--json", action="store_true", help="Emit the pruning comparison as JSON.")
    _add_manifest_argument(comparative_compare_pruning)

    ancestral = subparsers.add_parser(
        get_command_spec("ancestral").name,
        help=get_command_spec("ancestral").summary,
    )
    ancestral_subparsers = ancestral.add_subparsers(dest="ancestral_command", required=True)
    ancestral_continuous = ancestral_subparsers.add_parser(
        "continuous",
        help="Reconstruct ancestral states for a continuous trait.",
    )
    ancestral_continuous.add_argument("tree", type=Path)
    ancestral_continuous.add_argument("table", type=Path)
    ancestral_continuous.add_argument("--trait", required=True)
    ancestral_continuous.add_argument("--taxon-column")
    ancestral_continuous.add_argument("--model", choices=("brownian", "ou"), default="brownian")
    ancestral_continuous.add_argument("--alpha", type=float, default=1.0)
    ancestral_continuous.add_argument("--table-out", type=Path)
    ancestral_continuous.add_argument("--json", action="store_true", help="Emit the reconstruction as JSON.")
    _add_manifest_argument(ancestral_continuous)
    ancestral_discrete = ancestral_subparsers.add_parser(
        "discrete",
        help="Reconstruct ancestral states for a discrete trait.",
    )
    ancestral_discrete.add_argument("tree", type=Path)
    ancestral_discrete.add_argument("table", type=Path)
    ancestral_discrete.add_argument("--trait", required=True)
    ancestral_discrete.add_argument("--taxon-column")
    ancestral_discrete.add_argument("--model", choices=("fitch", "equal-rates", "symmetric", "all-rates-different"), default="fitch")
    ancestral_discrete.add_argument("--state-ordering", choices=("unordered", "ordered"), default="unordered")
    ancestral_discrete.add_argument("--ordered-states", help="Comma-delimited explicit ordered state vocabulary.")
    ancestral_discrete.add_argument("--table-out", type=Path)
    ancestral_discrete.add_argument("--json", action="store_true", help="Emit the reconstruction as JSON.")
    _add_manifest_argument(ancestral_discrete)
    ancestral_compare = ancestral_subparsers.add_parser(
        "compare",
        help="Compare two continuous ancestral-state models node by node.",
    )
    ancestral_compare.add_argument("tree", type=Path)
    ancestral_compare.add_argument("table", type=Path)
    ancestral_compare.add_argument("--trait", required=True)
    ancestral_compare.add_argument("--taxon-column")
    ancestral_compare.add_argument("--left-model", choices=("brownian", "ou"), default="brownian")
    ancestral_compare.add_argument("--right-model", choices=("brownian", "ou"), default="ou")
    ancestral_compare.add_argument("--left-alpha", type=float, default=1.0)
    ancestral_compare.add_argument("--right-alpha", type=float, default=1.0)
    ancestral_compare.add_argument("--json", action="store_true", help="Emit the comparison as JSON.")
    _add_manifest_argument(ancestral_compare)
    ancestral_sensitivity = ancestral_subparsers.add_parser(
        "sensitivity",
        help="Summarize how ancestral results change across model, tree, pruning, or coding choices.",
    )
    ancestral_sensitivity.add_argument("tree", type=Path)
    ancestral_sensitivity.add_argument("table", type=Path)
    ancestral_sensitivity.add_argument("--trait", required=True)
    ancestral_sensitivity.add_argument("--kind", choices=("continuous", "discrete"), required=True)
    ancestral_sensitivity.add_argument("--taxon-column")
    ancestral_sensitivity.add_argument("--model")
    ancestral_sensitivity.add_argument("--state-ordering", choices=("unordered", "ordered"), default="unordered")
    ancestral_sensitivity.add_argument("--ordered-states", help="Comma-delimited explicit ordered state vocabulary.")
    ancestral_sensitivity.add_argument("--alpha", type=float, default=1.0)
    ancestral_sensitivity.add_argument("--compare-model")
    ancestral_sensitivity.add_argument("--compare-tree", type=Path)
    ancestral_sensitivity.add_argument("--drop-taxa", nargs="+")
    ancestral_sensitivity.add_argument("--coding-map", help="Comma-delimited KEY=VALUE recoding map for discrete traits.")
    ancestral_sensitivity.add_argument("--json", action="store_true", help="Emit the sensitivity report as JSON.")
    _add_manifest_argument(ancestral_sensitivity)
    ancestral_render = ancestral_subparsers.add_parser(
        "render",
        help="Render a tree annotated with reconstructed ancestral states.",
    )
    ancestral_render.add_argument("tree", type=Path)
    ancestral_render.add_argument("table", type=Path)
    ancestral_render.add_argument("--trait", required=True)
    ancestral_render.add_argument("--kind", choices=("continuous", "discrete"), required=True)
    ancestral_render.add_argument("--taxon-column")
    ancestral_render.add_argument("--model")
    ancestral_render.add_argument("--state-ordering", choices=("unordered", "ordered"), default="unordered")
    ancestral_render.add_argument("--ordered-states", help="Comma-delimited explicit ordered state vocabulary.")
    ancestral_render.add_argument("--alpha", type=float, default=1.0)
    ancestral_render.add_argument("--layout", choices=("cladogram", "phylogram", "circular"), default="phylogram")
    ancestral_render.add_argument("--out", required=True, type=Path)
    ancestral_render.add_argument("--json", action="store_true", help="Emit the render result as JSON.")
    _add_manifest_argument(ancestral_render)
    ancestral_report = ancestral_subparsers.add_parser(
        "report",
        help="Render an HTML report for ancestral-state reconstruction.",
    )
    ancestral_report.add_argument("tree", type=Path)
    ancestral_report.add_argument("table", type=Path)
    ancestral_report.add_argument("--trait", required=True)
    ancestral_report.add_argument("--kind", choices=("continuous", "discrete"), required=True)
    ancestral_report.add_argument("--taxon-column")
    ancestral_report.add_argument("--model")
    ancestral_report.add_argument("--state-ordering", choices=("unordered", "ordered"), default="unordered")
    ancestral_report.add_argument("--ordered-states", help="Comma-delimited explicit ordered state vocabulary.")
    ancestral_report.add_argument("--alpha", type=float, default=1.0)
    ancestral_report.add_argument("--compare-model")
    ancestral_report.add_argument("--compare-tree", type=Path)
    ancestral_report.add_argument("--drop-taxa", nargs="+")
    ancestral_report.add_argument("--coding-map", help="Comma-delimited KEY=VALUE recoding map for discrete traits.")
    ancestral_report.add_argument("--out", required=True, type=Path)
    ancestral_report.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(ancestral_report)
    ancestral_package = ancestral_subparsers.add_parser(
        "package",
        help="Write a publication-ready ancestral-state figure package.",
    )
    ancestral_package.add_argument("tree", type=Path)
    ancestral_package.add_argument("table", type=Path)
    ancestral_package.add_argument("--trait", required=True)
    ancestral_package.add_argument("--kind", choices=("continuous", "discrete"), required=True)
    ancestral_package.add_argument("--taxon-column")
    ancestral_package.add_argument("--model")
    ancestral_package.add_argument("--state-ordering", choices=("unordered", "ordered"), default="unordered")
    ancestral_package.add_argument("--ordered-states", help="Comma-delimited explicit ordered state vocabulary.")
    ancestral_package.add_argument("--alpha", type=float, default=1.0)
    ancestral_package.add_argument("--layout", choices=("cladogram", "phylogram", "circular"), default="phylogram")
    ancestral_package.add_argument("--out-dir", required=True, type=Path)
    ancestral_package.add_argument("--json", action="store_true", help="Emit the package build result as JSON.")
    _add_manifest_argument(ancestral_package)

    discrete_evolution = subparsers.add_parser(
        get_command_spec("discrete-evolution").name,
        help=get_command_spec("discrete-evolution").summary,
    )
    discrete_evolution_subparsers = discrete_evolution.add_subparsers(
        dest="discrete_evolution_command",
        required=True,
    )
    discrete_validate = discrete_evolution_subparsers.add_parser(
        "validate-coding",
        help="Validate discrete-state labels against tree-overlapping taxa.",
    )
    discrete_validate.add_argument("tree", type=Path)
    discrete_validate.add_argument("table", type=Path)
    discrete_validate.add_argument("--trait", required=True)
    discrete_validate.add_argument("--taxon-column")
    discrete_validate.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, accept any single token state label.",
    )
    discrete_validate.add_argument("--state-ordering", choices=("unordered", "ordered"), default="unordered")
    discrete_validate.add_argument("--ordered-states", help="Comma-delimited explicit ordered state vocabulary.")
    discrete_validate.add_argument("--json", action="store_true", help="Emit the validation report as JSON.")
    _add_manifest_argument(discrete_validate)
    discrete_imbalance = discrete_evolution_subparsers.add_parser(
        "imbalance",
        help="Detect rare, dominant, or degenerate state balance problems.",
    )
    discrete_imbalance.add_argument("tree", type=Path)
    discrete_imbalance.add_argument("table", type=Path)
    discrete_imbalance.add_argument("--trait", required=True)
    discrete_imbalance.add_argument("--taxon-column")
    discrete_imbalance.add_argument("--json", action="store_true", help="Emit the imbalance report as JSON.")
    _add_manifest_argument(discrete_imbalance)
    discrete_reference = discrete_evolution_subparsers.add_parser(
        "reference",
        help="Validate deterministic discrete-state transition examples against built-in reference expectations.",
    )
    discrete_reference.add_argument("--json", action="store_true", help="Emit the reference-validation report as JSON.")
    _add_manifest_argument(discrete_reference)
    discrete_model = discrete_evolution_subparsers.add_parser(
        "model",
        help="Run one discrete-state transition model and export node or branch summaries.",
    )
    discrete_model.add_argument("tree", type=Path)
    discrete_model.add_argument("table", type=Path)
    discrete_model.add_argument("--trait", required=True)
    discrete_model.add_argument("--taxon-column")
    discrete_model.add_argument("--model", choices=("equal-rates", "symmetric", "all-rates-different"), default="equal-rates")
    discrete_model.add_argument("--state-ordering", choices=("unordered", "ordered"), default="unordered")
    discrete_model.add_argument("--ordered-states", help="Comma-delimited explicit ordered state vocabulary.")
    discrete_model.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_model.add_argument("--node-table-out", type=Path, help="Write node-state probabilities as TSV.")
    discrete_model.add_argument("--transitions-out", type=Path, help="Write branch transition summaries as TSV.")
    discrete_model.add_argument("--json", action="store_true", help="Emit the model report as JSON.")
    _add_manifest_argument(discrete_model)
    discrete_compare = discrete_evolution_subparsers.add_parser(
        "compare-models",
        help="Compare two supported discrete-state evolution models node by node.",
    )
    discrete_compare.add_argument("tree", type=Path)
    discrete_compare.add_argument("table", type=Path)
    discrete_compare.add_argument("--trait", required=True)
    discrete_compare.add_argument("--taxon-column")
    discrete_compare.add_argument("--left-model", choices=("equal-rates", "symmetric", "all-rates-different"), default="equal-rates")
    discrete_compare.add_argument(
        "--right-model",
        choices=("equal-rates", "symmetric", "all-rates-different"),
        default="all-rates-different",
    )
    discrete_compare.add_argument("--state-ordering", choices=("unordered", "ordered"), default="unordered")
    discrete_compare.add_argument("--ordered-states", help="Comma-delimited explicit ordered state vocabulary.")
    discrete_compare.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_compare.add_argument("--table-out", type=Path, help="Write node-wise model differences as TSV.")
    discrete_compare.add_argument("--json", action="store_true", help="Emit the comparison report as JSON.")
    _add_manifest_argument(discrete_compare)
    discrete_render = discrete_evolution_subparsers.add_parser(
        "render",
        help="Render a tree annotated with reconstructed geographic or other discrete states.",
    )
    discrete_render.add_argument("tree", type=Path)
    discrete_render.add_argument("table", type=Path)
    discrete_render.add_argument("--trait", required=True)
    discrete_render.add_argument("--taxon-column")
    discrete_render.add_argument("--model", choices=("equal-rates", "symmetric", "all-rates-different"), default="equal-rates")
    discrete_render.add_argument("--state-ordering", choices=("unordered", "ordered"), default="unordered")
    discrete_render.add_argument("--ordered-states", help="Comma-delimited explicit ordered state vocabulary.")
    discrete_render.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_render.add_argument("--layout", choices=("cladogram", "phylogram", "circular"), default="phylogram")
    discrete_render.add_argument("--out", required=True, type=Path)
    discrete_render.add_argument("--json", action="store_true", help="Emit the render result as JSON.")
    _add_manifest_argument(discrete_render)
    discrete_report = discrete_evolution_subparsers.add_parser(
        "report",
        help="Render an HTML report for one discrete-state evolution analysis.",
    )
    discrete_report.add_argument("tree", type=Path)
    discrete_report.add_argument("table", type=Path)
    discrete_report.add_argument("--trait", required=True)
    discrete_report.add_argument("--taxon-column")
    discrete_report.add_argument("--model", choices=("equal-rates", "symmetric", "all-rates-different"), default="equal-rates")
    discrete_report.add_argument("--state-ordering", choices=("unordered", "ordered"), default="unordered")
    discrete_report.add_argument("--ordered-states", help="Comma-delimited explicit ordered state vocabulary.")
    discrete_report.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_report.add_argument("--compare-model", choices=("equal-rates", "symmetric", "all-rates-different"))
    discrete_report.add_argument("--out", required=True, type=Path)
    discrete_report.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(discrete_report)

    diversification = subparsers.add_parser(
        get_command_spec("diversification").name,
        help=get_command_spec("diversification").summary,
    )
    diversification_subparsers = diversification.add_subparsers(dest="diversification_command", required=True)
    diversification_ltt = diversification_subparsers.add_parser(
        "ltt",
        help="Compute a lineage-through-time curve for one rooted ultrametric tree.",
    )
    diversification_ltt.add_argument("tree", type=Path)
    diversification_ltt.add_argument("--out", type=Path, help="Write the lineage-through-time table as TSV.")
    diversification_ltt.add_argument("--json", action="store_true", help="Emit the LTT report as JSON.")
    _add_manifest_argument(diversification_ltt)
    diversification_sampling = diversification_subparsers.add_parser(
        "sampling",
        help="Inspect taxon sampling-fraction metadata against the tree tips.",
    )
    diversification_sampling.add_argument("tree", type=Path)
    diversification_sampling.add_argument("table", type=Path)
    diversification_sampling.add_argument("--taxon-column")
    diversification_sampling.add_argument("--sampling-column")
    diversification_sampling.add_argument("--json", action="store_true", help="Emit the sampling report as JSON.")
    _add_manifest_argument(diversification_sampling)
    diversification_estimate = diversification_subparsers.add_parser(
        "estimate",
        help="Estimate a simple Yule or birth-death diversification model.",
    )
    diversification_estimate.add_argument("tree", type=Path)
    diversification_estimate.add_argument("--metadata", type=Path)
    diversification_estimate.add_argument("--taxon-column")
    diversification_estimate.add_argument("--sampling-column")
    diversification_estimate.add_argument("--model", choices=("yule", "birth-death"), default="birth-death")
    diversification_estimate.add_argument("--json", action="store_true", help="Emit the diversification estimate as JSON.")
    _add_manifest_argument(diversification_estimate)
    diversification_compare = diversification_subparsers.add_parser(
        "compare-models",
        help="Compare Yule and birth-death diversification fits.",
    )
    diversification_compare.add_argument("tree", type=Path)
    diversification_compare.add_argument("--metadata", type=Path)
    diversification_compare.add_argument("--taxon-column")
    diversification_compare.add_argument("--sampling-column")
    diversification_compare.add_argument("--json", action="store_true", help="Emit the model comparison as JSON.")
    _add_manifest_argument(diversification_compare)
    diversification_clades = diversification_subparsers.add_parser(
        "clades",
        help="Detect clades with unusually high or low diversification.",
    )
    diversification_clades.add_argument("tree", type=Path)
    diversification_clades.add_argument("--model", choices=("yule", "birth-death"), default="birth-death")
    diversification_clades.add_argument("--min-tip-count", type=int, default=2)
    diversification_clades.add_argument("--out", type=Path, help="Write the clade diversification table as TSV.")
    diversification_clades.add_argument("--json", action="store_true", help="Emit the clade scan report as JSON.")
    _add_manifest_argument(diversification_clades)
    diversification_trait = diversification_subparsers.add_parser(
        "trait-dependent",
        help="Summarize simple trait-linked diversification rates when states form interpretable clades.",
    )
    diversification_trait.add_argument("tree", type=Path)
    diversification_trait.add_argument("table", type=Path)
    diversification_trait.add_argument("--trait", required=True)
    diversification_trait.add_argument("--taxon-column")
    diversification_trait.add_argument("--out", type=Path, help="Write the trait-dependent diversification table as TSV.")
    diversification_trait.add_argument("--json", action="store_true", help="Emit the trait-dependent report as JSON.")
    _add_manifest_argument(diversification_trait)
    diversification_report = diversification_subparsers.add_parser(
        "report",
        help="Render an HTML diversification and macroevolution report.",
    )
    diversification_report.add_argument("tree", type=Path)
    diversification_report.add_argument("--metadata", type=Path)
    diversification_report.add_argument("--taxon-column")
    diversification_report.add_argument("--sampling-column")
    diversification_report.add_argument("--traits", type=Path)
    diversification_report.add_argument("--trait")
    diversification_report.add_argument("--out", required=True, type=Path)
    diversification_report.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(diversification_report)

    distance = subparsers.add_parser(get_command_spec("distance").name, help=get_command_spec("distance").summary)
    distance_subparsers = distance.add_subparsers(dest="distance_command", required=True)
    distance_validate = distance_subparsers.add_parser("validate", help="Validate an imported long-form distance matrix.")
    distance_validate.add_argument("matrix", type=Path)
    distance_validate.add_argument("--json", action="store_true", help="Emit the validation report as JSON.")
    _add_manifest_argument(distance_validate)
    distance_assumptions = distance_subparsers.add_parser(
        "assumptions",
        help="Audit NJ and UPGMA assumptions for an imported distance matrix.",
    )
    distance_assumptions.add_argument("matrix", type=Path)
    distance_assumptions.add_argument("--json", action="store_true", help="Emit the assumption audit as JSON.")
    _add_manifest_argument(distance_assumptions)
    distance_build_tree = distance_subparsers.add_parser(
        "build-tree",
        help="Build a Neighbor-Joining or UPGMA tree from an imported distance matrix.",
    )
    distance_build_tree.add_argument("matrix", type=Path)
    distance_build_tree.add_argument("--method", choices=("neighbor-joining", "upgma"), required=True)
    distance_build_tree.add_argument("--out", required=True, type=Path)
    distance_build_tree.add_argument("--json", action="store_true", help="Emit the build report as JSON.")
    _add_manifest_argument(distance_build_tree)
    distance_report = distance_subparsers.add_parser(
        "report",
        help="Render an HTML report for an imported distance matrix.",
    )
    distance_report.add_argument("matrix", type=Path)
    distance_report.add_argument("--out", required=True, type=Path)
    distance_report.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(distance_report)
    distance_explain = distance_subparsers.add_parser(
        "explain",
        help="Explain why distance-based tree building is approximate.",
    )
    distance_explain.add_argument("matrix", type=Path)
    distance_explain.add_argument("--json", action="store_true", help="Emit the explanation as JSON.")
    _add_manifest_argument(distance_explain)
    distance_reference = distance_subparsers.add_parser(
        "reference",
        help="Validate built-in reference examples for core distance calculations.",
    )
    distance_reference.add_argument("--json", action="store_true", help="Emit the validation report as JSON.")
    _add_manifest_argument(distance_reference)

    tree_set = subparsers.add_parser(get_command_spec("tree-set").name, help=get_command_spec("tree-set").summary)
    tree_set_subparsers = tree_set.add_subparsers(dest="tree_set_command", required=True)
    tree_set_inspect = tree_set_subparsers.add_parser(
        "inspect",
        help="Inspect a tree set for tree count and topology diversity.",
    )
    tree_set_inspect.add_argument("tree_set", type=Path)
    tree_set_inspect.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(tree_set_inspect)
    tree_set_consensus = tree_set_subparsers.add_parser(
        "consensus",
        help="Build a majority-rule consensus tree from a tree set.",
    )
    tree_set_consensus.add_argument("tree_set", type=Path)
    tree_set_consensus.add_argument("--out", required=True, type=Path)
    tree_set_consensus.add_argument("--json", action="store_true", help="Emit the consensus report as JSON.")
    _add_manifest_argument(tree_set_consensus)
    tree_set_clades = tree_set_subparsers.add_parser(
        "clade-frequencies",
        help="Compute clade support frequencies across a tree set.",
    )
    tree_set_clades.add_argument("tree_set", type=Path)
    tree_set_clades.add_argument("--out", type=Path, help="Write the clade-frequency table as TSV.")
    tree_set_clades.add_argument("--json", action="store_true", help="Emit the clade-frequency report as JSON.")
    _add_manifest_argument(tree_set_clades)
    tree_set_distances = tree_set_subparsers.add_parser(
        "distance-matrix",
        help="Compute pairwise RF distances across a tree set.",
    )
    tree_set_distances.add_argument("tree_set", type=Path)
    tree_set_distances.add_argument("--out", type=Path, help="Write the pairwise distance table as TSV.")
    tree_set_distances.add_argument("--json", action="store_true", help="Emit the distance report as JSON.")
    _add_manifest_argument(tree_set_distances)
    tree_set_clusters = tree_set_subparsers.add_parser(
        "cluster",
        help="Cluster trees by identical rooted topology signatures.",
    )
    tree_set_clusters.add_argument("tree_set", type=Path)
    tree_set_clusters.add_argument("--json", action="store_true", help="Emit the cluster report as JSON.")
    _add_manifest_argument(tree_set_clusters)
    tree_set_unstable_taxa = tree_set_subparsers.add_parser(
        "unstable-taxa",
        help="Detect taxa with inconsistent placements across a tree set.",
    )
    tree_set_unstable_taxa.add_argument("tree_set", type=Path)
    tree_set_unstable_taxa.add_argument("--json", action="store_true", help="Emit the instability report as JSON.")
    _add_manifest_argument(tree_set_unstable_taxa)
    tree_set_unstable_clades = tree_set_subparsers.add_parser(
        "unstable-clades",
        help="Detect non-unanimous and conflicting clades across a tree set.",
    )
    tree_set_unstable_clades.add_argument("tree_set", type=Path)
    tree_set_unstable_clades.add_argument("--json", action="store_true", help="Emit the instability report as JSON.")
    _add_manifest_argument(tree_set_unstable_clades)
    tree_set_compare = tree_set_subparsers.add_parser(
        "compare",
        help="Compare two posterior tree sets over clade support and topology distance.",
    )
    tree_set_compare.add_argument("left", type=Path)
    tree_set_compare.add_argument("right", type=Path)
    tree_set_compare.add_argument("--json", action="store_true", help="Emit the comparison report as JSON.")
    _add_manifest_argument(tree_set_compare)
    tree_set_report = tree_set_subparsers.add_parser(
        "report",
        help="Render an HTML uncertainty report for a tree set.",
    )
    tree_set_report.add_argument("tree_set", type=Path)
    tree_set_report.add_argument("--out", required=True, type=Path)
    tree_set_report.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(tree_set_report)

    simulate = subparsers.add_parser(get_command_spec("simulate").name, help=get_command_spec("simulate").summary)
    simulate_subparsers = simulate.add_subparsers(dest="simulate_command", required=True)
    simulate_birth_death = simulate_subparsers.add_parser(
        "tree-birth-death",
        help="Simulate one or more trees under a birth-death process.",
    )
    simulate_birth_death.add_argument("--tree-count", type=int, default=1)
    simulate_birth_death.add_argument("--tip-count", type=int, required=True)
    simulate_birth_death.add_argument("--birth-rate", type=float, default=1.0)
    simulate_birth_death.add_argument("--death-rate", type=float, default=0.25)
    simulate_birth_death.add_argument("--seed", type=int, default=1)
    simulate_birth_death.add_argument("--out", required=True, type=Path)
    simulate_birth_death.add_argument("--json", action="store_true", help="Emit the simulation report as JSON.")
    _add_manifest_argument(simulate_birth_death)
    simulate_coalescent = simulate_subparsers.add_parser(
        "tree-coalescent",
        help="Simulate one or more trees under a coalescent model.",
    )
    simulate_coalescent.add_argument("--tree-count", type=int, default=1)
    simulate_coalescent.add_argument("--tip-count", type=int, required=True)
    simulate_coalescent.add_argument("--population-size", type=float, default=1.0)
    simulate_coalescent.add_argument("--seed", type=int, default=1)
    simulate_coalescent.add_argument("--out", required=True, type=Path)
    simulate_coalescent.add_argument("--json", action="store_true", help="Emit the simulation report as JSON.")
    _add_manifest_argument(simulate_coalescent)
    simulate_brownian = simulate_subparsers.add_parser(
        "traits-brownian",
        help="Simulate a continuous tip trait under Brownian motion.",
    )
    simulate_brownian.add_argument("tree", type=Path)
    simulate_brownian.add_argument("--root-state", type=float, default=0.0)
    simulate_brownian.add_argument("--sigma", type=float, default=1.0)
    simulate_brownian.add_argument("--seed", type=int, default=1)
    simulate_brownian.add_argument("--out", required=True, type=Path)
    simulate_brownian.add_argument("--json", action="store_true", help="Emit the simulation report as JSON.")
    _add_manifest_argument(simulate_brownian)
    simulate_ou = simulate_subparsers.add_parser(
        "traits-ou",
        help="Simulate a continuous tip trait under an OU process.",
    )
    simulate_ou.add_argument("tree", type=Path)
    simulate_ou.add_argument("--root-state", type=float, default=0.0)
    simulate_ou.add_argument("--sigma", type=float, default=1.0)
    simulate_ou.add_argument("--alpha", type=float, default=1.0)
    simulate_ou.add_argument("--theta", type=float, default=0.0)
    simulate_ou.add_argument("--seed", type=int, default=1)
    simulate_ou.add_argument("--out", required=True, type=Path)
    simulate_ou.add_argument("--json", action="store_true", help="Emit the simulation report as JSON.")
    _add_manifest_argument(simulate_ou)
    simulate_discrete = simulate_subparsers.add_parser(
        "traits-discrete",
        help="Simulate a discrete tip trait under a symmetric jump process.",
    )
    simulate_discrete.add_argument("tree", type=Path)
    simulate_discrete.add_argument("--states", nargs="+", required=True)
    simulate_discrete.add_argument("--transition-rate", type=float, default=1.0)
    simulate_discrete.add_argument("--root-state")
    simulate_discrete.add_argument("--seed", type=int, default=1)
    simulate_discrete.add_argument("--out", required=True, type=Path)
    simulate_discrete.add_argument("--json", action="store_true", help="Emit the simulation report as JSON.")
    _add_manifest_argument(simulate_discrete)
    simulate_dna = simulate_subparsers.add_parser(
        "alignment-dna",
        help="Simulate a DNA alignment along a rooted tree.",
    )
    simulate_dna.add_argument("tree", type=Path)
    simulate_dna.add_argument("--sequence-length", type=int, required=True)
    simulate_dna.add_argument("--substitution-rate", type=float, default=1.0)
    simulate_dna.add_argument("--seed", type=int, default=1)
    simulate_dna.add_argument("--out", required=True, type=Path)
    simulate_dna.add_argument("--json", action="store_true", help="Emit the simulation report as JSON.")
    _add_manifest_argument(simulate_dna)
    simulate_protein = simulate_subparsers.add_parser(
        "alignment-protein",
        help="Simulate a protein alignment along a rooted tree.",
    )
    simulate_protein.add_argument("tree", type=Path)
    simulate_protein.add_argument("--sequence-length", type=int, required=True)
    simulate_protein.add_argument("--substitution-rate", type=float, default=1.0)
    simulate_protein.add_argument("--seed", type=int, default=1)
    simulate_protein.add_argument("--out", required=True, type=Path)
    simulate_protein.add_argument("--json", action="store_true", help="Emit the simulation report as JSON.")
    _add_manifest_argument(simulate_protein)

    benchmark = subparsers.add_parser(get_command_spec("benchmark").name, help=get_command_spec("benchmark").summary)
    benchmark_subparsers = benchmark.add_subparsers(dest="benchmark_command", required=True)
    benchmark_validate = benchmark_subparsers.add_parser(
        "tree-validation",
        help="Benchmark tree validation across size classes.",
    )
    benchmark_validate.add_argument("--replicates", type=int, default=3)
    benchmark_validate.add_argument("--json", action="store_true", help="Emit the benchmark report as JSON.")
    _add_manifest_argument(benchmark_validate)
    benchmark_compare = benchmark_subparsers.add_parser(
        "tree-comparison",
        help="Benchmark tree comparison across increasing taxon counts.",
    )
    benchmark_compare.add_argument("--replicates", type=int, default=3)
    benchmark_compare.add_argument("--json", action="store_true", help="Emit the benchmark report as JSON.")
    _add_manifest_argument(benchmark_compare)
    benchmark_alignment = benchmark_subparsers.add_parser(
        "alignment-diagnostics",
        help="Benchmark alignment diagnostics across increasing sequence counts.",
    )
    benchmark_alignment.add_argument("--replicates", type=int, default=3)
    benchmark_alignment.add_argument("--sequence-length", type=int, default=128)
    benchmark_alignment.add_argument("--json", action="store_true", help="Emit the benchmark report as JSON.")
    _add_manifest_argument(benchmark_alignment)

    validate = subparsers.add_parser(get_command_spec("validate").name, help=get_command_spec("validate").summary)
    validate.add_argument("tree", type=Path)
    validate.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    validate.add_argument("--allow-duplicates", action="store_true")
    validate.add_argument("--allow-negative-branches", action="store_true")
    validate.add_argument("--require-rooted", action="store_true")
    validate.add_argument("--require-ultrametric", action="store_true")
    validate.add_argument("--strict", action="store_true")
    validate.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(validate)

    inspect = subparsers.add_parser(get_command_spec("inspect").name, help=get_command_spec("inspect").summary)
    inspect.add_argument("tree", type=Path)
    inspect.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    inspect.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(inspect)

    normalize = subparsers.add_parser(get_command_spec("normalize").name, help=get_command_spec("normalize").summary)
    normalize.add_argument("tree", type=Path)
    normalize.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    normalize.add_argument("--out", required=True, type=Path)
    normalize.add_argument("--json", action="store_true", help="Emit the normalization result as JSON.")
    _add_manifest_argument(normalize)

    normalize_taxa = subparsers.add_parser(
        get_command_spec("normalize-taxa").name,
        help=get_command_spec("normalize-taxa").summary,
    )
    normalize_taxa.add_argument("tree", type=Path)
    normalize_taxa.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    normalize_taxa.add_argument("--policy", choices=("spaces-to-underscores",), required=True)
    normalize_taxa.add_argument("--out", required=True, type=Path)
    normalize_taxa.add_argument("--mapping-out", type=Path)
    normalize_taxa.add_argument("--json", action="store_true", help="Emit the normalization result as JSON.")
    _add_manifest_argument(normalize_taxa)

    taxonomy = subparsers.add_parser(get_command_spec("taxonomy").name, help=get_command_spec("taxonomy").summary)
    taxonomy_subparsers = taxonomy.add_subparsers(dest="taxonomy_command", required=True)
    taxonomy_synonyms = taxonomy_subparsers.add_parser(
        "synonyms",
        help="Audit a tree against a configurable taxon synonym table.",
    )
    taxonomy_synonyms.add_argument("tree", type=Path)
    taxonomy_synonyms.add_argument("--synonym-table", required=True, type=Path)
    taxonomy_synonyms.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    taxonomy_synonyms.add_argument("--json", action="store_true", help="Emit the audit as JSON.")
    _add_manifest_argument(taxonomy_synonyms)
    taxonomy_resolve = taxonomy_subparsers.add_parser(
        "resolve-synonyms",
        help="Resolve tree tip synonyms to accepted labels with a reversible mapping artifact.",
    )
    taxonomy_resolve.add_argument("tree", type=Path)
    taxonomy_resolve.add_argument("--synonym-table", required=True, type=Path)
    taxonomy_resolve.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    taxonomy_resolve.add_argument("--resolution-policy", choices=("reject-ambiguous",), default="reject-ambiguous")
    taxonomy_resolve.add_argument("--out", required=True, type=Path)
    taxonomy_resolve.add_argument("--mapping-out", type=Path)
    taxonomy_resolve.add_argument("--json", action="store_true", help="Emit the resolution report as JSON.")
    _add_manifest_argument(taxonomy_resolve)
    taxonomy_namespaces = taxonomy_subparsers.add_parser(
        "namespaces",
        help="Classify tree tip labels into accession, species, sample, isolate, or user-defined namespaces.",
    )
    taxonomy_namespaces.add_argument("tree", type=Path)
    taxonomy_namespaces.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    taxonomy_namespaces.add_argument("--json", action="store_true", help="Emit the namespace report as JSON.")
    _add_manifest_argument(taxonomy_namespaces)
    taxonomy_loss = taxonomy_subparsers.add_parser(
        "loss",
        help="Trace where taxa were lost across tree, alignment, metadata, traits, inference, and reporting stages.",
    )
    taxonomy_loss.add_argument("tree", type=Path)
    taxonomy_loss.add_argument("--metadata", required=True, type=Path)
    taxonomy_loss.add_argument("--traits", required=True, type=Path)
    taxonomy_loss.add_argument("--alignment", type=Path)
    taxonomy_loss.add_argument("--filtered-alignment", type=Path)
    taxonomy_loss.add_argument("--inference-tree", type=Path)
    taxonomy_loss.add_argument("--reported-taxa", type=Path)
    taxonomy_loss.add_argument("--json", action="store_true", help="Emit the workflow loss report as JSON.")
    _add_manifest_argument(taxonomy_loss)
    taxonomy_stability = taxonomy_subparsers.add_parser(
        "stability",
        help="Compare taxon retention across multiple named workflow artifacts.",
    )
    taxonomy_stability.add_argument(
        "--run",
        action="append",
        required=True,
        metavar="LABEL=PATH",
        help="Named run source in the form label=/path/to/tree_or_alignment_or_table",
    )
    taxonomy_stability.add_argument("--json", action="store_true", help="Emit the stability report as JSON.")
    _add_manifest_argument(taxonomy_stability)

    topology = subparsers.add_parser(get_command_spec("topology").name, help=get_command_spec("topology").summary)
    topology_subparsers = topology.add_subparsers(dest="topology_command", required=True)
    topology_outgroup = topology_subparsers.add_parser("root-outgroup", help="Root a tree on explicit outgroup taxa.")
    topology_outgroup.add_argument("tree", type=Path)
    topology_outgroup.add_argument("--taxa", nargs="+", required=True)
    topology_outgroup.add_argument("--out", required=True, type=Path)
    topology_outgroup.add_argument("--json", action="store_true", help="Emit the rooting report as JSON.")
    _add_manifest_argument(topology_outgroup)
    topology_midpoint = topology_subparsers.add_parser("reroot-midpoint", help="Reroot a tree by midpoint.")
    topology_midpoint.add_argument("tree", type=Path)
    topology_midpoint.add_argument("--out", required=True, type=Path)
    topology_midpoint.add_argument("--json", action="store_true", help="Emit the rerooting report as JSON.")
    _add_manifest_argument(topology_midpoint)
    topology_unroot = topology_subparsers.add_parser("unroot", help="Convert a rooted tree into an explicit unrooted trifurcation.")
    topology_unroot.add_argument("tree", type=Path)
    topology_unroot.add_argument("--out", required=True, type=Path)
    topology_unroot.add_argument("--json", action="store_true", help="Emit the unrooting report as JSON.")
    _add_manifest_argument(topology_unroot)

    compare = subparsers.add_parser(get_command_spec("compare").name, help=get_command_spec("compare").summary)
    compare.add_argument("left")
    compare.add_argument("right")
    compare.add_argument("third", nargs="?")
    compare.add_argument("--out", type=Path)
    compare.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(compare)

    annotate = subparsers.add_parser(get_command_spec("annotate").name, help=get_command_spec("annotate").summary)
    annotate.add_argument("tree", type=Path)
    annotate.add_argument("--metadata", required=True, type=Path)
    annotate.add_argument("--taxon-column")
    annotate.add_argument("--out", type=Path)
    annotate.add_argument("--joined-out", type=Path)
    annotate.add_argument("--json", action="store_true", help="Emit the linkage report as JSON.")
    _add_manifest_argument(annotate)

    diagnose = subparsers.add_parser(get_command_spec("diagnose").name, help=get_command_spec("diagnose").summary)
    diagnose.add_argument("target")
    diagnose.add_argument("tree", nargs="?", type=Path)
    diagnose.add_argument("--metadata", type=Path)
    diagnose.add_argument("--taxon-column")
    diagnose.add_argument("--out", type=Path)
    diagnose.add_argument("--tolerance", type=float, default=1e-6)
    diagnose.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(diagnose)

    render = subparsers.add_parser(get_command_spec("render").name, help=get_command_spec("render").summary)
    render.add_argument("tree", type=Path)
    render.add_argument("--metadata", type=Path)
    render.add_argument("--traits", type=Path)
    render.add_argument("--taxon-column")
    render.add_argument("--label-column")
    render.add_argument("--layout", choices=["cladogram", "phylogram", "circular"], default="cladogram")
    render.add_argument("--support-labels", action="store_true")
    render.add_argument("--categorical-column")
    render.add_argument("--continuous-column")
    render.add_argument("--metadata-strip-columns")
    render.add_argument("--heatmap-columns")
    render.add_argument("--collapse-clades")
    render.add_argument("--package-dir", type=Path)
    render.add_argument("--out", required=True, type=Path)
    render.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(render)

    evidence = subparsers.add_parser(get_command_spec("evidence").name, help=get_command_spec("evidence").summary)
    evidence_subparsers = evidence.add_subparsers(dest="evidence_command", required=True)
    evidence_bundle = evidence_subparsers.add_parser(
        "bundle",
        help="Bundle explicit phylogenetics inputs and outputs as evidence.",
    )
    evidence_bundle.add_argument("--inputs", nargs="+", required=True, type=Path)
    evidence_bundle.add_argument("--outputs", nargs="+", required=True, type=Path)
    evidence_bundle.add_argument("--out", required=True, type=Path)
    evidence_bundle.add_argument("--json", action="store_true", help="Emit the bundle report as JSON.")
    _add_manifest_argument(evidence_bundle)
    evidence_validate = evidence_subparsers.add_parser("validate", help="Validate an existing evidence bundle.")
    evidence_validate.add_argument("bundle_root", type=Path)
    evidence_validate.add_argument("--json", action="store_true", help="Emit the validation report as JSON.")
    _add_manifest_argument(evidence_validate)

    report = subparsers.add_parser(get_command_spec("report").name, help=get_command_spec("report").summary)
    report_subparsers = report.add_subparsers(dest="report_command", required=True)
    report_tree = report_subparsers.add_parser("tree", help="Render a deterministic single-tree HTML report.")
    report_tree.add_argument("tree", type=Path)
    report_tree.add_argument("--out", required=True, type=Path)
    report_tree.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(report_tree)
    report_dataset = report_subparsers.add_parser("dataset", help="Render a tree plus table dataset HTML report.")
    report_dataset.add_argument("--tree", required=True, type=Path)
    report_dataset.add_argument("--metadata", required=True, type=Path)
    report_dataset.add_argument("--traits", type=Path)
    report_dataset.add_argument("--alignment", type=Path)
    report_dataset.add_argument("--tip-dates", type=Path)
    report_dataset.add_argument("--calibrations", type=Path)
    report_dataset.add_argument("--out", required=True, type=Path)
    report_dataset.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(report_dataset)
    report_phylo_inputs = report_subparsers.add_parser(
        "phylo-inputs",
        help="Render a tree plus alignment HTML input report.",
    )
    report_phylo_inputs.add_argument("--tree", required=True, type=Path)
    report_phylo_inputs.add_argument("--alignment", required=True, type=Path)
    report_phylo_inputs.add_argument("--out", required=True, type=Path)
    report_phylo_inputs.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(report_phylo_inputs)

    demo = subparsers.add_parser(get_command_spec("demo").name, help=get_command_spec("demo").summary)
    demo_subparsers = demo.add_subparsers(dest="demo_command", required=True)
    demo_run = demo_subparsers.add_parser("run", help="Run the repository capability demo workflow.")
    demo_run.add_argument("--out", required=True, type=Path)
    demo_run.add_argument("--json", action="store_true", help="Emit the demo result as JSON.")
    _add_manifest_argument(demo_run)

    adapter = subparsers.add_parser(get_command_spec("adapter").name, help=get_command_spec("adapter").summary)
    adapter_subparsers = adapter.add_subparsers(dest="adapter_command", required=True)
    adapter_inspect = adapter_subparsers.add_parser("inspect", help="Report external engine version metadata.")
    adapter_inspect.add_argument("engine_name", choices=("mafft", "trimal", "iqtree", "FastTree", "MrBayes"))
    adapter_inspect.add_argument("--executable", type=str)
    adapter_inspect.add_argument("--json", action="store_true", help="Emit the adapter report as JSON.")
    _add_manifest_argument(adapter_inspect)
    adapter_report = adapter_subparsers.add_parser("report", help="Render an HTML report from an engine workflow manifest.")
    adapter_report.add_argument("manifest_path", type=Path)
    adapter_report.add_argument("--out", required=True, type=Path)
    adapter_report.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(adapter_report)
    adapter_align = adapter_subparsers.add_parser("align", help="Run multiple-sequence alignment on unaligned FASTA.")
    adapter_align.add_argument("input_path", type=Path)
    adapter_align.add_argument("--out", required=True, type=Path)
    adapter_align.add_argument("--executable", type=str)
    adapter_align.add_argument("--json", action="store_true", help="Emit the workflow report as JSON.")
    _add_manifest_argument(adapter_align)
    adapter_trim = adapter_subparsers.add_parser("trim", help="Run external alignment trimming.")
    adapter_trim.add_argument("input_path", type=Path)
    adapter_trim.add_argument("--out", required=True, type=Path)
    adapter_trim.add_argument("--gap-threshold", type=float, default=0.1)
    adapter_trim.add_argument("--executable", type=str)
    adapter_trim.add_argument("--json", action="store_true", help="Emit the workflow report as JSON.")
    _add_manifest_argument(adapter_trim)
    adapter_model = adapter_subparsers.add_parser("model-select", help="Run external sequence-model selection.")
    adapter_model.add_argument("input_path", type=Path)
    adapter_model.add_argument("--out-dir", required=True, type=Path)
    adapter_model.add_argument("--prefix", default="model-selection")
    adapter_model.add_argument("--sequence-type", choices=("dna", "rna", "protein", "unknown"))
    adapter_model.add_argument("--executable", type=str)
    adapter_model.add_argument("--json", action="store_true", help="Emit the workflow report as JSON.")
    _add_manifest_argument(adapter_model)
    adapter_ml = adapter_subparsers.add_parser("infer-ml", help="Run maximum-likelihood tree inference.")
    adapter_ml.add_argument("input_path", type=Path)
    adapter_ml.add_argument("--out-dir", required=True, type=Path)
    adapter_ml.add_argument("--model", required=True)
    adapter_ml.add_argument("--prefix", default="maximum-likelihood")
    adapter_ml.add_argument("--sequence-type", choices=("dna", "rna", "protein", "unknown"))
    adapter_ml.add_argument("--executable", type=str)
    adapter_ml.add_argument("--json", action="store_true", help="Emit the workflow report as JSON.")
    _add_manifest_argument(adapter_ml)
    adapter_bootstrap = adapter_subparsers.add_parser("bootstrap", help="Run bootstrap support estimation.")
    adapter_bootstrap.add_argument("input_path", type=Path)
    adapter_bootstrap.add_argument("--out-dir", required=True, type=Path)
    adapter_bootstrap.add_argument("--model", required=True)
    adapter_bootstrap.add_argument("--replicates", type=int, default=1000)
    adapter_bootstrap.add_argument("--prefix", default="bootstrap-support")
    adapter_bootstrap.add_argument("--sequence-type", choices=("dna", "rna", "protein", "unknown"))
    adapter_bootstrap.add_argument("--executable", type=str)
    adapter_bootstrap.add_argument("--json", action="store_true", help="Emit the workflow report as JSON.")
    _add_manifest_argument(adapter_bootstrap)
    adapter_consensus = adapter_subparsers.add_parser("consensus", help="Build a consensus tree from bootstrap trees.")
    adapter_consensus.add_argument("input_path", type=Path)
    adapter_consensus.add_argument("--out-dir", required=True, type=Path)
    adapter_consensus.add_argument("--prefix", default="bootstrap-consensus")
    adapter_consensus.add_argument("--minimum-support", type=float, default=0.5)
    adapter_consensus.add_argument("--executable", type=str)
    adapter_consensus.add_argument("--json", action="store_true", help="Emit the workflow report as JSON.")
    _add_manifest_argument(adapter_consensus)
    adapter_fast = adapter_subparsers.add_parser("infer-fast", help="Run fast approximate tree inference.")
    adapter_fast.add_argument("input_path", type=Path)
    adapter_fast.add_argument("--out", required=True, type=Path)
    adapter_fast.add_argument("--sequence-type", choices=("dna", "rna", "protein", "unknown"))
    adapter_fast.add_argument("--executable", type=str)
    adapter_fast.add_argument("--json", action="store_true", help="Emit the workflow report as JSON.")
    _add_manifest_argument(adapter_fast)
    adapter_compare = adapter_subparsers.add_parser("compare", help="Compare fast approximate and ML trees.")
    adapter_compare.add_argument("--fast-tree", required=True, type=Path)
    adapter_compare.add_argument("--ml-tree", required=True, type=Path)
    adapter_compare.add_argument("--out", required=True, type=Path)
    adapter_compare.add_argument("--json", action="store_true", help="Emit the comparison report as JSON.")
    _add_manifest_argument(adapter_compare)
    adapter_mrbayes_prepare = adapter_subparsers.add_parser(
        "mrbayes-prepare",
        help="Prepare a MrBayes NEXUS analysis from an aligned FASTA file.",
    )
    adapter_mrbayes_prepare.add_argument("input_path", type=Path)
    adapter_mrbayes_prepare.add_argument("--out", required=True, type=Path)
    adapter_mrbayes_prepare.add_argument("--model", default="gtr")
    adapter_mrbayes_prepare.add_argument("--rates", default="gamma")
    adapter_mrbayes_prepare.add_argument("--ngen", type=int, default=10000)
    adapter_mrbayes_prepare.add_argument("--nchains", type=int, default=4)
    adapter_mrbayes_prepare.add_argument("--samplefreq", type=int, default=100)
    adapter_mrbayes_prepare.add_argument("--printfreq", type=int, default=100)
    adapter_mrbayes_prepare.add_argument("--burnin-fraction", type=float, default=0.25)
    adapter_mrbayes_prepare.add_argument("--json", action="store_true", help="Emit the preparation report as JSON.")
    _add_manifest_argument(adapter_mrbayes_prepare)
    adapter_mrbayes_run = adapter_subparsers.add_parser(
        "mrbayes-run",
        help="Run a prepared MrBayes posterior inference workflow.",
    )
    adapter_mrbayes_run.add_argument("input_path", type=Path)
    adapter_mrbayes_run.add_argument("--executable", type=str)
    adapter_mrbayes_run.add_argument("--resume", action="store_true")
    adapter_mrbayes_run.add_argument("--json", action="store_true", help="Emit the workflow report as JSON.")
    _add_manifest_argument(adapter_mrbayes_run)
    adapter_mrbayes_summarize = adapter_subparsers.add_parser(
        "mrbayes-summarize",
        help="Summarize MrBayes posterior trees after burn-in removal.",
    )
    adapter_mrbayes_summarize.add_argument("input_path", type=Path)
    adapter_mrbayes_summarize.add_argument("--burnin-fraction", type=float, default=0.25)
    adapter_mrbayes_summarize.add_argument("--json", action="store_true", help="Emit the posterior summary as JSON.")
    _add_manifest_argument(adapter_mrbayes_summarize)
    adapter_mrbayes_traces = adapter_subparsers.add_parser(
        "mrbayes-traces",
        help="Parse a MrBayes parameter trace table.",
    )
    adapter_mrbayes_traces.add_argument("input_path", type=Path)
    adapter_mrbayes_traces.add_argument("--json", action="store_true", help="Emit the trace report as JSON.")
    _add_manifest_argument(adapter_mrbayes_traces)
    adapter_mrbayes_ess = adapter_subparsers.add_parser(
        "mrbayes-ess",
        help="Compute effective sample sizes from a MrBayes trace table.",
    )
    adapter_mrbayes_ess.add_argument("input_path", type=Path)
    adapter_mrbayes_ess.add_argument("--json", action="store_true", help="Emit the ESS report as JSON.")
    _add_manifest_argument(adapter_mrbayes_ess)
    adapter_mrbayes_convergence = adapter_subparsers.add_parser(
        "mrbayes-convergence",
        help="Assess MrBayes trace convergence from ESS and trace drift.",
    )
    adapter_mrbayes_convergence.add_argument("input_path", type=Path)
    adapter_mrbayes_convergence.add_argument("--ess-threshold", type=float, default=200.0)
    adapter_mrbayes_convergence.add_argument("--mean-shift-threshold", type=float, default=0.5)
    adapter_mrbayes_convergence.add_argument("--json", action="store_true", help="Emit the convergence report as JSON.")
    _add_manifest_argument(adapter_mrbayes_convergence)
    adapter_mrbayes_report = adapter_subparsers.add_parser(
        "mrbayes-report",
        help="Render an HTML Bayesian posterior report from posterior trees and traces.",
    )
    adapter_mrbayes_report.add_argument("posterior_trees", type=Path)
    adapter_mrbayes_report.add_argument("--traces", required=True, type=Path)
    adapter_mrbayes_report.add_argument("--out", required=True, type=Path)
    adapter_mrbayes_report.add_argument("--burnin-fraction", type=float, default=0.25)
    adapter_mrbayes_report.add_argument("--ess-threshold", type=float, default=200.0)
    adapter_mrbayes_report.add_argument("--mean-shift-threshold", type=float, default=0.5)
    adapter_mrbayes_report.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(adapter_mrbayes_report)
    adapter_beast_prepare = adapter_subparsers.add_parser(
        "beast-prepare",
        help="Prepare a BEAST-style time-tree XML analysis from aligned sequences and dating inputs.",
    )
    adapter_beast_prepare.add_argument("input_path", type=Path)
    adapter_beast_prepare.add_argument("--out", required=True, type=Path)
    adapter_beast_prepare.add_argument("--tree", type=Path)
    adapter_beast_prepare.add_argument("--calibrations", type=Path)
    adapter_beast_prepare.add_argument("--tip-dates", type=Path)
    adapter_beast_prepare.add_argument("--clock-model", default="strict")
    adapter_beast_prepare.add_argument("--tree-prior", default="yule")
    adapter_beast_prepare.add_argument("--chain-length", type=int, default=1000000)
    adapter_beast_prepare.add_argument("--log-every", type=int, default=1000)
    adapter_beast_prepare.add_argument("--json", action="store_true", help="Emit the preparation report as JSON.")
    _add_manifest_argument(adapter_beast_prepare)
    adapter_beast_calibrations = adapter_subparsers.add_parser(
        "beast-calibrations",
        help="Validate a fossil calibration table against a tree.",
    )
    adapter_beast_calibrations.add_argument("tree_path", type=Path)
    adapter_beast_calibrations.add_argument("calibration_path", type=Path)
    adapter_beast_calibrations.add_argument("--json", action="store_true", help="Emit the validation report as JSON.")
    _add_manifest_argument(adapter_beast_calibrations)
    adapter_beast_tip_dates = adapter_subparsers.add_parser(
        "beast-tip-dates",
        help="Validate tip-dating metadata against a tree and optional alignment.",
    )
    adapter_beast_tip_dates.add_argument("tree_path", type=Path)
    adapter_beast_tip_dates.add_argument("tip_dates_path", type=Path)
    adapter_beast_tip_dates.add_argument("--alignment", type=Path)
    adapter_beast_tip_dates.add_argument("--date-column", default="date")
    adapter_beast_tip_dates.add_argument("--json", action="store_true", help="Emit the validation report as JSON.")
    _add_manifest_argument(adapter_beast_tip_dates)
    adapter_beast_log = adapter_subparsers.add_parser(
        "beast-log",
        help="Parse a BEAST log file into a deterministic numeric trace table.",
    )
    adapter_beast_log.add_argument("input_path", type=Path)
    adapter_beast_log.add_argument("--json", action="store_true", help="Emit the parsed log report as JSON.")
    _add_manifest_argument(adapter_beast_log)
    adapter_beast_convergence = adapter_subparsers.add_parser(
        "beast-convergence",
        help="Assess BEAST log convergence from ESS and trace drift.",
    )
    adapter_beast_convergence.add_argument("input_path", type=Path)
    adapter_beast_convergence.add_argument("--ess-threshold", type=float, default=200.0)
    adapter_beast_convergence.add_argument("--mean-shift-threshold", type=float, default=0.5)
    adapter_beast_convergence.add_argument("--json", action="store_true", help="Emit the convergence report as JSON.")
    _add_manifest_argument(adapter_beast_convergence)
    adapter_beast_calibration_report = adapter_subparsers.add_parser(
        "beast-calibration-report",
        help="Render an HTML calibration audit report.",
    )
    adapter_beast_calibration_report.add_argument("tree_path", type=Path)
    adapter_beast_calibration_report.add_argument("calibration_path", type=Path)
    adapter_beast_calibration_report.add_argument("--out", required=True, type=Path)
    adapter_beast_calibration_report.add_argument("--tip-dates", type=Path)
    adapter_beast_calibration_report.add_argument("--alignment", type=Path)
    adapter_beast_calibration_report.add_argument("--date-column", default="date")
    adapter_beast_calibration_report.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(adapter_beast_calibration_report)
    adapter_bayesian_evidence = adapter_subparsers.add_parser(
        "bayesian-evidence",
        help="Bundle Bayesian configs, trees, logs, diagnostics, and reports into one evidence package.",
    )
    adapter_bayesian_evidence.add_argument("--out-dir", required=True, type=Path)
    adapter_bayesian_evidence.add_argument("--inputs", nargs="+", required=True, type=Path)
    adapter_bayesian_evidence.add_argument("--configs", nargs="+", required=True, type=Path)
    adapter_bayesian_evidence.add_argument("--trees", nargs="+", required=True, type=Path)
    adapter_bayesian_evidence.add_argument("--logs", nargs="+", required=True, type=Path)
    adapter_bayesian_evidence.add_argument("--diagnostics", nargs="+", required=True, type=Path)
    adapter_bayesian_evidence.add_argument("--reports", nargs="+", required=True, type=Path)
    adapter_bayesian_evidence.add_argument("--json", action="store_true", help="Emit the evidence-package report as JSON.")
    _add_manifest_argument(adapter_bayesian_evidence)

    return parser


def run_command(args: Any, *, parser: argparse.ArgumentParser) -> int:
    """Run the selected command."""
    try:
        if args.command == "commands":
            _print_commands(output_format=args.format)
            return 0
        if args.command == "env":
            report = inspect_environment()
            outputs = _finalize_outputs(args, command="env", inputs=[])
            _print_result(
                build_command_result(
                    command="env",
                    inputs=[],
                    outputs=outputs,
                    metrics={"dependency_count": len(report.dependencies)},
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "metadata":
            report = inspect_metadata_table(args.table, taxon_column=args.taxon_column)
            outputs = _finalize_outputs(args, command="metadata", inputs=[args.table])
            _print_result(
                build_command_result(
                    command="metadata",
                    inputs=[args.table],
                    outputs=outputs,
                    metrics={
                        "row_count": report.row_count,
                        "column_count": report.column_count,
                        "taxon_count": len(report.taxa),
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "traits":
            if args.traits_command == "validate":
                report = validate_traits_table(args.table, taxon_column=args.taxon_column)
                outputs = _finalize_outputs(args, command="traits", inputs=[args.table])
                _print_result(
                    build_command_result(
                        command="traits",
                        inputs=[args.table],
                        outputs=outputs,
                        metrics={
                            "row_count": report.row_count,
                            "trait_column_count": len(report.trait_columns),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.traits_command == "missing":
                report = detect_missing_trait_values(args.table, taxon_column=args.taxon_column)
                outputs = _finalize_outputs(args, command="traits", inputs=[args.table])
                _print_result(
                    build_command_result(
                        command="traits",
                        inputs=[args.table],
                        outputs=outputs,
                        metrics={"missing_value_count": len(report.missing_values)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.traits_command == "prune":
                rows, report = prune_traits_to_tree(
                    args.tree,
                    args.table,
                    taxon_column=args.taxon_column,
                )
                table = load_taxon_table(args.table, taxon_column=args.taxon_column)
                output_path = write_taxon_rows(args.out, columns=table.columns, rows=rows)
                outputs = _finalize_outputs(
                    args,
                    command="traits",
                    inputs=[args.tree, args.table],
                    outputs=[output_path],
                )
                _print_result(
                    build_command_result(
                        command="traits",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "original_row_count": report.original_row_count,
                            "kept_taxa": len(report.kept_taxa),
                            "removed_taxa": len(report.removed_taxa),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            report = link_tree_to_traits(
                args.tree,
                args.table,
                taxon_column=args.taxon_column,
                strict=args.strict,
            )
            outputs = _finalize_outputs(args, command="traits", inputs=[args.tree, args.table])
            _print_result(
                build_command_result(
                    command="traits",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                    metrics={
                        "tree_taxa": report.tree_taxa,
                        "trait_taxa": report.trait_taxa,
                        "linked_taxa": report.linked_taxa,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "validate":
            report = validate_tree_path(
                args.tree,
                source_format=args.format,
                allow_duplicates=args.allow_duplicates,
                strict=args.strict,
                allow_negative_branch_lengths=args.allow_negative_branches,
                require_rooted=args.require_rooted,
                require_ultrametric=args.require_ultrametric,
            )
            outputs = _finalize_outputs(args, command="validate", inputs=[args.tree])
            _print_result(
                build_command_result(
                    command="validate",
                    inputs=[args.tree],
                    outputs=outputs,
                    warnings=report.warnings,
                    metrics={
                        "tip_count": report.tip_count,
                        "internal_node_count": report.internal_node_count,
                        "validity_decision": report.validity_decision,
                        "syntax_valid": report.syntax_valid,
                        "biologically_safe": report.biologically_safe,
                        "polytomy_count": report.polytomy_count,
                        "missing_internal_branch_count": len(report.missing_internal_branch_nodes),
                        "missing_terminal_branch_count": len(report.missing_terminal_branch_taxa),
                        "singleton_internal_node_count": len(report.singleton_internal_nodes),
                        "integrity_issue_count": len(report.integrity_issues),
                        "unsafe_external_label_count": len(report.unsafe_external_labels),
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "prune":
            if args.keep_from is not None:
                tree, report = prune_tree_to_taxa(args.tree, args.keep_from, taxon_column=args.taxon_column)
                prune_inputs = [args.tree, args.keep_from]
            elif args.exclude_taxa is not None:
                tree, report = drop_tree_taxa(args.tree, list(args.exclude_taxa))
                prune_inputs = [args.tree]
            else:
                tree, report = prune_tree_to_requested_taxa(args.tree, list(args.taxa))
                prune_inputs = [args.tree]
            output_path = write_newick(args.out, tree)
            pruned_taxa_path = args.pruned_taxa_out or args.out.with_name("pruned_taxa.tsv")
            write_pruned_taxa(pruned_taxa_path, report.removed_taxa)
            outputs = _finalize_outputs(
                args,
                command="prune",
                inputs=prune_inputs,
                outputs=[output_path, pruned_taxa_path],
            )
            _print_result(
                build_command_result(
                    command="prune",
                    inputs=prune_inputs,
                    outputs=outputs,
                    metrics={
                        "kept_taxa": len(report.kept_taxa),
                        "removed_taxa": len(report.removed_taxa),
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "alignment":
            if args.alignment_command == "alphabet":
                records = load_fasta_alignment(args.alignment)
                alphabet = infer_alignment_alphabet(records)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={"alphabet": alphabet, "sequence_count": len(records)},
                        data={"alignment_path": args.alignment, "inferred_alphabet": alphabet},
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "profiles":
                report = list_alignment_filter_profiles()
                outputs = _finalize_outputs(args, command="alignment", inputs=[])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[],
                        outputs=outputs,
                        metrics={"profile_count": len(report)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "gc":
                report = summarise_fasta(args.alignment)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "alphabet": report.inferred_alphabet,
                            "gc_sequence_count": len(report.per_sequence_gc_content),
                        },
                        data={
                            "alignment_path": report.path,
                            "inferred_alphabet": report.inferred_alphabet,
                            "per_sequence_gc_content": report.per_sequence_gc_content,
                            "whole_alignment_gc_content": report.whole_alignment_gc_content,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "inspect":
                report = summarise_fasta(args.alignment)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "sequence_count": report.sequence_count,
                            "alignment_length": report.alignment_length,
                            "variable_site_count": report.variable_site_count,
                            "parsimony_informative_site_count": report.parsimony_informative_site_count,
                            "alphabet": report.inferred_alphabet,
                            "invalid_character_count": len(report.invalid_characters),
                            "composition_outlier_count": len(report.composition_outliers),
                            "duplicate_group_count": len(report.duplicate_sequence_groups),
                            "near_duplicate_count": len(report.near_duplicate_pairs),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "classify":
                report = classify_alignment_sequences(args.alignment)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "sequence_count": report.sequence_count,
                            "min_sequence_length": report.min_sequence_length,
                            "max_sequence_length": report.max_sequence_length,
                            "state": report.state,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "quality":
                report = build_alignment_quality_report(args.alignment)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "invalid_character_count": len(report.invalid_characters),
                            "composition_outlier_count": len(report.composition_outliers),
                            "sequence_length_outlier_count": len(report.sequence_length_outliers),
                            "duplicate_group_count": len(report.duplicate_sequence_groups),
                            "near_duplicate_count": len(report.near_duplicate_pairs),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "windows":
                windows = summarize_alignment_windows(
                    args.alignment,
                    window_size=args.window_size,
                    step_size=args.step_size,
                )
                over_aligned = detect_over_aligned_regions(
                    args.alignment,
                    window_size=args.window_size,
                    step_size=args.step_size,
                )
                under_aligned = detect_under_aligned_regions(
                    args.alignment,
                    window_size=args.window_size,
                    step_size=args.step_size,
                )
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        warnings=[region.note for region in over_aligned + under_aligned],
                        metrics={
                            "window_count": len(windows),
                            "over_aligned_region_count": len(over_aligned),
                            "under_aligned_region_count": len(under_aligned),
                        },
                        data={
                            "windows": windows,
                            "over_aligned_regions": over_aligned,
                            "under_aligned_regions": under_aligned,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "readiness":
                report = summarize_alignment_readiness(args.alignment)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "sequence_count": report.sequence_count,
                            "alignment_length": report.alignment_length,
                            "ready_method_count": sum(1 for method in report.methods if method.ready),
                            "blocked_method_count": sum(1 for method in report.methods if not method.ready),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "forensic":
                report = build_alignment_forensic_report(args.alignment)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "quality_score": report.quality.quality_score,
                            "safe_for_distance_analysis": report.safe_for_distance_analysis,
                            "safe_for_maximum_likelihood": report.safe_for_maximum_likelihood,
                            "safe_for_bayesian_inference": report.safe_for_bayesian_inference,
                            "safe_for_coding_analysis": report.safe_for_coding_analysis,
                            "safe_for_publication": report.safe_for_publication,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "composition":
                report = summarise_fasta(args.alignment)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "alphabet": report.inferred_alphabet,
                            "sequence_count": report.sequence_count,
                            "gc_sequence_count": len(report.per_sequence_gc_content),
                        },
                        data={
                            "path": report.path,
                            "inferred_alphabet": report.inferred_alphabet,
                            "nucleotide_composition": report.nucleotide_composition,
                            "amino_acid_composition": report.amino_acid_composition,
                            "per_sequence_gc_content": report.per_sequence_gc_content,
                            "whole_alignment_gc_content": report.whole_alignment_gc_content,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "invalid":
                report = detect_invalid_alignment_characters(args.alignment, alphabet=args.alphabet)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={"invalid_character_count": len(report), "alphabet": args.alphabet},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "duplicates":
                duplicates = detect_identical_duplicate_sequences(args.alignment)
                near_duplicates = detect_near_duplicate_sequences(
                    args.alignment,
                    identity_threshold=args.identity_threshold,
                )
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "duplicate_group_count": len(duplicates),
                            "near_duplicate_count": len(near_duplicates),
                            "identity_threshold": args.identity_threshold,
                        },
                        data={
                            "duplicate_sequence_groups": duplicates,
                            "near_duplicate_pairs": near_duplicates,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "outliers":
                report = detect_composition_outlier_sequences(
                    args.alignment,
                    deviation_threshold=args.deviation_threshold,
                )
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "composition_outlier_count": len(report),
                            "deviation_threshold": args.deviation_threshold,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "length-outliers":
                report = detect_sequence_length_outliers(args.alignment)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={"sequence_length_outlier_count": len(report)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "filter":
                group_columns = None if not args.group_columns else [column.strip() for column in args.group_columns.split(",") if column.strip()]
                records, report = clean_alignment_with_profile(
                    args.alignment,
                    profile_name=args.profile,
                    group_table_path=args.group_table,
                    group_columns=group_columns,
                )
                output_path = write_fasta_alignment(args.out, records)
                filter_inputs = [args.alignment, *([args.group_table] if args.group_table is not None else [])]
                outputs = _finalize_outputs(
                    args,
                    command="alignment",
                    inputs=filter_inputs,
                    outputs=[output_path],
                )
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=filter_inputs,
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "profile": report.profile.name,
                            "trimmed_sequence_count": report.trim.trimmed_sequence_count,
                            "trimmed_alignment_length": report.trim.trimmed_alignment_length,
                            "signal_warning_count": len(report.signal_warnings),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "compare":
                report = compare_alignment_versions(args.left_alignment, args.right_alignment)
                outputs = _finalize_outputs(
                    args,
                    command="alignment",
                    inputs=[args.left_alignment, args.right_alignment],
                )
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.left_alignment, args.right_alignment],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "shared_taxa": len(report.shared_taxa),
                            "left_only_taxa": len(report.left_only_taxa),
                            "right_only_taxa": len(report.right_only_taxa),
                            "left_alignment_length": report.left_alignment_length,
                            "right_alignment_length": report.right_alignment_length,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "trim":
                records, report = trim_alignment(
                    args.alignment,
                    remove_all_gap_sites=not args.keep_all_gap_sites,
                    remove_all_missing_sites=not args.keep_all_missing_sites,
                    site_missingness_threshold=args.site_missingness_threshold,
                    sequence_missingness_threshold=args.sequence_missingness_threshold,
                )
                output_path = write_fasta_alignment(args.out, records)
                outputs = _finalize_outputs(
                    args,
                    command="alignment",
                    inputs=[args.alignment],
                    outputs=[output_path],
                )
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "trimmed_sequence_count": report.trimmed_sequence_count,
                            "trimmed_alignment_length": report.trimmed_alignment_length,
                            "removed_column_count": len(report.removed_columns),
                            "removed_sequence_count": len(report.removed_sequences),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "identity-matrix":
                report = compute_pairwise_sequence_identity_matrix(args.alignment)
                outputs: list[Path | str] = []
                if args.out is not None:
                    outputs.append(write_sequence_identity_matrix(args.out, report))
                outputs = _finalize_outputs(
                    args,
                    command="alignment",
                    inputs=[args.alignment],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "sequence_count": len(report.identifiers),
                            "pair_count": len(report.pairs),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "distance-matrix":
                report = compute_pairwise_genetic_distance_matrix(
                    args.alignment,
                    model=args.model,
                    gap_handling=args.gap_handling,
                    ambiguity_policy=args.ambiguity_policy,
                )
                outputs: list[Path | str] = []
                if args.out is not None:
                    outputs.append(write_genetic_distance_matrix(args.out, report))
                outputs = _finalize_outputs(
                    args,
                    command="alignment",
                    inputs=[args.alignment],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "sequence_count": len(report.identifiers),
                            "pair_count": len(report.pairs),
                            "model": report.model,
                            "gap_handling": report.gap_handling,
                            "ambiguity_policy": report.ambiguity_policy,
                            "alphabet": report.inferred_alphabet,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "distance-quality":
                report = inspect_distance_matrix_quality(
                    args.alignment,
                    model=args.model,
                    gap_handling=args.gap_handling,
                    ambiguity_policy=args.ambiguity_policy,
                )
                outputs = _finalize_outputs(
                    args,
                    command="alignment",
                    inputs=[args.alignment],
                )
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "saturated_pair_count": len(report.saturated_pairs),
                            "low_information_pair_count": len(report.low_information_pairs),
                            "decision": report.method_assessment.decision,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "distance-assumptions":
                report = assess_distance_method_assumptions(
                    args.alignment,
                    model=args.model,
                    gap_handling=args.gap_handling,
                    ambiguity_policy=args.ambiguity_policy,
                )
                outputs = _finalize_outputs(
                    args,
                    command="alignment",
                    inputs=[args.alignment],
                )
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "ultrametric_compatible": report.ultrametric_compatible,
                            "upgma_violation_count": len(report.upgma_ultrametric_violations),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "build-tree":
                tree, report = build_distance_tree(
                    args.alignment,
                    method=args.method,
                    model=args.model,
                    gap_handling=args.gap_handling,
                    ambiguity_policy=args.ambiguity_policy,
                )
                output_path = write_newick(args.out, tree)
                outputs = _finalize_outputs(
                    args,
                    command="alignment",
                    inputs=[args.alignment],
                    outputs=[output_path],
                )
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "pair_count": report.pair_count,
                            "method": report.method,
                            "ambiguity_policy": report.ambiguity_policy,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "compare-distance-trees":
                report = compare_distance_tree_topologies(
                    args.alignment,
                    model=args.model,
                    gap_handling=args.gap_handling,
                    ambiguity_policy=args.ambiguity_policy,
                )
                outputs = _finalize_outputs(
                    args,
                    command="alignment",
                    inputs=[args.alignment],
                )
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "shared_taxa": len(report.shared_taxa),
                            "robinson_foulds_distance": report.robinson_foulds_distance,
                            "same_unrooted_topology": report.same_unrooted_topology,
                            "ambiguity_policy": report.ambiguity_policy,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "bootstrap-tree":
                trees, report = bootstrap_distance_trees(
                    args.alignment,
                    method=args.method,
                    model=args.model,
                    gap_handling=args.gap_handling,
                    ambiguity_policy=args.ambiguity_policy,
                    replicates=args.replicates,
                    seed=args.seed,
                )
                outputs: list[Path | str] = []
                if args.support_out is not None:
                    outputs.append(write_distance_bootstrap_support(args.support_out, report))
                if args.tree_set_out is not None:
                    outputs.append(write_tree_set(args.tree_set_out, trees))
                outputs = _finalize_outputs(
                    args,
                    command="alignment",
                    inputs=[args.alignment],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "replicate_count": report.tree_count,
                            "support_row_count": len(report.support),
                            "method": report.method,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "distance-bundle":
                report = write_distance_reproducibility_bundle(
                    args.out_dir,
                    alignment_path=args.alignment,
                    method=args.method,
                    model=args.model,
                    gap_handling=args.gap_handling,
                    ambiguity_policy=args.ambiguity_policy,
                    replicates=args.replicates,
                    seed=args.seed,
                )
                outputs = _finalize_outputs(
                    args,
                    command="alignment",
                    inputs=[args.alignment],
                    outputs=list(report.files),
                )
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "file_count": len(report.files),
                            "replicates": report.replicates,
                            "method": report.method,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "coding":
                report = inspect_coding_alignment(args.alignment)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "frameshift_like_sequence_count": len(report.frameshift_like_sequences),
                            "stop_codon_count": len(report.stop_codons),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "translate":
                records, report = translate_coding_alignment(args.alignment)
                output_path = write_fasta_alignment(args.out, records)
                outputs = _finalize_outputs(
                    args,
                    command="alignment",
                    inputs=[args.alignment],
                    outputs=[output_path],
                )
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "translated_sequence_count": report.translated_sequence_count,
                            "translated_alignment_length": report.translated_alignment_length,
                            "stop_codon_count": report.stop_codon_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            report = link_alignment_to_tree(args.tree, args.alignment, strict=args.strict)
            outputs = _finalize_outputs(args, command="alignment", inputs=[args.tree, args.alignment])
            _print_result(
                build_command_result(
                    command="alignment",
                    inputs=[args.tree, args.alignment],
                    outputs=outputs,
                    metrics={
                        "tree_taxa": report.tree_taxa,
                        "alignment_ids": report.alignment_ids,
                        "linked_taxa": report.linked_taxa,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "comparative":
            if args.comparative_command == "readiness":
                report = summarize_numeric_trait_readiness(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                )
                outputs = _finalize_outputs(args, command="comparative", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "tree_taxa": report.tree_taxa,
                            "analysis_taxa": len(report.analysis_taxa),
                            "ready": report.ready,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "summarize":
                report = summarize_numeric_trait(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                )
                outputs = _finalize_outputs(args, command="comparative", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "mean": report.mean,
                            "variance": report.variance,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "contrasts":
                report = compute_phylogenetic_independent_contrasts(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                )
                outputs = _finalize_outputs(args, command="comparative", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "contrast_count": len(report.contrasts),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "signal":
                blomberg = compute_blombergs_k(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                )
                lambda_report = estimate_pagels_lambda(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                )
                test_report = compute_phylogenetic_signal_test(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    permutations=args.permutations,
                    seed=args.seed,
                )
                outputs = _finalize_outputs(args, command="comparative", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "taxon_count": blomberg.taxon_count,
                            "blombergs_k": blomberg.k,
                            "pagels_lambda": lambda_report.lambda_value,
                            "signal_p_value": test_report.p_value,
                        },
                        data={
                            "blombergs_k": blomberg,
                            "pagels_lambda": lambda_report,
                            "signal_test": test_report,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "brownian":
                report = fit_brownian_motion_model(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                )
                outputs = _finalize_outputs(args, command="comparative", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.residual_diagnostics.warnings,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "root_state": report.root_state,
                            "rate": report.rate,
                            "log_likelihood": report.log_likelihood,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "ou":
                report = fit_ornstein_uhlenbeck_model(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                )
                outputs = _finalize_outputs(args, command="comparative", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=[
                            *report.residual_diagnostics.warnings,
                            *[warning.message for warning in report.identifiability_warnings],
                        ],
                        metrics={
                            "taxon_count": report.taxon_count,
                            "alpha": report.alpha,
                            "theta": report.theta,
                            "log_likelihood": report.log_likelihood,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "compare-models":
                report = compare_brownian_and_ou_models(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                )
                outputs = _finalize_outputs(args, command="comparative", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "better_model": report.better_model,
                            "model_count": len(report.rows),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "validate-reference":
                report = validate_comparative_reference_examples()
                outputs = _finalize_outputs(args, command="comparative", inputs=[])
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[],
                        outputs=outputs,
                        metrics={
                            "case_count": len(report.observations),
                            "all_passed": report.all_passed,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "sensitivity":
                report = run_comparative_sensitivity_analysis(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    model=args.model,
                    taxon_column=args.taxon_column,
                )
                outputs = _finalize_outputs(args, command="comparative", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "taxon_count": len(report.rows),
                            "model": report.model,
                            "influential_taxa": len(report.most_influential_taxa),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            lambda_value: float | str
            if args.lambda_value == "estimate":
                lambda_value = "estimate"
            else:
                lambda_value = float(args.lambda_value)
            if args.comparative_command == "multiple-testing":
                report = run_pgls_multiple_testing(
                    args.tree,
                    args.table,
                    responses=list(args.responses),
                    predictors=list(args.predictors),
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                )
                outputs = _finalize_outputs(args, command="comparative", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "response_count": len(report.responses),
                            "test_count": len(report.rows),
                            "significant_count": sum(1 for row in report.rows if row.significant),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "report":
                report = build_comparative_method_report(
                    args.tree,
                    args.table,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                )
                if args.out is not None:
                    write_comparative_method_report(args.out, report)
                outputs = _finalize_outputs(
                    args,
                    command="comparative",
                    inputs=[args.tree, args.table],
                    outputs=[args.out] if args.out else [],
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "taxon_count": report.snapshot.pgls_model.taxon_count,
                            "selected_model": report.snapshot.model_comparison.better_model,
                            "audit_row_count": len(report.snapshot.audit_rows),
                            "limitation_count": len(report.snapshot.limitations),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "influence":
                report = build_trait_influence_report(
                    args.tree,
                    args.table,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                )
                outputs = _finalize_outputs(args, command="comparative", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "predictor_count": len(report.predictor_rows),
                            "taxon_count": len(report.taxon_rows),
                            "selected_model": report.selected_model,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "compare-trees":
                report = compare_comparative_results_across_trees(
                    args.left_tree,
                    args.right_tree,
                    args.table,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                )
                outputs = _finalize_outputs(args, command="comparative", inputs=[args.left_tree, args.right_tree, args.table])
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.left_tree, args.right_tree, args.table],
                        outputs=outputs,
                        metrics={
                            "coefficient_delta_count": len(report.coefficient_deltas),
                            "left_selected_model": report.left_selected_model,
                            "right_selected_model": report.right_selected_model,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "compare-pruning":
                report = compare_comparative_results_across_pruning(
                    args.tree,
                    args.table,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    drop_taxa=list(args.drop_taxa or []),
                    keep_taxa=list(args.keep_taxa or []),
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                )
                outputs = _finalize_outputs(args, command="comparative", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "baseline_taxa": len(report.baseline_taxa),
                            "pruned_taxa": len(report.pruned_taxa),
                            "dropped_taxa": len(report.dropped_taxa),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            input_report = inspect_pgls_inputs(
                args.tree,
                args.table,
                response=args.response,
                predictors=list(args.predictors or []),
                formula=args.formula,
                taxon_column=args.taxon_column,
            )
            report = run_pgls(
                args.tree,
                args.table,
                response=args.response,
                predictors=list(args.predictors or []),
                formula=args.formula,
                taxon_column=args.taxon_column,
                lambda_value=lambda_value,
            )
            outputs = _finalize_outputs(args, command="comparative", inputs=[args.tree, args.table])
            _print_result(
                build_command_result(
                    command="comparative",
                    inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=input_report.warnings,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "predictor_count": len(report.predictors),
                            "encoded_predictor_count": len(report.encoded_columns) - 1,
                            "lambda_value": report.lambda_value,
                            "r_squared": report.r_squared,
                        },
                    data={
                        "inputs": input_report,
                        "model": report,
                    },
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "ancestral":
            if args.ancestral_command == "continuous":
                report = reconstruct_continuous_ancestral_states(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    alpha=args.alpha,
                )
                outputs: list[Path | str] = []
                if args.table_out is not None:
                    outputs.append(write_ancestral_state_table(args.table_out, report))
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "estimate_count": len(report.estimates),
                            "model": report.model,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "discrete":
                if args.state_ordering == "ordered" and args.model == "fitch":
                    parser.error("ordered ancestral discrete reconstruction requires a likelihood model")
                report = reconstruct_discrete_ancestral_states(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                )
                outputs: list[Path | str] = []
                if args.table_out is not None:
                    outputs.append(write_ancestral_state_table(args.table_out, report))
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "estimate_count": len(report.estimates),
                            "state_count": len(report.observed_states),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "compare":
                report = compare_continuous_ancestral_models(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    left_model=args.left_model,
                    right_model=args.right_model,
                    left_alpha=args.left_alpha,
                    right_alpha=args.right_alpha,
                )
                outputs = _finalize_outputs(args, command="ancestral", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "compared_node_count": len(report.rows),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "sensitivity":
                _validate_ancestral_discrete_model_arguments(args, parser)
                resolved_model = args.model or ("brownian" if args.kind == "continuous" else "fitch")
                report = build_ancestral_sensitivity_report(
                    tree_path=args.tree,
                    traits_path=args.table,
                    trait=args.trait,
                    reconstruction_kind=args.kind,
                    model=resolved_model,
                    taxon_column=args.taxon_column,
                    alpha=args.alpha,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                    compare_tree_path=args.compare_tree,
                    compare_model=args.compare_model,
                    drop_taxa=args.drop_taxa,
                    coding_map=_parse_assignment_map(args.coding_map) or None,
                )
                outputs = _finalize_outputs(args, command="ancestral", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "baseline_node_count": report.baseline_node_count,
                            "has_model_sensitivity": report.model_sensitivity is not None,
                            "has_tree_sensitivity": report.tree_sensitivity is not None,
                            "has_pruning_sensitivity": report.pruning_sensitivity is not None,
                            "has_trait_coding_sensitivity": report.trait_coding_sensitivity is not None,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "render":
                _validate_ancestral_discrete_model_arguments(args, parser)
                if args.kind == "continuous":
                    resolved_model = args.model or "brownian"
                    reconstruction = reconstruct_continuous_ancestral_states(
                        args.tree,
                        args.table,
                        trait=args.trait,
                        taxon_column=args.taxon_column,
                        model=resolved_model,
                        alpha=args.alpha,
                    )
                else:
                    resolved_model = args.model or "fitch"
                    reconstruction = reconstruct_discrete_ancestral_states(
                        args.tree,
                        args.table,
                        trait=args.trait,
                        taxon_column=args.taxon_column,
                        model=resolved_model,
                        state_ordering=args.state_ordering,
                        ordered_states=_split_csv_values(args.ordered_states) or None,
                    )
                result = render_ancestral_state_tree(
                    args.tree,
                    reconstruction,
                    out_path=args.out,
                    layout=args.layout,
                )
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=[result.output_path],
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=getattr(reconstruction, "warnings", []),
                        metrics={
                            "tip_count": result.tip_count,
                            "rendered_internal_annotation_count": result.rendered_internal_annotation_count,
                            "layout": result.layout,
                        },
                        data={
                            "reconstruction": reconstruction,
                            "render": result,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            resolved_model = args.model or ("brownian" if args.kind == "continuous" else "fitch")
            if args.ancestral_command == "package":
                _validate_ancestral_discrete_model_arguments(args, parser)
                result = build_ancestral_figure_package(
                    tree_path=args.tree,
                    traits_path=args.table,
                    trait=args.trait,
                    reconstruction_kind=args.kind,
                    out_dir=args.out_dir,
                    taxon_column=args.taxon_column,
                    model=resolved_model,
                    alpha=args.alpha,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                    layout=args.layout,
                )
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=[
                        result.figure_path,
                        result.node_table_path,
                        result.uncertainty_table_path,
                        result.legend_path,
                        result.model_description_path,
                        result.caption_path,
                        result.manifest_path,
                    ],
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "output_dir": str(result.output_dir),
                            "artifact_count": 7,
                        },
                        data=result,
                    ),
                    json_output=args.json,
                )
                return 0
            _validate_ancestral_discrete_model_arguments(args, parser)
            result = render_ancestral_state_report(
                tree_path=args.tree,
                traits_path=args.table,
                trait=args.trait,
                reconstruction_kind=args.kind,
                out_path=args.out,
                taxon_column=args.taxon_column,
                model=resolved_model,
                alpha=args.alpha,
                state_ordering=args.state_ordering,
                ordered_states=_split_csv_values(args.ordered_states) or None,
                compare_model=args.compare_model,
                compare_tree_path=args.compare_tree,
                drop_taxa=args.drop_taxa,
                coding_map=_parse_assignment_map(args.coding_map) or None,
            )
            outputs = _finalize_outputs(
                args,
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=[result.output_path, args.out.with_suffix(".svg")],
            )
            _print_result(
                build_command_result(
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                    metrics={
                        "report_kind": result.report_kind,
                        "reconstruction_kind": result.reconstruction_kind,
                    },
                    data=result,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "discrete-evolution":
            allowed_states = _split_csv_values(args.allowed_states) if hasattr(args, "allowed_states") else []
            if args.discrete_evolution_command == "validate-coding":
                report = validate_discrete_state_coding(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    allowed_states=allowed_states or None,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                )
                outputs = _finalize_outputs(args, command="discrete-evolution", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "valid": report.valid,
                            "issue_count": len(report.issues),
                            "observed_state_count": len(report.observed_states),
                            "state_ordering": report.state_ordering,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.discrete_evolution_command == "imbalance":
                report = detect_state_imbalance_problems(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                )
                outputs = _finalize_outputs(args, command="discrete-evolution", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=[warning.message for warning in report.warnings],
                        metrics={
                            "taxon_count": report.taxon_count,
                            "observed_state_count": len(report.observed_states),
                            "warning_count": len(report.warnings),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.discrete_evolution_command == "reference":
                report = validate_discrete_transition_reference_examples()
                outputs = _finalize_outputs(args, command="discrete-evolution", inputs=[])
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[],
                        outputs=outputs,
                        metrics={
                            "case_count": report.case_count,
                            "all_passed": report.all_passed,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.discrete_evolution_command == "model":
                report = estimate_ancestral_geographic_states(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    allowed_states=allowed_states or None,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                )
                outputs: list[Path | str] = []
                if args.node_table_out is not None:
                    outputs.append(write_node_state_probability_table(args.node_table_out, report))
                if args.transitions_out is not None:
                    outputs.append(write_transition_summary_table(args.transitions_out, report))
                outputs = _finalize_outputs(
                    args,
                    command="discrete-evolution",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "observed_state_count": len(report.observed_states),
                            "transition_count": report.transition_summary.transition_count,
                            "strongly_supported_transition_count": report.transition_summary.strongly_supported_transition_count,
                            "model": report.model,
                            "state_ordering": report.state_ordering,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.discrete_evolution_command == "render":
                report = estimate_ancestral_geographic_states(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    allowed_states=allowed_states or None,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                )
                result = render_tree_with_geographic_states(
                    args.tree,
                    report,
                    out_path=args.out,
                    layout=args.layout,
                )
                outputs = _finalize_outputs(
                    args,
                    command="discrete-evolution",
                    inputs=[args.tree, args.table],
                    outputs=[result.output_path],
                )
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "tip_count": result.tip_count,
                            "rendered_internal_annotation_count": result.rendered_internal_annotation_count,
                            "layout": result.layout,
                            "model": report.model,
                            "state_ordering": report.state_ordering,
                        },
                        data={
                            "reconstruction": report,
                            "render": result,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.discrete_evolution_command == "report":
                result = render_discrete_state_evolution_report(
                    tree_path=args.tree,
                    traits_path=args.table,
                    trait=args.trait,
                    out_path=args.out,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    allowed_states=allowed_states or None,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                    compare_model=args.compare_model,
                )
                outputs = _finalize_outputs(
                    args,
                    command="discrete-evolution",
                    inputs=[args.tree, args.table],
                    outputs=[result.output_path, args.out.with_suffix(".svg")],
                )
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "report_kind": result.report_kind,
                            "model": result.model,
                            "state_ordering": args.state_ordering,
                        },
                        data=result,
                    ),
                    json_output=args.json,
                )
                return 0
            comparison = compare_discrete_state_models(
                args.tree,
                args.table,
                trait=args.trait,
                taxon_column=args.taxon_column,
                left_model=args.left_model,
                right_model=args.right_model,
                allowed_states=allowed_states or None,
                state_ordering=args.state_ordering,
                ordered_states=_split_csv_values(args.ordered_states) or None,
            )
            outputs: list[Path | str] = []
            if args.table_out is not None:
                outputs.append(write_discrete_model_comparison_table(args.table_out, comparison))
            outputs = _finalize_outputs(
                args,
                command="discrete-evolution",
                inputs=[args.tree, args.table],
                outputs=outputs,
            )
            _print_result(
                build_command_result(
                    command="discrete-evolution",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                    metrics={
                        "better_model": comparison.better_model,
                        "model_count": len(comparison.rows),
                        "differing_node_count": sum(1 for row in comparison.node_differences if row.differs),
                        "state_ordering": args.state_ordering,
                    },
                    data=comparison,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "diversification":
            if args.diversification_command == "ltt":
                report = compute_lineage_through_time_curve(args.tree)
                outputs: list[Path | str] = []
                if args.out is not None:
                    outputs.append(write_lineage_through_time_table(args.out, report))
                outputs = _finalize_outputs(
                    args,
                    command="diversification",
                    inputs=[args.tree],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={
                            "tip_count": report.tip_count,
                            "root_age": report.root_age,
                            "point_count": len(report.points),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.diversification_command == "sampling":
                report = detect_incomplete_taxon_sampling_metadata(
                    args.tree,
                    args.table,
                    taxon_column=args.taxon_column,
                    sampling_column=args.sampling_column,
                )
                outputs = _finalize_outputs(args, command="diversification", inputs=[args.tree, args.table])
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "complete": report.complete,
                            "matched_taxon_count": len(report.matched_taxa),
                            "missing_taxon_count": len(report.missing_taxa),
                            "invalid_row_count": len(report.invalid_rows),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.diversification_command == "estimate":
                inputs = [args.tree]
                if args.metadata is not None:
                    inputs.append(args.metadata)
                report = estimate_diversification_rate(
                    args.tree,
                    metadata_path=args.metadata,
                    taxon_column=args.taxon_column,
                    sampling_column=args.sampling_column,
                    model=args.model,
                )
                outputs = _finalize_outputs(args, command="diversification", inputs=inputs)
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=inputs,
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "model": report.model,
                            "sampling_fraction": report.sampling_fraction,
                            "net_diversification_rate": report.net_diversification_rate,
                            "aic": report.aic,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.diversification_command == "compare-models":
                inputs = [args.tree]
                if args.metadata is not None:
                    inputs.append(args.metadata)
                report = compare_diversification_models(
                    args.tree,
                    metadata_path=args.metadata,
                    taxon_column=args.taxon_column,
                    sampling_column=args.sampling_column,
                )
                outputs = _finalize_outputs(args, command="diversification", inputs=inputs)
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=inputs,
                        outputs=outputs,
                        metrics={
                            "better_model": report.better_model,
                            "model_count": len(report.rows),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.diversification_command == "clades":
                report = detect_diversification_outlier_clades(
                    args.tree,
                    min_tip_count=args.min_tip_count,
                    model=args.model,
                )
                outputs: list[Path | str] = []
                if args.out is not None:
                    outputs.append(write_clade_diversification_table(args.out, report))
                outputs = _finalize_outputs(
                    args,
                    command="diversification",
                    inputs=[args.tree],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=[args.tree],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "global_rate": report.global_rate,
                            "high_clade_count": len(report.high_diversification_clades),
                            "low_clade_count": len(report.low_diversification_clades),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.diversification_command == "trait-dependent":
                report = run_trait_dependent_diversification_analysis(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                )
                outputs: list[Path | str] = []
                if args.out is not None:
                    outputs.append(write_trait_dependent_diversification_table(args.out, report))
                outputs = _finalize_outputs(
                    args,
                    command="diversification",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "state_count": len(report.states),
                            "monophyletic_state_count": sum(1 for row in report.states if row.monophyletic),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            inputs = [args.tree]
            if args.metadata is not None:
                inputs.append(args.metadata)
            if args.traits is not None:
                inputs.append(args.traits)
            result = render_diversification_report(
                tree_path=args.tree,
                out_path=args.out,
                metadata_path=args.metadata,
                taxon_column=args.taxon_column,
                sampling_column=args.sampling_column,
                traits_path=args.traits,
                trait=args.trait,
            )
            outputs = _finalize_outputs(
                args,
                command="diversification",
                inputs=inputs,
                outputs=[result.output_path],
            )
            _print_result(
                build_command_result(
                    command="diversification",
                    inputs=inputs,
                    outputs=outputs,
                    metrics={
                        "report_kind": result.report_kind,
                    },
                    data=result,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "distance":
            if args.distance_command == "reference":
                report = validate_distance_reference_examples()
                outputs = _finalize_outputs(args, command="distance", inputs=[])
                _print_result(
                    build_command_result(
                        command="distance",
                        inputs=[],
                        outputs=outputs,
                        metrics={
                            "case_count": len(report.observations),
                            "all_passed": report.all_passed,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.distance_command == "assumptions":
                report = assess_imported_distance_method_assumptions(args.matrix)
                outputs = _finalize_outputs(args, command="distance", inputs=[args.matrix])
                _print_result(
                    build_command_result(
                        command="distance",
                        inputs=[args.matrix],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "ultrametric_compatible": report.ultrametric_compatible,
                            "upgma_violation_count": len(report.upgma_ultrametric_violations),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.distance_command == "validate":
                report = validate_imported_distance_matrix(args.matrix)
                outputs = _finalize_outputs(args, command="distance", inputs=[args.matrix])
                _print_result(
                    build_command_result(
                        command="distance",
                        inputs=[args.matrix],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "taxon_count": len(report.identifiers),
                            "pair_count": report.pair_count,
                            "complete": report.complete,
                            "symmetric": report.symmetric,
                            "nonmetric_observation_count": len(report.nonmetric_observations),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.distance_command == "build-tree":
                tree, report = build_tree_from_imported_distance_matrix(args.matrix, method=args.method)
                output_path = write_newick(args.out, tree)
                outputs = _finalize_outputs(
                    args,
                    command="distance",
                    inputs=[args.matrix],
                    outputs=[output_path],
                )
                _print_result(
                    build_command_result(
                        command="distance",
                        inputs=[args.matrix],
                        outputs=outputs,
                        metrics={
                            "method": report.method,
                            "taxon_count": report.taxon_count,
                            "pair_count": report.pair_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.distance_command == "report":
                report = render_distance_report(out_path=args.out, matrix_path=args.matrix)
                outputs = _finalize_outputs(
                    args,
                    command="distance",
                    inputs=[args.matrix],
                    outputs=[args.out],
                )
                _print_result(
                    build_command_result(
                        command="distance",
                        inputs=[args.matrix],
                        outputs=outputs,
                        metrics={"section_count": len(report.machine_manifest["sections"])},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            limitations = distance_method_limitations()
            outputs = _finalize_outputs(args, command="distance", inputs=[args.matrix])
            _print_result(
                build_command_result(
                    command="distance",
                    inputs=[args.matrix],
                    outputs=outputs,
                    metrics={"limitation_count": len(limitations)},
                    data={"limitations": limitations},
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "tree-set":
            if args.tree_set_command == "inspect":
                report = load_tree_set(args.tree_set)
                outputs = _finalize_outputs(args, command="tree-set", inputs=[args.tree_set])
                _print_result(
                    build_command_result(
                        command="tree-set",
                        inputs=[args.tree_set],
                        outputs=outputs,
                        metrics={
                            "tree_count": report.tree_count,
                            "shared_taxon_count": len(report.shared_taxa),
                            "rooted_topology_count": report.rooted_topology_count,
                            "unrooted_topology_count": report.unrooted_topology_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.tree_set_command == "consensus":
                tree, report = compute_consensus_tree(args.tree_set)
                output_path = write_consensus_tree(args.out, tree)
                outputs = _finalize_outputs(
                    args,
                    command="tree-set",
                    inputs=[args.tree_set],
                    outputs=[output_path],
                )
                _print_result(
                    build_command_result(
                        command="tree-set",
                        inputs=[args.tree_set],
                        outputs=outputs,
                        metrics={"tree_count": report.tree_count, "shared_taxon_count": len(report.shared_taxa)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.tree_set_command == "clade-frequencies":
                report = compute_clade_frequency_table(args.tree_set)
                outputs = []
                if args.out is not None:
                    outputs.append(write_clade_frequency_table(args.out, report))
                outputs = _finalize_outputs(
                    args,
                    command="tree-set",
                    inputs=[args.tree_set],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="tree-set",
                        inputs=[args.tree_set],
                        outputs=outputs,
                        metrics={"tree_count": report.tree_count, "clade_count": len(report.clade_frequencies)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.tree_set_command == "distance-matrix":
                report = compute_tree_distance_matrix(args.tree_set)
                outputs = []
                if args.out is not None:
                    outputs.append(write_tree_distance_matrix(args.out, report))
                outputs = _finalize_outputs(
                    args,
                    command="tree-set",
                    inputs=[args.tree_set],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="tree-set",
                        inputs=[args.tree_set],
                        outputs=outputs,
                        metrics={"tree_count": report.tree_count, "pair_count": len(report.pairs)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.tree_set_command == "cluster":
                report = cluster_trees_by_topology(args.tree_set)
                outputs = _finalize_outputs(args, command="tree-set", inputs=[args.tree_set])
                _print_result(
                    build_command_result(
                        command="tree-set",
                        inputs=[args.tree_set],
                        outputs=outputs,
                        metrics={"tree_count": report.tree_count, "cluster_count": len(report.clusters)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.tree_set_command == "unstable-taxa":
                report = detect_unstable_taxa(args.tree_set)
                outputs = _finalize_outputs(args, command="tree-set", inputs=[args.tree_set])
                _print_result(
                    build_command_result(
                        command="tree-set",
                        inputs=[args.tree_set],
                        outputs=outputs,
                        metrics={"tree_count": report.tree_count, "unstable_taxon_count": len(report.taxa)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.tree_set_command == "unstable-clades":
                report = detect_unstable_clades(args.tree_set)
                outputs = _finalize_outputs(args, command="tree-set", inputs=[args.tree_set])
                _print_result(
                    build_command_result(
                        command="tree-set",
                        inputs=[args.tree_set],
                        outputs=outputs,
                        metrics={"tree_count": report.tree_count, "unstable_clade_count": len(report.clades)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.tree_set_command == "compare":
                report = compare_posterior_tree_sets(args.left, args.right)
                outputs = _finalize_outputs(args, command="tree-set", inputs=[args.left, args.right])
                _print_result(
                    build_command_result(
                        command="tree-set",
                        inputs=[args.left, args.right],
                        outputs=outputs,
                        metrics={
                            "left_tree_count": report.left_tree_count,
                            "right_tree_count": report.right_tree_count,
                            "shared_rooted_topology_count": report.shared_rooted_topology_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            report = render_tree_uncertainty_report(tree_set_path=args.tree_set, out_path=args.out)
            outputs = _finalize_outputs(
                args,
                command="tree-set",
                inputs=[args.tree_set],
                outputs=[args.out],
            )
            _print_result(
                build_command_result(
                    command="tree-set",
                    inputs=[args.tree_set],
                    outputs=outputs,
                    metrics={"tree_count": report.tree_count, "section_count": len(report.machine_manifest["sections"])},
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "simulate":
            if args.simulate_command == "tree-birth-death":
                trees, report = simulate_birth_death_trees(
                    tree_count=args.tree_count,
                    tip_count=args.tip_count,
                    birth_rate=args.birth_rate,
                    death_rate=args.death_rate,
                    seed=args.seed,
                )
                output_path = write_tree_set(args.out, trees)
                outputs = _finalize_outputs(args, command="simulate", inputs=[], outputs=[output_path])
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[],
                        outputs=outputs,
                        metrics={"tree_count": report.tree_count, "tip_count": report.tip_count},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "tree-coalescent":
                trees, report = simulate_coalescent_trees(
                    tree_count=args.tree_count,
                    tip_count=args.tip_count,
                    population_size=args.population_size,
                    seed=args.seed,
                )
                output_path = write_tree_set(args.out, trees)
                outputs = _finalize_outputs(args, command="simulate", inputs=[], outputs=[output_path])
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[],
                        outputs=outputs,
                        metrics={"tree_count": report.tree_count, "tip_count": report.tip_count},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "traits-brownian":
                report = simulate_brownian_traits(
                    args.tree,
                    root_state=args.root_state,
                    sigma=args.sigma,
                    seed=args.seed,
                )
                output_path = write_continuous_trait_table(args.out, report)
                outputs = _finalize_outputs(args, command="simulate", inputs=[args.tree], outputs=[output_path])
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={"tip_count": report.tip_count, "trait_count": len(report.traits)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "traits-ou":
                report = simulate_ou_traits(
                    args.tree,
                    root_state=args.root_state,
                    sigma=args.sigma,
                    alpha=args.alpha,
                    theta=args.theta,
                    seed=args.seed,
                )
                output_path = write_continuous_trait_table(args.out, report)
                outputs = _finalize_outputs(args, command="simulate", inputs=[args.tree], outputs=[output_path])
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={"tip_count": report.tip_count, "trait_count": len(report.traits)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "traits-discrete":
                report = simulate_discrete_traits(
                    args.tree,
                    states=args.states,
                    transition_rate=args.transition_rate,
                    root_state=args.root_state,
                    seed=args.seed,
                )
                output_path = write_discrete_trait_table(args.out, report)
                outputs = _finalize_outputs(args, command="simulate", inputs=[args.tree], outputs=[output_path])
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={"tip_count": report.tip_count, "trait_count": len(report.traits)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "alignment-dna":
                report = simulate_dna_alignment(
                    args.tree,
                    sequence_length=args.sequence_length,
                    substitution_rate=args.substitution_rate,
                    seed=args.seed,
                )
                output_path = write_simulated_alignment(args.out, report)
                outputs = _finalize_outputs(args, command="simulate", inputs=[args.tree], outputs=[output_path])
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={"tip_count": report.tip_count, "sequence_length": report.sequence_length},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            report = simulate_protein_alignment(
                args.tree,
                sequence_length=args.sequence_length,
                substitution_rate=args.substitution_rate,
                seed=args.seed,
            )
            output_path = write_simulated_alignment(args.out, report)
            outputs = _finalize_outputs(args, command="simulate", inputs=[args.tree], outputs=[output_path])
            _print_result(
                build_command_result(
                    command="simulate",
                    inputs=[args.tree],
                    outputs=outputs,
                    metrics={"tip_count": report.tip_count, "sequence_length": report.sequence_length},
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "benchmark":
            if args.benchmark_command == "tree-validation":
                report = benchmark_tree_validation(replicates=args.replicates)
            elif args.benchmark_command == "tree-comparison":
                report = benchmark_tree_comparison(replicates=args.replicates)
            else:
                report = benchmark_alignment_diagnostics(
                    replicates=args.replicates,
                    sequence_length=args.sequence_length,
                )
            outputs = _finalize_outputs(args, command="benchmark", inputs=[])
            _print_result(
                build_command_result(
                    command="benchmark",
                    inputs=[],
                    outputs=outputs,
                    metrics={"observation_count": len(report.observations), "replicates": report.replicates},
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "inspect":
            report = inspect_tree_path(args.tree, source_format=args.format)
            outputs = _finalize_outputs(args, command="inspect", inputs=[args.tree])
            _print_result(
                build_command_result(
                    command="inspect",
                    inputs=[args.tree],
                    outputs=outputs,
                    warnings=report.warnings,
                    metrics={
                        "tip_count": report.tip_count,
                        "node_count": report.node_count,
                        "internal_node_count": report.internal_node_count,
                        "edge_count": report.edge_count,
                        "clade_count": report.clade_count,
                        "is_binary": report.is_binary,
                        "polytomy_count": report.polytomy_count,
                        "branch_length_status": report.branch_length_status,
                        "is_ultrametric": report.is_ultrametric,
                        "tree_diameter": report.tree_diameter,
                        "colless_imbalance_index": report.colless_imbalance_index,
                        "sackin_imbalance_index": report.sackin_imbalance_index,
                        "tree_quality_score": report.tree_quality_score,
                        "zero_length_branch_count": report.zero_length_branch_count,
                        "cherry_count": report.cherry_count,
                        "missing_internal_branch_count": len(report.missing_internal_branch_nodes),
                        "missing_terminal_branch_count": len(report.missing_terminal_branch_taxa),
                        "singleton_internal_node_count": len(report.singleton_internal_nodes),
                        "long_branch_outlier_count": len(report.long_branch_outliers),
                        "short_branch_outlier_count": len(report.short_branch_outliers),
                        "likely_support_label_count": len(report.likely_support_labels),
                        "likely_named_internal_label_count": len(report.likely_named_internal_labels),
                        "suspicious_support_range_count": len(report.suspicious_support_value_ranges),
                        "root_classification": report.root_state_confidence.classification,
                        "internal_label_conflict_count": len(report.internal_label_conflicts),
                        "unsafe_external_label_count": len(report.unsafe_external_labels),
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "normalize":
            tree = load_tree(args.tree, source_format=args.format)
            output_path = write_newick(args.out, tree)
            outputs = _finalize_outputs(args, command="normalize", inputs=[args.tree], outputs=[output_path])
            if args.json:
                _print_result(
                    build_command_result(
                        command="normalize",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={"tip_count": tree.tip_count},
                        data={"source_format": tree.source_format, "output_format": "newick"},
                    ),
                    json_output=True,
                )
            else:
                print(output_path)
            return 0
        if args.command == "normalize-taxa":
            tree = load_tree(args.tree, source_format=args.format)
            normalized_tree, report = normalize_tree_taxa(tree, policy=args.policy)
            output_path = write_newick(args.out, normalized_tree)
            mapping_path = args.mapping_out or args.out.with_suffix(f"{args.out.suffix}.mapping.tsv")
            write_taxon_mapping(mapping_path, report.renamed_taxa)
            outputs = _finalize_outputs(
                args,
                command="normalize-taxa",
                inputs=[args.tree],
                outputs=[output_path, mapping_path],
            )
            if args.json:
                _print_result(
                    build_command_result(
                        command="normalize-taxa",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={"renamed_taxa": len(report.renamed_taxa)},
                        data=report,
                    ),
                    json_output=True,
                )
            else:
                print(output_path)
            return 0
        if args.command == "taxonomy":
            if args.taxonomy_command == "synonyms":
                tree = load_tree(args.tree, source_format=args.format)
                report = audit_tree_taxon_synonyms(tree, args.synonym_table)
                outputs = _finalize_outputs(args, command="taxonomy", inputs=[args.tree, args.synonym_table])
                _print_result(
                    build_command_result(
                        command="taxonomy",
                        inputs=[args.tree, args.synonym_table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "candidate_count": len(report.candidates),
                            "ambiguous_mapping_count": len(report.ambiguous_mappings),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.taxonomy_command == "resolve-synonyms":
                tree = load_tree(args.tree, source_format=args.format)
                resolved_tree, report = resolve_tree_taxon_synonyms(
                    tree,
                    synonym_table_path=args.synonym_table,
                    resolution_policy=args.resolution_policy,
                )
                output_path = write_newick(args.out, resolved_tree)
                mapping_path = args.mapping_out or args.out.with_suffix(f"{args.out.suffix}.synonyms.tsv")
                write_synonym_resolution_mapping(mapping_path, report)
                outputs = _finalize_outputs(
                    args,
                    command="taxonomy",
                    inputs=[args.tree, args.synonym_table],
                    outputs=[output_path, mapping_path],
                )
                _print_result(
                    build_command_result(
                        command="taxonomy",
                        inputs=[args.tree, args.synonym_table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "renamed_taxa": len(report.renamed_taxa),
                            "ambiguous_mapping_count": len(report.ambiguous_mappings),
                            "duplicate_resolved_label_count": len(report.duplicate_resolved_labels),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.taxonomy_command == "namespaces":
                tree = load_tree(args.tree, source_format=args.format)
                report = inspect_tree_taxon_namespaces(tree)
                outputs = _finalize_outputs(args, command="taxonomy", inputs=[args.tree])
                _print_result(
                    build_command_result(
                        command="taxonomy",
                        inputs=[args.tree],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "namespace_count": len(report.namespace_counts),
                            "mixed_namespaces": report.mixed_namespaces,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.taxonomy_command == "loss":
                report = build_taxon_workflow_loss_report(
                    args.tree,
                    args.metadata,
                    args.traits,
                    alignment_path=args.alignment,
                    filtered_alignment_path=args.filtered_alignment,
                    inference_tree_path=args.inference_tree,
                    reported_taxa_path=args.reported_taxa,
                )
                inputs = [args.tree, args.metadata, args.traits]
                for optional_path in (args.alignment, args.filtered_alignment, args.inference_tree, args.reported_taxa):
                    if optional_path is not None:
                        inputs.append(optional_path)
                outputs = _finalize_outputs(args, command="taxonomy", inputs=inputs)
                _print_result(
                    build_command_result(
                        command="taxonomy",
                        inputs=inputs,
                        outputs=outputs,
                        metrics={
                            "taxon_count": len(report.rows),
                            "loss_stage_count": len(report.loss_stage_counts),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            run_sources = [
                load_taxon_run_source(label=label, path=path)
                for label, path in (_parse_labelled_run(raw) for raw in args.run)
            ]
            report = build_taxon_stability_report(run_sources)
            inputs = [source.path for source in run_sources]
            outputs = _finalize_outputs(args, command="taxonomy", inputs=inputs)
            _print_result(
                build_command_result(
                    command="taxonomy",
                    inputs=inputs,
                    outputs=outputs,
                    metrics={
                        "source_count": len(report.sources),
                        "stable_taxa": len(report.stable_taxa),
                        "unstable_taxa": len(report.unstable_taxa),
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "topology":
            if args.topology_command == "root-outgroup":
                tree, report = root_tree_on_outgroup(args.tree, outgroup_taxa=list(args.taxa))
            elif args.topology_command == "reroot-midpoint":
                tree, report = reroot_tree_by_midpoint(args.tree)
            else:
                tree, report = unroot_tree(args.tree)
            output_path = write_newick(args.out, tree)
            outputs = _finalize_outputs(
                args,
                command="topology",
                inputs=[args.tree],
                outputs=[output_path],
            )
            _print_result(
                build_command_result(
                    command="topology",
                    inputs=[args.tree],
                    outputs=outputs,
                    metrics={
                        "tip_count": tree.tip_count,
                        "matched_taxa": len(report.matched_taxa),
                        "absent_taxa": len(report.absent_taxa),
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "diagnose":
            if args.target == "distances":
                if args.tree is None:
                    parser.exit(status=2, message="diagnose distances requires a tree path\n")
                report = compute_root_to_tip_distances(args.tree)
                outputs: list[Path | str] = []
                if args.out is not None:
                    output_path = write_root_to_tip_tsv(args.out, report)
                    outputs.append(output_path)
                outputs = _finalize_outputs(args, command="diagnose", inputs=[args.tree], outputs=outputs)
                _print_result(
                    build_command_result(
                        command="diagnose",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={"tip_count": len(report.distances)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.target == "ultrametric":
                if args.tree is None:
                    parser.exit(status=2, message="diagnose ultrametric requires a tree path\n")
                report = diagnose_ultrametricity(args.tree, tolerance=args.tolerance)
                outputs = _finalize_outputs(args, command="diagnose", inputs=[args.tree])
                _print_result(
                    build_command_result(
                        command="diagnose",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={
                            "tip_count": report.tip_count,
                            "tolerance": report.tolerance,
                            "max_deviation": report.max_deviation,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.target == "assumptions":
                if args.tree is None:
                    parser.exit(status=2, message="diagnose assumptions requires a tree path\n")
                report = assess_tree_assumptions(
                    args.tree,
                    metadata_path=args.metadata,
                    taxon_column=args.taxon_column,
                )
                diagnose_inputs: list[Path | str] = [args.tree]
                if args.metadata is not None:
                    diagnose_inputs.append(args.metadata)
                outputs = _finalize_outputs(args, command="diagnose", inputs=diagnose_inputs)
                _print_result(
                    build_command_result(
                        command="diagnose",
                        inputs=diagnose_inputs,
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "standardized_support_count": len(report.standardized_support_labels),
                            "time_tree_compatible": report.time_tree_compatible,
                            "substitution_tree_compatible": report.substitution_tree_compatible,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            tree_path = Path(args.target) if args.tree is None else args.tree
            report = diagnose_tree_path(tree_path)
            outputs = _finalize_outputs(args, command="diagnose", inputs=[tree_path])
            _print_result(
                build_command_result(
                    command="diagnose",
                    inputs=[tree_path],
                    outputs=outputs,
                    warnings=report.forensic.warnings,
                    metrics={
                        "tip_count": report.inspection.tip_count,
                        "validity_decision": report.validation.validity_decision,
                        "polytomy_count": report.validation.polytomy_count,
                        "cherry_count": report.inspection.cherry_count,
                        "tree_diameter": report.inspection.tree_diameter,
                        "tree_quality_score": report.inspection.tree_quality_score,
                        "safe_for_topology_comparison": report.forensic.safe_for_topology_comparison,
                        "safe_for_time_tree_analysis": report.forensic.safe_for_time_tree_analysis,
                        "safe_for_comparative_methods": report.forensic.safe_for_comparative_methods,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "compare":
            if args.left == "report":
                if args.third is None:
                    parser.exit(status=2, message="compare report requires two tree paths\n")
                if args.out is None:
                    parser.exit(status=2, message="compare report requires --out\n")
                left_path = Path(args.right)
                right_path = Path(args.third)
                report = build_tree_comparison_report(left_path, right_path, out_path=args.out)
                outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path], outputs=[args.out])
                _print_result(
                    build_command_result(
                        command="compare",
                        inputs=[left_path, right_path],
                        outputs=outputs,
                        metrics={"shared_taxa": len(report.topology.shared_taxa)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.left == "support":
                if args.third is None:
                    parser.exit(status=2, message="compare support requires two tree paths\n")
                left_path = Path(args.right)
                right_path = Path(args.third)
                report = compare_support_values(left_path, right_path)
                outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path])
                _print_result(
                    build_command_result(
                        command="compare",
                        inputs=[left_path, right_path],
                        outputs=outputs,
                        metrics={"shared_clades": len(report.shared_clades)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.left == "clades":
                if args.third is None:
                    parser.exit(status=2, message="compare clades requires two tree paths\n")
                left_path = Path(args.right)
                right_path = Path(args.third)
                report = compare_clade_sets(left_path, right_path)
                outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path])
                _print_result(
                    build_command_result(
                        command="compare",
                        inputs=[left_path, right_path],
                        outputs=outputs,
                        metrics={
                            "shared_clades": len(report.shared_clades),
                            "left_only_clades": len(report.left_only_clades),
                            "right_only_clades": len(report.right_only_clades),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.left == "prune":
                if args.third is None:
                    parser.exit(status=2, message="compare prune requires two tree paths\n")
                if args.out is None:
                    parser.exit(status=2, message="compare prune requires --out as an output directory\n")
                left_path = Path(args.right)
                right_path = Path(args.third)
                pruned_left, pruned_right, report = prune_trees_to_shared_taxa(left_path, right_path)
                args.out.mkdir(parents=True, exist_ok=True)
                left_out = write_newick(args.out / "left-shared.nwk", pruned_left)
                right_out = write_newick(args.out / "right-shared.nwk", pruned_right)
                outputs = _finalize_outputs(
                    args,
                    command="compare",
                    inputs=[left_path, right_path],
                    outputs=[left_out, right_out],
                )
                _print_result(
                    build_command_result(
                        command="compare",
                        inputs=[left_path, right_path],
                        outputs=outputs,
                        metrics={"shared_taxa": len(report.shared_taxa)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.left == "changes":
                if args.third is None:
                    parser.exit(status=2, message="compare changes requires two tree paths\n")
                left_path = Path(args.right)
                right_path = Path(args.third)
                report = detect_clade_changes(left_path, right_path)
                outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path])
                _print_result(
                    build_command_result(
                        command="compare",
                        inputs=[left_path, right_path],
                        outputs=outputs,
                        metrics={
                            "lost_clades": len(report.lost_clades),
                            "gained_clades": len(report.gained_clades),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.left == "branch-lengths":
                if args.third is None:
                    parser.exit(status=2, message="compare branch-lengths requires two tree paths\n")
                left_path = Path(args.right)
                right_path = Path(args.third)
                report = compare_branch_lengths(left_path, right_path)
                outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path])
                _print_result(
                    build_command_result(
                        command="compare",
                        inputs=[left_path, right_path],
                        outputs=outputs,
                        metrics={"shared_splits": len(report.shared_splits)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.left == "table":
                if args.third is None:
                    parser.exit(status=2, message="compare table requires two tree paths\n")
                if args.out is None:
                    parser.exit(status=2, message="compare table requires --out\n")
                left_path = Path(args.right)
                right_path = Path(args.third)
                output_path = write_tree_comparison_table(args.out, left_path, right_path)
                outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path], outputs=[output_path])
                _print_result(
                    build_command_result(
                        command="compare",
                        inputs=[left_path, right_path],
                        outputs=outputs,
                        metrics={"table_rows": sum(1 for _ in output_path.read_text(encoding='utf-8').splitlines()[1:])},
                        data={"table_path": output_path},
                    ),
                    json_output=args.json,
                )
                return 0
            left_path = Path(args.left)
            right_path = Path(args.right)
            report = compare_tree_paths(left_path, right_path)
            outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path])
            _print_result(
                build_command_result(
                    command="compare",
                    inputs=[left_path, right_path],
                    outputs=outputs,
                    metrics={
                        "shared_taxa": len(report.shared_taxa),
                        "robinson_foulds_distance": report.robinson_foulds_distance,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "annotate":
            report = annotate_tree_against_table(args.tree, args.metadata, taxon_column=args.taxon_column)
            outputs: list[Path | str] = []
            if args.out is not None:
                outputs.append(write_annotation_report(args.out, report))
            if args.joined_out is not None:
                table = load_taxon_table(args.metadata, taxon_column=args.taxon_column)
                outputs.append(
                    write_taxon_rows(
                        args.joined_out,
                        columns=["taxon", "matched", *[column for column in table.columns if column != table.taxon_column]],
                        rows=[
                            {
                                "taxon": row.taxon,
                                "matched": str(row.matched).lower(),
                                **{column: row.values.get(column, "") for column in table.columns if column != table.taxon_column},
                            }
                            for row in report.joined_rows
                        ],
                    )
                )
            outputs = _finalize_outputs(args, command="annotate", inputs=[args.tree, args.metadata], outputs=outputs)
            _print_result(
                build_command_result(
                    command="annotate",
                    inputs=[args.tree, args.metadata],
                    outputs=outputs,
                    metrics={
                        "tree_taxa": report.tree_taxa,
                        "table_rows": report.table_rows,
                        "linked_taxa": report.linked_taxa,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "render":
            metadata_table = load_taxon_table(args.metadata, taxon_column=args.taxon_column) if args.metadata is not None else None
            traits_table = load_taxon_table(args.traits, taxon_column=args.taxon_column) if args.traits is not None else None
            labels: dict[str, str] | None = None
            if metadata_table is not None and args.label_column is not None:
                if args.label_column not in metadata_table.columns:
                    raise MetadataJoinError(f"metadata table does not contain label column '{args.label_column}'")
                labels = {
                    row[metadata_table.taxon_column]: row[args.label_column]
                    for row in metadata_table.rows
                    if row[args.label_column]
                }
            categorical_traits = (
                _build_string_trait_map(traits_table, args.categorical_column)
                if traits_table is not None and args.categorical_column is not None
                else None
            )
            continuous_traits = (
                _build_numeric_trait_map(traits_table, args.continuous_column)
                if traits_table is not None and args.continuous_column is not None
                else None
            )
            metadata_strips = (
                _build_annotation_strips(metadata_table, _split_csv_values(args.metadata_strip_columns))
                if metadata_table is not None
                else []
            )
            heatmap_columns = (
                _build_annotation_strips(traits_table, _split_csv_values(args.heatmap_columns))
                if traits_table is not None
                else []
            )
            collapsed_clades = _split_csv_values(args.collapse_clades)
            result = render_tree_svg(
                args.tree,
                out_path=args.out,
                labels=labels,
                layout=args.layout,
                show_support_values=args.support_labels,
                categorical_traits=categorical_traits,
                continuous_traits=continuous_traits,
                metadata_strips=metadata_strips,
                heatmap_columns=heatmap_columns,
                collapsed_clades=collapsed_clades,
            )
            inputs = [args.tree]
            if args.metadata is not None:
                inputs.append(args.metadata)
            if args.traits is not None:
                inputs.append(args.traits)
            outputs = [result.output_path]
            package_result = None
            if args.package_dir is not None:
                package_result = build_tree_figure_package(
                    args.tree,
                    out_dir=args.package_dir,
                    labels=labels,
                    layout=args.layout,
                    show_support_values=args.support_labels,
                    categorical_traits=categorical_traits,
                    continuous_traits=continuous_traits,
                    metadata_strips=metadata_strips,
                    heatmap_columns=heatmap_columns,
                    collapsed_clades=collapsed_clades,
                )
                outputs.append(package_result.output_dir)
            outputs = _finalize_outputs(args, command="render", inputs=inputs, outputs=outputs)
            if args.json:
                _print_result(
                    build_command_result(
                        command="render",
                        inputs=inputs,
                        outputs=outputs,
                        warnings=result.missing_metadata_labels,
                        metrics={
                            "tip_count": result.tip_count,
                            "label_count": result.label_count,
                            "rendered_support_count": result.rendered_support_count,
                            "rendered_categorical_trait_count": result.rendered_categorical_trait_count,
                            "rendered_continuous_trait_count": result.rendered_continuous_trait_count,
                            "rendered_metadata_strip_count": result.rendered_metadata_strip_count,
                            "rendered_heatmap_column_count": result.rendered_heatmap_column_count,
                            "collapsed_clade_count": result.collapsed_clade_count,
                        },
                        data={
                            "render": result,
                            "figure_package_dir": package_result.output_dir if package_result is not None else None,
                        },
                    ),
                    json_output=True,
                )
                return 0
            print(result.output_path)
            return 0
        if args.command == "evidence":
            if args.evidence_command == "bundle":
                report = bundle_directory(args.inputs, args.outputs, args.out)
                inputs = [*args.inputs, *args.outputs]
                outputs = _finalize_outputs(args, command="evidence", inputs=inputs, outputs=[args.out])
                _print_result(
                    build_command_result(
                        command="evidence",
                        inputs=inputs,
                        outputs=outputs,
                        metrics={
                            "file_count": report.file_count,
                            "input_file_count": report.input_file_count,
                            "output_file_count": report.output_file_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            report = validate_bundle(args.bundle_root)
            if not report.valid:
                raise EvidenceContractError(f"evidence bundle validation failed with {len(report.mismatches)} mismatch(es)")
            outputs = _finalize_outputs(args, command="evidence", inputs=[args.bundle_root])
            _print_result(
                build_command_result(
                    command="evidence",
                    inputs=[args.bundle_root],
                    outputs=outputs,
                    metrics={"file_count": report.file_count},
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "demo":
            if args.demo_command == "run":
                result = run_capability_demo(args.out)
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.tree_report,
                        result.dataset_report,
                        result.phylo_inputs_report,
                        result.comparison_report,
                        result.evidence_bundle,
                        result.capability_summary,
                    ],
                )
                if args.json:
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={"artifact_count": 6},
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            raise NotImplementedError(f"unsupported demo command: {args.demo_command}")
        if args.command == "report":
            if args.report_command == "tree":
                result = render_tree_report(tree_path=args.tree, out_path=args.out)
                outputs = _finalize_outputs(args, command="report", inputs=[args.tree], outputs=[result.output_path])
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=[args.tree],
                            outputs=outputs,
                            warnings=result.validation.warnings + result.inspection.warnings,
                            metrics={"tip_count": result.inspection.tip_count},
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_path)
                return 0
            if args.report_command == "dataset":
                result = render_dataset_report(
                    tree_path=args.tree,
                    metadata_path=args.metadata,
                    traits_path=args.traits,
                    alignment_path=args.alignment,
                    tip_dates_path=args.tip_dates,
                    calibration_path=args.calibrations,
                    out_path=args.out,
                )
                inputs = [args.tree, args.metadata]
                if args.traits is not None:
                    inputs.append(args.traits)
                if args.alignment is not None:
                    inputs.append(args.alignment)
                if args.tip_dates is not None:
                    inputs.append(args.tip_dates)
                if args.calibrations is not None:
                    inputs.append(args.calibrations)
                outputs = _finalize_outputs(args, command="report", inputs=inputs, outputs=[result.output_path])
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=inputs,
                            outputs=outputs,
                            warnings=result.validation.warnings + result.inspection.warnings,
                            metrics={
                                "tip_count": result.inspection.tip_count,
                                "linked_taxa": result.metadata_linkage.linked_taxa,
                                "readiness_decision": None if result.dataset_audit is None else result.dataset_audit.readiness_decision,
                                "excluded_taxa": 0 if result.dataset_audit is None else len(result.dataset_audit.exclusion_table.rows),
                                "blocked_analysis_count": 0 if result.dataset_audit is None else len(result.dataset_audit.blocked_analyses),
                                "risky_analysis_count": 0 if result.dataset_audit is None else sum(
                                    1
                                    for row in result.dataset_audit.analysis_decisions
                                    if row.decision == "risky"
                                ),
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_path)
                return 0
            if args.report_command == "phylo-inputs":
                result = render_phylo_inputs_report(
                    tree_path=args.tree,
                    alignment_path=args.alignment,
                    out_path=args.out,
                )
                inputs = [args.tree, args.alignment]
                outputs = _finalize_outputs(args, command="report", inputs=inputs, outputs=[result.output_path])
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=inputs,
                            outputs=outputs,
                            warnings=result.validation.warnings + result.inspection.warnings,
                            metrics={
                                "tip_count": result.inspection.tip_count,
                                "alignment_length": result.alignment.alignment_length,
                                "linked_taxa": result.alignment_linkage.linked_taxa,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_path)
                return 0
            raise NotImplementedError(f"unsupported report command: {args.report_command}")
        if args.command == "adapter":
            if args.adapter_command == "inspect":
                executable = args.executable or args.engine_name
                report = read_engine_version(
                    args.engine_name,
                    executable,
                    version_args=_adapter_version_args(args.engine_name),
                )
                outputs = _finalize_outputs(args, command="adapter", inputs=[args.engine_name])
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.engine_name],
                        outputs=outputs,
                        metrics={"version_line_count": len(report.text.splitlines())},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "report":
                report = render_inference_workflow_report(manifest_path=args.manifest_path, out_path=args.out)
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.manifest_path],
                    outputs=[report.output_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.manifest_path],
                        outputs=outputs,
                        metrics={"warning_count": report.warning_count},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "align":
                report = run_multiple_sequence_alignment(
                    args.input_path,
                    args.out,
                    executable=args.executable or "mafft",
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={"warning_count": len(report.run.warning_lines)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "trim":
                report = run_alignment_trimming(
                    args.input_path,
                    args.out,
                    executable=args.executable or "trimal",
                    gap_threshold=args.gap_threshold,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={"warning_count": len(report.run.warning_lines)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "model-select":
                report = run_model_selection(
                    args.input_path,
                    out_dir=args.out_dir,
                    prefix=args.prefix,
                    executable=args.executable or "iqtree2",
                    sequence_type=args.sequence_type,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={"selected_model": report.selected_model},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "infer-ml":
                report = run_maximum_likelihood_tree_inference(
                    args.input_path,
                    out_dir=args.out_dir,
                    model=args.model,
                    prefix=args.prefix,
                    executable=args.executable or "iqtree2",
                    sequence_type=args.sequence_type,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={"selected_model": report.selected_model},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "bootstrap":
                report = run_bootstrap_support_estimation(
                    args.input_path,
                    out_dir=args.out_dir,
                    model=args.model,
                    replicates=args.replicates,
                    prefix=args.prefix,
                    executable=args.executable or "iqtree2",
                    sequence_type=args.sequence_type,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={"bootstrap_replicates": args.replicates},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "consensus":
                report = run_bootstrap_consensus_tree(
                    args.input_path,
                    out_dir=args.out_dir,
                    prefix=args.prefix,
                    executable=args.executable or "iqtree2",
                    minimum_support=args.minimum_support,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={"minimum_support": args.minimum_support},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "infer-fast":
                report = run_fast_tree_inference(
                    args.input_path,
                    args.out,
                    executable=args.executable or "FastTree",
                    sequence_type=args.sequence_type,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={"warning_count": len(report.run.warning_lines)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-prepare":
                report = prepare_mrbayes_analysis(
                    args.input_path,
                    args.out,
                    model=args.model,
                    rates=args.rates,
                    ngen=args.ngen,
                    nchains=args.nchains,
                    samplefreq=args.samplefreq,
                    printfreq=args.printfreq,
                    burnin_fraction=args.burnin_fraction,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[args.out],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={"taxon_count": report.taxon_count, "character_count": report.character_count},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-run":
                report = run_mrbayes_posterior_inference(
                    args.input_path,
                    executable=args.executable or "mb",
                    resume=args.resume,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={"warning_count": len(report.run.warning_lines), "resumed": report.resumed},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-summarize":
                consensus_tree, report = summarize_mrbayes_posterior_trees(
                    args.input_path,
                    burnin_fraction=args.burnin_fraction,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[report.filtered_tree_set_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "kept_tree_count": report.kept_tree_count,
                            "rooted_topology_count": report.rooted_topology_count,
                            "tip_count": consensus_tree.tip_count,
                        },
                        data={"summary": report, "consensus_newick": report.consensus_newick},
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-traces":
                report = parse_mrbayes_parameter_traces(args.input_path)
                outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={"row_count": report.row_count, "column_count": len(report.columns)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-ess":
                report = compute_mrbayes_effective_sample_sizes(args.input_path)
                outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={"parameter_count": len(report.effective_sample_sizes)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-convergence":
                report = assess_mrbayes_convergence(
                    args.input_path,
                    ess_threshold=args.ess_threshold,
                    mean_shift_threshold=args.mean_shift_threshold,
                )
                outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=[warning["message"] for warning in report.warnings],
                        metrics={"warning_count": len(report.warnings), "converged": report.converged},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-report":
                report = render_bayesian_posterior_report(
                    posterior_tree_path=args.posterior_trees,
                    trace_path=args.traces,
                    out_path=args.out,
                    burnin_fraction=args.burnin_fraction,
                    ess_threshold=args.ess_threshold,
                    mean_shift_threshold=args.mean_shift_threshold,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.posterior_trees, args.traces],
                    outputs=[report.output_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.posterior_trees, args.traces],
                        outputs=outputs,
                        metrics={"kept_tree_count": report.kept_tree_count, "warning_count": report.warning_count},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-prepare":
                report = prepare_beast_time_tree_analysis(
                    args.input_path,
                    args.out,
                    tree_path=args.tree,
                    calibration_path=args.calibrations,
                    tip_dates_path=args.tip_dates,
                    clock_model=args.clock_model,
                    tree_prior=args.tree_prior,
                    chain_length=args.chain_length,
                    log_every=args.log_every,
                )
                outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path], outputs=[args.out])
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={"taxon_count": report.taxon_count, "calibration_count": report.calibration_count, "tip_date_count": report.tip_date_count},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-calibrations":
                report = validate_fossil_calibration_table(args.tree_path, args.calibration_path)
                outputs = _finalize_outputs(args, command="adapter", inputs=[args.tree_path, args.calibration_path])
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.tree_path, args.calibration_path],
                        outputs=outputs,
                        metrics={"calibration_count": report.calibration_count, "invalid_calibration_count": report.invalid_calibration_count},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-tip-dates":
                report = validate_tip_dating_metadata(
                    args.tree_path,
                    args.tip_dates_path,
                    alignment_path=args.alignment,
                    date_column=args.date_column,
                )
                inputs = [args.tree_path, args.tip_dates_path, *([args.alignment] if args.alignment is not None else [])]
                outputs = _finalize_outputs(args, command="adapter", inputs=inputs)
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=inputs,
                        outputs=outputs,
                        metrics={"valid_tip_count": report.valid_tip_count, "invalid_tip_count": report.invalid_tip_count},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-log":
                report = parse_beast_log(args.input_path)
                outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={"row_count": report.row_count, "column_count": len(report.columns)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-convergence":
                report = assess_beast_convergence(
                    args.input_path,
                    ess_threshold=args.ess_threshold,
                    mean_shift_threshold=args.mean_shift_threshold,
                )
                outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=[warning["message"] for warning in report.warnings],
                        metrics={"warning_count": len(report.warnings), "converged": report.converged},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-calibration-report":
                report = render_calibration_audit_report(
                    tree_path=args.tree_path,
                    calibration_path=args.calibration_path,
                    out_path=args.out,
                    tip_dates_path=args.tip_dates,
                    alignment_path=args.alignment,
                    date_column=args.date_column,
                )
                inputs = [args.tree_path, args.calibration_path, *([args.tip_dates] if args.tip_dates is not None else []), *([args.alignment] if args.alignment is not None else [])]
                outputs = _finalize_outputs(args, command="adapter", inputs=inputs, outputs=[report.output_path])
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=inputs,
                        outputs=outputs,
                        metrics={"invalid_calibration_count": report.invalid_calibration_count, "warning_count": report.warning_count},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "bayesian-evidence":
                report = build_bayesian_evidence_package(
                    bundle_root=args.out_dir,
                    input_paths=args.inputs,
                    config_paths=args.configs,
                    tree_paths=args.trees,
                    log_paths=args.logs,
                    diagnostic_paths=args.diagnostics,
                    report_paths=args.reports,
                )
                inputs = [*args.inputs, *args.configs, *args.trees, *args.logs, *args.diagnostics, *args.reports]
                outputs = _finalize_outputs(args, command="adapter", inputs=inputs, outputs=[args.out_dir])
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=inputs,
                        outputs=outputs,
                        metrics={"file_count": report.file_count, "valid": report.valid},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "compare":
                report = compare_fast_and_ml_trees(args.fast_tree, args.ml_tree, out_path=args.out)
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.fast_tree, args.ml_tree],
                    outputs=[report.comparison_report.output_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.fast_tree, args.ml_tree],
                        outputs=outputs,
                        metrics={
                            "shared_taxa": len(report.comparison_report.topology.shared_taxa),
                            "robinson_foulds_distance": report.comparison_report.topology.robinson_foulds_distance,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            raise EngineUnavailableError(f"unsupported adapter command: {args.adapter_command}")
    except PhylogeneticsError as error:
        if _json_requested(args):
            _print_result(
                build_error_result(command=args.command, inputs=_command_inputs(args), error=error),
                json_output=True,
            )
            return 2
        parser.exit(status=2, message=f"{error.code}: {error.message}\n")
    except FileNotFoundError as error:
        parser.exit(status=2, message=f"{error}\n")
    except ValueError as error:
        parser.exit(status=2, message=f"{error}\n")
    except NotImplementedError as error:
        parser.exit(status=2, message=f"{error}\n")
    except Exception as error:  # pragma: no cover - defensive CLI guard
        parser.exit(status=1, message=f"unexpected error: {error}\n")

    parser.print_help(sys.stderr)
    return 2
