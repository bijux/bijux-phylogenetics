# ruff: noqa: F401
from __future__ import annotations

import hashlib
from pathlib import Path
import re

from bijux_phylogenetics.compare.presentation import build_tree_comparison_report
from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet,
    AlignmentRecord,
    AlignmentSummary,
    CodingSequenceExclusion,
)
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    build_partition_summary_report,
    normalize_partition_data_type,
    parse_locus_partitions,
    slice_partition_sequence,
    write_locus_partitions,
    write_partition_summary_table,
)
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.runtime.errors import EngineWorkflowError, PhylogeneticsError
from bijux_phylogenetics.io.fasta import (
    infer_alignment_alphabet,
    load_fasta_alignment,
    write_fasta_alignment,
)
from bijux_phylogenetics.io.fasta.coding import (
    back_translate_aligned_coding_sequences,
    classify_sequence_coding_behavior,
    prepare_coding_sequences_for_alignment,
    translate_prepared_coding_sequences,
)
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.trees import load_tree_set

from .models import (
    AlignmentTrimmingSummary,
    CodonAwareAlignmentWorkflowReport,
    EngineWorkflowReport,
    ExternalTreeComparisonReport,
    IqtreeSupportValue,
    IqtreeWorkflowSummary,
    PreparedIqtreePartitions as _PreparedIqtreePartitions,
)
from .alignment import (
    list_mafft_alignment_modes,
    list_trimal_trimming_modes,
    resolve_mafft_alignment_mode,
    resolve_trimal_trimming_mode,
    run_alignment_trimming,
    run_codon_aware_multiple_sequence_alignment,
    run_multiple_sequence_alignment,
)
from .state import (
    _build_alignment_trimming_summary,
    _ensure_inference_ready_alignment,
    _manifest_path_from_output,
    _partition_alignment_file_name,
    _partition_support_path,
    _persist_workflow_report,
    _prefix_path,
    _record_output_validation_failure,
    _require_nonempty_text_output,
    _resolve_incomplete_workflow_state,
    _restore_codon_aware_alignment_report,
    _restore_workflow_report,
    _resume_existing_codon_aware_alignment,
    _resume_existing_workflow,
    _resume_has_bootstrap_review_outputs,
    _resume_has_fasttree_review_outputs,
    _resume_has_sh_alrt_review_outputs,
    _sidecar,
    _validate_alignment_output,
    _validate_incomplete_run_policy,
    _validate_iqtree_required_artifacts,
    _validate_support_value_count,
    _validate_tree_output,
    _validate_tree_set_output,
    _write_coding_exclusion_table,
    _write_coding_summary_table,
)
from ..artifacts.bootstrap import (
    build_bootstrap_support_histogram_rows,
    build_bootstrap_support_rows,
    build_low_support_bootstrap_rows,
    write_bootstrap_support_histogram,
    write_bootstrap_support_table,
)
from ..common import (
    EngineRunReport,
    EngineVersionInfo,
    active_engine_run_is_live,
    build_engine_output_error,
    build_file_checksums,
    cleanup_incomplete_engine_run,
    clear_incomplete_engine_run,
    engine_active_marker_path,
    engine_incomplete_marker_path,
    execute_engine_command,
    load_active_engine_run,
    load_engine_manifest,
    load_incomplete_engine_run,
    load_unaligned_fasta,
    observe_engine_outputs,
    read_engine_version,
    resolve_engine_executable,
    update_incomplete_engine_run,
    validate_timeout_seconds,
    write_engine_manifest,
)
from ..artifacts.fasttree import (
    build_fasttree_low_support_rows,
    build_fasttree_support_histogram_rows,
    build_fasttree_support_rows,
    write_fasttree_support_histogram,
    write_fasttree_support_table,
)
from ..artifacts.iqtree import (
    IqtreeModelCandidate,
    IqtreeModelSelectionSummary,
    parse_best_model_file,
    parse_iqtree_model_selection_summary,
    parse_log_likelihood_file,
    resolve_iqtree_model_sidecar,
    write_iqtree_model_candidates_table,
)
from ..artifacts.sh_alrt import (
    build_conflicting_sh_alrt_support_rows,
    build_sh_alrt_support_rows,
    write_sh_alrt_support_table,
)
from ..validation import (
    BootstrapSupportNode,
    BootstrapSupportSummaryReport,
    FastTreeSupportNode,
    FastTreeSupportSummaryReport,
    ShAlrtSupportNode,
    ShAlrtSupportSummaryReport,
    WeakBackboneReport,
    detect_weakly_supported_backbone,
    summarize_bootstrap_support_distribution,
    summarize_fasttree_support_distribution,
    summarize_sh_alrt_support_distribution,
)

_MINIMUM_UFBOOT_REPLICATES = 1000


def _iqtree_partition_supports_fixed_model(
    *,
    model: str,
    mixed_data_types: bool,
) -> bool:
    if not mixed_data_types:
        return True
    normalized = model.strip().upper()
    return normalized in {
        "TEST",
        "TESTONLY",
        "TESTNEWONLY",
        "MF",
        "MFP",
        "TESTMERGE",
        "TESTMERGEONLY",
        "MF+MERGE",
        "MFP+MERGE",
    }

def _iqtree_sequence_type_flag(
    path: Path, sequence_type: AlignmentAlphabet | None
) -> list[str]:
    detected = sequence_type
    if detected is None:
        detected = infer_alignment_alphabet(load_fasta_alignment(path))
    if detected in {"dna", "rna"}:
        return ["-st", "DNA"]
    if detected == "protein":
        return ["-st", "AA"]
    return []


def _iqtree_execution_controls(*, seed: int, threads: int) -> list[str]:
    if seed < 1:
        raise ValueError(f"iqtree seed must be positive, got {seed}")
    if threads < 1:
        raise ValueError(f"iqtree threads must be positive, got {threads}")
    return ["-seed", str(seed), "-nt", str(threads)]


def _validate_ufboot_replicates(replicates: int) -> None:
    if replicates < _MINIMUM_UFBOOT_REPLICATES:
        raise EngineWorkflowError(
            "iqtree ultrafast bootstrap requires at least "
            f"{_MINIMUM_UFBOOT_REPLICATES} replicates, got {replicates}"
        )


def _validate_sh_alrt_replicates(replicates: int) -> None:
    if replicates < 1:
        raise ValueError(f"sh-alrt replicates must be positive, got {replicates}")


def _parse_best_model_artifact(prefix_path: Path) -> str | None:
    for candidate in (
        prefix_path.with_suffix(".iqtree"),
        resolve_iqtree_model_sidecar(prefix_path),
    ):
        if candidate is None:
            continue
        model = parse_best_model_file(candidate)
        if model is not None:
            return model
    return None


def _parse_log_likelihood_artifact(prefix_path: Path) -> float | None:
    for candidate in (
        prefix_path.with_suffix(".iqtree"),
        prefix_path.with_suffix(".log"),
    ):
        if not candidate.exists():
            continue
        log_likelihood = parse_log_likelihood_file(candidate)
        if log_likelihood is not None:
            return log_likelihood
    return None


def _validate_iqtree_model_result(
    prefix_path: Path,
    *,
    workflow: str,
    default_selected_model: str | None = None,
) -> str:
    selected_model = _parse_best_model_artifact(prefix_path)
    if selected_model is None and default_selected_model is not None:
        return default_selected_model
    if selected_model is None:
        raise build_engine_output_error(
            f"iqtree {workflow} did not expose a parsable best-fit model result",
            code="engine_model_result_missing",
            engine_name="iqtree",
            workflow=workflow,
            path=prefix_path.with_suffix(".iqtree"),
            output_name="iqtree_report",
            artifact_kind="iqtree-model-result",
            details={
                "model_sidecar_path": (
                    None
                    if resolve_iqtree_model_sidecar(prefix_path) is None
                    else str(resolve_iqtree_model_sidecar(prefix_path))
                )
            },
        )
    return selected_model


def _existing_iqtree_outputs(
    prefix_path: Path,
    *,
    include_tree: bool = False,
    include_bootstrap: bool = False,
    include_consensus: bool = False,
) -> dict[str, Path]:
    outputs: dict[str, Path] = {}
    for key, candidate in (
        ("iqtree_report", prefix_path.with_suffix(".iqtree")),
        ("iqtree_log", prefix_path.with_suffix(".log")),
    ):
        if candidate.exists():
            outputs[key] = candidate
    model_sidecar = resolve_iqtree_model_sidecar(prefix_path)
    if model_sidecar is not None:
        outputs["model_selection_sidecar"] = model_sidecar
    tree_candidate = prefix_path.with_suffix(".treefile")
    if include_tree and tree_candidate.exists():
        outputs["tree"] = tree_candidate
    bootstrap_candidate = prefix_path.with_suffix(".ufboot")
    if include_bootstrap and bootstrap_candidate.exists():
        outputs["bootstrap_trees"] = bootstrap_candidate
    consensus_candidate = prefix_path.with_suffix(".contree")
    if include_consensus and consensus_candidate.exists():
        outputs["consensus_tree"] = consensus_candidate
    return outputs


def _build_iqtree_summary(
    prefix_path: Path,
    *,
    default_selected_model: str | None,
    support_tree_path: Path | None = None,
) -> IqtreeWorkflowSummary:
    selected_model = _parse_best_model_artifact(prefix_path) or default_selected_model
    log_likelihood = _parse_log_likelihood_artifact(prefix_path)
    support_values: list[IqtreeSupportValue] = []
    minimum_support: float | None = None
    maximum_support: float | None = None
    if support_tree_path is not None and support_tree_path.exists():
        support_summary = summarize_bootstrap_support_distribution(support_tree_path)
        support_values = [
            IqtreeSupportValue(
                node=node.node,
                descendant_taxa=list(node.descendant_taxa),
                support=node.support,
                support_fraction=node.support_fraction,
                is_backbone=node.is_backbone,
            )
            for node in support_summary.nodes
        ]
        minimum_support = support_summary.minimum_support
        maximum_support = support_summary.maximum_support
    return IqtreeWorkflowSummary(
        selected_model=selected_model,
        log_likelihood=log_likelihood,
        support_value_count=len(support_values),
        minimum_support=minimum_support,
        maximum_support=maximum_support,
        support_values=support_values,
    )


def _build_iqtree_model_selection_summary(
    prefix_path: Path,
) -> IqtreeModelSelectionSummary | None:
    return parse_iqtree_model_selection_summary(
        iqtree_report_path=prefix_path.with_suffix(".iqtree"),
        model_sidecar_path=resolve_iqtree_model_sidecar(prefix_path),
    )


def _prepare_iqtree_partitions(
    input_path: Path,
    partition_path: Path,
    *,
    prefix_path: Path,
) -> _PreparedIqtreePartitions:
    records = load_fasta_alignment(input_path)
    alignment_summary = summarise_fasta(input_path)
    partitions = parse_locus_partitions(partition_path)
    summary = build_partition_summary_report(
        partitions,
        alignment_length=alignment_summary.alignment_length,
    )
    summary_path = _partition_support_path(prefix_path, "partition-summary.tsv")
    write_partition_summary_table(summary_path, summary)
    notes = [
        f"validated {summary.partition_count} partitions across {summary.assigned_site_count} assigned sites",
    ]
    output_paths: dict[str, Path] = {
        "partition_summary": summary_path,
    }

    declared_types = {
        normalize_partition_data_type(partition.data_type)
        for partition in partitions
        if partition.data_type is not None
    }
    if len(declared_types) <= 1:
        normalized_partition_path = _partition_support_path(
            prefix_path, "partition-scheme.partitions"
        )
        write_locus_partitions(normalized_partition_path, partitions)
        output_paths["partition_scheme"] = normalized_partition_path
        notes.append(
            "prepared a normalized partition scheme for single-alignment IQ-TREE analysis"
        )
        return _PreparedIqtreePartitions(
            command_args=[
                "-s",
                str(input_path.resolve()),
                "-p",
                str(normalized_partition_path.resolve()),
            ],
            summary=summary,
            output_paths=output_paths,
            notes=notes,
            mixed_data_types=False,
        )

    if any(partition.data_type is None for partition in partitions):
        raise EngineWorkflowError(
            "mixed partition analyses require every partition to declare a data_type"
        )
    unsupported_types = sorted(
        {
            data_type
            for data_type in declared_types
            if data_type is not None and data_type not in {"DNA", "RNA", "PROTEIN"}
        }
    )
    if unsupported_types:
        raise EngineWorkflowError(
            "mixed partition analyses currently support only DNA, RNA, and PROTEIN datatypes; "
            f"got: {', '.join(unsupported_types)}"
        )

    partition_alignment_dir = _partition_support_path(
        prefix_path, "partition-alignments"
    )
    partition_alignment_dir.mkdir(parents=True, exist_ok=True)
    lines = ["#nexus", "begin sets;"]
    for partition in partitions:
        partition_alignment_path = (
            partition_alignment_dir / _partition_alignment_file_name(partition)
        )
        write_fasta_alignment(
            partition_alignment_path,
            [
                AlignmentRecord(
                    identifier=record.identifier,
                    sequence=slice_partition_sequence(record.sequence, partition),
                )
                for record in records
            ],
        )
        output_paths[f"partition_alignment_{partition.name}"] = partition_alignment_path
        lines.append(
            f"    charset {partition.name} = {partition_alignment_path.name}: *;"
        )
    lines.append("end;")
    mixed_partition_path = _partition_support_path(prefix_path, "partition-scheme.nex")
    mixed_partition_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    output_paths["partition_scheme"] = mixed_partition_path
    notes.append(
        "prepared a mixed-datatype NEXUS partition scheme with one extracted alignment per partition"
    )
    return _PreparedIqtreePartitions(
        command_args=["-p", str(mixed_partition_path.resolve())],
        summary=summary,
        output_paths=output_paths,
        notes=notes,
        mixed_data_types=True,
    )


def run_model_selection(
    input_path: Path,
    *,
    out_dir: Path,
    prefix: str = "model-selection",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
    partition_path: Path | None = None,
    resume: bool = False,
    seed: int = 1,
    threads: int = 1,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run a model-selection workflow on an aligned FASTA file."""
    _ensure_inference_ready_alignment(input_path)
    validate_timeout_seconds(timeout_seconds)
    prefix_path = _prefix_path(out_dir, prefix)
    manifest_path = prefix_path.with_suffix(".manifest.json")
    version = read_engine_version(
        "iqtree",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    iqtree_report_path = prefix_path.with_suffix(".iqtree")
    iqtree_log_path = prefix_path.with_suffix(".log")
    prepared_partitions = (
        None
        if partition_path is None
        else _prepare_iqtree_partitions(
            input_path,
            partition_path,
            prefix_path=prefix_path,
        )
    )
    command = [
        resolved,
        *(
            prepared_partitions.command_args
            if prepared_partitions is not None
            else ["-s", str(input_path.resolve())]
        ),
        *(
            []
            if prepared_partitions is not None and prepared_partitions.mixed_data_types
            else _iqtree_sequence_type_flag(input_path, sequence_type)
        ),
        *_iqtree_execution_controls(seed=seed, threads=threads),
        "-m",
        "MF",
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=(
                [input_path] if partition_path is None else [input_path, partition_path]
            ),
            expected_command=command,
            expected_version=version,
        )
        if resumed is not None:
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="model-selection",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "iqtree_report": iqtree_report_path,
            "iqtree_log": iqtree_log_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _validate_iqtree_required_artifacts(prefix_path, workflow="model-selection")
        selected_model = _validate_iqtree_model_result(
            prefix_path,
            workflow="model-selection",
        )
        iqtree_summary = _build_iqtree_summary(
            prefix_path,
            default_selected_model=selected_model,
        )
        model_selection_summary = _build_iqtree_model_selection_summary(prefix_path)
        if (
            model_selection_summary is None
            or model_selection_summary.candidate_count < 1
        ):
            raise build_engine_output_error(
                "iqtree model-selection did not expose a parsable candidate-model table",
                code="iqtree_model_candidates_missing",
                engine_name="iqtree",
                workflow="model-selection",
                path=iqtree_report_path,
                output_name="iqtree_report",
                artifact_kind="iqtree-model-candidates",
            )
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    selected_model_path = prefix_path.with_suffix(".selected-model.txt")
    selected_model_path.write_text(selected_model + "\n", encoding="utf-8")
    model_candidates_path = prefix_path.with_suffix(".model-candidates.tsv")
    write_iqtree_model_candidates_table(model_candidates_path, model_selection_summary)
    report = EngineWorkflowReport(
        workflow="model-selection",
        engine_name="iqtree",
        input_paths=(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_paths={
            **({} if prepared_partitions is None else prepared_partitions.output_paths),
            **_existing_iqtree_outputs(prefix_path),
            "selected_model": selected_model_path,
            "model_candidates": model_candidates_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_checksums={},
        config={
            "sequence_type": sequence_type,
            "partition_path": None if partition_path is None else str(partition_path),
            "seed": seed,
            "threads": threads,
            "timeout_seconds": timeout_seconds,
        },
        selected_model=selected_model,
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        model_selection_summary=model_selection_summary,
        notes=[
            *([] if prepared_partitions is None else prepared_partitions.notes),
            f"iqtree random seed: {seed}",
            f"iqtree threads: {threads}",
            "best-fit substitution model parsed from engine output",
            f"parsed {model_selection_summary.candidate_count} candidate substitution models from iqtree output",
            *(
                []
                if model_selection_summary.selected_criterion is None
                else [
                    "model-selection workflow exposed the governing information criterion"
                ]
            ),
            *(
                []
                if iqtree_summary.log_likelihood is None
                else [
                    "model-selection workflow exposed a parsable log-likelihood score"
                ]
            ),
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)


def run_maximum_likelihood_tree_inference(
    input_path: Path,
    *,
    out_dir: Path,
    model: str,
    prefix: str = "maximum-likelihood",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
    partition_path: Path | None = None,
    resume: bool = False,
    seed: int = 1,
    threads: int = 1,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run an external maximum-likelihood tree inference workflow."""
    _ensure_inference_ready_alignment(input_path)
    validate_timeout_seconds(timeout_seconds)
    prefix_path = _prefix_path(out_dir, prefix)
    manifest_path = prefix_path.with_suffix(".manifest.json")
    version = read_engine_version(
        "iqtree",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    prepared_partitions = (
        None
        if partition_path is None
        else _prepare_iqtree_partitions(
            input_path,
            partition_path,
            prefix_path=prefix_path,
        )
    )
    if prepared_partitions is not None and not _iqtree_partition_supports_fixed_model(
        model=model,
        mixed_data_types=prepared_partitions.mixed_data_types,
    ):
        raise EngineWorkflowError(
            "mixed DNA/protein partition analyses require a model-selection keyword such as MF, MFP, TEST, or TESTMERGE"
        )
    tree_path = prefix_path.with_suffix(".treefile")
    report_path = prefix_path.with_suffix(".iqtree")
    log_path = prefix_path.with_suffix(".log")
    command = [
        resolved,
        *(
            prepared_partitions.command_args
            if prepared_partitions is not None
            else ["-s", str(input_path.resolve())]
        ),
        *(
            []
            if prepared_partitions is not None and prepared_partitions.mixed_data_types
            else _iqtree_sequence_type_flag(input_path, sequence_type)
        ),
        *_iqtree_execution_controls(seed=seed, threads=threads),
        "-m",
        model,
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=(
                [input_path] if partition_path is None else [input_path, partition_path]
            ),
            expected_command=command,
            expected_version=version,
        )
        if resumed is not None:
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="maximum-likelihood-tree",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "tree": tree_path,
            "iqtree_report": report_path,
            "iqtree_log": log_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _validate_iqtree_required_artifacts(
            prefix_path,
            workflow="maximum-likelihood-tree",
        )
        _validate_tree_output(
            tree_path,
            engine_name="iqtree",
            workflow="maximum-likelihood-tree",
            output_name="tree",
            artifact_kind="maximum-likelihood-tree",
        )
        selected_model = _validate_iqtree_model_result(
            prefix_path,
            workflow="maximum-likelihood-tree",
            default_selected_model=model,
        )
        iqtree_summary = _build_iqtree_summary(
            prefix_path,
            default_selected_model=selected_model,
            support_tree_path=tree_path,
        )
        model_selection_summary = _build_iqtree_model_selection_summary(prefix_path)
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    report = EngineWorkflowReport(
        workflow="maximum-likelihood-tree",
        engine_name="iqtree",
        input_paths=(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_paths={
            **({} if prepared_partitions is None else prepared_partitions.output_paths),
            **_existing_iqtree_outputs(prefix_path, include_tree=True),
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_checksums={},
        config={
            "model": model,
            "sequence_type": sequence_type,
            "partition_path": None if partition_path is None else str(partition_path),
            "seed": seed,
            "threads": threads,
            "timeout_seconds": timeout_seconds,
        },
        selected_model=iqtree_summary.selected_model,
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        model_selection_summary=model_selection_summary,
        notes=[
            *([] if prepared_partitions is None else prepared_partitions.notes),
            f"iqtree random seed: {seed}",
            f"iqtree threads: {threads}",
            "maximum-likelihood tree validated as parseable Newick output",
            *(
                []
                if iqtree_summary.log_likelihood is None
                else ["log-likelihood parsed from iqtree inference artifacts"]
            ),
            *(
                []
                if iqtree_summary.support_value_count == 0
                else ["support values parsed from the inferred maximum-likelihood tree"]
            ),
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)


def run_bootstrap_support_estimation(
    input_path: Path,
    *,
    out_dir: Path,
    model: str,
    replicates: int = 1000,
    prefix: str = "bootstrap-support",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
    partition_path: Path | None = None,
    resume: bool = False,
    seed: int = 1,
    threads: int = 1,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run external bootstrap support estimation and retain bootstrap trees."""
    if replicates < 1:
        raise ValueError(f"replicates must be positive, got {replicates}")
    _validate_ufboot_replicates(replicates)
    _ensure_inference_ready_alignment(input_path)
    validate_timeout_seconds(timeout_seconds)
    prefix_path = _prefix_path(out_dir, prefix)
    manifest_path = prefix_path.with_suffix(".manifest.json")
    version = read_engine_version(
        "iqtree",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    prepared_partitions = (
        None
        if partition_path is None
        else _prepare_iqtree_partitions(
            input_path,
            partition_path,
            prefix_path=prefix_path,
        )
    )
    if prepared_partitions is not None and not _iqtree_partition_supports_fixed_model(
        model=model,
        mixed_data_types=prepared_partitions.mixed_data_types,
    ):
        raise EngineWorkflowError(
            "mixed DNA/protein partition analyses require a model-selection keyword such as MF, MFP, TEST, or TESTMERGE"
        )
    support_tree_path = prefix_path.with_suffix(".treefile")
    bootstrap_tree_path = prefix_path.with_suffix(".ufboot")
    support_table_path = prefix_path.with_suffix(".support.tsv")
    low_support_branches_path = prefix_path.with_suffix(".low-support.tsv")
    support_histogram_path = prefix_path.with_suffix(".support-histogram.tsv")
    report_path = prefix_path.with_suffix(".iqtree")
    log_path = prefix_path.with_suffix(".log")
    command = [
        resolved,
        *(
            prepared_partitions.command_args
            if prepared_partitions is not None
            else ["-s", str(input_path.resolve())]
        ),
        *(
            []
            if prepared_partitions is not None and prepared_partitions.mixed_data_types
            else _iqtree_sequence_type_flag(input_path, sequence_type)
        ),
        *_iqtree_execution_controls(seed=seed, threads=threads),
        "-m",
        model,
        "-bb",
        str(replicates),
        "-wbt",
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=(
                [input_path] if partition_path is None else [input_path, partition_path]
            ),
            expected_command=command,
            expected_version=version,
        )
        if resumed is not None and _resume_has_bootstrap_review_outputs(resumed):
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="bootstrap-support",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "support_tree": support_tree_path,
            "bootstrap_trees": bootstrap_tree_path,
            "iqtree_report": report_path,
            "iqtree_log": log_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _validate_iqtree_required_artifacts(prefix_path, workflow="bootstrap-support")
        _validate_tree_output(
            support_tree_path,
            engine_name="iqtree",
            workflow="bootstrap-support",
            output_name="support_tree",
            artifact_kind="bootstrap-supported-tree",
        )
        _validate_tree_set_output(
            bootstrap_tree_path,
            engine_name="iqtree",
            workflow="bootstrap-support",
            output_name="bootstrap_trees",
            artifact_kind="bootstrap-tree-set",
        )
        selected_model = _validate_iqtree_model_result(
            prefix_path,
            workflow="bootstrap-support",
            default_selected_model=model,
        )
        iqtree_summary = _build_iqtree_summary(
            prefix_path,
            default_selected_model=selected_model,
            support_tree_path=support_tree_path,
        )
        bootstrap_support_summary = summarize_bootstrap_support_distribution(
            support_tree_path
        )
        _validate_support_value_count(
            engine_name="iqtree",
            workflow="bootstrap-support",
            path=support_tree_path,
            output_name="support_tree",
            artifact_kind="bootstrap-supported-tree",
            support_value_count=bootstrap_support_summary.supported_node_count,
            support_kind="bootstrap support",
        )
        weak_backbone_report = detect_weakly_supported_backbone(support_tree_path)
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    write_bootstrap_support_table(
        support_table_path,
        build_bootstrap_support_rows(bootstrap_support_summary),
    )
    write_bootstrap_support_table(
        low_support_branches_path,
        build_low_support_bootstrap_rows(bootstrap_support_summary),
    )
    write_bootstrap_support_histogram(
        support_histogram_path,
        build_bootstrap_support_histogram_rows(bootstrap_support_summary),
    )
    model_selection_summary = _build_iqtree_model_selection_summary(prefix_path)
    report = EngineWorkflowReport(
        workflow="bootstrap-support",
        engine_name="iqtree",
        input_paths=(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_paths={
            **({} if prepared_partitions is None else prepared_partitions.output_paths),
            **_existing_iqtree_outputs(
                prefix_path,
                include_tree=False,
                include_bootstrap=True,
                include_consensus=True,
            ),
            "support_tree": support_tree_path,
            "support_table": support_table_path,
            "low_support_branches": low_support_branches_path,
            "support_histogram": support_histogram_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_checksums={},
        config={
            "model": model,
            "replicates": replicates,
            "sequence_type": sequence_type,
            "partition_path": None if partition_path is None else str(partition_path),
            "seed": seed,
            "threads": threads,
            "timeout_seconds": timeout_seconds,
        },
        selected_model=iqtree_summary.selected_model,
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        model_selection_summary=model_selection_summary,
        bootstrap_support_summary=bootstrap_support_summary,
        weak_backbone_report=weak_backbone_report,
        notes=[
            *([] if prepared_partitions is None else prepared_partitions.notes),
            f"iqtree random seed: {seed}",
            f"iqtree threads: {threads}",
            "bootstrap tree set retained for downstream consensus construction",
            *(
                []
                if iqtree_summary.log_likelihood is None
                else ["log-likelihood parsed from iqtree bootstrap inference artifacts"]
            ),
            *(
                []
                if iqtree_summary.support_value_count == 0
                else [
                    "support values parsed from the bootstrap-supported tree artifact"
                ]
            ),
            "branch-level support table exported for bootstrap review",
            "low-support branch ledger exported for weak-clade review",
            "support histogram exported for reviewer-facing support distribution checks",
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)


def run_sh_alrt_support_estimation(
    input_path: Path,
    *,
    out_dir: Path,
    model: str,
    sh_alrt_replicates: int = 1000,
    bootstrap_replicates: int = 1000,
    prefix: str = "sh-alrt-support",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
    partition_path: Path | None = None,
    resume: bool = False,
    seed: int = 1,
    threads: int = 1,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run combined SH-aLRT and UFBoot branch-support estimation."""
    _validate_sh_alrt_replicates(sh_alrt_replicates)
    _validate_ufboot_replicates(bootstrap_replicates)
    _ensure_inference_ready_alignment(input_path)
    validate_timeout_seconds(timeout_seconds)
    prefix_path = _prefix_path(out_dir, prefix)
    manifest_path = prefix_path.with_suffix(".manifest.json")
    version = read_engine_version(
        "iqtree",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    prepared_partitions = (
        None
        if partition_path is None
        else _prepare_iqtree_partitions(
            input_path,
            partition_path,
            prefix_path=prefix_path,
        )
    )
    if prepared_partitions is not None and not _iqtree_partition_supports_fixed_model(
        model=model,
        mixed_data_types=prepared_partitions.mixed_data_types,
    ):
        raise EngineWorkflowError(
            "mixed DNA/protein partition analyses require a model-selection keyword such as MF, MFP, TEST, or TESTMERGE"
        )
    support_tree_path = prefix_path.with_suffix(".treefile")
    bootstrap_tree_path = prefix_path.with_suffix(".ufboot")
    support_table_path = prefix_path.with_suffix(".support.tsv")
    conflicting_support_branches_path = prefix_path.with_suffix(
        ".conflicting-support.tsv"
    )
    report_path = prefix_path.with_suffix(".iqtree")
    log_path = prefix_path.with_suffix(".log")
    command = [
        resolved,
        *(
            prepared_partitions.command_args
            if prepared_partitions is not None
            else ["-s", str(input_path.resolve())]
        ),
        *(
            []
            if prepared_partitions is not None and prepared_partitions.mixed_data_types
            else _iqtree_sequence_type_flag(input_path, sequence_type)
        ),
        *_iqtree_execution_controls(seed=seed, threads=threads),
        "-m",
        model,
        "-alrt",
        str(sh_alrt_replicates),
        "-bb",
        str(bootstrap_replicates),
        "-wbt",
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=(
                [input_path] if partition_path is None else [input_path, partition_path]
            ),
            expected_command=command,
            expected_version=version,
        )
        if resumed is not None and _resume_has_sh_alrt_review_outputs(resumed):
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="sh-alrt-support",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "support_tree": support_tree_path,
            "bootstrap_trees": bootstrap_tree_path,
            "iqtree_report": report_path,
            "iqtree_log": log_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _validate_iqtree_required_artifacts(prefix_path, workflow="sh-alrt-support")
        _validate_tree_output(
            support_tree_path,
            engine_name="iqtree",
            workflow="sh-alrt-support",
            output_name="support_tree",
            artifact_kind="sh-alrt-supported-tree",
        )
        _validate_tree_set_output(
            bootstrap_tree_path,
            engine_name="iqtree",
            workflow="sh-alrt-support",
            output_name="bootstrap_trees",
            artifact_kind="bootstrap-tree-set",
        )
        selected_model = _validate_iqtree_model_result(
            prefix_path,
            workflow="sh-alrt-support",
            default_selected_model=model,
        )
        iqtree_summary = _build_iqtree_summary(
            prefix_path,
            default_selected_model=selected_model,
            support_tree_path=support_tree_path,
        )
        bootstrap_support_summary = summarize_bootstrap_support_distribution(
            support_tree_path
        )
        _validate_support_value_count(
            engine_name="iqtree",
            workflow="sh-alrt-support",
            path=support_tree_path,
            output_name="support_tree",
            artifact_kind="sh-alrt-supported-tree",
            support_value_count=bootstrap_support_summary.supported_node_count,
            support_kind="ultrafast bootstrap support",
        )
        sh_alrt_support_summary = summarize_sh_alrt_support_distribution(
            support_tree_path
        )
        _validate_support_value_count(
            engine_name="iqtree",
            workflow="sh-alrt-support",
            path=support_tree_path,
            output_name="support_tree",
            artifact_kind="sh-alrt-supported-tree",
            support_value_count=sh_alrt_support_summary.annotated_node_count,
            support_kind="sh-alrt support",
        )
        _validate_support_value_count(
            engine_name="iqtree",
            workflow="sh-alrt-support",
            path=support_tree_path,
            output_name="support_tree",
            artifact_kind="sh-alrt-supported-tree",
            support_value_count=sh_alrt_support_summary.fully_scored_node_count,
            support_kind="joint sh-alrt and ultrafast bootstrap support",
        )
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    write_sh_alrt_support_table(
        support_table_path,
        build_sh_alrt_support_rows(sh_alrt_support_summary),
    )
    write_sh_alrt_support_table(
        conflicting_support_branches_path,
        build_conflicting_sh_alrt_support_rows(sh_alrt_support_summary),
    )
    model_selection_summary = _build_iqtree_model_selection_summary(prefix_path)
    report = EngineWorkflowReport(
        workflow="sh-alrt-support",
        engine_name="iqtree",
        input_paths=(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_paths={
            **({} if prepared_partitions is None else prepared_partitions.output_paths),
            **_existing_iqtree_outputs(
                prefix_path,
                include_tree=False,
                include_bootstrap=True,
                include_consensus=True,
            ),
            "support_tree": support_tree_path,
            "support_table": support_table_path,
            "conflicting_support_branches": conflicting_support_branches_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_checksums={},
        config={
            "model": model,
            "sh_alrt_replicates": sh_alrt_replicates,
            "bootstrap_replicates": bootstrap_replicates,
            "sequence_type": sequence_type,
            "partition_path": None if partition_path is None else str(partition_path),
            "seed": seed,
            "threads": threads,
            "timeout_seconds": timeout_seconds,
        },
        selected_model=iqtree_summary.selected_model,
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        model_selection_summary=model_selection_summary,
        bootstrap_support_summary=bootstrap_support_summary,
        sh_alrt_support_summary=sh_alrt_support_summary,
        notes=[
            *([] if prepared_partitions is None else prepared_partitions.notes),
            f"iqtree random seed: {seed}",
            f"iqtree threads: {threads}",
            f"sh-alrt replicates: {sh_alrt_replicates}",
            f"ultrafast bootstrap replicates: {bootstrap_replicates}",
            "combined sh-alrt and ufboot branch-support table exported for review",
            "conflicting sh-alrt versus ufboot support signals exported for review",
            *(
                []
                if iqtree_summary.log_likelihood is None
                else ["log-likelihood parsed from iqtree support inference artifacts"]
            ),
            *(
                []
                if sh_alrt_support_summary.annotated_node_count == 0
                else [
                    "combined sh-alrt and ufboot support values parsed from the supported tree artifact"
                ]
            ),
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)


def run_bootstrap_consensus_tree(
    bootstrap_trees_path: Path,
    *,
    out_dir: Path,
    prefix: str = "bootstrap-consensus",
    executable: str | Path = "iqtree2",
    minimum_support: float = 0.5,
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Construct a consensus tree from bootstrap trees."""
    if not 0.0 <= minimum_support <= 1.0:
        raise ValueError(
            f"minimum_support must be between 0 and 1 inclusive, got {minimum_support}"
        )
    if not bootstrap_trees_path.exists():
        raise FileNotFoundError(bootstrap_trees_path)
    validate_timeout_seconds(timeout_seconds)
    prefix_path = _prefix_path(out_dir, prefix)
    manifest_path = prefix_path.with_suffix(".manifest.json")
    version = read_engine_version(
        "iqtree",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    consensus_tree_path = prefix_path.with_suffix(".contree")
    log_path = prefix_path.with_suffix(".log")
    command = [
        resolved,
        "-t",
        str(bootstrap_trees_path.resolve()),
        "-con",
        "-minsup",
        str(minimum_support),
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=[bootstrap_trees_path],
            expected_command=command,
            expected_version=version,
        )
        if resumed is not None:
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="bootstrap-consensus",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "consensus_tree": consensus_tree_path,
            "iqtree_log": log_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _require_nonempty_text_output(
            log_path,
            engine_name="iqtree",
            workflow="bootstrap-consensus",
            output_name="iqtree_log",
            artifact_kind="iqtree-log",
        )
        _validate_tree_output(
            consensus_tree_path,
            engine_name="iqtree",
            workflow="bootstrap-consensus",
            output_name="consensus_tree",
            artifact_kind="bootstrap-consensus-tree",
        )
        iqtree_summary = _build_iqtree_summary(
            prefix_path,
            default_selected_model=None,
            support_tree_path=consensus_tree_path,
        )
        _validate_support_value_count(
            engine_name="iqtree",
            workflow="bootstrap-consensus",
            path=consensus_tree_path,
            output_name="consensus_tree",
            artifact_kind="bootstrap-consensus-tree",
            support_value_count=iqtree_summary.support_value_count,
            support_kind="bootstrap consensus support",
        )
        model_selection_summary = _build_iqtree_model_selection_summary(prefix_path)
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    report = EngineWorkflowReport(
        workflow="bootstrap-consensus",
        engine_name="iqtree",
        input_paths=[bootstrap_trees_path],
        output_paths=_existing_iqtree_outputs(prefix_path, include_consensus=True),
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([bootstrap_trees_path]),
        output_checksums={},
        config={
            "minimum_support": minimum_support,
            "timeout_seconds": timeout_seconds,
        },
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        model_selection_summary=model_selection_summary,
        notes=[
            "consensus tree validated as parseable Newick output",
            *(
                []
                if iqtree_summary.support_value_count == 0
                else ["support values parsed from the bootstrap consensus tree"]
            ),
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)
