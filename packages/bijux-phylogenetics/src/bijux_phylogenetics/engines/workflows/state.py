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
from ..fasttree_artifacts import (
    build_fasttree_low_support_rows,
    build_fasttree_support_histogram_rows,
    build_fasttree_support_rows,
    write_fasttree_support_histogram,
    write_fasttree_support_table,
)
from ..iqtree_artifacts import (
    IqtreeModelCandidate,
    IqtreeModelSelectionSummary,
    parse_best_model_file,
    parse_iqtree_model_selection_summary,
    parse_log_likelihood_file,
    resolve_iqtree_model_sidecar,
    write_iqtree_model_candidates_table,
)
from ..sh_alrt_artifacts import (
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

_INCOMPLETE_RUN_POLICIES = {"reject", "clean"}

def _sidecar(path: Path, label: str) -> Path:
    return path.parent / f"{path.name}.{label}"


def _prefix_path(out_dir: Path, prefix: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / prefix


def _manifest_path_from_output(path: Path) -> Path:
    return _sidecar(path, "manifest.json")


def _validate_incomplete_run_policy(policy: str) -> str:
    if policy not in _INCOMPLETE_RUN_POLICIES:
        available = ", ".join(sorted(_INCOMPLETE_RUN_POLICIES))
        raise ValueError(
            f"incomplete_run_policy must be one of: {available}; got {policy}"
        )
    return policy


def _resolve_incomplete_workflow_state(
    *,
    manifest_path: Path,
    incomplete_run_policy: str,
) -> list[str]:
    _validate_incomplete_run_policy(incomplete_run_policy)
    active_record = load_active_engine_run(manifest_path)
    if active_record is not None and active_engine_run_is_live(active_record):
        raise EngineWorkflowError(
            "engine workflow is already running for the requested output manifest",
            code="engine_workflow_already_running",
            details={
                "manifest_path": str(manifest_path),
                "marker_path": str(engine_active_marker_path(manifest_path)),
                "running_process_id": active_record.process_id,
                "running_workflow": active_record.workflow,
                "running_engine_name": active_record.engine_name,
            },
        )
    record = load_incomplete_engine_run(manifest_path)
    if record is None:
        return []
    if incomplete_run_policy == "clean":
        cleanup_incomplete_engine_run(manifest_path)
        return [
            "removed outputs from a previously incomplete engine run before restarting"
        ]
    marker_path = engine_incomplete_marker_path(manifest_path)
    observed_outputs = [
        {
            "output_name": observation.output_name,
            "path": str(observation.path),
            "exists": observation.exists,
            "path_kind": observation.path_kind,
            "size_bytes": observation.size_bytes,
            "sha256": observation.sha256,
        }
        for observation in record.observed_outputs
    ]
    raise EngineWorkflowError(
        "a previous engine run left incomplete outputs and resume could not safely "
        f"reuse them; marker: {marker_path}",
        code="engine_incomplete_outputs_present",
        details={
            "manifest_path": str(manifest_path),
            "marker_path": str(marker_path),
            "engine_name": record.engine_name,
            "workflow": record.workflow,
            "failure_reason": (
                record.failure_reason
                if record.failure_reason is not None
                else "engine_run_incomplete"
            ),
            "failure_message": record.failure_message,
            "timed_out": record.timed_out,
            "exit_code": record.exit_code,
            "timeout_seconds": record.timeout_seconds,
            "missing_output_names": list(record.missing_output_names),
            "observed_outputs": observed_outputs,
            "incomplete_run_policy": incomplete_run_policy,
            "available_actions": ["resume", "clean"],
        },
    )


def _partition_support_path(prefix_path: Path, suffix: str) -> Path:
    return prefix_path.parent / f"{prefix_path.name}.{suffix}"


def _record_output_validation_failure(
    manifest_path: Path,
    run: EngineRunReport,
    error: PhylogeneticsError,
) -> None:
    observations = observe_engine_outputs(run.output_paths)
    update_incomplete_engine_run(
        manifest_path,
        ended_at_utc=run.ended_at_utc,
        timed_out=run.timed_out,
        exit_code=run.exit_code,
        failure_reason=error.code,
        failure_message=(
            f"{run.engine_name} {run.workflow} produced outputs that failed "
            f"validation: {error.code}"
        ),
        missing_output_names=[
            observation.output_name for observation in observations if not observation.exists
        ],
        observed_outputs=observations,
    )


def _require_nonempty_text_output(
    path: Path,
    *,
    engine_name: str,
    workflow: str,
    output_name: str,
    artifact_kind: str,
) -> str:
    if not path.exists():
        raise build_engine_output_error(
            f"{engine_name} {workflow} did not produce required output '{output_name}': {path}",
            code="engine_required_output_missing",
            engine_name=engine_name,
            workflow=workflow,
            path=path,
            output_name=output_name,
            artifact_kind=artifact_kind,
        )
    text = path.read_text(encoding="utf-8", errors="replace")
    if not text.strip():
        raise build_engine_output_error(
            f"{engine_name} {workflow} produced an empty required output '{output_name}': {path}",
            code="engine_output_empty",
            engine_name=engine_name,
            workflow=workflow,
            path=path,
            output_name=output_name,
            artifact_kind=artifact_kind,
        )
    return text


def _validate_alignment_output(
    path: Path,
    *,
    engine_name: str,
    workflow: str,
    output_name: str,
    artifact_kind: str,
) -> AlignmentSummary:
    _require_nonempty_text_output(
        path,
        engine_name=engine_name,
        workflow=workflow,
        output_name=output_name,
        artifact_kind=artifact_kind,
    )
    records = load_fasta_alignment(path)
    if not records or len(records[0].sequence) < 1:
        raise build_engine_output_error(
            f"{engine_name} {workflow} produced an empty alignment after validation: {path}",
            code="engine_output_empty",
            engine_name=engine_name,
            workflow=workflow,
            path=path,
            output_name=output_name,
            artifact_kind=artifact_kind,
        )
    return summarise_fasta(path)


def _write_coding_exclusion_table(
    path: Path, exclusions: list[CodingSequenceExclusion]
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "identifier",
                "comparable_length",
                "reason",
                "invalid_codon_count",
                "premature_stop_count",
                "terminal_stop_count",
                "trailing_bases",
                "note",
            ]
        )
    ]
    lines.extend(
        "\t".join(
            [
                row.identifier,
                str(row.comparable_length),
                row.reason,
                str(row.invalid_codon_count),
                str(row.premature_stop_count),
                str(row.terminal_stop_count),
                str(row.trailing_bases),
                row.note,
            ]
        )
        for row in exclusions
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_coding_summary_table(
    path: Path,
    *,
    input_path: Path,
    genetic_code: int,
    exclusions: list[CodingSequenceExclusion],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    exclusion_by_identifier = {row.identifier: row for row in exclusions}
    behaviors = classify_sequence_coding_behavior(
        input_path,
        genetic_code=genetic_code,
    )
    header = "\t".join(
        [
            "identifier",
            "status",
            "comparable_length",
            "divisible_by_three",
            "invalid_codon_count",
            "premature_stop_count",
            "terminal_stop_count",
            "exclusion_reason",
            "note",
        ]
    )
    rows = [header]
    for behavior in behaviors:
        exclusion = exclusion_by_identifier.get(behavior.identifier)
        rows.append(
            "\t".join(
                [
                    behavior.identifier,
                    "excluded" if exclusion is not None else "accepted",
                    str(behavior.comparable_length),
                    "yes" if behavior.divisible_by_three else "no",
                    str(behavior.invalid_codon_count),
                    str(behavior.premature_stop_count),
                    str(behavior.terminal_stop_count),
                    "" if exclusion is None else exclusion.reason,
                    behavior.note,
                ]
            )
        )
    ordered_rows = [rows[0], *sorted(rows[1:])]
    path.write_text("\n".join(ordered_rows) + "\n", encoding="utf-8")
    return path


def _build_alignment_trimming_summary(
    *,
    mode: str,
    gap_threshold: float,
    input_summary: AlignmentSummary,
    trimmed_summary: AlignmentSummary,
) -> AlignmentTrimmingSummary:
    if trimmed_summary.alignment_length > input_summary.alignment_length:
        raise EngineWorkflowError(
            "trimmed alignment is longer than the input alignment, which is not a valid trimAl result"
        )
    retained_site_count = trimmed_summary.alignment_length
    removed_site_count = (
        input_summary.alignment_length - trimmed_summary.alignment_length
    )
    retained_site_fraction = retained_site_count / input_summary.alignment_length
    removed_site_fraction = removed_site_count / input_summary.alignment_length
    return AlignmentTrimmingSummary(
        mode=mode,
        gap_threshold=gap_threshold if mode == "gap-threshold" else None,
        input_alignment_length=input_summary.alignment_length,
        trimmed_alignment_length=trimmed_summary.alignment_length,
        retained_site_count=retained_site_count,
        removed_site_count=removed_site_count,
        retained_site_fraction=retained_site_fraction,
        removed_site_fraction=removed_site_fraction,
        input_gap_fraction=input_summary.gap_fraction,
        trimmed_gap_fraction=trimmed_summary.gap_fraction,
        input_gap_percentage=input_summary.gap_fraction * 100.0,
        trimmed_gap_percentage=trimmed_summary.gap_fraction * 100.0,
    )


def _validate_tree_output(
    path: Path,
    *,
    engine_name: str,
    workflow: str,
    output_name: str,
    artifact_kind: str,
) -> None:
    tree_text = _require_nonempty_text_output(
        path,
        engine_name=engine_name,
        workflow=workflow,
        output_name=output_name,
        artifact_kind=artifact_kind,
    )
    loads_newick(tree_text)
    validate_tree_path(path)


def _validate_tree_set_output(
    path: Path,
    *,
    engine_name: str,
    workflow: str,
    output_name: str,
    artifact_kind: str,
) -> None:
    _require_nonempty_text_output(
        path,
        engine_name=engine_name,
        workflow=workflow,
        output_name=output_name,
        artifact_kind=artifact_kind,
    )
    load_tree_set(path)


def _validate_iqtree_required_artifacts(
    prefix_path: Path,
    *,
    workflow: str,
) -> None:
    _require_nonempty_text_output(
        prefix_path.with_suffix(".iqtree"),
        engine_name="iqtree",
        workflow=workflow,
        output_name="iqtree_report",
        artifact_kind="iqtree-report",
    )
    _require_nonempty_text_output(
        prefix_path.with_suffix(".log"),
        engine_name="iqtree",
        workflow=workflow,
        output_name="iqtree_log",
        artifact_kind="iqtree-log",
    )


def _validate_support_value_count(
    *,
    engine_name: str,
    workflow: str,
    path: Path,
    output_name: str,
    artifact_kind: str,
    support_value_count: int,
    support_kind: str,
) -> None:
    if support_value_count > 0:
        return
    raise build_engine_output_error(
        f"{engine_name} {workflow} did not expose any parsable {support_kind} values in '{output_name}': {path}",
        code="engine_support_values_missing",
        engine_name=engine_name,
        workflow=workflow,
        path=path,
        output_name=output_name,
        artifact_kind=artifact_kind,
        details={"support_kind": support_kind},
    )


def _ensure_inference_ready_alignment(path: Path) -> None:
    records = load_fasta_alignment(path)
    if not records or len(records[0].sequence) < 1:
        raise EngineWorkflowError(
            f"inference alignment is empty after filtering: {path}"
        )
    summarise_fasta(path)


def _partition_alignment_file_name(partition: LocusPartition) -> str:
    normalized = re.sub(r"[^A-Za-z0-9._-]+", "-", partition.name.strip().lower())
    normalized = normalized.strip("-._") or "partition"
    digest = hashlib.sha1(partition.name.encode("utf-8")).hexdigest()[:8]  # nosec B324
    return f"{normalized}-{digest}.fasta"


def _restore_workflow_report(payload: dict[str, object]) -> EngineWorkflowReport:
    run_payload = payload["run"]
    if not isinstance(run_payload, dict):
        raise EngineWorkflowError("engine manifest contains an invalid run payload")
    version_payload = run_payload["version"]
    if not isinstance(version_payload, dict):
        raise EngineWorkflowError("engine manifest contains an invalid version payload")
    version = EngineVersionInfo(
        engine_name=str(version_payload["engine_name"]),
        executable=str(version_payload["executable"]),
        command=[str(item) for item in version_payload["command"]],
        text=str(version_payload["text"]),
    )
    run = EngineRunReport(
        engine_name=str(run_payload["engine_name"]),
        workflow=str(run_payload["workflow"]),
        executable=str(run_payload["executable"]),
        working_directory=Path(run_payload["working_directory"]),
        version=version,
        command=[str(item) for item in run_payload["command"]],
        exit_code=int(run_payload["exit_code"]),
        stdout_path=Path(run_payload["stdout_path"]),
        stderr_path=Path(run_payload["stderr_path"]),
        output_paths={
            str(key): Path(value)
            for key, value in dict(run_payload["output_paths"]).items()
        },
        warning_lines=[str(item) for item in run_payload["warning_lines"]],
        started_at_utc=str(run_payload.get("started_at_utc", "")),
        ended_at_utc=str(run_payload.get("ended_at_utc", "")),
        runtime_seconds=float(run_payload.get("runtime_seconds", 0.0)),
        timeout_seconds=(
            None
            if run_payload.get("timeout_seconds") is None
            else float(run_payload["timeout_seconds"])
        ),
        timed_out=bool(run_payload.get("timed_out", False)),
    )
    return EngineWorkflowReport(
        workflow=str(payload["workflow"]),
        engine_name=str(payload["engine_name"]),
        input_paths=[Path(item) for item in payload["input_paths"]],
        output_paths={
            str(key): Path(value)
            for key, value in dict(payload["output_paths"]).items()
        },
        run=run,
        manifest_path=Path(payload["manifest_path"]),
        input_checksums={
            str(key): str(value)
            for key, value in dict(payload.get("input_checksums", {})).items()
        },
        output_checksums={
            str(key): str(value)
            for key, value in dict(payload.get("output_checksums", {})).items()
        },
        config={
            str(key): value for key, value in dict(payload.get("config", {})).items()
        },
        selected_model=None
        if payload.get("selected_model") is None
        else str(payload["selected_model"]),
        log_likelihood=(
            None
            if payload.get("log_likelihood") is None
            else float(payload["log_likelihood"])
        ),
        iqtree_summary=(
            None
            if payload.get("iqtree_summary") is None
            else IqtreeWorkflowSummary(
                selected_model=(
                    None
                    if dict(payload["iqtree_summary"]).get("selected_model") is None
                    else str(dict(payload["iqtree_summary"])["selected_model"])
                ),
                log_likelihood=(
                    None
                    if dict(payload["iqtree_summary"]).get("log_likelihood") is None
                    else float(dict(payload["iqtree_summary"])["log_likelihood"])
                ),
                support_value_count=int(
                    dict(payload["iqtree_summary"])["support_value_count"]
                ),
                minimum_support=(
                    None
                    if dict(payload["iqtree_summary"]).get("minimum_support") is None
                    else float(dict(payload["iqtree_summary"])["minimum_support"])
                ),
                maximum_support=(
                    None
                    if dict(payload["iqtree_summary"]).get("maximum_support") is None
                    else float(dict(payload["iqtree_summary"])["maximum_support"])
                ),
                support_values=[
                    IqtreeSupportValue(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        support=float(dict(item)["support"]),
                        support_fraction=float(dict(item)["support_fraction"]),
                        is_backbone=bool(dict(item)["is_backbone"]),
                    )
                    for item in list(dict(payload["iqtree_summary"])["support_values"])
                ],
            )
        ),
        model_selection_summary=(
            None
            if payload.get("model_selection_summary") is None
            else IqtreeModelSelectionSummary(
                selected_model=(
                    None
                    if dict(payload["model_selection_summary"]).get("selected_model")
                    is None
                    else str(dict(payload["model_selection_summary"])["selected_model"])
                ),
                selected_criterion=(
                    None
                    if dict(payload["model_selection_summary"]).get(
                        "selected_criterion"
                    )
                    is None
                    else str(
                        dict(payload["model_selection_summary"])["selected_criterion"]
                    )
                ),
                best_model_aic=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_model_aic")
                    is None
                    else str(dict(payload["model_selection_summary"])["best_model_aic"])
                ),
                best_model_aicc=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_model_aicc")
                    is None
                    else str(
                        dict(payload["model_selection_summary"])["best_model_aicc"]
                    )
                ),
                best_model_bic=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_model_bic")
                    is None
                    else str(dict(payload["model_selection_summary"])["best_model_bic"])
                ),
                best_score_aic=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_score_aic")
                    is None
                    else float(
                        dict(payload["model_selection_summary"])["best_score_aic"]
                    )
                ),
                best_score_aicc=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_score_aicc")
                    is None
                    else float(
                        dict(payload["model_selection_summary"])["best_score_aicc"]
                    )
                ),
                best_score_bic=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_score_bic")
                    is None
                    else float(
                        dict(payload["model_selection_summary"])["best_score_bic"]
                    )
                ),
                candidate_count=int(
                    dict(payload["model_selection_summary"])["candidate_count"]
                ),
                candidates=[
                    IqtreeModelCandidate(
                        rank=int(dict(item)["rank"]),
                        model=str(dict(item)["model"]),
                        log_likelihood=float(dict(item)["log_likelihood"]),
                        parameter_count=(
                            None
                            if dict(item).get("parameter_count") is None
                            else int(dict(item)["parameter_count"])
                        ),
                        aic=float(dict(item)["aic"]),
                        aicc=float(dict(item)["aicc"]),
                        bic=float(dict(item)["bic"]),
                    )
                    for item in list(
                        dict(payload["model_selection_summary"])["candidates"]
                    )
                ],
                bic_near_best_models=[
                    str(item)
                    for item in list(
                        dict(payload["model_selection_summary"])["bic_near_best_models"]
                    )
                ],
            )
        ),
        bootstrap_support_summary=(
            None
            if payload.get("bootstrap_support_summary") is None
            else BootstrapSupportSummaryReport(
                tree_path=Path(dict(payload["bootstrap_support_summary"])["tree_path"]),
                internal_node_count=int(
                    dict(payload["bootstrap_support_summary"])["internal_node_count"]
                ),
                supported_node_count=int(
                    dict(payload["bootstrap_support_summary"])["supported_node_count"]
                ),
                minimum_support=(
                    None
                    if dict(payload["bootstrap_support_summary"]).get("minimum_support")
                    is None
                    else float(
                        dict(payload["bootstrap_support_summary"])["minimum_support"]
                    )
                ),
                maximum_support=(
                    None
                    if dict(payload["bootstrap_support_summary"]).get("maximum_support")
                    is None
                    else float(
                        dict(payload["bootstrap_support_summary"])["maximum_support"]
                    )
                ),
                median_support=(
                    None
                    if dict(payload["bootstrap_support_summary"]).get("median_support")
                    is None
                    else float(
                        dict(payload["bootstrap_support_summary"])["median_support"]
                    )
                ),
                weakly_supported_clade_count=int(
                    dict(payload["bootstrap_support_summary"])[
                        "weakly_supported_clade_count"
                    ]
                ),
                support_histogram={
                    str(key): int(value)
                    for key, value in dict(
                        dict(payload["bootstrap_support_summary"])["support_histogram"]
                    ).items()
                },
                nodes=[
                    BootstrapSupportNode(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        support=float(dict(item)["support"]),
                        support_fraction=float(dict(item)["support_fraction"]),
                        is_backbone=bool(dict(item)["is_backbone"]),
                    )
                    for item in list(
                        dict(payload["bootstrap_support_summary"])["nodes"]
                    )
                ],
                warnings=[
                    str(item)
                    for item in list(
                        dict(payload["bootstrap_support_summary"])["warnings"]
                    )
                ],
            )
        ),
        fasttree_support_summary=(
            None
            if payload.get("fasttree_support_summary") is None
            else FastTreeSupportSummaryReport(
                tree_path=Path(dict(payload["fasttree_support_summary"])["tree_path"]),
                internal_node_count=int(
                    dict(payload["fasttree_support_summary"])["internal_node_count"]
                ),
                annotated_node_count=int(
                    dict(payload["fasttree_support_summary"])["annotated_node_count"]
                ),
                minimum_local_support=(
                    None
                    if dict(payload["fasttree_support_summary"]).get(
                        "minimum_local_support"
                    )
                    is None
                    else float(
                        dict(payload["fasttree_support_summary"])[
                            "minimum_local_support"
                        ]
                    )
                ),
                maximum_local_support=(
                    None
                    if dict(payload["fasttree_support_summary"]).get(
                        "maximum_local_support"
                    )
                    is None
                    else float(
                        dict(payload["fasttree_support_summary"])[
                            "maximum_local_support"
                        ]
                    )
                ),
                median_local_support=(
                    None
                    if dict(payload["fasttree_support_summary"]).get(
                        "median_local_support"
                    )
                    is None
                    else float(
                        dict(payload["fasttree_support_summary"])[
                            "median_local_support"
                        ]
                    )
                ),
                weakly_supported_clade_count=int(
                    dict(payload["fasttree_support_summary"])[
                        "weakly_supported_clade_count"
                    ]
                ),
                support_histogram={
                    str(key): int(value)
                    for key, value in dict(
                        dict(payload["fasttree_support_summary"])["support_histogram"]
                    ).items()
                },
                approximate_method=bool(
                    dict(payload["fasttree_support_summary"])["approximate_method"]
                ),
                support_label_kind=str(
                    dict(payload["fasttree_support_summary"])["support_label_kind"]
                ),
                support_scale=str(
                    dict(payload["fasttree_support_summary"])["support_scale"]
                ),
                nodes=[
                    FastTreeSupportNode(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        local_support=float(dict(item)["local_support"]),
                        support_fraction=float(dict(item)["support_fraction"]),
                        is_backbone=bool(dict(item)["is_backbone"]),
                    )
                    for item in list(dict(payload["fasttree_support_summary"])["nodes"])
                ],
                warnings=[
                    str(item)
                    for item in list(
                        dict(payload["fasttree_support_summary"])["warnings"]
                    )
                ],
            )
        ),
        sh_alrt_support_summary=(
            None
            if payload.get("sh_alrt_support_summary") is None
            else ShAlrtSupportSummaryReport(
                tree_path=Path(dict(payload["sh_alrt_support_summary"])["tree_path"]),
                internal_node_count=int(
                    dict(payload["sh_alrt_support_summary"])["internal_node_count"]
                ),
                annotated_node_count=int(
                    dict(payload["sh_alrt_support_summary"])["annotated_node_count"]
                ),
                fully_scored_node_count=int(
                    dict(payload["sh_alrt_support_summary"])["fully_scored_node_count"]
                ),
                minimum_sh_alrt_support=(
                    None
                    if dict(payload["sh_alrt_support_summary"]).get(
                        "minimum_sh_alrt_support"
                    )
                    is None
                    else float(
                        dict(payload["sh_alrt_support_summary"])[
                            "minimum_sh_alrt_support"
                        ]
                    )
                ),
                maximum_sh_alrt_support=(
                    None
                    if dict(payload["sh_alrt_support_summary"]).get(
                        "maximum_sh_alrt_support"
                    )
                    is None
                    else float(
                        dict(payload["sh_alrt_support_summary"])[
                            "maximum_sh_alrt_support"
                        ]
                    )
                ),
                minimum_ufboot_support=(
                    None
                    if dict(payload["sh_alrt_support_summary"]).get(
                        "minimum_ufboot_support"
                    )
                    is None
                    else float(
                        dict(payload["sh_alrt_support_summary"])[
                            "minimum_ufboot_support"
                        ]
                    )
                ),
                maximum_ufboot_support=(
                    None
                    if dict(payload["sh_alrt_support_summary"]).get(
                        "maximum_ufboot_support"
                    )
                    is None
                    else float(
                        dict(payload["sh_alrt_support_summary"])[
                            "maximum_ufboot_support"
                        ]
                    )
                ),
                weak_sh_alrt_clade_count=int(
                    dict(payload["sh_alrt_support_summary"])["weak_sh_alrt_clade_count"]
                ),
                weak_ufboot_clade_count=int(
                    dict(payload["sh_alrt_support_summary"])["weak_ufboot_clade_count"]
                ),
                conflicting_support_signal_count=int(
                    dict(payload["sh_alrt_support_summary"])[
                        "conflicting_support_signal_count"
                    ]
                ),
                nodes=[
                    ShAlrtSupportNode(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        sh_alrt_support=(
                            None
                            if dict(item).get("sh_alrt_support") is None
                            else float(dict(item)["sh_alrt_support"])
                        ),
                        sh_alrt_support_fraction=(
                            None
                            if dict(item).get("sh_alrt_support_fraction") is None
                            else float(dict(item)["sh_alrt_support_fraction"])
                        ),
                        ufboot_support=(
                            None
                            if dict(item).get("ufboot_support") is None
                            else float(dict(item)["ufboot_support"])
                        ),
                        ufboot_support_fraction=(
                            None
                            if dict(item).get("ufboot_support_fraction") is None
                            else float(dict(item)["ufboot_support_fraction"])
                        ),
                        is_backbone=bool(dict(item)["is_backbone"]),
                        sh_alrt_strong=bool(dict(item)["sh_alrt_strong"]),
                        ufboot_strong=bool(dict(item)["ufboot_strong"]),
                        conflicting_support_signal=bool(
                            dict(item)["conflicting_support_signal"]
                        ),
                        support_agreement=str(dict(item)["support_agreement"]),
                    )
                    for item in list(dict(payload["sh_alrt_support_summary"])["nodes"])
                ],
                warnings=[
                    str(item)
                    for item in list(
                        dict(payload["sh_alrt_support_summary"])["warnings"]
                    )
                ],
            )
        ),
        weak_backbone_report=(
            None
            if payload.get("weak_backbone_report") is None
            else WeakBackboneReport(
                tree_path=Path(dict(payload["weak_backbone_report"])["tree_path"]),
                threshold=float(dict(payload["weak_backbone_report"])["threshold"]),
                evaluated_backbone_node_count=int(
                    dict(payload["weak_backbone_report"])[
                        "evaluated_backbone_node_count"
                    ]
                ),
                weak_backbone_node_count=int(
                    dict(payload["weak_backbone_report"])["weak_backbone_node_count"]
                ),
                weak_nodes=[
                    BootstrapSupportNode(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        support=float(dict(item)["support"]),
                        support_fraction=float(dict(item)["support_fraction"]),
                        is_backbone=bool(dict(item)["is_backbone"]),
                    )
                    for item in list(
                        dict(payload["weak_backbone_report"])["weak_nodes"]
                    )
                ],
                warnings=[
                    str(item)
                    for item in list(dict(payload["weak_backbone_report"])["warnings"])
                ],
            )
        ),
        trimming_summary=None
        if payload.get("trimming_summary") is None
        else AlignmentTrimmingSummary(
            mode=str(dict(payload["trimming_summary"])["mode"]),
            gap_threshold=(
                None
                if dict(payload["trimming_summary"]).get("gap_threshold") is None
                else float(dict(payload["trimming_summary"])["gap_threshold"])
            ),
            input_alignment_length=int(
                dict(payload["trimming_summary"])["input_alignment_length"]
            ),
            trimmed_alignment_length=int(
                dict(payload["trimming_summary"])["trimmed_alignment_length"]
            ),
            retained_site_count=int(
                dict(payload["trimming_summary"])["retained_site_count"]
            ),
            removed_site_count=int(
                dict(payload["trimming_summary"])["removed_site_count"]
            ),
            retained_site_fraction=float(
                dict(payload["trimming_summary"])["retained_site_fraction"]
            ),
            removed_site_fraction=float(
                dict(payload["trimming_summary"])["removed_site_fraction"]
            ),
            input_gap_fraction=float(
                dict(payload["trimming_summary"])["input_gap_fraction"]
            ),
            trimmed_gap_fraction=float(
                dict(payload["trimming_summary"])["trimmed_gap_fraction"]
            ),
            input_gap_percentage=float(
                dict(payload["trimming_summary"])["input_gap_percentage"]
            ),
            trimmed_gap_percentage=float(
                dict(payload["trimming_summary"])["trimmed_gap_percentage"]
            ),
        ),
        resumed=bool(payload.get("resumed", False)),
        notes=[str(item) for item in payload.get("notes", [])],
    )


def _resume_existing_workflow(
    *,
    manifest_path: Path,
    input_paths: list[Path],
    expected_command: list[str],
    expected_version: EngineVersionInfo,
) -> EngineWorkflowReport | None:
    if not manifest_path.exists():
        return None
    payload = load_engine_manifest(manifest_path)
    report = _restore_workflow_report(payload)
    if report.run.command != expected_command:
        return None
    if report.run.version.text != expected_version.text:
        return None
    current_input_checksums = build_file_checksums(input_paths)
    if report.input_checksums != current_input_checksums:
        return None
    if any(not path.exists() for path in report.output_paths.values()):
        return None
    current_output_checksums = build_file_checksums(list(report.output_paths.values()))
    if report.output_checksums != current_output_checksums:
        return None
    clear_incomplete_engine_run(manifest_path)
    report.resumed = True
    return report


def _persist_workflow_report(report: EngineWorkflowReport) -> EngineWorkflowReport:
    report.output_checksums = build_file_checksums(list(report.output_paths.values()))
    clear_incomplete_engine_run(report.manifest_path)
    write_engine_manifest(report.manifest_path, report)
    return report


def _restore_codon_aware_alignment_report(
    payload: dict[str, object],
) -> CodonAwareAlignmentWorkflowReport:
    run_payload = payload["run"]
    if not isinstance(run_payload, dict):
        raise EngineWorkflowError("engine manifest contains an invalid run payload")
    version_payload = run_payload["version"]
    if not isinstance(version_payload, dict):
        raise EngineWorkflowError("engine manifest contains an invalid version payload")
    version = EngineVersionInfo(
        engine_name=str(version_payload["engine_name"]),
        executable=str(version_payload["executable"]),
        command=[str(item) for item in version_payload["command"]],
        text=str(version_payload["text"]),
    )
    run = EngineRunReport(
        engine_name=str(run_payload["engine_name"]),
        workflow=str(run_payload["workflow"]),
        executable=str(run_payload["executable"]),
        working_directory=Path(run_payload["working_directory"]),
        version=version,
        command=[str(item) for item in run_payload["command"]],
        exit_code=int(run_payload["exit_code"]),
        stdout_path=Path(run_payload["stdout_path"]),
        stderr_path=Path(run_payload["stderr_path"]),
        output_paths={
            str(key): Path(value)
            for key, value in dict(run_payload["output_paths"]).items()
        },
        warning_lines=[str(item) for item in run_payload["warning_lines"]],
        started_at_utc=str(run_payload.get("started_at_utc", "")),
        ended_at_utc=str(run_payload.get("ended_at_utc", "")),
        runtime_seconds=float(run_payload.get("runtime_seconds", 0.0)),
        timeout_seconds=(
            None
            if run_payload.get("timeout_seconds") is None
            else float(run_payload["timeout_seconds"])
        ),
        timed_out=bool(run_payload.get("timed_out", False)),
    )
    return CodonAwareAlignmentWorkflowReport(
        workflow=str(payload["workflow"]),
        engine_name=str(payload["engine_name"]),
        input_path=Path(payload["input_path"]),
        output_paths={
            str(key): Path(value)
            for key, value in dict(payload["output_paths"]).items()
        },
        run=run,
        manifest_path=Path(payload["manifest_path"]),
        input_checksums={
            str(key): str(value)
            for key, value in dict(payload.get("input_checksums", {})).items()
        },
        output_checksums={
            str(key): str(value)
            for key, value in dict(payload.get("output_checksums", {})).items()
        },
        config={
            str(key): value for key, value in dict(payload.get("config", {})).items()
        },
        sequence_type=str(payload["sequence_type"]),
        genetic_code_id=int(payload.get("genetic_code_id", 1)),
        genetic_code_name=str(payload.get("genetic_code_name", "Standard")),
        input_sequence_count=int(payload.get("input_sequence_count", 0)),
        accepted_sequence_count=int(payload["accepted_sequence_count"]),
        invalid_codon_sequence_count=int(
            payload.get("invalid_codon_sequence_count", 0)
        ),
        excluded_sequences=[
            CodingSequenceExclusion(
                identifier=str(item["identifier"]),
                comparable_length=int(item["comparable_length"]),
                reason=str(item["reason"]),
                invalid_codon_count=int(item.get("invalid_codon_count", 0)),
                premature_stop_count=int(item["premature_stop_count"]),
                terminal_stop_count=int(item["terminal_stop_count"]),
                trailing_bases=int(item["trailing_bases"]),
                note=str(item["note"]),
            )
            for item in payload.get("excluded_sequences", [])
        ],
        terminal_stop_sequence_count=int(payload["terminal_stop_sequence_count"]),
        notes=[str(item) for item in payload.get("notes", [])],
        warnings=[str(item) for item in payload.get("warnings", [])],
        resumed=bool(payload.get("resumed", False)),
    )


def _resume_existing_codon_aware_alignment(
    *,
    manifest_path: Path,
    input_path: Path,
    expected_command: list[str],
    expected_version: EngineVersionInfo,
    expected_sequence_type: AlignmentAlphabet,
    expected_genetic_code_id: int,
) -> CodonAwareAlignmentWorkflowReport | None:
    if not manifest_path.exists():
        return None
    payload = load_engine_manifest(manifest_path)
    report = _restore_codon_aware_alignment_report(payload)
    if report.run.command != expected_command:
        return None
    if report.run.version.text != expected_version.text:
        return None
    if report.sequence_type != expected_sequence_type:
        return None
    if report.genetic_code_id != expected_genetic_code_id:
        return None
    current_input_checksums = build_file_checksums([input_path])
    if report.input_checksums != current_input_checksums:
        return None
    if any(not path.exists() for path in report.output_paths.values()):
        return None
    current_output_checksums = build_file_checksums(list(report.output_paths.values()))
    if report.output_checksums != current_output_checksums:
        return None
    clear_incomplete_engine_run(manifest_path)
    report.resumed = True
    return report


def _resume_has_bootstrap_review_outputs(report: EngineWorkflowReport) -> bool:
    return (
        report.bootstrap_support_summary is not None
        and report.weak_backbone_report is not None
        and report.output_paths.get("support_table") is not None
        and report.output_paths.get("low_support_branches") is not None
        and report.output_paths.get("support_histogram") is not None
    )


def _resume_has_fasttree_review_outputs(report: EngineWorkflowReport) -> bool:
    return (
        report.fasttree_support_summary is not None
        and report.output_paths.get("support_table") is not None
        and report.output_paths.get("low_support_branches") is not None
        and report.output_paths.get("support_histogram") is not None
    )


def _resume_has_sh_alrt_review_outputs(report: EngineWorkflowReport) -> bool:
    return (
        report.sh_alrt_support_summary is not None
        and report.output_paths.get("support_table") is not None
        and report.output_paths.get("conflicting_support_branches") is not None
    )
