from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet
from bijux_phylogenetics.runtime.errors import EngineWorkflowError, PhylogeneticsError

from ...artifacts.sh_alrt import (
    build_conflicting_sh_alrt_support_rows,
    build_sh_alrt_support_rows,
    write_sh_alrt_support_table,
)
from ...common import (
    build_file_checksums,
    execute_engine_command,
    read_engine_version,
    resolve_engine_executable,
    validate_timeout_seconds,
)
from ...validation.audits import (
    summarize_bootstrap_support_distribution,
    summarize_sh_alrt_support_distribution,
)
from ...validation.preflight import require_external_engine_surface
from ..models import EngineWorkflowReport
from ..state import (
    _ensure_inference_ready_alignment,
    _persist_workflow_report,
    _prefix_path,
    _record_output_validation_failure,
    _resolve_incomplete_workflow_state,
    _resume_existing_workflow,
    _resume_has_sh_alrt_review_outputs,
    _validate_complete_support_coverage,
    _validate_iqtree_required_artifacts,
    _validate_matching_tree_taxa,
    _validate_support_value_count,
    _validate_tree_output,
    _validate_tree_set_output,
)
from .partitions import _prepare_iqtree_partitions
from .shared import (
    _build_iqtree_model_selection_summary,
    _build_iqtree_summary,
    _existing_iqtree_outputs,
    _iqtree_execution_controls,
    _iqtree_partition_supports_fixed_model,
    _iqtree_sequence_type_flag,
    _validate_iqtree_model_result,
    _validate_sh_alrt_replicates,
    _validate_ufboot_replicates,
)


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
    require_external_engine_surface(
        workflow_id="sh-alrt-support",
        summary="IQ-TREE SH-aLRT and UFBoot support workflow.",
        required_engines=("iqtree",),
        executables={"iqtree": executable},
        preserve_missing_error=True,
    )
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
        _validate_matching_tree_taxa(
            engine_name="iqtree",
            workflow="sh-alrt-support",
            reference_tree_path=support_tree_path,
            comparison_tree_set_path=bootstrap_tree_path,
            reference_output_name="support_tree",
            comparison_output_name="bootstrap_trees",
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
        _validate_complete_support_coverage(
            engine_name="iqtree",
            workflow="sh-alrt-support",
            path=support_tree_path,
            output_name="support_tree",
            artifact_kind="sh-alrt-supported-tree",
            annotated_branch_count=bootstrap_support_summary.supported_node_count,
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
        _validate_complete_support_coverage(
            engine_name="iqtree",
            workflow="sh-alrt-support",
            path=support_tree_path,
            output_name="support_tree",
            artifact_kind="sh-alrt-supported-tree",
            annotated_branch_count=sh_alrt_support_summary.annotated_node_count,
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
        _validate_complete_support_coverage(
            engine_name="iqtree",
            workflow="sh-alrt-support",
            path=support_tree_path,
            output_name="support_tree",
            artifact_kind="sh-alrt-supported-tree",
            annotated_branch_count=sh_alrt_support_summary.fully_scored_node_count,
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
