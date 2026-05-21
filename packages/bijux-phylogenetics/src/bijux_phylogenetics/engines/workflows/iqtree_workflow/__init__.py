# ruff: noqa: F401
from __future__ import annotations

import hashlib
from pathlib import Path
import re

from bijux_phylogenetics.compare.presentation import build_tree_comparison_report
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.fasta import write_fasta_alignment
from bijux_phylogenetics.io.fasta.coding import (
    back_translate_aligned_coding_sequences,
    classify_sequence_coding_behavior,
    prepare_coding_sequences_for_alignment,
    translate_prepared_coding_sequences,
)
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet,
    AlignmentSummary,
    CodingSequenceExclusion,
)
from bijux_phylogenetics.phylo.alignment.partitions import LocusPartition
from bijux_phylogenetics.runtime.errors import EngineWorkflowError, PhylogeneticsError
from bijux_phylogenetics.trees import load_tree_set

from ...artifacts.bootstrap import (
    build_bootstrap_support_histogram_rows,
    build_bootstrap_support_rows,
    build_low_support_bootstrap_rows,
    write_bootstrap_support_histogram,
    write_bootstrap_support_table,
)
from ...artifacts.fasttree import (
    build_fasttree_low_support_rows,
    build_fasttree_support_histogram_rows,
    build_fasttree_support_rows,
    write_fasttree_support_histogram,
    write_fasttree_support_table,
)
from ...artifacts.iqtree import (
    IqtreeModelCandidate,
    IqtreeModelSelectionSummary,
    write_iqtree_model_candidates_table,
)
from ...artifacts.sh_alrt import (
    build_conflicting_sh_alrt_support_rows,
    build_sh_alrt_support_rows,
    write_sh_alrt_support_table,
)
from ...artifacts.support import (
    BootstrapSupportNode,
    BootstrapSupportSummaryReport,
    FastTreeSupportNode,
    FastTreeSupportSummaryReport,
    ShAlrtSupportNode,
    ShAlrtSupportSummaryReport,
    WeakBackboneReport,
)
from ...common import (
    EngineRunReport,
    EngineVersionInfo,
    active_engine_run_is_live,
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
from ...validation import (
    detect_weakly_supported_backbone,
    summarize_bootstrap_support_distribution,
    summarize_fasttree_support_distribution,
    summarize_sh_alrt_support_distribution,
)
from ..alignment import (
    list_mafft_alignment_modes,
    list_trimal_trimming_modes,
    resolve_mafft_alignment_mode,
    resolve_trimal_trimming_mode,
    run_alignment_trimming,
    run_codon_aware_multiple_sequence_alignment,
    run_multiple_sequence_alignment,
)
from ..models import (
    AlignmentTrimmingSummary,
    CodonAwareAlignmentWorkflowReport,
    EngineWorkflowReport,
    ExternalTreeComparisonReport,
    IqtreeSupportValue,
    IqtreeWorkflowSummary,
)
from ..models import (
    PreparedIqtreePartitions as _PreparedIqtreePartitions,
)
from ..state import (
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
from .partitions import _prepare_iqtree_partitions
from .bootstrap_support import run_bootstrap_support_estimation
from .maximum_likelihood import run_maximum_likelihood_tree_inference
from .model_selection import run_model_selection
from .sh_alrt_support import run_sh_alrt_support_estimation
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
