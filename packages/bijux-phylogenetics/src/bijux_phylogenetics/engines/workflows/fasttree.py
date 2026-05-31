from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.presentation import build_tree_comparison_report
from bijux_phylogenetics.io.fasta import (
    infer_alignment_alphabet,
    load_fasta_alignment,
)
from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from ..artifacts.fasttree import (
    build_fasttree_low_support_rows,
    build_fasttree_support_histogram_rows,
    build_fasttree_support_rows,
    write_fasttree_support_histogram,
    write_fasttree_support_table,
)
from ..common import (
    build_file_checksums,
    execute_engine_command,
    read_engine_version,
    resolve_engine_executable,
    validate_timeout_seconds,
)
from ..validation.audits import summarize_fasttree_support_distribution
from ..validation.preflight import require_external_engine_surface
from .models import EngineWorkflowReport, ExternalTreeComparisonReport
from .state import (
    _ensure_inference_ready_alignment,
    _manifest_path_from_output,
    _persist_workflow_report,
    _record_output_validation_failure,
    _resolve_incomplete_workflow_state,
    _resume_existing_workflow,
    _resume_has_fasttree_review_outputs,
    _sidecar,
    _validate_support_value_count,
    _validate_tree_output,
)


def _fasttree_args(path: Path, sequence_type: AlignmentAlphabet | None) -> list[str]:
    detected = sequence_type
    if detected is None:
        detected = infer_alignment_alphabet(load_fasta_alignment(path))
    if detected in {"dna", "rna"}:
        return ["-gtr", "-nt", str(path)]
    if detected == "protein":
        return ["-lg", str(path)]
    return [str(path)]


def run_fast_tree_inference(
    input_path: Path,
    out_path: Path,
    *,
    executable: str | Path = "FastTree",
    sequence_type: AlignmentAlphabet | None = None,
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run a fast approximate tree inference engine against an aligned FASTA file."""
    _ensure_inference_ready_alignment(input_path)
    validate_timeout_seconds(timeout_seconds)
    require_external_engine_surface(
        workflow_id="fast-approximate-tree",
        summary="FastTree approximate tree-inference workflow.",
        required_engines=("fasttree",),
        executables={"fasttree": executable},
        preserve_missing_error=True,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    support_table_path = _sidecar(out_path, "support.tsv")
    low_support_branches_path = _sidecar(out_path, "low-support.tsv")
    support_histogram_path = _sidecar(out_path, "support-histogram.tsv")
    version = read_engine_version(
        "FastTree",
        executable,
        version_args=("-help",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    manifest_path = _manifest_path_from_output(out_path)
    command = [resolved, *_fasttree_args(input_path.resolve(), sequence_type)]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=[input_path],
            expected_command=command,
            expected_version=version,
        )
        if resumed is not None and _resume_has_fasttree_review_outputs(resumed):
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )
    run = execute_engine_command(
        engine_name="FastTree",
        workflow="fast-approximate-tree",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_path.parent,
        stdout_path=out_path,
        stderr_path=_sidecar(out_path, "stderr.log"),
        output_paths={"tree": out_path},
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _validate_tree_output(
            out_path,
            engine_name="FastTree",
            workflow="fast-approximate-tree",
            output_name="tree",
            artifact_kind="fast-approximate-tree",
        )
        fasttree_support_summary = summarize_fasttree_support_distribution(out_path)
        _validate_support_value_count(
            engine_name="FastTree",
            workflow="fast-approximate-tree",
            path=out_path,
            output_name="tree",
            artifact_kind="fast-approximate-tree",
            support_value_count=fasttree_support_summary.annotated_node_count,
            support_kind="FastTree local support",
        )
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    write_fasttree_support_table(
        support_table_path,
        build_fasttree_support_rows(fasttree_support_summary),
    )
    write_fasttree_support_table(
        low_support_branches_path,
        build_fasttree_low_support_rows(fasttree_support_summary),
    )
    write_fasttree_support_histogram(
        support_histogram_path,
        build_fasttree_support_histogram_rows(fasttree_support_summary),
    )
    report = EngineWorkflowReport(
        workflow="fast-approximate-tree",
        engine_name="FastTree",
        input_paths=[input_path],
        output_paths={
            "tree": out_path,
            "support_table": support_table_path,
            "low_support_branches": low_support_branches_path,
            "support_histogram": support_histogram_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        config={
            "sequence_type": sequence_type,
            "timeout_seconds": timeout_seconds,
        },
        fasttree_support_summary=fasttree_support_summary,
        notes=[
            "fast approximate tree validated as parseable Newick output",
            "FastTree is an approximately maximum-likelihood engine and should be reviewed as an exploratory or large-alignment inference method",
            "FastTree local support values are SH-like support proportions on the documented 0-1 scale",
            "branch-level FastTree local support table exported for review",
            "low-support FastTree branch ledger exported for weak-clade review",
            "FastTree local support histogram exported for reviewer-facing distribution checks",
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)


def compare_fast_and_ml_trees(
    fast_tree_path: Path,
    ml_tree_path: Path,
    *,
    out_path: Path,
) -> ExternalTreeComparisonReport:
    """Compare a fast approximate tree against a maximum-likelihood tree."""
    comparison_report = build_tree_comparison_report(
        fast_tree_path, ml_tree_path, out_path=out_path
    )
    return ExternalTreeComparisonReport(
        fast_tree_path=fast_tree_path,
        ml_tree_path=ml_tree_path,
        comparison_report=comparison_report,
    )
