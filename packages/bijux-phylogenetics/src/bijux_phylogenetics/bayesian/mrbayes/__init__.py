from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines.common import (
    build_file_checksums,
    execute_engine_command,
    read_engine_version,
    resolve_engine_executable,
    validate_timeout_seconds,
)
from bijux_phylogenetics.engines.workflows.models import EngineWorkflowReport
from bijux_phylogenetics.engines.workflows.state import (
    _persist_workflow_report,
    _record_output_validation_failure,
    _resolve_incomplete_workflow_state,
    _resume_existing_workflow,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError, PhylogeneticsError
from .artifacts import _mrbayes_artifact_error
from .diagnostics import (
    assess_mrbayes_burnin_sensitivity,
    assess_mrbayes_convergence,
    compute_mrbayes_effective_sample_sizes,
    summarize_mrbayes_parameter_diagnostics,
    write_mrbayes_burnin_sensitivity_slice_table,
    write_mrbayes_parameter_summary_table,
)
from .models import (
    EffectiveSampleSize,
    MrBayesBurninSensitivityReport,
    MrBayesBurninSensitivitySlice,
    MrBayesConsensusTreeReport,
    MrBayesConvergenceReport,
    MrBayesESSReport,
    MrBayesMcmcReport,
    MrBayesMcmcRow,
    MrBayesParameterDiagnosticsReport,
    MrBayesParameterSummary,
    MrBayesPosteriorSummaryReport,
    MrBayesPosteriorTreeSample,
    MrBayesPosteriorTreeSetReport,
    MrBayesPreparationReport,
    MrBayesTraceReport,
    MrBayesTraceRow,
)
from .preparation import prepare_mrbayes_analysis
from .posterior_trees import (
    parse_mrbayes_consensus_tree,
    parse_mrbayes_posterior_tree_samples,
    summarize_mrbayes_posterior_trees,
)
from .tabular import (
    parse_mrbayes_mcmc_diagnostics,
    parse_mrbayes_parameter_traces,
)


def run_mrbayes_posterior_inference(
    nexus_path: Path,
    *,
    executable: str | Path = "mb",
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run a MrBayes posterior tree inference workflow from a prepared NEXUS file."""
    if not nexus_path.exists():
        raise _mrbayes_artifact_error(
            f"MrBayes analysis NEXUS file was not found: {nexus_path}",
            code="mrbayes_analysis_missing_file",
            path=nexus_path,
            artifact_kind="mrbayes-analysis-nexus",
            details={"expected_section": "analysis nexus file"},
        )
    validate_timeout_seconds(timeout_seconds)
    resolved = resolve_engine_executable(executable)
    prefix_path = nexus_path.with_suffix("")
    trace_path = Path(f"{nexus_path}.run1.p")
    tree_path = Path(f"{nexus_path}.run1.t")
    mcmc_path = Path(f"{nexus_path}.mcmc")
    consensus_path = Path(f"{nexus_path}.con.tre")
    manifest_path = prefix_path.with_suffix(".manifest.json")
    version = read_engine_version(
        "MrBayes",
        executable,
        version_args=("-v",),
        timeout_seconds=timeout_seconds,
    )
    command = [resolved, nexus_path.name]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=[nexus_path],
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
        engine_name="MrBayes",
        workflow="posterior-tree-inference",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=nexus_path.parent,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "posterior_trees": tree_path,
            "parameter_traces": trace_path,
            "mcmc_diagnostics": mcmc_path,
            "consensus_tree": consensus_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        parse_mrbayes_parameter_traces(trace_path)
        parse_mrbayes_mcmc_diagnostics(mcmc_path)
        parse_mrbayes_consensus_tree(consensus_path)
        summarize_mrbayes_posterior_trees(tree_path, burnin_fraction=0.25)
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    report = EngineWorkflowReport(
        workflow="posterior-tree-inference",
        engine_name="MrBayes",
        input_paths=[nexus_path],
        output_paths={
            "posterior_trees": tree_path,
            "parameter_traces": trace_path,
            "mcmc_diagnostics": mcmc_path,
            "consensus_tree": consensus_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([nexus_path]),
        output_checksums={},
        config={
            "timeout_seconds": timeout_seconds,
        },
        notes=[
            "MrBayes posterior trees, parameter traces, consensus tree, and MCMC diagnostics validated after engine execution",
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)
