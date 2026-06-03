from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from bijux_phylogenetics.bayesian.posterior_execution import (
    run_bayesian_posterior_execution,
)

from .artifacts import _mrbayes_artifact_error
from .posterior_trees import (
    parse_mrbayes_consensus_tree,
    summarize_mrbayes_posterior_trees,
)
from .tabular import (
    parse_mrbayes_mcmc_diagnostics,
    parse_mrbayes_parameter_traces,
)

if TYPE_CHECKING:
    from bijux_phylogenetics.engines.workflows.models import EngineWorkflowReport


def run_mrbayes_posterior_inference(
    nexus_path: Path,
    *,
    executable: str | Path = "mb",
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run a MrBayes posterior tree inference workflow from a prepared NEXUS file."""
    from bijux_phylogenetics.engines.common import (
        read_engine_version,
        resolve_engine_executable,
        validate_timeout_seconds,
    )
    from bijux_phylogenetics.engines.validation.preflight import (
        require_external_engine_surface,
    )

    if not nexus_path.exists():
        raise _mrbayes_artifact_error(
            f"MrBayes analysis NEXUS file was not found: {nexus_path}",
            code="mrbayes_analysis_missing_file",
            path=nexus_path,
            artifact_kind="mrbayes-analysis-nexus",
            details={"expected_section": "analysis nexus file"},
        )
    validate_timeout_seconds(timeout_seconds)
    require_external_engine_surface(
        workflow_id="mrbayes-posterior",
        summary="MrBayes posterior inference workflow.",
        required_engines=("mrbayes",),
        executables={"mrbayes": executable},
        preserve_missing_error=True,
    )
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

    def _validate_consistent_posterior_taxa() -> None:
        _consensus_tree, consensus_report = parse_mrbayes_consensus_tree(consensus_path)
        _summary_tree, tree_report = summarize_mrbayes_posterior_trees(
            tree_path,
            burnin_fraction=0.25,
        )
        if sorted(consensus_report.tip_names) == sorted(tree_report.shared_taxa):
            return
        raise _mrbayes_artifact_error(
            f"MrBayes posterior outputs disagree on sampled taxa: {consensus_path}",
            code="mrbayes_outputs_inconsistent_taxa",
            path=consensus_path,
            artifact_kind="mrbayes-consensus-tree",
            details={
                "posterior_trees_path": str(tree_path),
                "consensus_tree_taxa": sorted(consensus_report.tip_names),
                "posterior_tree_taxa": sorted(tree_report.shared_taxa),
            },
        )

    def validate_outputs() -> None:
        parse_mrbayes_parameter_traces(trace_path)
        parse_mrbayes_mcmc_diagnostics(mcmc_path)
        _validate_consistent_posterior_taxa()

    return run_bayesian_posterior_execution(
        engine_name="MrBayes",
        executable=resolved,
        version=version,
        command=command,
        input_paths=[nexus_path],
        output_paths={
            "posterior_trees": tree_path,
            "parameter_traces": trace_path,
            "mcmc_diagnostics": mcmc_path,
            "consensus_tree": consensus_path,
        },
        manifest_path=manifest_path,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        work_dir=nexus_path.parent,
        timeout_seconds=timeout_seconds,
        config={
            "timeout_seconds": timeout_seconds,
        },
        notes=[
            "MrBayes posterior trees, parameter traces, consensus tree, and MCMC diagnostics validated after engine execution",
        ],
        resume=resume,
        incomplete_run_policy=incomplete_run_policy,
        validate_outputs=validate_outputs,
    )
