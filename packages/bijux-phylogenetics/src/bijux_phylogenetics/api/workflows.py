from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.comparative.pgls import run_pgls
from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.engines.inference import (
    run_fasta_to_tree_workflow,
    run_phylo_workflow_config,
)
from bijux_phylogenetics.engines.workflows.alignment import (
    run_alignment_trimming,
    run_multiple_sequence_alignment,
)
from bijux_phylogenetics.engines.workflows.iqtree import (
    run_bootstrap_support_estimation,
    run_maximum_likelihood_tree_inference,
)
from bijux_phylogenetics.io.fasta.records import validate_fasta_input
from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet
from bijux_phylogenetics.reports.service import render_phylogenetics_report

from .workflow_results import (
    AlignmentWorkflowResult,
    AncestralReconstructionWorkflowResult,
    ComparativeModelWorkflowResult,
    ConfiguredPhyloWorkflowResult,
    FastaValidationResult,
    InferenceWorkflowResult,
    ReportWorkflowResult,
    SequenceToTreeWorkflowResult,
    SupportWorkflowResult,
    TreeComparisonWorkflowResult,
    TrimmingWorkflowResult,
)

__all__ = [
    "AlignmentWorkflowResult",
    "AncestralReconstructionWorkflowResult",
    "ComparativeModelWorkflowResult",
    "ConfiguredPhyloWorkflowResult",
    "FastaValidationResult",
    "InferenceWorkflowResult",
    "ReportWorkflowResult",
    "SequenceToTreeWorkflowResult",
    "SupportWorkflowResult",
    "TreeComparisonWorkflowResult",
    "TrimmingWorkflowResult",
    "render_report_workflow",
    "run_alignment_workflow",
    "run_ancestral_reconstruction_workflow",
    "run_comparative_model_workflow",
    "run_configured_phylo_workflow",
    "run_fasta_validation_workflow",
    "run_sequence_to_tree_workflow",
    "run_support_workflow",
    "run_tree_comparison_workflow",
    "run_tree_inference_workflow",
    "run_trimming_workflow",
]


def run_fasta_validation_workflow(
    input_path: Path,
) -> FastaValidationResult:
    """Run the CLI-grade FASTA validation surface from Python."""
    return FastaValidationResult(report=validate_fasta_input(input_path))


def run_alignment_workflow(
    input_path: Path,
    out_path: Path,
    *,
    executable: str | Path = "mafft",
    mode: str = "auto",
    extra_args: tuple[str, ...] = (),
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> AlignmentWorkflowResult:
    """Run the serious multiple-sequence-alignment workflow from Python."""
    return AlignmentWorkflowResult(
        report=run_multiple_sequence_alignment(
            input_path,
            out_path,
            executable=executable,
            mode=mode,
            extra_args=extra_args,
            resume=resume,
            timeout_seconds=timeout_seconds,
            incomplete_run_policy=incomplete_run_policy,
        )
    )


def run_trimming_workflow(
    input_path: Path,
    out_path: Path,
    *,
    executable: str | Path = "trimal",
    mode: str = "gap-threshold",
    gap_threshold: float = 0.1,
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> TrimmingWorkflowResult:
    """Run the serious alignment-trimming workflow from Python."""
    return TrimmingWorkflowResult(
        report=run_alignment_trimming(
            input_path,
            out_path,
            executable=executable,
            mode=mode,
            gap_threshold=gap_threshold,
            resume=resume,
            timeout_seconds=timeout_seconds,
            incomplete_run_policy=incomplete_run_policy,
        )
    )


def run_tree_inference_workflow(
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
) -> InferenceWorkflowResult:
    """Run the serious maximum-likelihood tree-inference workflow from Python."""
    return InferenceWorkflowResult(
        report=run_maximum_likelihood_tree_inference(
            input_path,
            out_dir=out_dir,
            model=model,
            prefix=prefix,
            executable=executable,
            sequence_type=sequence_type,
            partition_path=partition_path,
            resume=resume,
            seed=seed,
            threads=threads,
            timeout_seconds=timeout_seconds,
            incomplete_run_policy=incomplete_run_policy,
        )
    )


def run_support_workflow(
    input_path: Path,
    *,
    out_dir: Path,
    model: str,
    prefix: str = "bootstrap-support",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
    partition_path: Path | None = None,
    replicates: int = 1000,
    seed: int = 1,
    threads: int = 1,
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> SupportWorkflowResult:
    """Run the serious bootstrap-support workflow from Python."""
    return SupportWorkflowResult(
        report=run_bootstrap_support_estimation(
            input_path,
            out_dir=out_dir,
            model=model,
            prefix=prefix,
            executable=executable,
            sequence_type=sequence_type,
            partition_path=partition_path,
            replicates=replicates,
            seed=seed,
            threads=threads,
            resume=resume,
            timeout_seconds=timeout_seconds,
            incomplete_run_policy=incomplete_run_policy,
        )
    )


def run_sequence_to_tree_workflow(
    input_path: Path,
    *,
    out_dir: Path,
    prefix: str = "fasta-to-tree",
    mafft_executable: str | Path = "mafft",
    mafft_mode: str = "auto",
    trimal_executable: str | Path = "trimal",
    trimal_mode: str = "gap-threshold",
    trim_gap_threshold: float = 0.1,
    iqtree_executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
    bootstrap_replicates: int = 1000,
    seed: int = 1,
    threads: int = 1,
    normalize_identifiers: bool = False,
    remove_invalid_records: bool = False,
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> SequenceToTreeWorkflowResult:
    """Run the full owned FASTA-to-tree workflow from Python."""
    return SequenceToTreeWorkflowResult(
        report=run_fasta_to_tree_workflow(
            input_path,
            out_dir=out_dir,
            prefix=prefix,
            mafft_executable=mafft_executable,
            alignment_mode=mafft_mode,
            trimal_executable=trimal_executable,
            trimming_mode=trimal_mode,
            trim_gap_threshold=trim_gap_threshold,
            iqtree_executable=iqtree_executable,
            sequence_type=sequence_type,
            bootstrap_replicates=bootstrap_replicates,
            iqtree_seed=seed,
            iqtree_threads=threads,
            normalize_identifiers=normalize_identifiers,
            remove_invalid_records=remove_invalid_records,
            resume=resume,
            timeout_seconds=timeout_seconds,
            incomplete_run_policy=incomplete_run_policy,
        )
    )


def run_tree_comparison_workflow(
    left_path: Path,
    right_path: Path,
    *,
    rf_mode: str = "rooted",
    taxon_overlap_policy: str = "prune-to-shared",
) -> TreeComparisonWorkflowResult:
    """Run the serious topology-comparison workflow from Python."""
    return TreeComparisonWorkflowResult(
        report=compare_tree_paths(
            left_path,
            right_path,
            rf_mode=rf_mode,
            taxon_overlap_policy=taxon_overlap_policy,
        )
    )


def run_comparative_model_workflow(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativeModelWorkflowResult:
    """Run the serious PGLS comparative-model workflow from Python."""
    return ComparativeModelWorkflowResult(
        report=run_pgls(
            tree_path,
            traits_path,
            response=response,
            predictors=predictors,
            formula=formula,
            taxon_column=taxon_column,
            lambda_value=lambda_value,
        )
    )


def run_ancestral_reconstruction_workflow(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "fitch",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    root_prior_mode: str = "equal",
    fixed_root_state: str | None = None,
    allowed_transition_pairs: list[tuple[str, str]] | None = None,
) -> AncestralReconstructionWorkflowResult:
    """Run the serious discrete ancestral-reconstruction workflow from Python."""
    return AncestralReconstructionWorkflowResult(
        report=reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
            root_prior_mode=root_prior_mode,
            fixed_root_state=fixed_root_state,
            allowed_transition_pairs=allowed_transition_pairs,
        )
    )


def render_report_workflow(
    *,
    tree_path: Path,
    out_path: Path,
    alignment_path: Path | None = None,
    traits_path: Path | None = None,
    metadata_path: Path | None = None,
) -> ReportWorkflowResult:
    """Render the governed phylogenetics report surface from Python."""
    return ReportWorkflowResult(
        report=render_phylogenetics_report(
            tree_path=tree_path,
            out_path=out_path,
            alignment_path=alignment_path,
            traits_path=traits_path,
            metadata_path=metadata_path,
        )
    )


def run_configured_phylo_workflow(config_path: Path) -> ConfiguredPhyloWorkflowResult:
    """Run the one-command governed workflow config surface from Python."""
    return ConfiguredPhyloWorkflowResult(report=run_phylo_workflow_config(config_path))
