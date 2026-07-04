"""Stable Python workflow entrypoints for notebook and pipeline use."""

from bijux_phylogenetics.engines import (
    NATIVE_MAXIMUM_LIKELIHOOD_ONLY_ENV_VAR,
    is_native_maximum_likelihood_only_mode_enabled,
    native_maximum_likelihood_only_mode,
)

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
from .workflows import (
    render_report_workflow,
    run_alignment_workflow,
    run_ancestral_reconstruction_workflow,
    run_comparative_model_workflow,
    run_configured_phylo_workflow,
    run_fasta_validation_workflow,
    run_sequence_to_tree_workflow,
    run_support_workflow,
    run_tree_comparison_workflow,
    run_tree_inference_workflow,
    run_trimming_workflow,
)

__all__ = [
    "AlignmentWorkflowResult",
    "AncestralReconstructionWorkflowResult",
    "ComparativeModelWorkflowResult",
    "ConfiguredPhyloWorkflowResult",
    "FastaValidationResult",
    "InferenceWorkflowResult",
    "NATIVE_MAXIMUM_LIKELIHOOD_ONLY_ENV_VAR",
    "ReportWorkflowResult",
    "SequenceToTreeWorkflowResult",
    "SupportWorkflowResult",
    "TreeComparisonWorkflowResult",
    "TrimmingWorkflowResult",
    "is_native_maximum_likelihood_only_mode_enabled",
    "native_maximum_likelihood_only_mode",
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
