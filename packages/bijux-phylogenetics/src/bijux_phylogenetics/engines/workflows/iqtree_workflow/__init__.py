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
from .bootstrap_consensus import run_bootstrap_consensus_tree
from .bootstrap_support import run_bootstrap_support_estimation
from .maximum_likelihood import run_maximum_likelihood_tree_inference
from .model_selection import run_model_selection
from .partitions import _prepare_iqtree_partitions
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
