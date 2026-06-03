# ruff: noqa: F401
from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.presentation import build_tree_comparison_report
from bijux_phylogenetics.io.fasta import (
    infer_alignment_alphabet,
    load_fasta_alignment,
    write_fasta_alignment,
)
from bijux_phylogenetics.io.fasta.coding import (
    back_translate_aligned_coding_sequences,
    prepare_coding_sequences_for_alignment,
    translate_prepared_coding_sequences,
)
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet,
    AlignmentRecord,
    AlignmentSummary,
    CodingSequenceExclusion,
)
from bijux_phylogenetics.phylo.alignment.partitions import LocusPartition
from bijux_phylogenetics.runtime.errors import EngineWorkflowError, PhylogeneticsError

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
    parse_best_model_file,
    parse_iqtree_model_selection_summary,
    parse_log_likelihood_file,
    resolve_iqtree_model_sidecar,
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
    build_engine_output_error,
    build_file_checksums,
    clear_incomplete_engine_run,
    execute_engine_command,
    load_engine_manifest,
    load_unaligned_fasta,
    read_engine_version,
    resolve_engine_executable,
    validate_timeout_seconds,
    write_engine_manifest,
)
from ...validation import (
    detect_weakly_supported_backbone,
    summarize_bootstrap_support_distribution,
    summarize_fasttree_support_distribution,
    summarize_sh_alrt_support_distribution,
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
from .artifact_validation import (
    _ensure_inference_ready_alignment,
    _require_nonempty_text_output,
    _validate_alignment_output,
    _validate_complete_support_coverage,
    _validate_iqtree_required_artifacts,
    _validate_matching_tree_taxa,
    _validate_support_value_count,
    _validate_tree_output,
    _validate_tree_set_output,
)
from .coding_alignment import (
    _build_alignment_trimming_summary,
    _write_alignment_trimming_summary_table,
    _write_coding_exclusion_table,
    _write_coding_summary_table,
)
from .incomplete_runs import (
    _record_output_validation_failure,
    _resolve_incomplete_workflow_state,
    _validate_incomplete_run_policy,
)
from .paths import (
    _manifest_path_from_output,
    _partition_alignment_file_name,
    _partition_support_path,
    _prefix_path,
    _sidecar,
)
from .report_restore import (
    _restore_codon_aware_alignment_report,
    _restore_workflow_report,
)
from .resume_runtime import (
    _persist_workflow_report,
    _resume_existing_codon_aware_alignment,
    _resume_existing_workflow,
    _resume_has_bootstrap_review_outputs,
    _resume_has_fasttree_review_outputs,
    _resume_has_sh_alrt_review_outputs,
)
