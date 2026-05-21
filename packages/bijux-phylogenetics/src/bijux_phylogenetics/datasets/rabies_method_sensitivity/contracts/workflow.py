from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.compare.topology import TreeComparisonReport
from bijux_phylogenetics.engines.inference import (
    InferenceComparisonWorkflowReport,
)
from bijux_phylogenetics.engines.workflows.models import EngineWorkflowReport
from bijux_phylogenetics.phylo.topology import TreeRootingReport

from .dataset import (
    RabiesMethodSensitivityPanelDataset,
    RabiesMethodSensitivityVariant,
)


@dataclass(slots=True)
class RabiesMethodSensitivityVariantRun:
    """One executed variant in the governed method-sensitivity matrix."""

    config: RabiesMethodSensitivityVariant
    alignment_workflow: EngineWorkflowReport
    trimming_workflow: EngineWorkflowReport
    inference_comparison: InferenceComparisonWorkflowReport
    rooted_fasttree_path: Path
    rooted_iqtree_path: Path
    fasttree_rooting: TreeRootingReport
    iqtree_rooting: TreeRootingReport
    rooted_engine_comparison: TreeComparisonReport
    rooted_engine_comparison_table_path: Path
    alignment_length: int
    trimmed_alignment_length: int


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityTaskRecord:
    """One isolated variant execution record within the workflow batch."""

    variant_id: str
    label: str
    execution_mode: str
    status: str
    log_path: Path
    output_root: Path
    error_code: str | None
    error_message: str | None


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityPreprocessingComparisonRow:
    """One rooted IQ-TREE comparison across two preprocessing variants."""

    left_variant_id: str
    right_variant_id: str
    comparison_axis: str
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    same_taxa_different_rooting: bool


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityCladeRow:
    """One aggregated stable or changed clade-level conclusion across variants."""

    split_id: str
    conclusion_class: str
    evidence_class: str
    occurrence_count: int
    variant_count: int
    detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityConclusionRow:
    """One high-level biological or analytical conclusion from the workflow."""

    conclusion_id: str
    method_axis: str
    stability_status: str
    claim: str
    evidence: str
    caution: str


@dataclass(slots=True)
class RabiesMethodSensitivityPanelWorkflowReport:
    """One governed method-sensitivity workflow run over the packaged rabies panel."""

    dataset: RabiesMethodSensitivityPanelDataset
    execution_record_path: Path
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    parallel_workers: int
    execution_mode: str
    task_records: tuple[RabiesMethodSensitivityTaskRecord, ...]
    variant_runs: tuple[RabiesMethodSensitivityVariantRun, ...]
    preprocessing_comparison_rows: tuple[
        RabiesMethodSensitivityPreprocessingComparisonRow, ...
    ]
    stable_clade_rows: tuple[RabiesMethodSensitivityCladeRow, ...]
    changed_clade_rows: tuple[RabiesMethodSensitivityCladeRow, ...]
    conclusion_rows: tuple[RabiesMethodSensitivityConclusionRow, ...]
