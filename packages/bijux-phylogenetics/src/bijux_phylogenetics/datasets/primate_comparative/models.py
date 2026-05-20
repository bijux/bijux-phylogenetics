from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral import (
    ContinuousAncestralReport,
    DiscreteAncestralReport,
)
from bijux_phylogenetics.comparative import (
    BrownianTraitEvolutionSummaryReport,
    OUTraitEvolutionSummaryReport,
    PGLSResult,
    PhylogeneticSignalSummaryReport,
)


@dataclass(slots=True)
class PrimateComparativeDataset:
    """Packaged real mammal dataset for comparative workflow review."""

    dataset_id: str
    label: str
    dataset_root: Path
    tree_path: Path
    traits_path: Path
    reference_output_root: Path
    taxon_column: str
    taxon_count: int
    continuous_traits: tuple[str, ...]
    categorical_traits: tuple[str, ...]
    workflow_continuous_trait: str
    workflow_pgls_predictor: str
    workflow_discrete_trait: str
    source_locator: str


@dataclass(slots=True)
class PrimateComparativeDatasetExportResult:
    """Materialized copy of the packaged primate dataset."""

    output_root: Path
    readme_path: Path
    tree_path: Path
    traits_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class PrimateComparativeWorkflowReport:
    """One governed comparative workflow run over the packaged primate dataset."""

    dataset: PrimateComparativeDataset
    pgls: PGLSResult
    brownian: BrownianTraitEvolutionSummaryReport
    ou: OUTraitEvolutionSummaryReport
    signal: PhylogeneticSignalSummaryReport
    continuous_ancestral: ContinuousAncestralReport
    discrete_ancestral: DiscreteAncestralReport


@dataclass(slots=True)
class PrimateComparativeWorkflowBundle:
    """Written comparative workflow outputs for the packaged primate dataset."""

    output_root: Path
    summary_path: Path
    pgls_lambda_profile_path: Path
    brownian_summary_path: Path
    brownian_exclusion_path: Path
    ou_summary_path: Path
    ou_exclusion_path: Path
    signal_summary_path: Path
    signal_permutations_path: Path
    continuous_ancestral_summary_path: Path
    continuous_ancestral_uncertainty_path: Path
    continuous_ancestral_exclusion_path: Path
    discrete_ancestral_summary_path: Path
    discrete_ancestral_probability_path: Path
    discrete_ancestral_exclusion_path: Path


@dataclass(slots=True)
class PrimateComparativeDemoResult:
    """Dataset export plus workflow outputs for the public primate demo."""

    output_root: Path
    dataset: PrimateComparativeDataset
    dataset_export: PrimateComparativeDatasetExportResult
    workflow_bundle: PrimateComparativeWorkflowBundle
    overview_path: Path
