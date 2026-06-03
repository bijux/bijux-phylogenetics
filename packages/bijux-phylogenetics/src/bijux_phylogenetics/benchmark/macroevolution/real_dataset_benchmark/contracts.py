from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.datasets.central_european_seashore_flora import (
    CentralEuropeanSeashoreFloraDataset,
    CentralEuropeanSeashoreFloraDatasetExportResult,
)


@dataclass(slots=True)
class RealDatasetMacroevolutionSummaryRow:
    """One benchmark-level conclusion row for a native or review surface."""

    surface_id: str
    trait: str
    trait_kind: str
    review_scope: str
    bijux_selected_model: str
    geiger_selected_model: str
    selection_matches_geiger: bool
    bijux_selected_model_akaike_weight: float | None
    geiger_selected_model_akaike_weight: float | None
    stable_conclusion_supported: bool
    aligned_taxa_count: int
    dropped_tree_taxon_count: int
    dropped_trait_taxon_count: int
    dropped_missing_value_taxon_count: int
    biological_interpretation: str
    notes: list[str]


@dataclass(slots=True)
class RealDatasetMacroevolutionModelRow:
    """One native-model table row comparing Bijux and local geiger on one model."""

    surface_id: str
    trait: str
    trait_kind: str
    model: str
    bijux_rank: int
    geiger_rank: int
    bijux_selected: bool
    geiger_selected: bool
    bijux_parameter_count: int
    geiger_parameter_count: int
    bijux_log_likelihood: float
    geiger_log_likelihood: float
    bijux_aic: float
    geiger_aic: float
    bijux_aicc: float
    geiger_aicc: float
    bijux_akaike_weight: float
    geiger_akaike_weight: float
    bijux_parameter_name: str | None
    geiger_parameter_name: str | None
    bijux_parameter_value: float | None
    geiger_parameter_value: float | None
    bijux_rate: float | None
    geiger_rate: float | None
    bijux_root_state: float | None
    geiger_root_state: float | None
    notes: list[str]


@dataclass(slots=True)
class RealDatasetMacroevolutionAlignmentReviewRow:
    """One missing-or-mismatched-taxon handling row over the real dataset review input."""

    surface_id: str
    trait: str
    model: str
    original_tree_taxa: int
    original_trait_taxa: int
    aligned_taxa_count: int
    dropped_tree_taxa: list[str]
    dropped_trait_taxa: list[str]
    dropped_missing_value_taxa: list[str]
    geiger_overlap_taxa: int
    geiger_usable_taxa: int
    notes: list[str]


@dataclass(slots=True)
class RealDatasetMacroevolutionParityRow:
    """One parity comparison row between Bijux and local geiger benchmark evidence."""

    surface_id: str
    trait: str
    model: str
    comparison_scope: str
    bijux_log_likelihood: float
    geiger_log_likelihood: float
    absolute_log_likelihood_delta: float
    bijux_aicc: float
    geiger_aicc: float
    absolute_aicc_delta: float
    bijux_parameter_name: str | None
    geiger_parameter_name: str | None
    bijux_parameter_value: float | None
    geiger_parameter_value: float | None
    absolute_parameter_delta: float | None
    within_log_likelihood_tolerance: bool
    within_aicc_tolerance: bool
    within_parameter_tolerance: bool | None
    notes: list[str]


@dataclass(slots=True)
class RealDatasetMacroevolutionBenchmarkReport:
    """Benchmark report over a real published dataset with governed geiger comparison."""

    dataset: CentralEuropeanSeashoreFloraDataset
    provenance_citation: str
    provenance_doi: str
    summary_rows: list[RealDatasetMacroevolutionSummaryRow]
    model_rows: list[RealDatasetMacroevolutionModelRow]
    alignment_review_rows: list[RealDatasetMacroevolutionAlignmentReviewRow]
    parity_rows: list[RealDatasetMacroevolutionParityRow]
    limitations: list[str]


@dataclass(slots=True)
class RealDatasetMacroevolutionBenchmarkBundle:
    """Written bundle for the real-dataset macroevolution benchmark."""

    output_root: Path
    review_traits_path: Path
    summary_path: Path
    model_table_path: Path
    alignment_review_path: Path
    parity_table_path: Path
    geiger_reference_path: Path


@dataclass(slots=True)
class RealDatasetMacroevolutionBenchmarkDemoResult:
    """Dataset export plus benchmark outputs for the release-style demo surface."""

    output_root: Path
    dataset: CentralEuropeanSeashoreFloraDataset
    dataset_export: CentralEuropeanSeashoreFloraDatasetExportResult
    benchmark_bundle: RealDatasetMacroevolutionBenchmarkBundle
    overview_path: Path
