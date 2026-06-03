from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from ..geiger_reference import (
    GEIGER_REAL_DATASET_MACROEVOLUTION_REFERENCE_PAYLOADS,
)
from .contracts import RealDatasetMacroevolutionBenchmarkReport
from .shared import format_float, format_optional_float


def write_real_dataset_macroevolution_summary_table(
    path: Path,
    report: RealDatasetMacroevolutionBenchmarkReport,
) -> Path:
    """Write one summary row per native or review benchmark surface."""
    return write_taxon_rows(
        path,
        columns=[
            "surface_id",
            "trait",
            "trait_kind",
            "review_scope",
            "bijux_selected_model",
            "geiger_selected_model",
            "selection_matches_geiger",
            "bijux_selected_model_akaike_weight",
            "geiger_selected_model_akaike_weight",
            "stable_conclusion_supported",
            "aligned_taxa_count",
            "dropped_tree_taxon_count",
            "dropped_trait_taxon_count",
            "dropped_missing_value_taxon_count",
            "biological_interpretation",
            "notes",
        ],
        rows=[
            {
                "surface_id": row.surface_id,
                "trait": row.trait,
                "trait_kind": row.trait_kind,
                "review_scope": row.review_scope,
                "bijux_selected_model": row.bijux_selected_model,
                "geiger_selected_model": row.geiger_selected_model,
                "selection_matches_geiger": row.selection_matches_geiger,
                "bijux_selected_model_akaike_weight": format_optional_float(
                    row.bijux_selected_model_akaike_weight
                ),
                "geiger_selected_model_akaike_weight": format_optional_float(
                    row.geiger_selected_model_akaike_weight
                ),
                "stable_conclusion_supported": row.stable_conclusion_supported,
                "aligned_taxa_count": row.aligned_taxa_count,
                "dropped_tree_taxon_count": row.dropped_tree_taxon_count,
                "dropped_trait_taxon_count": row.dropped_trait_taxon_count,
                "dropped_missing_value_taxon_count": row.dropped_missing_value_taxon_count,
                "biological_interpretation": row.biological_interpretation,
                "notes": " | ".join(row.notes),
            }
            for row in report.summary_rows
        ],
    )


def write_real_dataset_macroevolution_model_table(
    path: Path,
    report: RealDatasetMacroevolutionBenchmarkReport,
) -> Path:
    """Write the native continuous and discrete model tables against geiger."""
    return write_taxon_rows(
        path,
        columns=[
            "surface_id",
            "trait",
            "trait_kind",
            "model",
            "bijux_rank",
            "geiger_rank",
            "bijux_selected",
            "geiger_selected",
            "bijux_parameter_count",
            "geiger_parameter_count",
            "bijux_log_likelihood",
            "geiger_log_likelihood",
            "bijux_aic",
            "geiger_aic",
            "bijux_aicc",
            "geiger_aicc",
            "bijux_akaike_weight",
            "geiger_akaike_weight",
            "bijux_parameter_name",
            "geiger_parameter_name",
            "bijux_parameter_value",
            "geiger_parameter_value",
            "bijux_rate",
            "geiger_rate",
            "bijux_root_state",
            "geiger_root_state",
            "notes",
        ],
        rows=[
            {
                "surface_id": row.surface_id,
                "trait": row.trait,
                "trait_kind": row.trait_kind,
                "model": row.model,
                "bijux_rank": row.bijux_rank,
                "geiger_rank": row.geiger_rank,
                "bijux_selected": row.bijux_selected,
                "geiger_selected": row.geiger_selected,
                "bijux_parameter_count": row.bijux_parameter_count,
                "geiger_parameter_count": row.geiger_parameter_count,
                "bijux_log_likelihood": format_float(row.bijux_log_likelihood),
                "geiger_log_likelihood": format_float(row.geiger_log_likelihood),
                "bijux_aic": format_float(row.bijux_aic),
                "geiger_aic": format_float(row.geiger_aic),
                "bijux_aicc": format_float(row.bijux_aicc),
                "geiger_aicc": format_float(row.geiger_aicc),
                "bijux_akaike_weight": format_float(row.bijux_akaike_weight),
                "geiger_akaike_weight": format_float(row.geiger_akaike_weight),
                "bijux_parameter_name": row.bijux_parameter_name or "",
                "geiger_parameter_name": row.geiger_parameter_name or "",
                "bijux_parameter_value": format_optional_float(
                    row.bijux_parameter_value
                ),
                "geiger_parameter_value": format_optional_float(
                    row.geiger_parameter_value
                ),
                "bijux_rate": format_optional_float(row.bijux_rate),
                "geiger_rate": format_optional_float(row.geiger_rate),
                "bijux_root_state": format_optional_float(row.bijux_root_state),
                "geiger_root_state": format_optional_float(row.geiger_root_state),
                "notes": " | ".join(row.notes),
            }
            for row in report.model_rows
        ],
    )


def write_real_dataset_macroevolution_alignment_review_table(
    path: Path,
    report: RealDatasetMacroevolutionBenchmarkReport,
) -> Path:
    """Write missing and mismatched taxon handling rows for the review input."""
    return write_taxon_rows(
        path,
        columns=[
            "surface_id",
            "trait",
            "model",
            "original_tree_taxa",
            "original_trait_taxa",
            "aligned_taxa_count",
            "dropped_tree_taxa",
            "dropped_trait_taxa",
            "dropped_missing_value_taxa",
            "geiger_overlap_taxa",
            "geiger_usable_taxa",
            "notes",
        ],
        rows=[
            {
                "surface_id": row.surface_id,
                "trait": row.trait,
                "model": row.model,
                "original_tree_taxa": row.original_tree_taxa,
                "original_trait_taxa": row.original_trait_taxa,
                "aligned_taxa_count": row.aligned_taxa_count,
                "dropped_tree_taxa": ",".join(row.dropped_tree_taxa),
                "dropped_trait_taxa": ",".join(row.dropped_trait_taxa),
                "dropped_missing_value_taxa": ",".join(row.dropped_missing_value_taxa),
                "geiger_overlap_taxa": row.geiger_overlap_taxa,
                "geiger_usable_taxa": row.geiger_usable_taxa,
                "notes": " | ".join(row.notes),
            }
            for row in report.alignment_review_rows
        ],
    )


def write_real_dataset_macroevolution_parity_table(
    path: Path,
    report: RealDatasetMacroevolutionBenchmarkReport,
) -> Path:
    """Write model-table and review-fit parity deltas against stored geiger evidence."""
    return write_taxon_rows(
        path,
        columns=[
            "surface_id",
            "trait",
            "model",
            "comparison_scope",
            "bijux_log_likelihood",
            "geiger_log_likelihood",
            "absolute_log_likelihood_delta",
            "bijux_aicc",
            "geiger_aicc",
            "absolute_aicc_delta",
            "bijux_parameter_name",
            "geiger_parameter_name",
            "bijux_parameter_value",
            "geiger_parameter_value",
            "absolute_parameter_delta",
            "within_log_likelihood_tolerance",
            "within_aicc_tolerance",
            "within_parameter_tolerance",
            "notes",
        ],
        rows=[
            {
                "surface_id": row.surface_id,
                "trait": row.trait,
                "model": row.model,
                "comparison_scope": row.comparison_scope,
                "bijux_log_likelihood": format_float(row.bijux_log_likelihood),
                "geiger_log_likelihood": format_float(row.geiger_log_likelihood),
                "absolute_log_likelihood_delta": format_float(
                    row.absolute_log_likelihood_delta
                ),
                "bijux_aicc": format_float(row.bijux_aicc),
                "geiger_aicc": format_float(row.geiger_aicc),
                "absolute_aicc_delta": format_float(row.absolute_aicc_delta),
                "bijux_parameter_name": row.bijux_parameter_name or "",
                "geiger_parameter_name": row.geiger_parameter_name or "",
                "bijux_parameter_value": format_optional_float(
                    row.bijux_parameter_value
                ),
                "geiger_parameter_value": format_optional_float(
                    row.geiger_parameter_value
                ),
                "absolute_parameter_delta": format_optional_float(
                    row.absolute_parameter_delta
                ),
                "within_log_likelihood_tolerance": row.within_log_likelihood_tolerance,
                "within_aicc_tolerance": row.within_aicc_tolerance,
                "within_parameter_tolerance": (
                    ""
                    if row.within_parameter_tolerance is None
                    else row.within_parameter_tolerance
                ),
                "notes": " | ".join(row.notes),
            }
            for row in report.parity_rows
        ],
    )


def write_geiger_real_dataset_reference_payload_table(
    path: Path,
    report: RealDatasetMacroevolutionBenchmarkReport,
) -> Path:
    """Write the stored local geiger payloads used by the real-dataset benchmark."""
    surface_ids = {row.surface_id for row in report.summary_rows}
    return write_taxon_rows(
        path,
        columns=["surface_id", "reference_payload_json"],
        rows=[
            {
                "surface_id": surface_id,
                "reference_payload_json": json.dumps(
                    GEIGER_REAL_DATASET_MACROEVOLUTION_REFERENCE_PAYLOADS[surface_id],
                    sort_keys=True,
                ),
            }
            for surface_id in sorted(surface_ids)
        ],
    )
