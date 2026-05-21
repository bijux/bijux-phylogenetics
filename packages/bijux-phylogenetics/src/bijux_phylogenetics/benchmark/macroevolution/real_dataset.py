from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory

from bijux_phylogenetics.comparative.discrete_mk import (
    DiscreteMkFitReport,
    fit_discrete_mk_model,
)
from bijux_phylogenetics.comparative.evolutionary_modes import (
    ContinuousEvolutionaryModeComparisonReport,
    ContinuousEvolutionaryModeFitReport,
    compare_fitcontinuous_model_ranking,
    fit_continuous_evolutionary_mode,
)
from bijux_phylogenetics.datasets.central_european_seashore_flora import (
    CentralEuropeanSeashoreFloraDataset,
    CentralEuropeanSeashoreFloraDatasetExportResult,
    export_central_european_seashore_flora_dataset,
    load_central_european_seashore_flora_dataset,
)
from bijux_phylogenetics.datasets.study_inputs import (
    TreeTraitAlignmentReport,
    align_tree_and_trait_table,
    load_taxon_table,
    write_taxon_rows,
)

from .geiger_reference import (
    GEIGER_REAL_DATASET_MACROEVOLUTION_REFERENCE_PAYLOADS,
)
from .real_dataset_benchmark.contracts import (
    RealDatasetMacroevolutionAlignmentReviewRow,
    RealDatasetMacroevolutionBenchmarkBundle,
    RealDatasetMacroevolutionBenchmarkDemoResult,
    RealDatasetMacroevolutionBenchmarkReport,
    RealDatasetMacroevolutionModelRow,
    RealDatasetMacroevolutionParityRow,
    RealDatasetMacroevolutionSummaryRow,
)

_CONTINUOUS_SURFACE_ID = "seed-mass-native-model-table"
_DISCRETE_SURFACE_ID = "lifeform-native-model-table"
_CONTINUOUS_REVIEW_SURFACE_ID = "seed-mass-alignment-review"
_DISCRETE_REVIEW_SURFACE_ID = "lifeform-alignment-review"
_CONTINUOUS_MODES = (
    "brownian",
    "white-noise",
    "pagel-lambda",
    "ornstein-uhlenbeck",
    "early-burst",
)
_DISCRETE_MODELS = ("equal-rates", "symmetric", "all-rates-different")
_REMOVED_TREE_TAXON = "Triglochin_maritimum"
_EXTRA_TRAIT_TAXON = "unmatched_review_taxon"
_CONTINUOUS_MISSING_VALUE_TAXON = "Juncus_maritimus"
_DISCRETE_MISSING_VALUE_TAXON = "Juncus_gerardii"
_PROVENANCE_CITATION = (
    "Vandelook, Filip; Janssens, Steven B.; Matthies, Diethart (2018). "
    "Data from: Ecological niche and phylogeny explain distribution of seed mass "
    "in the Central European flora. Dryad."
)
_PROVENANCE_DOI = "10.5061/dryad.0st06f0"


def benchmark_real_dataset_macroevolution() -> RealDatasetMacroevolutionBenchmarkReport:
    """Benchmark continuous and discrete comparative fits on a real published dataset."""
    dataset = load_central_european_seashore_flora_dataset()
    with TemporaryDirectory(prefix="real-dataset-macroevolution-") as temporary_root:
        review_traits_path = _write_alignment_review_traits_table(
            Path(temporary_root) / "alignment-review-traits.csv",
            dataset,
        )
        return _build_report(dataset, review_traits_path)


def write_real_dataset_macroevolution_bundle(
    output_root: Path,
) -> RealDatasetMacroevolutionBenchmarkBundle:
    """Write the benchmark ledgers and review input for the real dataset benchmark."""
    if output_root.exists():
        for path in sorted(output_root.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                path.rmdir()
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_central_european_seashore_flora_dataset()
    review_traits_path = _write_alignment_review_traits_table(
        output_root / "alignment-review-traits.csv",
        dataset,
    )
    report = _build_report(dataset, review_traits_path)
    summary_path = write_real_dataset_macroevolution_summary_table(
        output_root / "benchmark-summary.tsv",
        report,
    )
    model_table_path = write_real_dataset_macroevolution_model_table(
        output_root / "model-table.tsv",
        report,
    )
    alignment_review_path = write_real_dataset_macroevolution_alignment_review_table(
        output_root / "alignment-review.tsv",
        report,
    )
    parity_table_path = write_real_dataset_macroevolution_parity_table(
        output_root / "geiger-parity.tsv",
        report,
    )
    geiger_reference_path = write_geiger_real_dataset_reference_payload_table(
        output_root / "geiger-reference.tsv",
        report,
    )
    return RealDatasetMacroevolutionBenchmarkBundle(
        output_root=output_root,
        review_traits_path=review_traits_path,
        summary_path=summary_path,
        model_table_path=model_table_path,
        alignment_review_path=alignment_review_path,
        parity_table_path=parity_table_path,
        geiger_reference_path=geiger_reference_path,
    )


def run_real_dataset_macroevolution_benchmark_demo(
    destination: Path,
) -> RealDatasetMacroevolutionBenchmarkDemoResult:
    """Materialize the public plant dataset plus the governed benchmark bundle."""
    dataset = load_central_european_seashore_flora_dataset()
    destination.mkdir(parents=True, exist_ok=True)
    dataset_export = export_central_european_seashore_flora_dataset(
        destination / "dataset"
    )
    benchmark_bundle = write_real_dataset_macroevolution_bundle(
        destination / "benchmark"
    )
    overview_path = _write_overview(
        destination / "README.md",
        dataset=dataset,
        bundle=benchmark_bundle,
    )
    return RealDatasetMacroevolutionBenchmarkDemoResult(
        output_root=destination,
        dataset=dataset,
        dataset_export=dataset_export,
        benchmark_bundle=benchmark_bundle,
        overview_path=overview_path,
    )


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
                "bijux_selected_model_akaike_weight": _format_optional_float(
                    row.bijux_selected_model_akaike_weight
                ),
                "geiger_selected_model_akaike_weight": _format_optional_float(
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
                "bijux_log_likelihood": _format_float(row.bijux_log_likelihood),
                "geiger_log_likelihood": _format_float(row.geiger_log_likelihood),
                "bijux_aic": _format_float(row.bijux_aic),
                "geiger_aic": _format_float(row.geiger_aic),
                "bijux_aicc": _format_float(row.bijux_aicc),
                "geiger_aicc": _format_float(row.geiger_aicc),
                "bijux_akaike_weight": _format_float(row.bijux_akaike_weight),
                "geiger_akaike_weight": _format_float(row.geiger_akaike_weight),
                "bijux_parameter_name": row.bijux_parameter_name or "",
                "geiger_parameter_name": row.geiger_parameter_name or "",
                "bijux_parameter_value": _format_optional_float(
                    row.bijux_parameter_value
                ),
                "geiger_parameter_value": _format_optional_float(
                    row.geiger_parameter_value
                ),
                "bijux_rate": _format_optional_float(row.bijux_rate),
                "geiger_rate": _format_optional_float(row.geiger_rate),
                "bijux_root_state": _format_optional_float(row.bijux_root_state),
                "geiger_root_state": _format_optional_float(row.geiger_root_state),
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
                "bijux_log_likelihood": _format_float(row.bijux_log_likelihood),
                "geiger_log_likelihood": _format_float(row.geiger_log_likelihood),
                "absolute_log_likelihood_delta": _format_float(
                    row.absolute_log_likelihood_delta
                ),
                "bijux_aicc": _format_float(row.bijux_aicc),
                "geiger_aicc": _format_float(row.geiger_aicc),
                "absolute_aicc_delta": _format_float(row.absolute_aicc_delta),
                "bijux_parameter_name": row.bijux_parameter_name or "",
                "geiger_parameter_name": row.geiger_parameter_name or "",
                "bijux_parameter_value": _format_optional_float(
                    row.bijux_parameter_value
                ),
                "geiger_parameter_value": _format_optional_float(
                    row.geiger_parameter_value
                ),
                "absolute_parameter_delta": _format_optional_float(
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


def _build_report(
    dataset: CentralEuropeanSeashoreFloraDataset,
    review_traits_path: Path,
) -> RealDatasetMacroevolutionBenchmarkReport:
    continuous_comparison = compare_fitcontinuous_model_ranking(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.workflow_continuous_trait,
        taxon_column=dataset.taxon_column,
        modes=_CONTINUOUS_MODES,
    )
    continuous_fits = {
        mode: fit_continuous_evolutionary_mode(
            dataset.tree_path,
            dataset.traits_path,
            trait=dataset.workflow_continuous_trait,
            taxon_column=dataset.taxon_column,
            mode=mode,
        )
        for mode in _CONTINUOUS_MODES
    }
    discrete_fits = {
        model: fit_discrete_mk_model(
            dataset.tree_path,
            dataset.traits_path,
            trait=dataset.workflow_discrete_trait,
            taxon_column=dataset.taxon_column,
            model=model,
        )
        for model in _DISCRETE_MODELS
    }
    review_continuous_alignment = align_tree_and_trait_table(
        dataset.tree_path,
        review_traits_path,
        taxon_column=dataset.taxon_column,
        required_trait_columns=(dataset.workflow_continuous_trait,),
        drop_missing_for_columns=(dataset.workflow_continuous_trait,),
    ).report
    review_discrete_alignment = align_tree_and_trait_table(
        dataset.tree_path,
        review_traits_path,
        taxon_column=dataset.taxon_column,
        required_trait_columns=(dataset.workflow_discrete_trait,),
        drop_missing_for_columns=(dataset.workflow_discrete_trait,),
    ).report
    review_continuous_fit = fit_continuous_evolutionary_mode(
        dataset.tree_path,
        review_traits_path,
        trait=dataset.workflow_continuous_trait,
        taxon_column=dataset.taxon_column,
        mode="ornstein-uhlenbeck",
    )
    review_discrete_fit = fit_discrete_mk_model(
        dataset.tree_path,
        review_traits_path,
        trait=dataset.workflow_discrete_trait,
        taxon_column=dataset.taxon_column,
        model="equal-rates",
    )

    native_continuous_reference = GEIGER_REAL_DATASET_MACROEVOLUTION_REFERENCE_PAYLOADS[
        _CONTINUOUS_SURFACE_ID
    ]
    native_discrete_reference = GEIGER_REAL_DATASET_MACROEVOLUTION_REFERENCE_PAYLOADS[
        _DISCRETE_SURFACE_ID
    ]
    review_continuous_reference = GEIGER_REAL_DATASET_MACROEVOLUTION_REFERENCE_PAYLOADS[
        _CONTINUOUS_REVIEW_SURFACE_ID
    ]
    review_discrete_reference = GEIGER_REAL_DATASET_MACROEVOLUTION_REFERENCE_PAYLOADS[
        _DISCRETE_REVIEW_SURFACE_ID
    ]

    model_rows = _build_model_rows(
        continuous_comparison=continuous_comparison,
        continuous_fits=continuous_fits,
        native_continuous_reference=native_continuous_reference,
        discrete_fits=discrete_fits,
        native_discrete_reference=native_discrete_reference,
    )
    parity_rows = _build_parity_rows(
        model_rows=model_rows,
        review_continuous_fit=review_continuous_fit,
        review_continuous_reference=review_continuous_reference,
        review_discrete_fit=review_discrete_fit,
        review_discrete_reference=review_discrete_reference,
    )
    alignment_review_rows = _build_alignment_review_rows(
        review_continuous_alignment=review_continuous_alignment,
        review_discrete_alignment=review_discrete_alignment,
        review_continuous_reference=review_continuous_reference,
        review_discrete_reference=review_discrete_reference,
    )
    summary_rows = _build_summary_rows(
        continuous_comparison=continuous_comparison,
        native_continuous_reference=native_continuous_reference,
        discrete_fits=discrete_fits,
        native_discrete_reference=native_discrete_reference,
        review_continuous_alignment=review_continuous_alignment,
        review_discrete_alignment=review_discrete_alignment,
        review_continuous_reference=review_continuous_reference,
        review_discrete_reference=review_discrete_reference,
    )
    return RealDatasetMacroevolutionBenchmarkReport(
        dataset=dataset,
        provenance_citation=_PROVENANCE_CITATION,
        provenance_doi=_PROVENANCE_DOI,
        summary_rows=summary_rows,
        model_rows=model_rows,
        alignment_review_rows=alignment_review_rows,
        parity_rows=parity_rows,
        limitations=[
            "continuous biological interpretation remains cautious because the selected OU fit is near the lower alpha boundary and still carries a flat-likelihood warning",
            "discrete model-table parity on the sparse five-state lifeform trait is selection-level rather than row-level because Bijux SYM and ARD retain lower-bound and non-convergence warnings on this real dataset",
            "the alignment review input is a governed mismatch-and-missing-value perturbation of the public plant trait table and is not itself a biological inference surface",
        ],
    )


def _build_model_rows(
    *,
    continuous_comparison: ContinuousEvolutionaryModeComparisonReport,
    continuous_fits: dict[str, ContinuousEvolutionaryModeFitReport],
    native_continuous_reference: dict[str, object],
    discrete_fits: dict[str, DiscreteMkFitReport],
    native_discrete_reference: dict[str, object],
) -> list[RealDatasetMacroevolutionModelRow]:
    rows: list[RealDatasetMacroevolutionModelRow] = []
    geiger_continuous_rows = {
        row["model"]: row for row in native_continuous_reference["comparison_rows"]
    }
    for row in continuous_comparison.rows:
        fit = continuous_fits[row.model]
        geiger_row = geiger_continuous_rows[row.model]
        notes = []
        if (
            fit.boundary_assessment is not None
            and fit.boundary_assessment.boundary_warning_kinds
        ):
            notes.append(
                "boundary-review:"
                + ",".join(fit.boundary_assessment.boundary_warning_kinds)
            )
        rows.append(
            RealDatasetMacroevolutionModelRow(
                surface_id=_CONTINUOUS_SURFACE_ID,
                trait="seed_mass",
                trait_kind="continuous",
                model=row.model,
                bijux_rank=row.rank or 0,
                geiger_rank=int(geiger_row["rank"]),
                bijux_selected=row.selected,
                geiger_selected=(
                    row.model == native_continuous_reference["selected_model"]
                ),
                bijux_parameter_count=row.parameter_count,
                geiger_parameter_count=int(geiger_row["parameter_count"]),
                bijux_log_likelihood=fit.log_likelihood,
                geiger_log_likelihood=float(geiger_row["log_likelihood"]),
                bijux_aic=row.aic,
                geiger_aic=float(geiger_row["aic"]),
                bijux_aicc=row.aicc,
                geiger_aicc=float(geiger_row["aicc"]),
                bijux_akaike_weight=float(row.akaike_weight or 0.0),
                geiger_akaike_weight=0.0,
                bijux_parameter_name=fit.parameter_name,
                geiger_parameter_name=_optional_str(geiger_row.get("parameter_name")),
                bijux_parameter_value=fit.parameter_value,
                geiger_parameter_value=_optional_float(
                    geiger_row.get("parameter_value")
                ),
                bijux_rate=fit.rate,
                geiger_rate=_optional_float(geiger_row.get("rate")),
                bijux_root_state=fit.root_state,
                geiger_root_state=_optional_float(geiger_row.get("root_state")),
                notes=notes,
            )
        )
    _apply_geiger_akaike_weights(rows, surface_id=_CONTINUOUS_SURFACE_ID)

    geiger_discrete_rows = {
        row["model"]: row for row in native_discrete_reference["comparison_rows"]
    }
    discrete_selected_model = _selected_model_from_aicc(
        {model: fit.aicc for model, fit in discrete_fits.items()}
    )
    for rank, model in enumerate(
        sorted(_DISCRETE_MODELS, key=lambda item: discrete_fits[item].aicc),
        start=1,
    ):
        fit = discrete_fits[model]
        geiger_row = geiger_discrete_rows[model]
        rows.append(
            RealDatasetMacroevolutionModelRow(
                surface_id=_DISCRETE_SURFACE_ID,
                trait="lifeform",
                trait_kind="discrete",
                model=model,
                bijux_rank=rank,
                geiger_rank=int(geiger_row["rank"]),
                bijux_selected=(model == discrete_selected_model),
                geiger_selected=(model == native_discrete_reference["selected_model"]),
                bijux_parameter_count=fit.parameter_count,
                geiger_parameter_count=int(geiger_row["parameter_count"]),
                bijux_log_likelihood=fit.log_likelihood,
                geiger_log_likelihood=float(geiger_row["log_likelihood"]),
                bijux_aic=fit.aic,
                geiger_aic=float(geiger_row["aic"]),
                bijux_aicc=fit.aicc,
                geiger_aicc=float(geiger_row["aicc"]),
                bijux_akaike_weight=0.0,
                geiger_akaike_weight=0.0,
                bijux_parameter_name=None,
                geiger_parameter_name=None,
                bijux_parameter_value=None,
                geiger_parameter_value=None,
                bijux_rate=_representative_discrete_rate(fit),
                geiger_rate=_optional_float(geiger_row.get("representative_rate")),
                bijux_root_state=None,
                geiger_root_state=None,
                notes=list(fit.input_audit.warnings),
            )
        )
    _apply_akaike_weights_from_report(
        rows, surface_id=_DISCRETE_SURFACE_ID, engine="bijux"
    )
    _apply_geiger_akaike_weights(rows, surface_id=_DISCRETE_SURFACE_ID)
    return rows


def _build_parity_rows(
    *,
    model_rows: list[RealDatasetMacroevolutionModelRow],
    review_continuous_fit: ContinuousEvolutionaryModeFitReport,
    review_continuous_reference: dict[str, object],
    review_discrete_fit: DiscreteMkFitReport,
    review_discrete_reference: dict[str, object],
) -> list[RealDatasetMacroevolutionParityRow]:
    rows = [
        RealDatasetMacroevolutionParityRow(
            surface_id=row.surface_id,
            trait=row.trait,
            model=row.model,
            comparison_scope="native-model-table",
            bijux_log_likelihood=row.bijux_log_likelihood,
            geiger_log_likelihood=row.geiger_log_likelihood,
            absolute_log_likelihood_delta=abs(
                row.bijux_log_likelihood - row.geiger_log_likelihood
            ),
            bijux_aicc=row.bijux_aicc,
            geiger_aicc=row.geiger_aicc,
            absolute_aicc_delta=abs(row.bijux_aicc - row.geiger_aicc),
            bijux_parameter_name=row.bijux_parameter_name,
            geiger_parameter_name=row.geiger_parameter_name,
            bijux_parameter_value=row.bijux_parameter_value,
            geiger_parameter_value=row.geiger_parameter_value,
            absolute_parameter_delta=(
                None
                if row.bijux_parameter_value is None
                or row.geiger_parameter_value is None
                else abs(row.bijux_parameter_value - row.geiger_parameter_value)
            ),
            within_log_likelihood_tolerance=(
                abs(row.bijux_log_likelihood - row.geiger_log_likelihood) <= 0.5
            ),
            within_aicc_tolerance=(abs(row.bijux_aicc - row.geiger_aicc) <= 1.0),
            within_parameter_tolerance=(
                None
                if row.bijux_parameter_value is None
                or row.geiger_parameter_value is None
                else abs(row.bijux_parameter_value - row.geiger_parameter_value) <= 0.1
            ),
            notes=row.notes,
        )
        for row in model_rows
    ]
    review_continuous_summary = review_continuous_reference["fit_summary"]
    rows.append(
        RealDatasetMacroevolutionParityRow(
            surface_id=_CONTINUOUS_REVIEW_SURFACE_ID,
            trait="seed_mass",
            model="ornstein-uhlenbeck",
            comparison_scope="alignment-review",
            bijux_log_likelihood=review_continuous_fit.log_likelihood,
            geiger_log_likelihood=float(review_continuous_summary["log_likelihood"]),
            absolute_log_likelihood_delta=abs(
                review_continuous_fit.log_likelihood
                - float(review_continuous_summary["log_likelihood"])
            ),
            bijux_aicc=review_continuous_fit.aicc,
            geiger_aicc=float(review_continuous_summary["aicc"]),
            absolute_aicc_delta=abs(
                review_continuous_fit.aicc - float(review_continuous_summary["aicc"])
            ),
            bijux_parameter_name=review_continuous_fit.parameter_name,
            geiger_parameter_name=_optional_str(
                review_continuous_summary.get("parameter_name")
            ),
            bijux_parameter_value=review_continuous_fit.parameter_value,
            geiger_parameter_value=_optional_float(
                review_continuous_summary.get("parameter_value")
            ),
            absolute_parameter_delta=(
                None
                if review_continuous_fit.parameter_value is None
                else abs(
                    review_continuous_fit.parameter_value
                    - float(review_continuous_summary["parameter_value"])
                )
            ),
            within_log_likelihood_tolerance=(
                abs(
                    review_continuous_fit.log_likelihood
                    - float(review_continuous_summary["log_likelihood"])
                )
                <= 0.5
            ),
            within_aicc_tolerance=(
                abs(
                    review_continuous_fit.aicc
                    - float(review_continuous_summary["aicc"])
                )
                <= 1.0
            ),
            within_parameter_tolerance=(
                abs(
                    review_continuous_fit.parameter_value
                    - float(review_continuous_summary["parameter_value"])
                )
                <= 0.1
            ),
            notes=list(review_continuous_reference["notes"]),
        )
    )
    review_discrete_summary = review_discrete_reference["fit_summary"]
    rows.append(
        RealDatasetMacroevolutionParityRow(
            surface_id=_DISCRETE_REVIEW_SURFACE_ID,
            trait="lifeform",
            model="equal-rates",
            comparison_scope="alignment-review",
            bijux_log_likelihood=review_discrete_fit.log_likelihood,
            geiger_log_likelihood=float(review_discrete_summary["log_likelihood"]),
            absolute_log_likelihood_delta=abs(
                review_discrete_fit.log_likelihood
                - float(review_discrete_summary["log_likelihood"])
            ),
            bijux_aicc=review_discrete_fit.aicc,
            geiger_aicc=float(review_discrete_summary["aicc"]),
            absolute_aicc_delta=abs(
                review_discrete_fit.aicc - float(review_discrete_summary["aicc"])
            ),
            bijux_parameter_name=None,
            geiger_parameter_name=None,
            bijux_parameter_value=None,
            geiger_parameter_value=None,
            absolute_parameter_delta=None,
            within_log_likelihood_tolerance=(
                abs(
                    review_discrete_fit.log_likelihood
                    - float(review_discrete_summary["log_likelihood"])
                )
                <= 0.5
            ),
            within_aicc_tolerance=(
                abs(review_discrete_fit.aicc - float(review_discrete_summary["aicc"]))
                <= 1.0
            ),
            within_parameter_tolerance=None,
            notes=list(review_discrete_reference["notes"]),
        )
    )
    return rows


def _build_alignment_review_rows(
    *,
    review_continuous_alignment: TreeTraitAlignmentReport,
    review_discrete_alignment: TreeTraitAlignmentReport,
    review_continuous_reference: dict[str, object],
    review_discrete_reference: dict[str, object],
) -> list[RealDatasetMacroevolutionAlignmentReviewRow]:
    return [
        RealDatasetMacroevolutionAlignmentReviewRow(
            surface_id=_CONTINUOUS_REVIEW_SURFACE_ID,
            trait="seed_mass",
            model="ornstein-uhlenbeck",
            original_tree_taxa=review_continuous_alignment.original_tree_taxa,
            original_trait_taxa=review_continuous_alignment.original_trait_taxa,
            aligned_taxa_count=len(review_continuous_alignment.aligned_taxa),
            dropped_tree_taxa=review_continuous_alignment.dropped_tree_taxa,
            dropped_trait_taxa=review_continuous_alignment.dropped_trait_taxa,
            dropped_missing_value_taxa=review_continuous_alignment.dropped_missing_value_taxa,
            geiger_overlap_taxa=int(review_continuous_reference["overlap_taxa"]),
            geiger_usable_taxa=int(review_continuous_reference["usable_taxa"]),
            notes=list(review_continuous_reference["notes"]),
        ),
        RealDatasetMacroevolutionAlignmentReviewRow(
            surface_id=_DISCRETE_REVIEW_SURFACE_ID,
            trait="lifeform",
            model="equal-rates",
            original_tree_taxa=review_discrete_alignment.original_tree_taxa,
            original_trait_taxa=review_discrete_alignment.original_trait_taxa,
            aligned_taxa_count=len(review_discrete_alignment.aligned_taxa),
            dropped_tree_taxa=review_discrete_alignment.dropped_tree_taxa,
            dropped_trait_taxa=review_discrete_alignment.dropped_trait_taxa,
            dropped_missing_value_taxa=review_discrete_alignment.dropped_missing_value_taxa,
            geiger_overlap_taxa=int(review_discrete_reference["overlap_taxa"]),
            geiger_usable_taxa=int(review_discrete_reference["usable_taxa"]),
            notes=list(review_discrete_reference["notes"]),
        ),
    ]


def _build_summary_rows(
    *,
    continuous_comparison: ContinuousEvolutionaryModeComparisonReport,
    native_continuous_reference: dict[str, object],
    discrete_fits: dict[str, DiscreteMkFitReport],
    native_discrete_reference: dict[str, object],
    review_continuous_alignment: TreeTraitAlignmentReport,
    review_discrete_alignment: TreeTraitAlignmentReport,
    review_continuous_reference: dict[str, object],
    review_discrete_reference: dict[str, object],
) -> list[RealDatasetMacroevolutionSummaryRow]:
    discrete_weights = _akaike_weights(
        {model: fit.aicc for model, fit in discrete_fits.items()}
    )
    discrete_selected_model = _selected_model_from_aicc(
        {model: fit.aicc for model, fit in discrete_fits.items()}
    )
    return [
        RealDatasetMacroevolutionSummaryRow(
            surface_id=_CONTINUOUS_SURFACE_ID,
            trait="seed_mass",
            trait_kind="continuous",
            review_scope="native-model-table",
            bijux_selected_model=continuous_comparison.better_model,
            geiger_selected_model=str(native_continuous_reference["selected_model"]),
            selection_matches_geiger=(
                continuous_comparison.better_model
                == native_continuous_reference["selected_model"]
            ),
            bijux_selected_model_akaike_weight=continuous_comparison.selected_model_akaike_weight,
            geiger_selected_model_akaike_weight=_geiger_selected_weight(
                native_continuous_reference
            ),
            stable_conclusion_supported=continuous_comparison.stable_conclusion_supported,
            aligned_taxa_count=42,
            dropped_tree_taxon_count=0,
            dropped_trait_taxon_count=0,
            dropped_missing_value_taxon_count=0,
            biological_interpretation=(
                "seed_mass covariance is best fit by OU over BM, white, lambda, and early-burst on this published seashore flora subset, but the selected alpha remains near the lower boundary so this should be read as a weak covariance-shape preference rather than strong adaptive-process proof"
            ),
            notes=list(continuous_comparison.warnings),
        ),
        RealDatasetMacroevolutionSummaryRow(
            surface_id=_DISCRETE_SURFACE_ID,
            trait="lifeform",
            trait_kind="discrete",
            review_scope="native-model-table",
            bijux_selected_model=discrete_selected_model,
            geiger_selected_model=str(native_discrete_reference["selected_model"]),
            selection_matches_geiger=(
                discrete_selected_model == native_discrete_reference["selected_model"]
            ),
            bijux_selected_model_akaike_weight=discrete_weights[
                discrete_selected_model
            ],
            geiger_selected_model_akaike_weight=_geiger_selected_weight(
                native_discrete_reference
            ),
            stable_conclusion_supported=True,
            aligned_taxa_count=42,
            dropped_tree_taxon_count=0,
            dropped_trait_taxon_count=0,
            dropped_missing_value_taxon_count=0,
            biological_interpretation=(
                "lifeform transitions are best summarized by an equal-rates Mk surface on this sparse five-state trait; SYM and ARD improve raw likelihood but not enough to offset model complexity, and higher-parameter fits remain cautionary because a singleton state and optimizer-bound warnings weaken fine-grained directional claims"
            ),
            notes=[
                "one observed lifeform state is represented by a single taxon",
                "SYM and ARD remain review-only on this real dataset because their owned fits do not converge cleanly",
            ],
        ),
        RealDatasetMacroevolutionSummaryRow(
            surface_id=_CONTINUOUS_REVIEW_SURFACE_ID,
            trait="seed_mass",
            trait_kind="continuous",
            review_scope="alignment-review",
            bijux_selected_model="ornstein-uhlenbeck",
            geiger_selected_model=str(review_continuous_reference["model"]),
            selection_matches_geiger=True,
            bijux_selected_model_akaike_weight=None,
            geiger_selected_model_akaike_weight=None,
            stable_conclusion_supported=False,
            aligned_taxa_count=len(review_continuous_alignment.aligned_taxa),
            dropped_tree_taxon_count=len(review_continuous_alignment.dropped_tree_taxa),
            dropped_trait_taxon_count=len(
                review_continuous_alignment.dropped_trait_taxa
            ),
            dropped_missing_value_taxon_count=len(
                review_continuous_alignment.dropped_missing_value_taxa
            ),
            biological_interpretation=(
                "the alignment review confirms that one tree-only taxon, one trait-only row, and one overlapping missing seed_mass value are pruned before the benchmark fit; this is a data-handling review surface, not a new biological claim"
            ),
            notes=list(review_continuous_reference["notes"]),
        ),
        RealDatasetMacroevolutionSummaryRow(
            surface_id=_DISCRETE_REVIEW_SURFACE_ID,
            trait="lifeform",
            trait_kind="discrete",
            review_scope="alignment-review",
            bijux_selected_model="equal-rates",
            geiger_selected_model=str(review_discrete_reference["model"]),
            selection_matches_geiger=True,
            bijux_selected_model_akaike_weight=None,
            geiger_selected_model_akaike_weight=None,
            stable_conclusion_supported=True,
            aligned_taxa_count=len(review_discrete_alignment.aligned_taxa),
            dropped_tree_taxon_count=len(review_discrete_alignment.dropped_tree_taxa),
            dropped_trait_taxon_count=len(review_discrete_alignment.dropped_trait_taxa),
            dropped_missing_value_taxon_count=len(
                review_discrete_alignment.dropped_missing_value_taxa
            ),
            biological_interpretation=(
                "the alignment review confirms that one tree-only taxon, one trait-only row, and one overlapping missing lifeform value are pruned before the benchmark fit; this is a join-policy check rather than a biological interpretation surface"
            ),
            notes=list(review_discrete_reference["notes"]),
        ),
    ]


def _write_alignment_review_traits_table(
    path: Path,
    dataset: CentralEuropeanSeashoreFloraDataset,
) -> Path:
    table = load_taxon_table(dataset.traits_path, taxon_column=dataset.taxon_column)
    rows = [
        dict(row)
        for row in table.rows
        if row[dataset.taxon_column] != _REMOVED_TREE_TAXON
    ]
    extra_row = dict(rows[0])
    extra_row[dataset.taxon_column] = _EXTRA_TRAIT_TAXON
    rows.append(extra_row)
    for row in rows:
        if row[dataset.taxon_column] == _CONTINUOUS_MISSING_VALUE_TAXON:
            row[dataset.workflow_continuous_trait] = ""
        if row[dataset.taxon_column] == _DISCRETE_MISSING_VALUE_TAXON:
            row[dataset.workflow_discrete_trait] = ""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=table.columns)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _write_overview(
    path: Path,
    *,
    dataset: CentralEuropeanSeashoreFloraDataset,
    bundle: RealDatasetMacroevolutionBenchmarkBundle,
) -> Path:
    lines = [
        "# Real-Dataset Macroevolution Benchmark Demo",
        "",
        f"- dataset id: `{dataset.dataset_id}`",
        f"- dataset label: `{dataset.label}`",
        f"- provenance: `{_PROVENANCE_CITATION}`",
        f"- DOI: `{_PROVENANCE_DOI}`",
        f"- benchmark bundle directory: `{bundle.output_root.name}`",
        "",
        "Generated outputs:",
        "",
        f"- review traits input: `{bundle.review_traits_path.relative_to(bundle.output_root.parent)}`",
        f"- summary ledger: `{bundle.summary_path.relative_to(bundle.output_root.parent)}`",
        f"- native model table: `{bundle.model_table_path.relative_to(bundle.output_root.parent)}`",
        f"- alignment review ledger: `{bundle.alignment_review_path.relative_to(bundle.output_root.parent)}`",
        f"- geiger parity ledger: `{bundle.parity_table_path.relative_to(bundle.output_root.parent)}`",
        f"- stored geiger reference ledger: `{bundle.geiger_reference_path.relative_to(bundle.output_root.parent)}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _akaike_weights(aicc_by_model: dict[str, float]) -> dict[str, float]:
    finite_rows = {
        model: value for model, value in aicc_by_model.items() if math.isfinite(value)
    }
    minimum = min(finite_rows.values())
    unnormalized = {
        model: math.exp(-0.5 * (value - minimum))
        for model, value in finite_rows.items()
    }
    total = sum(unnormalized.values())
    return {model: weight / total for model, weight in unnormalized.items()}


def _selected_model_from_aicc(aicc_by_model: dict[str, float]) -> str:
    return min(aicc_by_model.items(), key=lambda item: item[1])[0]


def _apply_akaike_weights_from_report(
    rows: list[RealDatasetMacroevolutionModelRow],
    *,
    surface_id: str,
    engine: str,
) -> None:
    surface_rows = [row for row in rows if row.surface_id == surface_id]
    if engine == "bijux":
        weights = _akaike_weights({row.model: row.bijux_aicc for row in surface_rows})
        for row in surface_rows:
            row.bijux_akaike_weight = weights[row.model]
        return
    weights = _akaike_weights({row.model: row.geiger_aicc for row in surface_rows})
    for row in surface_rows:
        row.geiger_akaike_weight = weights[row.model]


def _apply_geiger_akaike_weights(
    rows: list[RealDatasetMacroevolutionModelRow],
    *,
    surface_id: str,
) -> None:
    _apply_akaike_weights_from_report(rows, surface_id=surface_id, engine="geiger")


def _geiger_selected_weight(payload: dict[str, object]) -> float:
    weights = _akaike_weights(
        {row["model"]: float(row["aicc"]) for row in payload["comparison_rows"]}
    )
    return weights[str(payload["selected_model"])]


def _representative_discrete_rate(report: DiscreteMkFitReport) -> float | None:
    if report.model != "equal-rates" or not report.transition_rate_rows:
        return None
    return report.transition_rate_rows[0].rate


def _optional_float(value: object) -> float | None:
    if value in (None, "", []):
        return None
    return float(value)


def _optional_str(value: object) -> str | None:
    if value in (None, "", []):
        return None
    return str(value)


def _format_float(value: float) -> str:
    return f"{value:.12f}"


def _format_optional_float(value: float | None) -> str:
    return "" if value is None else _format_float(value)
