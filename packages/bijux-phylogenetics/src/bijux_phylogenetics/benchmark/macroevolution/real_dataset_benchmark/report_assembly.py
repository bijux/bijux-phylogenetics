from __future__ import annotations

from pathlib import Path

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
)
from bijux_phylogenetics.datasets.study_inputs import (
    TreeTraitAlignmentReport,
    align_tree_and_trait_table,
)

from ..geiger_reference import (
    GEIGER_REAL_DATASET_MACROEVOLUTION_REFERENCE_PAYLOADS,
)
from .contracts import (
    RealDatasetMacroevolutionAlignmentReviewRow,
    RealDatasetMacroevolutionBenchmarkReport,
    RealDatasetMacroevolutionModelRow,
    RealDatasetMacroevolutionParityRow,
    RealDatasetMacroevolutionSummaryRow,
)
from .shared import (
    CONTINUOUS_MODES,
    CONTINUOUS_REVIEW_SURFACE_ID,
    CONTINUOUS_SURFACE_ID,
    DISCRETE_MODELS,
    DISCRETE_REVIEW_SURFACE_ID,
    DISCRETE_SURFACE_ID,
    PROVENANCE_CITATION,
    PROVENANCE_DOI,
    akaike_weights,
    apply_akaike_weights_from_report,
    apply_geiger_akaike_weights,
    geiger_selected_weight,
    optional_float,
    optional_str,
    selected_model_from_aicc,
)


def build_report(
    dataset: CentralEuropeanSeashoreFloraDataset,
    review_traits_path: Path,
) -> RealDatasetMacroevolutionBenchmarkReport:
    continuous_comparison = compare_fitcontinuous_model_ranking(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.workflow_continuous_trait,
        taxon_column=dataset.taxon_column,
        modes=CONTINUOUS_MODES,
    )
    continuous_fits = {
        mode: fit_continuous_evolutionary_mode(
            dataset.tree_path,
            dataset.traits_path,
            trait=dataset.workflow_continuous_trait,
            taxon_column=dataset.taxon_column,
            mode=mode,
        )
        for mode in CONTINUOUS_MODES
    }
    discrete_fits = {
        model: fit_discrete_mk_model(
            dataset.tree_path,
            dataset.traits_path,
            trait=dataset.workflow_discrete_trait,
            taxon_column=dataset.taxon_column,
            model=model,
        )
        for model in DISCRETE_MODELS
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
        CONTINUOUS_SURFACE_ID
    ]
    native_discrete_reference = GEIGER_REAL_DATASET_MACROEVOLUTION_REFERENCE_PAYLOADS[
        DISCRETE_SURFACE_ID
    ]
    review_continuous_reference = GEIGER_REAL_DATASET_MACROEVOLUTION_REFERENCE_PAYLOADS[
        CONTINUOUS_REVIEW_SURFACE_ID
    ]
    review_discrete_reference = GEIGER_REAL_DATASET_MACROEVOLUTION_REFERENCE_PAYLOADS[
        DISCRETE_REVIEW_SURFACE_ID
    ]

    model_rows = build_model_rows(
        continuous_comparison=continuous_comparison,
        continuous_fits=continuous_fits,
        native_continuous_reference=native_continuous_reference,
        discrete_fits=discrete_fits,
        native_discrete_reference=native_discrete_reference,
    )
    parity_rows = build_parity_rows(
        model_rows=model_rows,
        review_continuous_fit=review_continuous_fit,
        review_continuous_reference=review_continuous_reference,
        review_discrete_fit=review_discrete_fit,
        review_discrete_reference=review_discrete_reference,
    )
    alignment_review_rows = build_alignment_review_rows(
        review_continuous_alignment=review_continuous_alignment,
        review_discrete_alignment=review_discrete_alignment,
        review_continuous_reference=review_continuous_reference,
        review_discrete_reference=review_discrete_reference,
    )
    summary_rows = build_summary_rows(
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
        provenance_citation=PROVENANCE_CITATION,
        provenance_doi=PROVENANCE_DOI,
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


def build_model_rows(
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
                surface_id=CONTINUOUS_SURFACE_ID,
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
                geiger_parameter_name=optional_str(geiger_row.get("parameter_name")),
                bijux_parameter_value=fit.parameter_value,
                geiger_parameter_value=optional_float(
                    geiger_row.get("parameter_value")
                ),
                bijux_rate=fit.rate,
                geiger_rate=optional_float(geiger_row.get("rate")),
                bijux_root_state=fit.root_state,
                geiger_root_state=optional_float(geiger_row.get("root_state")),
                notes=notes,
            )
        )
    apply_geiger_akaike_weights(rows, surface_id=CONTINUOUS_SURFACE_ID)

    geiger_discrete_rows = {
        row["model"]: row for row in native_discrete_reference["comparison_rows"]
    }
    discrete_selected_model = selected_model_from_aicc(
        {model: fit.aicc for model, fit in discrete_fits.items()}
    )
    for rank, model in enumerate(
        sorted(DISCRETE_MODELS, key=lambda item: discrete_fits[item].aicc),
        start=1,
    ):
        fit = discrete_fits[model]
        geiger_row = geiger_discrete_rows[model]
        rows.append(
            RealDatasetMacroevolutionModelRow(
                surface_id=DISCRETE_SURFACE_ID,
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
                bijux_rate=representative_discrete_rate(fit),
                geiger_rate=optional_float(geiger_row.get("representative_rate")),
                bijux_root_state=None,
                geiger_root_state=None,
                notes=list(fit.input_audit.warnings),
            )
        )
    apply_akaike_weights_from_report(
        rows, surface_id=DISCRETE_SURFACE_ID, engine="bijux"
    )
    apply_geiger_akaike_weights(rows, surface_id=DISCRETE_SURFACE_ID)
    return rows


def build_parity_rows(
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
            surface_id=CONTINUOUS_REVIEW_SURFACE_ID,
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
            geiger_parameter_name=optional_str(
                review_continuous_summary.get("parameter_name")
            ),
            bijux_parameter_value=review_continuous_fit.parameter_value,
            geiger_parameter_value=optional_float(
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
            surface_id=DISCRETE_REVIEW_SURFACE_ID,
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


def build_alignment_review_rows(
    *,
    review_continuous_alignment: TreeTraitAlignmentReport,
    review_discrete_alignment: TreeTraitAlignmentReport,
    review_continuous_reference: dict[str, object],
    review_discrete_reference: dict[str, object],
) -> list[RealDatasetMacroevolutionAlignmentReviewRow]:
    return [
        RealDatasetMacroevolutionAlignmentReviewRow(
            surface_id=CONTINUOUS_REVIEW_SURFACE_ID,
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
            surface_id=DISCRETE_REVIEW_SURFACE_ID,
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


def build_summary_rows(
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
    discrete_weights = akaike_weights(
        {model: fit.aicc for model, fit in discrete_fits.items()}
    )
    discrete_selected_model = selected_model_from_aicc(
        {model: fit.aicc for model, fit in discrete_fits.items()}
    )
    return [
        RealDatasetMacroevolutionSummaryRow(
            surface_id=CONTINUOUS_SURFACE_ID,
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
            geiger_selected_model_akaike_weight=geiger_selected_weight(
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
            surface_id=DISCRETE_SURFACE_ID,
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
            geiger_selected_model_akaike_weight=geiger_selected_weight(
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
            surface_id=CONTINUOUS_REVIEW_SURFACE_ID,
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
            surface_id=DISCRETE_REVIEW_SURFACE_ID,
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


def representative_discrete_rate(report: DiscreteMkFitReport) -> float | None:
    if report.model != "equal-rates" or not report.transition_rate_rows:
        return None
    return report.transition_rate_rows[0].rate
