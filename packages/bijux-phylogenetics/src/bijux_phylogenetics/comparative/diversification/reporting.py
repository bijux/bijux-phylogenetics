from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.render.html import write_html_report

from .clades import (
    detect_diversification_outlier_clades,
)
from .lineage import (
    compute_lineage_through_time_curve,
)
from .models import (
    DiversificationGammaStatisticReport,
    DiversificationMethodReport,
    DiversificationMethodsSummaryTextResult,
    DiversificationModelComparisonReport,
    DiversificationReportBuildResult,
    GeigerBirthDeathExclusionReport,
    MedusaExclusionReport,
)
from .rates import (
    compare_diversification_models,
    compute_diversification_gamma_statistic,
    estimate_diversification_rate,
)
from .rates import (
    geiger_birth_death_exclusion_reason as _geiger_birth_death_exclusion_reason,
)
from .sampling import detect_incomplete_taxon_sampling_metadata
from .traits import (
    run_trait_dependent_diversification_analysis,
)
from .trees import (
    validate_time_tree_for_diversification,
)


def _bullet_list(values: list[str]) -> str:
    if not values:
        return "none"
    return ", ".join(f"`{value}`" for value in values)


def _deduplicate_text(values: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        normalized = value.strip()
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        ordered.append(normalized)
    return ordered


def _medusa_supported_surfaces() -> list[str]:
    return [
        "time-tree validation",
        "lineage-through-time curve",
        "sampling-fraction review",
        "yule and birth-death rate estimation",
        "gamma-statistic summary",
        "two-model diversification AIC comparison",
        "descriptive clade diversification outlier scan",
        "trait-linked diversification summary",
    ]


def _medusa_missing_surfaces() -> list[str]:
    return [
        "stepwise branch-specific rate-shift search",
        "shift-count model growth with information-criterion stopping",
        "best shift-placement ranking across candidate branch partitions",
        "background-versus-shift diversification parameterization",
        "weak-support review for near-tied shift configurations",
    ]


def _medusa_exclusion_reason() -> str:
    return (
        "geiger::medusa parity is explicitly excluded in this round because bijux "
        "does not yet implement MEDUSA's stepwise diversification rate-shift search, "
        "shift-count model-growth selection, or branch-placement ranking; existing "
        "clade diversification summaries are descriptive reviews and are not claimed "
        "as MEDUSA-equivalent rate-shift detection"
    )


def build_diversification_method_report(
    tree_path: Path,
    *,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
    sampling_column: str | None = None,
    traits_path: Path | None = None,
    trait: str | None = None,
    estimate_model: str = "birth-death",
    clade_model: str = "birth-death",
    clade_min_tip_count: int = 2,
) -> DiversificationMethodReport:
    """Build one integrated diversification-method report from governed evidence."""
    validation = validate_time_tree_for_diversification(tree_path)
    lineage = compute_lineage_through_time_curve(tree_path)
    gamma_statistic = compute_diversification_gamma_statistic(
        tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
    )
    primary_estimate = estimate_diversification_rate(
        tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
        model=estimate_model,
    )
    model_comparison = compare_diversification_models(
        tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
    )
    clade_scan = detect_diversification_outlier_clades(
        tree_path,
        min_tip_count=clade_min_tip_count,
        model=clade_model,
    )
    sampling_report = (
        detect_incomplete_taxon_sampling_metadata(
            tree_path,
            metadata_path,
            taxon_column=taxon_column,
            sampling_column=sampling_column,
        )
        if metadata_path is not None
        else None
    )
    trait_report = (
        run_trait_dependent_diversification_analysis(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
        )
        if traits_path is not None and trait is not None
        else None
    )
    return DiversificationMethodReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
        estimate_model=estimate_model,
        clade_model=clade_model,
        clade_min_tip_count=clade_min_tip_count,
        validation=validation,
        lineage=lineage,
        gamma_statistic=gamma_statistic,
        primary_estimate=primary_estimate,
        model_comparison=model_comparison,
        clade_scan=clade_scan,
        sampling_report=sampling_report,
        trait_report=trait_report,
    )


def _diversification_methods_summary_assumptions(
    report: DiversificationMethodReport,
) -> list[str]:
    assumptions = [
        *report.gamma_statistic.assumptions,
        *report.primary_estimate.assumptions,
        _medusa_exclusion_reason(),
    ]
    if report.sampling_report is None:
        assumptions.append(
            "no sampling metadata was supplied, so all rate summaries default to the complete-sampling assumption"
        )
    elif report.sampling_report.heterogeneous_values:
        assumptions.append(
            "heterogeneous taxon sampling fractions are collapsed to the mean sampling fraction before correction is applied"
        )
    if report.trait_report is not None:
        assumptions.append(
            "trait-linked diversification summaries remain descriptive and only convert state-specific crown ages into rates when the state taxa are monophyletic"
        )
    return _deduplicate_text(assumptions)


def _diversification_methods_summary_warnings(
    report: DiversificationMethodReport,
) -> list[str]:
    warnings = [
        *report.validation.warnings,
        *report.gamma_statistic.warnings,
        *report.primary_estimate.warnings,
        *report.clade_scan.warnings,
    ]
    if report.sampling_report is not None:
        warnings.extend(report.sampling_report.warnings)
    if report.trait_report is not None:
        warnings.extend(report.trait_report.warnings)
    return _deduplicate_text(warnings)


def _diversification_report_limitations(
    report: DiversificationMethodReport,
) -> list[str]:
    limitations = [
        "diversification rate, gamma-statistic, and clade outlier summaries depend on the supplied ultrametric tree and should not be treated as direct proof of diversification mechanism",
        "model-ranking differences do not by themselves establish a biological process without checking sampling completeness, model adequacy, and tree uncertainty",
        _medusa_exclusion_reason(),
        *report.validation.warnings,
        *report.gamma_statistic.assumptions,
        *report.primary_estimate.assumptions,
        *_diversification_methods_summary_warnings(report),
    ]
    return _deduplicate_text(limitations)


def build_diversification_methods_summary_text(
    report: DiversificationMethodReport,
) -> str:
    """Build reviewer-facing Markdown methods text for one diversification analysis."""
    comparison_rows = sorted(report.model_comparison.rows, key=lambda row: row.aic)
    better_row = comparison_rows[0]
    runner_up_delta = (
        comparison_rows[1].aic - comparison_rows[0].aic
        if len(comparison_rows) > 1
        else 0.0
    )
    sampling_report = report.sampling_report
    trait_report = report.trait_report
    warnings = _diversification_methods_summary_warnings(report)
    assumptions = _diversification_methods_summary_assumptions(report)
    sampling_fraction_text = (
        "not available"
        if sampling_report is None or sampling_report.sampling_fraction is None
        else format(sampling_report.sampling_fraction, ".15g")
    )
    trait_state_count = 0 if trait_report is None else len(trait_report.states)
    trait_monophyletic_count = (
        0
        if trait_report is None
        else sum(1 for row in trait_report.states if row.monophyletic)
    )
    return (
        "# Diversification Analysis Methods Summary\n\n"
        f"This diversification analysis evaluated rooted ultrametric time tree `{report.tree_path.name}`"
        + (
            f" with sampling metadata `{report.metadata_path.name}`"
            if report.metadata_path is not None
            else " without an external sampling-metadata table"
        )
        + (
            f" and trait table `{report.traits_path.name}` for state column `{trait_report.trait}`."
            if trait_report is not None
            else "."
        )
        + "\n\n## Time Tree And Inputs\n\n"
        f"- rooted time tree required: `{'yes' if report.validation.rooted else 'no'}`\n"
        f"- ultrametric time tree required: `{'yes' if report.validation.ultrametric else 'no'}`\n"
        f"- branch-length completeness: `{report.validation.branch_length_status}`\n"
        f"- analyzed tip count: `{report.validation.tip_count}`\n"
        f"- crown age: `{format(report.validation.root_age, '.15g')}`\n"
        f"- lineage-through-time points retained: `{len(report.lineage.points)}`\n"
        f"- tree validation warnings: {_bullet_list(report.validation.warnings)}\n\n"
        "## Sampling Correction\n\n"
        + (
            "- sampling metadata: not supplied, so diversification rates follow the complete-sampling assumption\n"
            if sampling_report is None
            else (
                f"- taxon column: `{sampling_report.taxon_column}`\n"
                f"- sampling column: `{sampling_report.sampling_column or 'missing'}`\n"
                f"- matched taxa with sampling rows: `{len(sampling_report.matched_taxa)}`\n"
                f"- tree tips missing sampling rows: `{len(sampling_report.missing_taxa)}`\n"
                f"- invalid sampling rows: `{len(sampling_report.invalid_rows)}`\n"
                f"- sampling metadata complete: `{'yes' if sampling_report.complete else 'no'}`\n"
                f"- mean sampling fraction used for correction: `{sampling_fraction_text}`\n"
                f"- heterogeneous sampling fractions: `{'yes' if sampling_report.heterogeneous_values else 'no'}`\n"
                f"- sampling warnings: {_bullet_list(sampling_report.warnings)}\n"
            )
        )
        + "\n## Diversification Models And Rates\n\n"
        f"- primary reported rate model: `{report.estimate_model}`\n"
        f"- compared candidate models: {_bullet_list([row.model for row in comparison_rows])}\n"
        f"- better-supported model by AIC: `{report.model_comparison.better_model}`\n"
        f"- better-model AIC: `{format(better_row.aic, '.15g')}`\n"
        f"- runner-up delta AIC: `{format(runner_up_delta, '.15g')}`\n"
        f"- reported net diversification rate: `{format(report.primary_estimate.net_diversification_rate, '.15g')}`\n"
        f"- reported relative extinction: `{format(report.primary_estimate.relative_extinction, '.15g')}`\n"
        f"- corrected tip count under the reported model: `{format(report.primary_estimate.corrected_tip_count, '.15g')}`\n"
        f"- Pybus-Harvey gamma statistic: `{format(report.gamma_statistic.gamma_statistic, '.15g')}`\n\n"
        "## Clade And Trait Review\n\n"
        f"- clade outlier scan model: `{report.clade_model}`\n"
        f"- minimum clade size for outlier review: `{report.clade_min_tip_count}`\n"
        f"- evaluated clades: `{len(report.clade_scan.observations)}`\n"
        f"- high-diversification outliers: `{len(report.clade_scan.high_diversification_clades)}`\n"
        f"- low-diversification outliers: `{len(report.clade_scan.low_diversification_clades)}`\n"
        + (
            "- trait-linked diversification surface: not requested\n"
            if trait_report is None
            else (
                f"- trait-linked diversification trait: `{trait_report.trait}`\n"
                f"- observed trait states reviewed: `{trait_state_count}`\n"
                f"- monophyletic states with interpretable crown ages: `{trait_monophyletic_count}`\n"
                f"- trait-linked warnings: {_bullet_list(trait_report.warnings)}\n"
            )
        )
        + f"- MEDUSA parity claim: {_medusa_exclusion_reason()}\n"
        + "\n## Assumptions And Caveats\n\n"
        + "\n".join(f"- {item}" for item in assumptions)
        + "\n\n## Reviewer Warnings\n\n"
        + f"- combined warning count: `{len(warnings)}`\n"
        + f"- warning details: {_bullet_list(warnings)}\n"
    )


def summarize_medusa_exclusion(
    tree_path: Path,
    *,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
    sampling_column: str | None = None,
) -> MedusaExclusionReport:
    """Explain the current exclusion boundary for geiger::medusa parity."""
    validation = validate_time_tree_for_diversification(tree_path)
    sampling_report = (
        None
        if metadata_path is None
        else detect_incomplete_taxon_sampling_metadata(
            tree_path,
            metadata_path,
            taxon_column=taxon_column,
            sampling_column=sampling_column,
        )
    )
    exclusion_reason = _medusa_exclusion_reason()
    warnings = list(validation.warnings)
    if sampling_report is not None:
        warnings.extend(sampling_report.warnings)
    warnings.append(exclusion_reason)
    return MedusaExclusionReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        validation=validation,
        sampling_report=sampling_report,
        supported_surfaces=_medusa_supported_surfaces(),
        missing_surfaces=_medusa_missing_surfaces(),
        exclusion_code="geiger_medusa_explicitly_excluded_this_round",
        exclusion_reason=exclusion_reason,
        warnings=_deduplicate_text(warnings),
    )


def summarize_geiger_birth_death_exclusion(
    tree_path: Path,
    *,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
    sampling_column: str | None = None,
) -> GeigerBirthDeathExclusionReport:
    """Explain why current diversification summaries do not claim geiger bd.ms parity."""
    validation = validate_time_tree_for_diversification(tree_path)
    sampling_report = (
        None
        if metadata_path is None
        else detect_incomplete_taxon_sampling_metadata(
            tree_path,
            metadata_path,
            taxon_column=taxon_column,
            sampling_column=sampling_column,
        )
    )
    exclusion_reason = _geiger_birth_death_exclusion_reason()
    warnings = list(validation.warnings)
    if sampling_report is not None:
        warnings.extend(sampling_report.warnings)
    warnings.append(exclusion_reason)
    return GeigerBirthDeathExclusionReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        validation=validation,
        sampling_report=sampling_report,
        geiger_reference_surface="geiger::bd.ms",
        geiger_reference_arguments=["phy", "time", "n", "missing", "crown", "epsilon"],
        owned_surface="heuristic-yule-and-birth-death-diversification-summary",
        exclusion_code="geiger_birth_death_explicitly_excluded_this_round",
        exclusion_reason=exclusion_reason,
        warnings=_deduplicate_text(warnings),
    )


def write_diversification_methods_summary_text(
    path: Path,
    report: DiversificationMethodReport,
) -> DiversificationMethodsSummaryTextResult:
    """Write reviewer-facing Markdown methods text for one diversification analysis."""
    text = build_diversification_methods_summary_text(report)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return DiversificationMethodsSummaryTextResult(
        output_path=path,
        title="Diversification Analysis Methods Summary",
        warning_count=len(_diversification_methods_summary_warnings(report)),
        better_model=report.model_comparison.better_model,
        sampling_metadata_complete=(
            None if report.sampling_report is None else report.sampling_report.complete
        ),
        clade_observation_count=len(report.clade_scan.observations),
        text=text,
        report=report,
    )


def write_diversification_gamma_statistic_table(
    path: Path, report: DiversificationGammaStatisticReport
) -> Path:
    """Export one deterministic diversification gamma-statistic ledger."""
    rows = [
        {
            "tip_count": str(report.tip_count),
            "rooted": str(report.rooted).lower(),
            "ultrametric": str(report.ultrametric).lower(),
            "bifurcating": str(report.bifurcating).lower(),
            "root_age": format(report.root_age, ".15g"),
            "branching_time_count": str(report.branching_time_count),
            "interval_count": str(report.interval_count),
            "minimum_branching_time": format(report.minimum_branching_time, ".15g"),
            "maximum_branching_time": format(report.maximum_branching_time, ".15g"),
            "gamma_statistic": format(report.gamma_statistic, ".15g"),
            "sampling_fraction": ""
            if report.sampling_fraction is None
            else format(report.sampling_fraction, ".15g"),
            "assumptions": "; ".join(report.assumptions),
            "warnings": "; ".join(report.warnings),
        }
    ]
    return write_taxon_rows(
        path,
        columns=[
            "tip_count",
            "rooted",
            "ultrametric",
            "bifurcating",
            "root_age",
            "branching_time_count",
            "interval_count",
            "minimum_branching_time",
            "maximum_branching_time",
            "gamma_statistic",
            "sampling_fraction",
            "assumptions",
            "warnings",
        ],
        rows=rows,
    )


def write_diversification_model_comparison_table(
    path: Path, report: DiversificationModelComparisonReport
) -> Path:
    """Export one deterministic diversification model-comparison ledger."""
    rows = [
        {
            "model": row.model,
            "parameter_count": str(row.parameter_count),
            "log_likelihood": format(row.log_likelihood, ".15g"),
            "aic": format(row.aic, ".15g"),
            "sampling_fraction": format(row.sampling_fraction, ".15g"),
            "net_diversification_rate": format(row.net_diversification_rate, ".15g"),
            "relative_extinction": format(row.relative_extinction, ".15g"),
            "better_model": str(row.model == report.better_model).lower(),
        }
        for row in report.rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "model",
            "parameter_count",
            "log_likelihood",
            "aic",
            "sampling_fraction",
            "net_diversification_rate",
            "relative_extinction",
            "better_model",
        ],
        rows=rows,
    )


def render_diversification_report(
    *,
    tree_path: Path,
    out_path: Path,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
    sampling_column: str | None = None,
    traits_path: Path | None = None,
    trait: str | None = None,
    methods_summary_path: Path | None = None,
) -> DiversificationReportBuildResult:
    """Render a deterministic HTML report for diversification and macroevolution summaries."""
    report = build_diversification_method_report(
        tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
        traits_path=traits_path,
        trait=trait,
        estimate_model="birth-death",
        clade_model="birth-death",
        clade_min_tip_count=2,
    )
    methods_summary_text = build_diversification_methods_summary_text(report)
    methods_summary = (
        None
        if methods_summary_path is None
        else write_diversification_methods_summary_text(methods_summary_path, report)
    )
    sections = [
        ("methods-summary-text", methods_summary_text),
        ("lineage-through-time", json.dumps(report.lineage, default=str, indent=2)),
        (
            "diversification-gamma-statistic",
            json.dumps(report.gamma_statistic, default=str, indent=2),
        ),
        (
            "diversification-estimate",
            json.dumps(report.primary_estimate, default=str, indent=2),
        ),
        (
            "diversification-model-comparison",
            json.dumps(report.model_comparison, default=str, indent=2),
        ),
        (
            "clade-diversification-scan",
            json.dumps(report.clade_scan, default=str, indent=2),
        ),
    ]
    if report.trait_report is not None:
        sections.append(
            (
                "trait-dependent-diversification",
                json.dumps(report.trait_report, default=str, indent=2),
            )
        )
    limitations = _diversification_report_limitations(report)
    sections.append(("limitations", json.dumps(limitations, indent=2)))
    title = "Bijux Diversification Report"
    manifest = {
        "report_kind": "diversification",
        "tree_path": str(tree_path),
        "metadata_path": None if metadata_path is None else str(metadata_path),
        "traits_path": None if traits_path is None else str(traits_path),
        "trait": trait,
        "sections": [name for name, _value in sections],
        "limitations": limitations,
        "outputs": {
            "methods_summary_path": None
            if methods_summary is None
            else str(methods_summary.output_path)
        },
        "metrics": {
            "methods_summary_warning_count": len(
                _diversification_methods_summary_warnings(report)
            ),
            "better_model": report.model_comparison.better_model,
        },
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=manifest,
    )
    return DiversificationReportBuildResult(
        output_path=out_path,
        report_kind="diversification",
        title=title,
        tree_path=tree_path,
        machine_manifest=manifest,
        methods_summary_text=methods_summary_text,
        methods_summary_warning_count=len(
            _diversification_methods_summary_warnings(report)
        ),
        methods_summary_path=None
        if methods_summary is None
        else methods_summary.output_path,
        report=report,
    )
