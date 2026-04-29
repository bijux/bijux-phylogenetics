from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import cast

from bijux_phylogenetics.ancestral.common import reconstruction_manifest, write_ancestral_rows
from bijux_phylogenetics.ancestral.continuous import ContinuousAncestralReport, reconstruct_continuous_ancestral_states
from bijux_phylogenetics.ancestral.discrete import DiscreteAncestralReport, reconstruct_discrete_ancestral_states
from bijux_phylogenetics.discrete_evolution import run_discrete_state_transition_model
from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.render.svg import TreeRenderResult, render_tree_svg


@dataclass(slots=True)
class ContinuousAncestralComparisonRow:
    """Node-by-node comparison between two continuous ancestral reconstructions."""

    node: str
    descendant_taxa: list[str]
    left_estimate: float
    right_estimate: float
    estimate_delta: float
    intervals_overlap: bool


@dataclass(slots=True)
class ContinuousAncestralComparisonReport:
    """Comparison between two continuous ancestral-state models."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    left_model: str
    right_model: str
    rows: list[ContinuousAncestralComparisonRow]


@dataclass(slots=True)
class DiscreteAncestralModelComparisonRow:
    """Model-comparison summary for one discrete ancestral-state model."""

    model: str
    parameter_count: int
    pseudo_log_likelihood: float
    aic: float
    selected: bool


@dataclass(slots=True)
class DiscreteAncestralModelDifference:
    """Node-level state difference between the selected and another discrete model."""

    comparison_model: str
    node: str
    descendant_taxa: list[str]
    selected_state: str
    comparison_state: str
    differs: bool


@dataclass(slots=True)
class DiscreteAncestralModelComparisonReport:
    """Comparison of ER, SYM, and ARD discrete ancestral reconstructions."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    selected_model: str
    rows: list[DiscreteAncestralModelComparisonRow]
    node_differences: list[DiscreteAncestralModelDifference]


@dataclass(slots=True)
class AncestralStateReportBuildResult:
    """HTML report artifact for ancestral-state reconstruction."""

    output_path: Path
    report_kind: str
    title: str
    tree_path: Path
    traits_path: Path
    trait: str
    reconstruction_kind: str
    model: str
    machine_manifest: dict[str, object]


def compare_continuous_ancestral_models(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    left_model: str = "brownian",
    right_model: str = "ou",
    left_alpha: float = 1.0,
    right_alpha: float = 1.0,
) -> ContinuousAncestralComparisonReport:
    """Compare two continuous ancestral reconstructions node by node."""
    left = reconstruct_continuous_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=left_model,
        alpha=left_alpha,
    )
    right = reconstruct_continuous_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=right_model,
        alpha=right_alpha,
    )
    right_by_node = {estimate.node: estimate for estimate in right.estimates}
    rows: list[ContinuousAncestralComparisonRow] = []
    for left_estimate in left.estimates:
        if left_estimate.node not in right_by_node:
            continue
        right_estimate = right_by_node[left_estimate.node]
        rows.append(
            ContinuousAncestralComparisonRow(
                node=left_estimate.node,
                descendant_taxa=left_estimate.descendant_taxa,
                left_estimate=left_estimate.estimate,
                right_estimate=right_estimate.estimate,
                estimate_delta=right_estimate.estimate - left_estimate.estimate,
                intervals_overlap=not (
                    left_estimate.upper_95_interval < right_estimate.lower_95_interval
                    or right_estimate.upper_95_interval < left_estimate.lower_95_interval
                ),
            )
        )
    return ContinuousAncestralComparisonReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=left.taxon_count,
        left_model=left_model,
        right_model=right_model,
        rows=rows,
    )


def compare_discrete_ancestral_models(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    models: tuple[str, str, str] = ("equal-rates", "symmetric", "all-rates-different"),
) -> DiscreteAncestralModelComparisonReport:
    """Compare supported discrete ancestral-state likelihood models and select the lowest-AIC fit."""
    reconstructions = {
        model: reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
        )
        for model in models
    }
    transition_reports = {
        model: run_discrete_state_transition_model(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
        )
        for model in models
    }
    rows = [
        DiscreteAncestralModelComparisonRow(
            model=model,
            parameter_count=transition_reports[model].transition_model.parameter_count,
            pseudo_log_likelihood=transition_reports[model].transition_model.pseudo_log_likelihood,
            aic=transition_reports[model].transition_model.aic,
            selected=False,
        )
        for model in models
    ]
    selected_model = min(rows, key=lambda row: row.aic).model
    for row in rows:
        row.selected = row.model == selected_model
    selected_report = reconstructions[selected_model]
    selected_by_node = {estimate.node: estimate for estimate in selected_report.estimates}
    node_differences: list[DiscreteAncestralModelDifference] = []
    for model, report in reconstructions.items():
        if model == selected_model:
            continue
        for estimate in report.estimates:
            selected_estimate = selected_by_node[estimate.node]
            node_differences.append(
                DiscreteAncestralModelDifference(
                    comparison_model=model,
                    node=estimate.node,
                    descendant_taxa=estimate.descendant_taxa,
                    selected_state=selected_estimate.most_likely_state,
                    comparison_state=estimate.most_likely_state,
                    differs=selected_estimate.most_likely_state != estimate.most_likely_state,
                )
            )
    return DiscreteAncestralModelComparisonReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=selected_report.taxon_count,
        selected_model=selected_model,
        rows=rows,
        node_differences=node_differences,
    )


def write_ancestral_state_table(path: Path, report: ContinuousAncestralReport | DiscreteAncestralReport) -> Path:
    """Export an ancestral-state report as a deterministic TSV table."""
    if isinstance(report, ContinuousAncestralReport):
        rows = [
            {
                "node": estimate.node,
                "node_name": estimate.node_name or "",
                "is_tip": str(estimate.is_tip).lower(),
                "descendant_taxa": ",".join(estimate.descendant_taxa),
                "estimate": str(estimate.estimate),
                "standard_error": str(estimate.standard_error),
                "lower_95_interval": str(estimate.lower_95_interval),
                "upper_95_interval": str(estimate.upper_95_interval),
            }
            for estimate in report.estimates
        ]
        return write_ancestral_rows(
            path,
            columns=[
                "node",
                "node_name",
                "is_tip",
                "descendant_taxa",
                "estimate",
                "standard_error",
                "lower_95_interval",
                "upper_95_interval",
            ],
            rows=rows,
        )
    rows = [
        {
            "node": estimate.node,
            "node_name": estimate.node_name or "",
            "is_tip": str(estimate.is_tip).lower(),
            "descendant_taxa": ",".join(estimate.descendant_taxa),
            "most_likely_state": estimate.most_likely_state,
            "state_set": ",".join(estimate.state_set),
            "state_probabilities": json.dumps(estimate.state_probabilities, sort_keys=True),
            "ambiguous": str(estimate.ambiguous).lower(),
        }
        for estimate in report.estimates
    ]
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "node_name",
            "is_tip",
            "descendant_taxa",
            "most_likely_state",
            "state_set",
            "state_probabilities",
            "ambiguous",
        ],
        rows=rows,
    )


def render_ancestral_state_tree(
    tree_path: Path,
    report: ContinuousAncestralReport | DiscreteAncestralReport,
    *,
    out_path: Path,
    layout: str = "cladogram",
) -> TreeRenderResult:
    """Render one tree with internal ancestral-state annotations."""
    internal_annotations: dict[str, str]
    internal_annotation_colors: dict[str, str]
    if isinstance(report, ContinuousAncestralReport):
        internal_annotations = {
            estimate.node: format(estimate.estimate, ".3g")
            for estimate in report.estimates
            if not estimate.is_tip
        }
        internal_annotation_colors = {estimate.node: "#6d28d9" for estimate in report.estimates if not estimate.is_tip}
        continuous_traits = {
            estimate.node_name: estimate.estimate
            for estimate in report.estimates
            if estimate.is_tip and estimate.node_name is not None
        }
        return render_tree_svg(
            tree_path,
            out_path=out_path,
            layout=layout,
            continuous_traits=continuous_traits,
            internal_annotations=internal_annotations,
            internal_annotation_colors=internal_annotation_colors,
        )

    internal_annotations = {
        estimate.node: estimate.most_likely_state if not estimate.ambiguous else "/".join(estimate.state_set)
        for estimate in report.estimates
        if not estimate.is_tip
    }
    discrete_tip_states = {
        estimate.node_name: estimate.most_likely_state
        for estimate in report.estimates
        if estimate.is_tip and estimate.node_name is not None
    }
    palette = {
        state: color
        for state, color in zip(
            sorted(report.observed_states),
            ("#0f766e", "#1d4ed8", "#c2410c", "#7c3aed", "#b91c1c", "#047857"),
            strict=False,
        )
    }
    internal_annotation_colors = {
        estimate.node: palette.get(estimate.most_likely_state, "#6d28d9")
        for estimate in report.estimates
        if not estimate.is_tip
    }
    return render_tree_svg(
        tree_path,
        out_path=out_path,
        layout=layout,
        categorical_traits=discrete_tip_states,
        internal_annotations=internal_annotations,
        internal_annotation_colors=internal_annotation_colors,
    )


def render_ancestral_state_report(
    *,
    tree_path: Path,
    traits_path: Path,
    trait: str,
    reconstruction_kind: str,
    out_path: Path,
    taxon_column: str | None = None,
    model: str = "brownian",
    alpha: float = 1.0,
    compare_model: str | None = None,
) -> AncestralStateReportBuildResult:
    """Build a deterministic HTML report for ancestral-state reconstruction."""
    if reconstruction_kind == "continuous":
        report = reconstruct_continuous_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
            alpha=alpha,
        )
        comparison = (
            compare_continuous_ancestral_models(
                tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
                left_model=model,
                right_model=compare_model,
                left_alpha=alpha,
                right_alpha=alpha,
            )
            if compare_model is not None
            else None
        )
    else:
        report = reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
        )
        comparison = None
    render_path = out_path.with_suffix(".svg")
    render_result = render_ancestral_state_tree(tree_path, report, out_path=render_path, layout="phylogram")
    sections = [
        ("ancestral-reconstruction", json.dumps(asdict(report), indent=2, sort_keys=True, default=str)),
        ("ancestral-render", json.dumps(asdict(render_result), indent=2, sort_keys=True, default=str)),
    ]
    if comparison is not None:
        sections.append(("ancestral-comparison", json.dumps(asdict(comparison), indent=2, sort_keys=True, default=str)))
    title = f"Bijux Ancestral State Report: {trait}"
    machine_manifest = reconstruction_manifest(
        report_kind="ancestral-state",
        title=title,
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        model=model,
        rendered_tree=str(render_path),
    )
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return AncestralStateReportBuildResult(
        output_path=out_path,
        report_kind="ancestral-state",
        title=title,
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        reconstruction_kind=reconstruction_kind,
        model=model,
        machine_manifest=machine_manifest,
    )
