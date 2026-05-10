from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    load_continuous_dataset,
    reconstruction_manifest,
    write_ancestral_rows,
)
from bijux_phylogenetics.ancestral.continuous import (
    ContinuousAncestralReport,
    _reconstruct_continuous_from_dataset,
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.comparative.evolutionary_modes import (
    transform_tree_for_evolutionary_mode,
)
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
class ContinuousAncestralTreeComparisonRow:
    """Node-wise comparison of continuous ancestral estimates across two trees."""

    node: str
    descendant_taxa: list[str]
    left_estimate: float
    right_estimate: float
    estimate_delta: float


@dataclass(slots=True)
class ContinuousAncestralTreeComparisonReport:
    """Continuous ancestral reconstruction differences across alternative trees."""

    left_tree_path: Path
    right_tree_path: Path
    traits_path: Path
    trait: str
    model: str
    rows: list[ContinuousAncestralTreeComparisonRow]


@dataclass(slots=True)
class DiscreteAncestralTreeComparisonRow:
    """Node-wise comparison of discrete ancestral states across two trees."""

    node: str
    descendant_taxa: list[str]
    left_state: str
    right_state: str
    differs: bool
    left_confidence: float
    right_confidence: float


@dataclass(slots=True)
class DiscreteAncestralTreeComparisonReport:
    """Discrete ancestral reconstruction differences across alternative trees."""

    left_tree_path: Path
    right_tree_path: Path
    traits_path: Path
    trait: str
    model: str
    rows: list[DiscreteAncestralTreeComparisonRow]


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
    supplement_sections: list[str]
    sensitivity: object | None
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class ContinuousEvolutionaryModeAncestralReport:
    """Continuous ancestral reconstruction over a governed transformed-tree mode."""

    tree_path: Path
    traits_path: Path
    trait: str
    mode: str
    parameter_name: str | None
    parameter_value: float | None
    transformed_tree_newick: str
    reconstruction: ContinuousAncestralReport


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
                    or right_estimate.upper_95_interval
                    < left_estimate.lower_95_interval
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


def reconstruct_continuous_evolutionary_mode_states(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    mode: str = "brownian",
    rate_change: float = 0.0,
) -> ContinuousEvolutionaryModeAncestralReport:
    """Reconstruct continuous ancestral states under Brownian or early-burst tree modes."""
    if mode == "brownian":
        reconstruction = reconstruct_continuous_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model="brownian",
        )
        return ContinuousEvolutionaryModeAncestralReport(
            tree_path=tree_path,
            traits_path=traits_path,
            trait=trait,
            mode="brownian",
            parameter_name=None,
            parameter_value=None,
            transformed_tree_newick=reconstruction.analysis_tree_newick,
            reconstruction=reconstruction,
        )
    if mode != "early-burst":
        raise ValueError(
            "unsupported evolutionary mode for continuous ancestral reconstruction"
        )
    dataset = load_continuous_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    transformed_tree = transform_tree_for_evolutionary_mode(
        dataset.tree,
        mode="early-burst",
        parameter_value=rate_change,
    )
    reconstruction = _reconstruct_continuous_from_dataset(
        dataset,
        working_tree=transformed_tree,
        model="brownian",
        alpha=1.0,
    )
    return ContinuousEvolutionaryModeAncestralReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        mode="early-burst",
        parameter_name="rate_change",
        parameter_value=rate_change,
        transformed_tree_newick=reconstruction.analysis_tree_newick,
        reconstruction=reconstruction,
    )


def compare_discrete_ancestral_models(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    models: tuple[str, str, str] = ("equal-rates", "symmetric", "all-rates-different"),
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> DiscreteAncestralModelComparisonReport:
    """Compare supported discrete ancestral-state likelihood models and select the lowest-AIC fit."""
    reconstructions = {
        model: reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
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
            state_ordering=state_ordering,
            ordered_states=ordered_states,
        )
        for model in models
    }
    rows = [
        DiscreteAncestralModelComparisonRow(
            model=model,
            parameter_count=transition_reports[model].transition_model.parameter_count,
            pseudo_log_likelihood=transition_reports[
                model
            ].transition_model.pseudo_log_likelihood,
            aic=transition_reports[model].transition_model.aic,
            selected=False,
        )
        for model in models
    ]
    selected_model = min(rows, key=lambda row: row.aic).model
    for row in rows:
        row.selected = row.model == selected_model
    selected_report = reconstructions[selected_model]
    selected_by_node = {
        estimate.node: estimate for estimate in selected_report.estimates
    }
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
                    differs=selected_estimate.most_likely_state
                    != estimate.most_likely_state,
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


def compare_continuous_ancestral_trees(
    left_tree_path: Path,
    right_tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "brownian",
    alpha: float = 1.0,
) -> ContinuousAncestralTreeComparisonReport:
    """Compare continuous ancestral estimates across two alternative trees."""
    left = reconstruct_continuous_ancestral_states(
        left_tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        alpha=alpha,
    )
    right = reconstruct_continuous_ancestral_states(
        right_tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        alpha=alpha,
    )
    right_by_node = {
        estimate.node: estimate for estimate in right.estimates if not estimate.is_tip
    }
    rows = [
        ContinuousAncestralTreeComparisonRow(
            node=estimate.node,
            descendant_taxa=estimate.descendant_taxa,
            left_estimate=estimate.estimate,
            right_estimate=right_by_node[estimate.node].estimate,
            estimate_delta=right_by_node[estimate.node].estimate - estimate.estimate,
        )
        for estimate in left.estimates
        if not estimate.is_tip and estimate.node in right_by_node
    ]
    return ContinuousAncestralTreeComparisonReport(
        left_tree_path=left_tree_path,
        right_tree_path=right_tree_path,
        traits_path=traits_path,
        trait=trait,
        model=model,
        rows=rows,
    )


def compare_discrete_ancestral_trees(
    left_tree_path: Path,
    right_tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "fitch",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> DiscreteAncestralTreeComparisonReport:
    """Compare discrete ancestral states across two alternative trees."""
    left = reconstruct_discrete_ancestral_states(
        left_tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    right = reconstruct_discrete_ancestral_states(
        right_tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    right_by_node = {
        estimate.node: estimate for estimate in right.estimates if not estimate.is_tip
    }
    rows = [
        DiscreteAncestralTreeComparisonRow(
            node=estimate.node,
            descendant_taxa=estimate.descendant_taxa,
            left_state=estimate.most_likely_state,
            right_state=right_by_node[estimate.node].most_likely_state,
            differs=estimate.most_likely_state
            != right_by_node[estimate.node].most_likely_state,
            left_confidence=estimate.confidence,
            right_confidence=right_by_node[estimate.node].confidence,
        )
        for estimate in left.estimates
        if not estimate.is_tip and estimate.node in right_by_node
    ]
    return DiscreteAncestralTreeComparisonReport(
        left_tree_path=left_tree_path,
        right_tree_path=right_tree_path,
        traits_path=traits_path,
        trait=trait,
        model=model,
        rows=rows,
    )


def write_ancestral_state_table(
    path: Path, report: ContinuousAncestralReport | DiscreteAncestralReport
) -> Path:
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
            "state_probabilities": json.dumps(
                estimate.state_probabilities, sort_keys=True
            ),
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
        internal_annotation_colors = {
            estimate.node: "#6d28d9"
            for estimate in report.estimates
            if not estimate.is_tip
        }
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
        estimate.node: estimate.most_likely_state
        if not estimate.ambiguous
        else "/".join(estimate.state_set)
        for estimate in report.estimates
        if not estimate.is_tip
    }
    discrete_tip_states = {
        estimate.node_name: estimate.most_likely_state
        for estimate in report.estimates
        if estimate.is_tip and estimate.node_name is not None
    }
    palette = dict(
        zip(
            sorted(report.observed_states),
            ("#0f766e", "#1d4ed8", "#c2410c", "#7c3aed", "#b91c1c", "#047857"),
            strict=False,
        )
    )
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
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    compare_tree_path: Path | None = None,
    drop_taxa: list[str] | None = None,
    coding_map: dict[str, str] | None = None,
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
            state_ordering=state_ordering,
            ordered_states=ordered_states,
        )
        comparison = (
            compare_discrete_ancestral_models(
                tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
                models=(model, compare_model, "all-rates-different"),
                state_ordering=state_ordering,
                ordered_states=ordered_states,
            )
            if compare_model is not None and model != "fitch"
            else None
        )
    from bijux_phylogenetics.ancestral.sensitivity import (
        build_ancestral_sensitivity_report,
    )

    sensitivity = build_ancestral_sensitivity_report(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        reconstruction_kind=reconstruction_kind,
        model=model,
        taxon_column=taxon_column,
        alpha=alpha,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
        compare_tree_path=compare_tree_path,
        compare_model=compare_model,
        drop_taxa=drop_taxa,
        coding_map=coding_map,
    )
    render_path = out_path.with_suffix(".svg")
    render_result = render_ancestral_state_tree(
        tree_path, report, out_path=render_path, layout="phylogram"
    )
    supplement_sections = [
        "ancestral-methods",
        "ancestral-exclusions",
        "ancestral-node-table",
        "ancestral-uncertainty",
        "ancestral-sensitivity",
    ]
    sections = [
        (
            "ancestral-reconstruction",
            json.dumps(asdict(report), indent=2, sort_keys=True, default=str),
        ),
        (
            "ancestral-render",
            json.dumps(asdict(render_result), indent=2, sort_keys=True, default=str),
        ),
        (
            "ancestral-methods",
            json.dumps(
                _report_methods(report, reconstruction_kind=reconstruction_kind),
                indent=2,
                sort_keys=True,
                default=str,
            ),
        ),
        (
            "ancestral-exclusions",
            json.dumps(
                _report_exclusions(report), indent=2, sort_keys=True, default=str
            ),
        ),
        (
            "ancestral-node-table",
            json.dumps(
                _report_node_table(report), indent=2, sort_keys=True, default=str
            ),
        ),
        (
            "ancestral-uncertainty",
            json.dumps(
                _report_uncertainty(report), indent=2, sort_keys=True, default=str
            ),
        ),
        (
            "ancestral-sensitivity",
            json.dumps(asdict(sensitivity), indent=2, sort_keys=True, default=str),
        ),
    ]
    if comparison is not None:
        sections.append(
            (
                "ancestral-comparison",
                json.dumps(asdict(comparison), indent=2, sort_keys=True, default=str),
            )
        )
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
    machine_manifest["supplement_sections"] = supplement_sections
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return AncestralStateReportBuildResult(
        output_path=out_path,
        report_kind="ancestral-state",
        title=title,
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        reconstruction_kind=reconstruction_kind,
        model=model,
        supplement_sections=supplement_sections,
        sensitivity=sensitivity,
        machine_manifest=machine_manifest,
    )


def _report_methods(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
    *,
    reconstruction_kind: str,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "reconstruction_kind": reconstruction_kind,
        "model": report.model,
        "trait": report.trait,
        "taxon_count": report.taxon_count,
        "warning_count": len(report.warnings),
    }
    if isinstance(report, ContinuousAncestralReport):
        payload["alpha"] = report.alpha
        payload["supports_explicit_models"] = ["brownian", "ou"]
    else:
        payload["observed_states"] = report.observed_states
        payload["supports_explicit_models"] = [
            "fitch",
            "equal-rates",
            "symmetric",
            "all-rates-different",
        ]
    return payload


def _report_exclusions(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "dropped_missing_taxa": report.dropped_missing_taxa,
        "warnings": report.warnings,
    }
    if isinstance(report, ContinuousAncestralReport):
        payload["dropped_non_numeric_taxa"] = report.dropped_non_numeric_taxa
    return payload


def _report_node_table(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> list[dict[str, object]]:
    if isinstance(report, ContinuousAncestralReport):
        return [
            {
                "node": estimate.node,
                "descendant_taxa": estimate.descendant_taxa,
                "estimate": estimate.estimate,
                "confidence": estimate.confidence,
                "interpretation": estimate.interpretation,
            }
            for estimate in report.estimates
            if not estimate.is_tip
        ]
    return [
        {
            "node": estimate.node,
            "descendant_taxa": estimate.descendant_taxa,
            "state": estimate.most_likely_state,
            "confidence": estimate.confidence,
            "unstable": estimate.unstable,
            "interpretation": estimate.interpretation,
        }
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def _report_uncertainty(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> dict[str, object]:
    if isinstance(report, ContinuousAncestralReport):
        return {
            "unstable_nodes": report.unstable_nodes,
            "weak_support_nodes": report.weak_support_nodes,
            "interval_rows": [
                {
                    "node": estimate.node,
                    "lower_95_interval": estimate.lower_95_interval,
                    "upper_95_interval": estimate.upper_95_interval,
                    "downstream_risks": estimate.downstream_risks,
                }
                for estimate in report.estimates
                if not estimate.is_tip
            ],
        }
    return {
        "unstable_nodes": report.unstable_nodes,
        "weak_support_nodes": report.weak_support_nodes,
        "probability_rows": [
            {
                "node": estimate.node,
                "state_probabilities": estimate.state_probabilities,
                "downstream_risks": estimate.downstream_risks,
            }
            for estimate in report.estimates
            if not estimate.is_tip
        ],
    }
