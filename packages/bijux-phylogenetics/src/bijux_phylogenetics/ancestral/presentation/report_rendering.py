from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    reconstruction_manifest,
    write_ancestral_rows,
)
from bijux_phylogenetics.ancestral.comparison import (
    compare_continuous_ancestral_models,
    compare_discrete_ancestral_reconstructions,
)
from bijux_phylogenetics.ancestral.continuous import (
    ContinuousAncestralReport,
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.ancestral.presentation.visualization import (
    render_ancestral_state_visualization,
)
from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.render.tree_svg import TreeRenderResult


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
    visualization = render_ancestral_state_visualization(
        tree_path,
        report,
        out_path=out_path,
        layout=layout,
        discrete_node_style="labels",
        branch_coloring="none",
    )
    return visualization.tree_render


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
            compare_discrete_ancestral_reconstructions(
                tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
                left_model=model,
                right_model=compare_model,
                state_ordering=state_ordering,
                ordered_states=ordered_states,
            )
            if compare_model is not None
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
        "limitations",
    ]
    limitations = _report_limitations(
        report,
        reconstruction_kind=reconstruction_kind,
        sensitivity=sensitivity,
        comparison=comparison,
    )
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
        ("limitations", json.dumps(limitations, indent=2, sort_keys=True)),
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
    machine_manifest["limitations"] = limitations
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


def _report_limitations(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
    *,
    reconstruction_kind: str,
    sensitivity: object,
    comparison: object | None,
) -> list[str]:
    limitations = list(report.warnings)
    if reconstruction_kind == "continuous":
        limitations.append(
            "continuous ancestral estimates are model-based internal-node expectations and should not be treated as direct ancestral measurements"
        )
    else:
        limitations.append(
            "discrete ancestral state labels and probabilities summarize model support and should not be treated as direct proof of historical state occupancy"
        )
    if report.dropped_missing_taxa:
        limitations.append(
            f"{len(report.dropped_missing_taxa)} taxa were excluded because the requested trait is missing"
        )
    if (
        isinstance(report, ContinuousAncestralReport)
        and report.dropped_non_numeric_taxa
    ):
        limitations.append(
            f"{len(report.dropped_non_numeric_taxa)} taxa were excluded because the requested continuous trait is not numeric"
        )
    if comparison is not None:
        limitations.append(
            "alternative model comparison is included because internal-node reconstructions can change under different evolutionary assumptions"
        )
    for summary in (
        sensitivity.model_sensitivity,
        sensitivity.tree_sensitivity,
        sensitivity.pruning_sensitivity,
        sensitivity.trait_coding_sensitivity,
    ):
        if summary is None:
            continue
        limitations.extend(summary.notes)
    return sorted(dict.fromkeys(item.strip() for item in limitations if item.strip()))


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
