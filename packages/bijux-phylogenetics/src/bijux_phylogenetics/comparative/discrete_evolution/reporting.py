from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import tempfile

from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.render.tree_svg import TreeRenderResult, render_tree_svg

from .analysis import (
    assess_geographic_state_analysis_readiness,
    estimate_ancestral_geographic_states,
    run_discrete_state_transition_model,
)
from .comparison import compare_discrete_state_models
from .models import (
    BiogeographicComputedResult,
    BiogeographicInterpretationReport,
    DiscreteEvolutionNarrative,
    DiscreteEvolutionReportBuildResult,
    DiscreteModelComparisonReport,
    DiscreteStateEvolutionReport,
    DiscreteTransitionReferenceObservation,
    DiscreteTransitionReferenceRate,
    DiscreteTransitionReferenceValidationReport,
)
from .palette import _DEFAULT_STATE_COLORS
from .state_coding import audit_discrete_state_coding


def _build_discrete_evolution_narrative(
    report: DiscreteStateEvolutionReport,
    *,
    comparison: DiscreteModelComparisonReport | None = None,
) -> DiscreteEvolutionNarrative:
    root_state = next(
        estimate.most_likely_state
        for estimate in report.estimates
        if estimate.node == node_signature(load_tree(report.tree_path).root)
    )
    summary = (
        f"reconstructed {report.trait} across {report.taxon_count} taxa under the {report.model} model"
        f" with root state '{root_state}' and {report.transition_summary.transition_count} inferred branch transitions"
    )
    strong = report.transition_summary.strongly_supported_transition_count
    transition_summary = (
        f"{strong} of {report.transition_summary.transition_count} changed branches remain strongly supported after"
        " weighting by parent-child state probabilities"
    )
    caveats = list(report.warnings)
    if comparison is not None and comparison.sensitive_region_count > 0:
        caveats.append(
            f"{comparison.sensitive_region_count} internal regions change reconstructed state between {comparison.left_model}"
            f" and {comparison.right_model}"
        )
    caveats.extend(
        [
            "deterministic node probabilities are an approximation and do not replace full Bayesian or likelihood-marginal ancestral mapping",
            "transition uncertainty intervals summarize effective event counts in the fitted deterministic model and are not posterior credible intervals",
        ]
    )
    return DiscreteEvolutionNarrative(
        summary=summary,
        transition_summary=transition_summary,
        interpretation_boundary=(
            "treat these outputs as computational evidence about state histories, not as direct proof of dispersal timing, mechanism, or causal biogeography"
        ),
        caveats=caveats,
    )


def _report_limitations(
    narrative: DiscreteEvolutionNarrative,
    interpretation: BiogeographicInterpretationReport,
) -> list[str]:
    limitations = [
        narrative.interpretation_boundary,
        *narrative.caveats,
        *interpretation.readiness_blockers,
        *interpretation.caveats,
    ]
    return sorted(dict.fromkeys(item.strip() for item in limitations if item.strip()))


def build_biogeographic_interpretation_report(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "equal-rates",
    compare_model: str | None = None,
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    coding_map: dict[str, str] | None = None,
) -> BiogeographicInterpretationReport:
    """Separate computed biogeographic outputs from downstream interpretation guidance."""
    readiness = assess_geographic_state_analysis_readiness(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    report = estimate_ancestral_geographic_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    comparison = (
        compare_discrete_state_models(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            left_model=model,
            right_model=compare_model,
            allowed_states=allowed_states,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
        )
        if compare_model is not None
        else None
    )
    coding_audit = audit_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
        coding_map=coding_map,
    )
    root_estimate = next(
        estimate
        for estimate in report.estimates
        if estimate.node == node_signature(load_tree(tree_path).root)
    )
    computed_results = [
        BiogeographicComputedResult(
            label="root_state", value=root_estimate.most_likely_state
        ),
        BiogeographicComputedResult(
            label="transition_count",
            value=str(report.transition_summary.transition_count),
        ),
        BiogeographicComputedResult(
            label="strongly_supported_transition_count",
            value=str(report.transition_summary.strongly_supported_transition_count),
        ),
        BiogeographicComputedResult(
            label="state_count", value=str(len(report.observed_states))
        ),
        BiogeographicComputedResult(label="model", value=report.model),
    ]
    return BiogeographicInterpretationReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        model=model,
        compare_model=compare_model,
        computed_results=computed_results,
        model_sensitive_regions=[]
        if comparison is None
        else comparison.sensitive_regions,
        coding_audit_summary={
            "row_count": coding_audit.row_count,
            "included_row_count": coding_audit.included_row_count,
            "excluded_row_count": coding_audit.excluded_row_count,
        },
        readiness_blockers=readiness.blockers,
        caveats=_build_discrete_evolution_narrative(
            report, comparison=comparison
        ).caveats,
        interpretation_guidance=[
            "computed ancestral regions summarize model-conditioned state histories, not direct evidence of dispersal mechanism",
            "biological interpretation should be restricted to patterns that remain stable across supported models and coding assumptions",
            "sampling gaps, sparse states, and dominant-state bias should be discussed before turning transitions into historical claims",
        ],
    )


def validate_discrete_transition_reference_examples(
    *,
    tolerance: float = 1e-9,
) -> DiscreteTransitionReferenceValidationReport:
    """Validate deterministic discrete-state transition outputs against built-in small references."""
    cases = (
        {
            "label": "toy-geography-er",
            "model": "equal-rates",
            "tree_newick": "((A:0.1,B:0.1):0.2,(C:0.2,D:0.2):0.1);",
            "trait_rows": [
                ("A", "north"),
                ("B", "north"),
                ("C", "south"),
                ("D", "island"),
            ],
            "expected_parameter_count": 1,
            "expected_transition_count": 2,
            "expected_root_state": "north",
            "expected_pseudo_log_likelihood": -7.28950386047299,
            "expected_rates": {
                ("island", "north"): 0.1,
                ("north", "south"): 0.1,
                ("south", "south"): 0.8,
            },
        },
        {
            "label": "toy-geography-sym",
            "model": "symmetric",
            "tree_newick": "((A:0.1,B:0.1):0.2,(C:0.2,D:0.2):0.1);",
            "trait_rows": [
                ("A", "north"),
                ("B", "north"),
                ("C", "south"),
                ("D", "island"),
            ],
            "expected_parameter_count": 3,
            "expected_transition_count": 2,
            "expected_root_state": "north",
            "expected_pseudo_log_likelihood": -6.5047894874944,
            "expected_rates": {
                ("island", "south"): 0.133333333333333,
                ("north", "south"): 0.133333333333333,
                ("south", "island"): 0.1,
            },
        },
        {
            "label": "toy-geography-ard",
            "model": "all-rates-different",
            "tree_newick": "((A:0.1,B:0.1):0.2,(C:0.2,D:0.2):0.1);",
            "trait_rows": [
                ("A", "north"),
                ("B", "north"),
                ("C", "south"),
                ("D", "island"),
            ],
            "expected_parameter_count": 6,
            "expected_transition_count": 2,
            "expected_root_state": "north",
            "expected_pseudo_log_likelihood": -6.2424252230423,
            "expected_rates": {
                ("island", "north"): 0.1,
                ("north", "south"): 0.133333333333333,
                ("south", "island"): 0.133333333333333,
            },
        },
    )
    observations: list[DiscreteTransitionReferenceObservation] = []
    for case in cases:
        with tempfile.TemporaryDirectory(prefix="bijux-discrete-reference-") as tmp_dir:
            tree_path = Path(tmp_dir) / "reference-tree.nwk"
            traits_path = Path(tmp_dir) / "reference-traits.tsv"
            tree_path.write_text(f"{case['tree_newick']}\n", encoding="utf-8")
            traits_path.write_text(
                "taxon\tregion\n"
                + "".join(f"{taxon}\t{state}\n" for taxon, state in case["trait_rows"]),
                encoding="utf-8",
            )
            report = run_discrete_state_transition_model(
                tree_path,
                traits_path,
                trait="region",
                model=str(case["model"]),
            )
        rate_lookup = {
            (row.source_state, target): value
            for row in report.transition_model.transition_matrix
            for target, value in row.target_rates.items()
        }
        rate_rows = [
            DiscreteTransitionReferenceRate(
                source_state=source_state,
                target_state=target_state,
                expected_rate=expected_rate,
                observed_rate=rate_lookup[(source_state, target_state)],
                absolute_delta=abs(
                    rate_lookup[(source_state, target_state)] - expected_rate
                ),
            )
            for (source_state, target_state), expected_rate in sorted(
                case["expected_rates"].items()
            )
        ]
        max_rate_delta = max((row.absolute_delta for row in rate_rows), default=0.0)
        root_state = next(
            estimate.most_likely_state
            for estimate in report.estimates
            if estimate.node == "A|B|C|D"
        )
        passed = (
            report.transition_model.parameter_count == case["expected_parameter_count"]
            and report.transition_summary.transition_count
            == case["expected_transition_count"]
            and root_state == case["expected_root_state"]
            and abs(
                report.transition_model.pseudo_log_likelihood
                - case["expected_pseudo_log_likelihood"]
            )
            <= tolerance
            and max_rate_delta <= tolerance
        )
        observations.append(
            DiscreteTransitionReferenceObservation(
                label=str(case["label"]),
                model=str(case["model"]),
                expected_parameter_count=int(case["expected_parameter_count"]),
                observed_parameter_count=report.transition_model.parameter_count,
                expected_transition_count=int(case["expected_transition_count"]),
                observed_transition_count=report.transition_summary.transition_count,
                expected_root_state=str(case["expected_root_state"]),
                observed_root_state=root_state,
                expected_pseudo_log_likelihood=float(
                    case["expected_pseudo_log_likelihood"]
                ),
                observed_pseudo_log_likelihood=report.transition_model.pseudo_log_likelihood,
                max_rate_delta=max_rate_delta,
                rate_rows=rate_rows,
                passed=passed,
            )
        )
    return DiscreteTransitionReferenceValidationReport(
        case_count=len(observations),
        all_passed=all(observation.passed for observation in observations),
        tolerance=tolerance,
        observations=observations,
    )


def write_node_state_probability_table(
    path: Path, report: DiscreteStateEvolutionReport
) -> Path:
    """Export one deterministic node-probability table for a discrete-state reconstruction."""
    rows = [
        {
            "node": estimate.node,
            "node_name": estimate.node_name or "",
            "is_tip": str(estimate.is_tip).lower(),
            "descendant_taxa": ",".join(estimate.descendant_taxa),
            "most_likely_state": estimate.most_likely_state,
            "state_probabilities": json.dumps(
                estimate.state_probabilities, sort_keys=True
            ),
            "ambiguous": str(estimate.ambiguous).lower(),
        }
        for estimate in report.estimates
    ]
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "is_tip",
            "descendant_taxa",
            "most_likely_state",
            "state_probabilities",
            "ambiguous",
        ],
        rows=rows,
    )


def write_transition_summary_table(
    path: Path, report: DiscreteStateEvolutionReport
) -> Path:
    """Export one branch-by-branch transition summary table."""
    support_by_branch = {
        (row.parent_node, row.child_node): row
        for row in report.transition_summary.support_rows
    }
    rows = [
        {
            "parent_node": event.parent_node,
            "child_node": event.child_node,
            "source_state": event.source_state,
            "target_state": event.target_state,
            "changed": str(event.changed).lower(),
            "support": support_by_branch[(event.parent_node, event.child_node)].support,
            "strongly_supported": str(
                support_by_branch[
                    (event.parent_node, event.child_node)
                ].strongly_supported
            ).lower(),
        }
        for event in report.transition_summary.events
    ]
    return write_taxon_rows(
        path,
        columns=[
            "parent_node",
            "child_node",
            "source_state",
            "target_state",
            "changed",
            "support",
            "strongly_supported",
        ],
        rows=rows,
    )


def write_discrete_model_comparison_table(
    path: Path, report: DiscreteModelComparisonReport
) -> Path:
    """Export one node-wise comparison table across two discrete-state models."""
    rows = [
        {
            "node": difference.node,
            "descendant_taxa": ",".join(difference.descendant_taxa),
            "left_state": difference.left_state,
            "right_state": difference.right_state,
            "differs": str(difference.differs).lower(),
            "left_probabilities": json.dumps(
                difference.left_probabilities, sort_keys=True
            ),
            "right_probabilities": json.dumps(
                difference.right_probabilities, sort_keys=True
            ),
        }
        for difference in report.node_differences
    ]
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "descendant_taxa",
            "left_state",
            "right_state",
            "differs",
            "left_probabilities",
            "right_probabilities",
        ],
        rows=rows,
    )


def render_tree_with_geographic_states(
    tree_path: Path,
    report: DiscreteStateEvolutionReport,
    *,
    out_path: Path,
    layout: str = "phylogram",
    state_colors: dict[str, str] | None = None,
) -> TreeRenderResult:
    """Render tip and internal discrete states onto a tree SVG."""
    palette = (
        dict(state_colors)
        if state_colors is not None
        else dict(
            zip(sorted(report.observed_states), _DEFAULT_STATE_COLORS, strict=False)
        )
    )
    internal_pies = {
        estimate.node: dict(estimate.state_probabilities)
        for estimate in report.estimates
        if not estimate.is_tip
    }
    internal_annotations = {
        estimate.node: f"{estimate.state_probabilities.get(estimate.most_likely_state, 0.0):.2f}"
        for estimate in report.estimates
        if not estimate.is_tip
    }
    internal_annotation_colors = {
        estimate.node: palette.get(estimate.most_likely_state, "#6d28d9")
        for estimate in report.estimates
        if not estimate.is_tip
    }
    categorical_traits = {
        estimate.node_name: estimate.most_likely_state
        for estimate in report.estimates
        if estimate.is_tip and estimate.node_name is not None
    }
    return render_tree_svg(
        tree_path,
        out_path=out_path,
        layout=layout,
        categorical_traits=categorical_traits,
        internal_annotations=internal_annotations,
        internal_annotation_colors=internal_annotation_colors,
        internal_pies=internal_pies,
        internal_pie_colors=palette,
    )


def render_discrete_state_evolution_report(
    *,
    tree_path: Path,
    traits_path: Path,
    trait: str,
    out_path: Path,
    taxon_column: str | None = None,
    model: str = "equal-rates",
    allowed_states: list[str] | None = None,
    compare_model: str | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> DiscreteEvolutionReportBuildResult:
    """Build a deterministic HTML report for one discrete-state evolution analysis."""
    report = run_discrete_state_transition_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    comparison = (
        compare_discrete_state_models(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            left_model=model,
            right_model=compare_model,
            allowed_states=allowed_states,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
        )
        if compare_model is not None
        else None
    )
    render_path = out_path.with_suffix(".svg")
    render_result = render_tree_with_geographic_states(
        tree_path, report, out_path=render_path, layout="phylogram"
    )
    narrative = _build_discrete_evolution_narrative(report, comparison=comparison)
    interpretation = build_biogeographic_interpretation_report(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        compare_model=compare_model,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    limitations = _report_limitations(narrative, interpretation)
    sections = [
        (
            "discrete-state-summary",
            json.dumps(asdict(narrative), default=str, indent=2, sort_keys=True),
        ),
        (
            "discrete-state-evolution",
            json.dumps(asdict(report), default=str, indent=2, sort_keys=True),
        ),
        (
            "biogeographic-interpretation",
            json.dumps(asdict(interpretation), default=str, indent=2, sort_keys=True),
        ),
        (
            "discrete-state-render",
            json.dumps(asdict(render_result), default=str, indent=2, sort_keys=True),
        ),
        ("limitations", json.dumps(limitations, indent=2, sort_keys=True)),
    ]
    if comparison is not None:
        sections.append(
            (
                "discrete-state-comparison",
                json.dumps(asdict(comparison), default=str, indent=2, sort_keys=True),
            )
        )
    title = f"Bijux Discrete-State Evolution Report: {trait}"
    machine_manifest = {
        "report_kind": "discrete-state-evolution",
        "title": title,
        "tree_path": str(tree_path),
        "traits_path": str(traits_path),
        "trait": trait,
        "model": model,
        "likelihood_method": report.likelihood_method,
        "state_ordering": report.state_ordering,
        "ordered_states": report.ordered_states,
        "caveat_count": len(narrative.caveats),
        "limitations": limitations,
        "interpretation_sections": [
            "computed_results",
            "caveats",
            "interpretation_guidance",
        ],
        "rendered_tree": str(render_path),
        "sections": [name for name, _ in sections],
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return DiscreteEvolutionReportBuildResult(
        output_path=out_path,
        report_kind="discrete-state-evolution",
        title=title,
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        model=model,
        machine_manifest=machine_manifest,
    )
