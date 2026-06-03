from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.continuous import ContinuousAncestralReport
from bijux_phylogenetics.ancestral.discrete import DiscreteAncestralReport


@dataclass(slots=True)
class AncestralMethodsSummaryTextResult:
    """Reviewer-facing Markdown methods summary for one ancestral reconstruction."""

    output_path: Path
    title: str
    reconstruction_kind: str
    model: str
    analyzed_taxon_count: int
    unstable_node_count: int
    warning_count: int
    text: str
    reconstruction: ContinuousAncestralReport | DiscreteAncestralReport


def write_ancestral_methods_summary_text(
    path: Path,
    *,
    reconstruction_kind: str,
    reconstruction: ContinuousAncestralReport | DiscreteAncestralReport,
) -> AncestralMethodsSummaryTextResult:
    """Write reviewer-facing methods text for one ancestral reconstruction."""
    text = build_ancestral_methods_summary_text(
        reconstruction_kind=reconstruction_kind,
        reconstruction=reconstruction,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return AncestralMethodsSummaryTextResult(
        output_path=path,
        title="Ancestral Reconstruction Methods Summary",
        reconstruction_kind=reconstruction_kind,
        model=reconstruction.model,
        analyzed_taxon_count=reconstruction.taxon_count,
        unstable_node_count=len(reconstruction.unstable_nodes),
        warning_count=len(reconstruction.warnings),
        text=text,
        reconstruction=reconstruction,
    )


def build_ancestral_methods_summary_text(
    *,
    reconstruction_kind: str,
    reconstruction: ContinuousAncestralReport | DiscreteAncestralReport,
) -> str:
    """Build reviewer-facing methods text for one ancestral reconstruction."""
    if reconstruction_kind == "continuous":
        if not isinstance(reconstruction, ContinuousAncestralReport):
            raise ValueError(
                "continuous ancestral methods summary requires a continuous reconstruction"
            )
        return _build_continuous_methods_summary_text(reconstruction)
    if reconstruction_kind == "discrete":
        if not isinstance(reconstruction, DiscreteAncestralReport):
            raise ValueError(
                "discrete ancestral methods summary requires a discrete reconstruction"
            )
        return _build_discrete_methods_summary_text(reconstruction)
    raise ValueError(
        f"unsupported ancestral methods-summary reconstruction kind: {reconstruction_kind}"
    )


def _build_continuous_methods_summary_text(
    reconstruction: ContinuousAncestralReport,
) -> str:
    internal_estimates = [
        estimate for estimate in reconstruction.estimates if not estimate.is_tip
    ]
    weak_support_count = len(reconstruction.weak_support_nodes)
    root_estimate = _root_internal_estimate(internal_estimates)
    top_interpretation = (
        internal_estimates[0].interpretation if internal_estimates else ""
    )
    downstream_risks = sorted(
        {
            risk
            for estimate in internal_estimates
            for risk in estimate.downstream_risks
            if risk
        }
    )
    diagnostics = reconstruction.brownian_fit_diagnostics
    optimizer = reconstruction.optimizer_diagnostics
    return (
        "# Ancestral Reconstruction Methods Summary\n\n"
        f"This continuous ancestral reconstruction evaluated trait `{reconstruction.trait}` on tree "
        f"`{reconstruction.tree_path.name}` and trait table `{reconstruction.traits_path.name}` using "
        f"model `{reconstruction.model}` with estimator `{reconstruction.estimator}`. The reconstruction retained "
        f"`{reconstruction.taxon_count}` taxa after pruning missing or non-numeric tips, estimated internal-node values "
        f"on the pruned analysis tree, and propagated node uncertainty through standard errors and 95% intervals.\n\n"
        "## Trait And Tree Preparation\n\n"
        f"- analyzed taxon count: `{reconstruction.taxon_count}`\n"
        f"- dropped tips with missing trait values: `{len(reconstruction.dropped_missing_taxa)}`\n"
        f"- dropped tips with non-numeric trait values: `{len(reconstruction.dropped_non_numeric_taxa)}`\n"
        f"- reconstructed internal node count: `{len(internal_estimates)}`\n"
        f"- model alpha setting: `{format(reconstruction.alpha, '.15g')}`\n\n"
        "## Continuous Reconstruction Model\n\n"
        f"- continuous model: `{reconstruction.model}`\n"
        f"- estimator: `{reconstruction.estimator}`\n"
        f"- root estimate: `{format(root_estimate.estimate, '.15g')}`\n"
        f"- root standard error: `{format(root_estimate.standard_error, '.15g')}`\n"
        f"- weak-support internal nodes: `{weak_support_count}`\n"
        f"- unstable internal nodes: `{len(reconstruction.unstable_nodes)}`\n"
        + (
            ""
            if diagnostics is None
            else (
                f"- covariance model: `{diagnostics.covariance_model}`\n"
                f"- tree ultrametric: `{'yes' if diagnostics.tree_is_ultrametric else 'no'}`\n"
                f"- covariance near singular: `{'yes' if diagnostics.covariance_near_singular else 'no'}`\n"
                f"- covariance condition number: `{format(diagnostics.covariance_condition_number, '.15g')}`\n"
                f"- log likelihood: `{format(diagnostics.log_likelihood, '.15g')}`\n"
                f"- residual sigma squared: `{format(diagnostics.residual_sigma_squared, '.15g')}`\n"
            )
        )
        + (
            ""
            if optimizer is None
            else (
                f"- optimizer: `{optimizer.optimizer_name}`\n"
                f"- optimizer converged: `{'yes' if optimizer.converged else 'no'}`\n"
                f"- optimizer iteration count: `{optimizer.iteration_count}`\n"
                f"- optimizer function evaluations: `{optimizer.function_evaluation_count}`\n"
            )
        )
        + "\n## Uncertainty And Node Interpretation\n\n"
        f"- node uncertainty is reported as standard error plus 95% interval for each internal node\n"
        f"- example internal interpretation: {top_interpretation or 'no internal-node interpretation text was recorded'}\n"
        f"- downstream interpretation risks: {_format_item_list(downstream_risks)}\n"
        f"- reconstruction warnings: {_format_item_list(reconstruction.warnings)}\n"
    )


def _build_discrete_methods_summary_text(
    reconstruction: DiscreteAncestralReport,
) -> str:
    internal_estimates = [
        estimate for estimate in reconstruction.estimates if not estimate.is_tip
    ]
    root_estimate = _root_internal_estimate(internal_estimates)
    top_interpretation = (
        internal_estimates[0].interpretation if internal_estimates else ""
    )
    downstream_risks = sorted(
        {
            risk
            for estimate in internal_estimates
            for risk in estimate.downstream_risks
            if risk
        }
    )
    observed_state_text = ", ".join(
        f"`{state}` ({reconstruction.state_counts[state]})"
        for state in reconstruction.observed_states
    )
    return (
        "# Ancestral Reconstruction Methods Summary\n\n"
        f"This discrete ancestral reconstruction evaluated trait `{reconstruction.trait}` on tree "
        f"`{reconstruction.tree_path.name}` and trait table `{reconstruction.traits_path.name}` using "
        f"model `{reconstruction.model}` with state-ordering policy `{reconstruction.state_ordering}`. The reconstruction retained "
        f"`{reconstruction.taxon_count}` taxa after pruning missing tips, inferred marginal internal-node state probabilities on the "
        f"pruned analysis tree, and recorded branchwise transition review evidence from the same fitted model.\n\n"
        "## Trait States And Pruning\n\n"
        f"- analyzed taxon count: `{reconstruction.taxon_count}`\n"
        f"- dropped tips with missing trait values: `{len(reconstruction.dropped_missing_taxa)}`\n"
        f"- observed states: {observed_state_text or 'none'}\n"
        f"- sparse observed states: `{len(reconstruction.sparse_states)}`\n"
        f"- reconstructed internal node count: `{len(internal_estimates)}`\n\n"
        "## Discrete Reconstruction Model\n\n"
        f"- discrete model: `{reconstruction.model}`\n"
        f"- root prior mode: `{reconstruction.root_prior_mode}`\n"
        f"- fixed root state: `{reconstruction.fixed_root_state or 'none'}`\n"
        f"- phytools rerooting-method comparable: `{'yes' if reconstruction.rerooting_method_compatibility.comparable else 'no'}`\n"
        f"- parameter count: `{_format_optional_number(reconstruction.parameter_count)}`\n"
        f"- log likelihood: `{_format_optional_number(reconstruction.log_likelihood)}`\n"
        f"- AIC: `{_format_optional_number(reconstruction.aic)}`\n"
        f"- transition rate count: `{len(reconstruction.transition_rate_rows)}`\n"
        f"- ambiguous internal nodes: `{sum(1 for estimate in internal_estimates if estimate.ambiguous)}`\n"
        f"- unstable internal nodes: `{len(reconstruction.unstable_nodes)}`\n"
        f"- root most-likely state: `{root_estimate.most_likely_state}`\n"
        f"- root confidence: `{format(root_estimate.confidence, '.15g')}`\n"
        + (
            ""
            if reconstruction.optimizer_diagnostics is None
            else (
                f"- optimizer: `{reconstruction.optimizer_diagnostics.optimizer_name}`\n"
                f"- optimizer converged: `{'yes' if reconstruction.optimizer_diagnostics.converged else 'no'}`\n"
                f"- optimizer iteration count: `{reconstruction.optimizer_diagnostics.iteration_count}`\n"
                f"- optimizer function evaluations: `{reconstruction.optimizer_diagnostics.function_evaluation_count}`\n"
            )
        )
        + (
            ""
            if reconstruction.baseline_comparison is None
            else (
                f"- baseline comparison model: `{reconstruction.baseline_comparison.baseline_model}`\n"
                f"- preferred model by AIC: `{reconstruction.baseline_comparison.preferred_model_by_aic}`\n"
                f"- baseline delta AIC: `{format(reconstruction.baseline_comparison.delta_aic, '.15g')}`\n"
            )
        )
        + "\n## Uncertainty And Node Interpretation\n\n"
        "- node uncertainty is reported as marginal state probabilities, confidence, and ambiguity flags for each internal node\n"
        f"- example internal interpretation: {top_interpretation or 'no internal-node interpretation text was recorded'}\n"
        f"- downstream interpretation risks: {_format_item_list(downstream_risks)}\n"
        f"- reconstruction warnings: {_format_item_list(reconstruction.warnings)}\n"
    )


def _format_item_list(items: list[str]) -> str:
    if not items:
        return "none"
    return ", ".join(f"`{item}`" for item in items)


def _format_optional_number(value: float | int | None) -> str:
    if value is None:
        return "none"
    return format(value, ".15g")


def _root_internal_estimate(internal_estimates: list[object]) -> object:
    if not internal_estimates:
        raise ValueError(
            "ancestral methods summary requires at least one internal-node estimate"
        )
    return max(
        internal_estimates,
        key=lambda estimate: (len(estimate.descendant_taxa), estimate.node),
    )
