from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.benchmark import (
    WorkflowPracticalLimitEntry,
    benchmark_workflow_practical_limits,
)


@dataclass(frozen=True, slots=True)
class ProductionScaleThreshold:
    scale: str
    minimum_taxa: int | None
    minimum_sites: int | None
    minimum_tree_count: int | None
    minimum_posterior_size: int | None
    description: str


@dataclass(frozen=True, slots=True)
class WorkflowScaleDecision:
    scale: str
    ready: bool
    limiting_dimensions: list[str]
    rationale: str


@dataclass(slots=True)
class WorkflowProductionScaleReadinessEntry:
    workflow: str
    evidence_source: str
    tested_taxon_limit: int | None
    tested_site_limit: int | None
    tested_tree_limit: int | None
    tested_posterior_size: int | None
    highest_ready_scale: str
    scale_decisions: list[WorkflowScaleDecision]
    notes: list[str]


@dataclass(slots=True)
class ProductionScaleReadinessReport:
    goal_id: int
    replicates: int
    stress_tiers: list[str]
    scale_definitions: list[ProductionScaleThreshold]
    entries: list[WorkflowProductionScaleReadinessEntry]
    limitations: list[str]


_DEFAULT_SCALE_DEFINITIONS: tuple[ProductionScaleThreshold, ...] = (
    ProductionScaleThreshold(
        scale="small",
        minimum_taxa=16,
        minimum_sites=128,
        minimum_tree_count=32,
        minimum_posterior_size=32,
        description="reviewer-facing small datasets or tree sets that exceed toy-example size",
    ),
    ProductionScaleThreshold(
        scale="medium",
        minimum_taxa=128,
        minimum_sites=512,
        minimum_tree_count=128,
        minimum_posterior_size=128,
        description="governed medium-scale inputs that start to reflect realistic review loads",
    ),
    ProductionScaleThreshold(
        scale="large",
        minimum_taxa=512,
        minimum_sites=1024,
        minimum_tree_count=256,
        minimum_posterior_size=256,
        description="large reviewer-facing workloads that should not be treated as incidental local toy runs",
    ),
    ProductionScaleThreshold(
        scale="hpc",
        minimum_taxa=1024,
        minimum_sites=2048,
        minimum_tree_count=1024,
        minimum_posterior_size=1024,
        description="governed high-scale workloads that approach HPC-oriented planning territory",
    ),
)


def _dimension_decision(
    *,
    observed: int | None,
    required: int | None,
    label: str,
) -> str | None:
    if observed is None or required is None:
        return None
    if observed >= required:
        return None
    return f"{label} tested to {observed}, below {required}"


def _evaluate_scale(
    entry: WorkflowPracticalLimitEntry,
    threshold: ProductionScaleThreshold,
) -> WorkflowScaleDecision:
    limiting_dimensions = [
        reason
        for reason in (
            _dimension_decision(
                observed=entry.tested_taxon_limit,
                required=threshold.minimum_taxa,
                label="taxa",
            ),
            _dimension_decision(
                observed=entry.tested_site_limit,
                required=threshold.minimum_sites,
                label="sites",
            ),
            _dimension_decision(
                observed=entry.tested_tree_limit,
                required=threshold.minimum_tree_count,
                label="tree count",
            ),
            _dimension_decision(
                observed=entry.tested_posterior_size,
                required=threshold.minimum_posterior_size,
                label="posterior size",
            ),
        )
        if reason is not None
    ]
    ready = not limiting_dimensions
    if ready:
        rationale = (
            f"{entry.workflow} meets every applicable {threshold.scale} threshold"
        )
    else:
        rationale = "; ".join(limiting_dimensions)
    return WorkflowScaleDecision(
        scale=threshold.scale,
        ready=ready,
        limiting_dimensions=limiting_dimensions,
        rationale=rationale,
    )


def _highest_ready_scale(decisions: list[WorkflowScaleDecision]) -> str:
    ready_scales = [decision.scale for decision in decisions if decision.ready]
    if not ready_scales:
        return "below-small"
    return ready_scales[-1]


def build_production_scale_readiness_report(
    *,
    replicates: int = 1,
    tree_tip_counts: list[int] | None = None,
    alignment_size_classes: list[tuple[str, int, int]] | None = None,
    tree_set_size_classes: list[tuple[str, int, int]] | None = None,
    stress_tiers: list[str] | None = None,
    scale_definitions: list[ProductionScaleThreshold] | None = None,
) -> ProductionScaleReadinessReport:
    """Classify which owned workflows are currently supported at small, medium, large, and HPC scale."""
    thresholds = list(scale_definitions or _DEFAULT_SCALE_DEFINITIONS)
    if not thresholds:
        raise ValueError("scale_definitions must contain at least one scale threshold")
    practical_limits = benchmark_workflow_practical_limits(
        replicates=replicates,
        tree_tip_counts=tree_tip_counts,
        alignment_size_classes=alignment_size_classes,
        tree_set_size_classes=tree_set_size_classes,
        stress_tiers=stress_tiers,
    )
    entries: list[WorkflowProductionScaleReadinessEntry] = []
    for entry in practical_limits.entries:
        decisions = [_evaluate_scale(entry, threshold) for threshold in thresholds]
        entries.append(
            WorkflowProductionScaleReadinessEntry(
                workflow=entry.workflow,
                evidence_source=entry.evidence_source,
                tested_taxon_limit=entry.tested_taxon_limit,
                tested_site_limit=entry.tested_site_limit,
                tested_tree_limit=entry.tested_tree_limit,
                tested_posterior_size=entry.tested_posterior_size,
                highest_ready_scale=_highest_ready_scale(decisions),
                scale_decisions=decisions,
                notes=entry.notes,
            )
        )
    limitations = list(practical_limits.limitations)
    limitations.append(
        "production-scale readiness is derived from governed tested maxima and should be read as evidence-backed support levels, not guarantees for larger or different hardware"
    )
    return ProductionScaleReadinessReport(
        goal_id=225,
        replicates=practical_limits.replicates,
        stress_tiers=practical_limits.stress_tiers,
        scale_definitions=thresholds,
        entries=entries,
        limitations=limitations,
    )
