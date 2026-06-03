from __future__ import annotations

from pathlib import Path

from .benchmark_dashboards import build_method_accuracy_dashboard
from .contracts import (
    MethodLimitationEntry,
    MethodLimitationRegistry,
    ScientificValidationClaim,
    ScientificValidationReport,
)
from .dataset_corpora import default_fixtures_root


def build_method_limitation_registry() -> MethodLimitationRegistry:
    """Enumerate major method families with explicit assumptions and trust boundaries."""
    return MethodLimitationRegistry(
        goal_id=250,
        entries=[
            MethodLimitationEntry(
                method="tree-validation",
                status="validated",
                validated_by=["level1-reference-validation", "broken-benchmark-corpus"],
                assumptions=[
                    "input tree syntax is parseable",
                    "tip labels are intended to identify taxa",
                ],
                invalid_inputs=[
                    "duplicate tip labels",
                    "negative branch lengths",
                    "malformed rootedness assumptions",
                ],
                limitations=[
                    "validation does not by itself prove biological interpretation or engine compatibility beyond audited surfaces"
                ],
            ),
            MethodLimitationEntry(
                method="dataset-audit",
                status="validated",
                validated_by=[
                    "clean-benchmark-corpus",
                    "messy-benchmark-corpus",
                    "regression-dataset-corpus",
                ],
                assumptions=[
                    "metadata and trait tables use stable taxon keys",
                    "caller-provided files correspond to one biological dataset",
                ],
                invalid_inputs=[
                    "missing tree taxa in metadata or traits",
                    "unsafe alignments for inference",
                    "invalid tip dates or calibrations",
                ],
                limitations=[
                    "dataset readiness is reviewer-facing triage, not a substitute for downstream method-specific validation"
                ],
            ),
            MethodLimitationEntry(
                method="distance-methods",
                status="validated",
                validated_by=[
                    "distance reference fixtures",
                    "runtime benchmark dashboard",
                ],
                assumptions=[
                    "aligned homologous sequences",
                    "distance model matches alphabet and saturation regime",
                ],
                invalid_inputs=[
                    "too few comparable sites",
                    "severe saturation",
                    "unsafe ambiguity handling for the chosen policy",
                ],
                limitations=[
                    "distance trees remain approximations and can disagree with likelihood or Bayesian inference"
                ],
            ),
            MethodLimitationEntry(
                method="comparative-models",
                status="validated",
                validated_by=[
                    "comparative reference fixtures",
                    "level1-reference-validation",
                ],
                assumptions=[
                    "tree/taxon linkage is correct",
                    "traits satisfy declared model assumptions",
                ],
                invalid_inputs=[
                    "overfit models",
                    "non-identifiable OU settings",
                    "missing trait coverage after pruning",
                ],
                limitations=[
                    "validated examples do not remove the need for biological judgment about causality or model adequacy"
                ],
            ),
            MethodLimitationEntry(
                method="ancestral-reconstruction",
                status="experimental",
                validated_by=[
                    "simulation validation surfaces",
                    "ancestral reference examples",
                ],
                assumptions=[
                    "chosen transition or continuous model is defensible",
                    "tree uncertainty has been considered",
                ],
                invalid_inputs=[
                    "impossible discrete coding",
                    "low-information or unstable internal nodes",
                    "unsafe pruning sensitivity",
                ],
                limitations=[
                    "external tool comparison and maturity gates are still incomplete for every reconstruction mode"
                ],
            ),
            MethodLimitationEntry(
                method="bayesian-time-tree",
                status="experimental",
                validated_by=[
                    "tip-date validation",
                    "calibration validation",
                    "BEAST and MrBayes workflow surfaces",
                ],
                assumptions=[
                    "valid dates, calibrations, priors, and convergence diagnostics",
                    "posterior sampling has mixed adequately",
                ],
                invalid_inputs=[
                    "invalid tip dates",
                    "invalid calibrations",
                    "low ESS or conflicting independent runs",
                ],
                limitations=[
                    "workflow support exists, but cross-environment reproducibility and broader benchmark validation remain incomplete"
                ],
            ),
        ],
        limitations=[
            "registry statuses summarize the current checked-in evidence and should move only when new validation surfaces are added or removed",
        ],
    )


def build_scientific_validation_report(
    *, fixtures_root: Path | None = None
) -> ScientificValidationReport:
    """Separate validated, unvalidated, experimental, and unsafe claims for reviewers."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    accuracy = build_method_accuracy_dashboard(fixtures_root=root)
    registry = build_method_limitation_registry()
    claims = [
        ScientificValidationClaim(
            status="validated",
            claim="checked-in Level 1 reference fixtures and benchmark corpora currently pass as expected",
            evidence=[
                f"{row.surface}:{row.passed_count}/{row.coverage_count}"
                for row in accuracy.rows
            ],
        ),
        ScientificValidationClaim(
            status="validated",
            claim="clean, broken, messy, and regression dataset corpora all preserve their expected trust signals",
            evidence=[
                "clean-benchmark-corpus",
                "broken-benchmark-corpus",
                "messy-benchmark-corpus",
                "regression-dataset-corpus",
            ],
        ),
        ScientificValidationClaim(
            status="experimental",
            claim="ancestral reconstruction and Bayesian time-tree workflows have substantive support but still carry incomplete maturity evidence",
            evidence=[
                entry.method
                for entry in registry.entries
                if entry.status == "experimental"
            ],
        ),
        ScientificValidationClaim(
            status="unvalidated",
            claim="cross-environment reproducibility is not yet claimed, and external R ecosystem comparisons remain narrower than the broader validation roadmap",
            evidence=[
                "goal-252 remains outside this iteration",
                "evidence-book/studies/primate-longevity-signal/evidence-001 demonstrates one checked-in R comparison bundle",
            ],
        ),
        ScientificValidationClaim(
            status="unsafe",
            claim="publication-grade time-tree conclusions remain unsafe when tip dates, calibrations, or convergence diagnostics fail",
            evidence=[
                "dataset-audit blockers include invalid tip dates and invalid calibrations",
                "bayesian-time-tree is still marked experimental in the limitation registry",
            ],
        ),
    ]
    return ScientificValidationReport(
        goal_id=249,
        claims=claims,
        limitations=[
            "the report summarizes current checked-in evidence; it does not replace method-specific diagnostics on a new biological dataset",
        ],
    )
