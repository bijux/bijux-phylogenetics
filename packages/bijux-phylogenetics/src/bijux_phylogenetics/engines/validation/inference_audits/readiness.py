from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta.quality import summarize_alignment_readiness

from .contracts import InferenceReadinessAuditReport, InferenceReadinessDecision


def audit_alignment_inference_readiness(path: Path) -> InferenceReadinessAuditReport:
    """Classify whether one alignment is suitable for ML, fast approximate, Bayesian, or none."""
    readiness = summarize_alignment_readiness(path)
    method_by_name = {method.analysis: method for method in readiness.methods}
    maximum_likelihood = method_by_name["maximum_likelihood"]
    bayesian = method_by_name["bayesian"]
    fast_approximate_ready = (
        method_by_name["distance"].ready or maximum_likelihood.ready
    )
    fast_approximate_blockers = (
        []
        if fast_approximate_ready
        else sorted(
            dict.fromkeys(
                method_by_name["distance"].blockers + maximum_likelihood.blockers
            )
        )
    )
    generic_warnings = sorted(
        dict.fromkeys(
            readiness.warnings
            + method_by_name["distance"].warnings
            + maximum_likelihood.warnings
            + bayesian.warnings
        )
    )
    decisions = [
        InferenceReadinessDecision(
            workflow="maximum_likelihood",
            ready=maximum_likelihood.ready,
            blockers=maximum_likelihood.blockers,
            warnings=maximum_likelihood.warnings,
        ),
        InferenceReadinessDecision(
            workflow="fast_approximate",
            ready=fast_approximate_ready,
            blockers=fast_approximate_blockers,
            warnings=method_by_name["distance"].warnings,
        ),
        InferenceReadinessDecision(
            workflow="bayesian",
            ready=bayesian.ready,
            blockers=bayesian.blockers,
            warnings=bayesian.warnings,
        ),
    ]
    if maximum_likelihood.ready:
        overall_decision = "ready"
        recommended_workflow = "maximum_likelihood"
    elif fast_approximate_ready:
        overall_decision = "ready_with_limits"
        recommended_workflow = "fast_approximate"
    elif bayesian.ready:
        overall_decision = "ready_with_limits"
        recommended_workflow = "bayesian"
    else:
        overall_decision = "blocked"
        recommended_workflow = "unsuitable"
    return InferenceReadinessAuditReport(
        alignment_path=path,
        sequence_count=readiness.sequence_count,
        alignment_length=readiness.alignment_length,
        inferred_alphabet=readiness.inferred_alphabet,
        overall_decision=overall_decision,
        recommended_workflow=recommended_workflow,
        decisions=decisions,
        warnings=generic_warnings,
    )
