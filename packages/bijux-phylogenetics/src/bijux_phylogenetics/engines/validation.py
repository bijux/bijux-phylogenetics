from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path
import re

from bijux_phylogenetics.engines.common import load_engine_manifest
from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.io.fasta import summarize_alignment_readiness
from bijux_phylogenetics.io.trees import load_tree

_BEST_MODEL_PATTERN = re.compile(
    r"(?:best-fit model(?: according to [A-Z0-9]+)?|best model)\s*[:=]\s*(?P<model>[A-Za-z0-9+._-]+)",
    re.IGNORECASE,
)


@dataclass(slots=True)
class InferenceReadinessDecision:
    workflow: str
    ready: bool
    blockers: list[str]
    warnings: list[str]


@dataclass(slots=True)
class InferenceReadinessAuditReport:
    alignment_path: Path
    sequence_count: int
    alignment_length: int | None
    inferred_alphabet: str
    overall_decision: str
    recommended_workflow: str
    decisions: list[InferenceReadinessDecision]
    warnings: list[str]


@dataclass(slots=True)
class ModelSelectionValidationReport:
    manifest_path: Path
    manifest_selected_model: str | None
    report_selected_model: str | None
    artifact_selected_model: str | None
    valid: bool
    issues: list[str]


@dataclass(slots=True)
class MLTreeTaxonValidationReport:
    manifest_path: Path
    expected_taxa: list[str]
    observed_taxa: list[str]
    missing_taxa: list[str]
    unexpected_taxa: list[str]
    valid: bool
    issues: list[str]


def audit_alignment_inference_readiness(path: Path) -> InferenceReadinessAuditReport:
    """Classify whether one alignment is suitable for ML, fast approximate, Bayesian, or none."""
    readiness = summarize_alignment_readiness(path)
    method_by_name = {method.analysis: method for method in readiness.methods}
    maximum_likelihood = method_by_name["maximum_likelihood"]
    bayesian = method_by_name["bayesian"]
    fast_approximate_ready = method_by_name["distance"].ready or maximum_likelihood.ready
    fast_approximate_blockers = [] if fast_approximate_ready else sorted(
        dict.fromkeys(method_by_name["distance"].blockers + maximum_likelihood.blockers)
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


def _parse_best_model_text(text: str) -> str | None:
    match = _BEST_MODEL_PATTERN.search(text)
    return None if match is None else match.group("model")


def _parse_best_model_file(path: Path) -> str | None:
    if not path.exists():
        return None
    return _parse_best_model_text(path.read_text(encoding="utf-8"))


def validate_model_selection_against_engine_outputs(manifest_path: Path) -> ModelSelectionValidationReport:
    """Verify that the selected model exposed in workflow artifacts matches the engine outputs exactly."""
    manifest = load_engine_manifest(manifest_path)
    output_paths = {key: Path(value) for key, value in dict(manifest["output_paths"]).items()}
    iqtree_report = output_paths.get("iqtree_report")
    selected_model_file = output_paths.get("selected_model")
    report_selected_model = None if iqtree_report is None else _parse_best_model_file(iqtree_report)
    artifact_selected_model = None
    if selected_model_file is not None and selected_model_file.exists():
        artifact_selected_model = selected_model_file.read_text(encoding="utf-8").strip() or None
    if report_selected_model is None and iqtree_report is not None:
        model_sidecar = iqtree_report.with_suffix(".model")
        if model_sidecar.exists():
            report_selected_model = _parse_best_model_file(model_sidecar)
        else:
            model_gz_sidecar = iqtree_report.with_suffix(".model.gz")
            if model_gz_sidecar.exists():
                report_selected_model = _parse_best_model_text(
                    gzip.decompress(model_gz_sidecar.read_bytes()).decode("utf-8", errors="replace")
                )
    manifest_selected_model = manifest.get("selected_model")
    issues: list[str] = []
    if manifest.get("workflow") != "model-selection":
        issues.append("manifest does not describe a model-selection workflow")
    if manifest_selected_model is None:
        issues.append("manifest selected_model field is missing")
    if report_selected_model is None:
        issues.append("engine report does not expose a parsable best-fit model")
    if artifact_selected_model is None:
        issues.append("selected-model artifact is missing or blank")
    comparable = [value for value in (manifest_selected_model, report_selected_model, artifact_selected_model) if value is not None]
    if comparable and len(set(comparable)) > 1:
        issues.append("selected model disagrees across manifest, report, and exported artifact")
    return ModelSelectionValidationReport(
        manifest_path=manifest_path,
        manifest_selected_model=None if manifest_selected_model is None else str(manifest_selected_model),
        report_selected_model=report_selected_model,
        artifact_selected_model=artifact_selected_model,
        valid=not issues,
        issues=issues,
    )


def validate_ml_tree_contains_expected_taxa(manifest_path: Path) -> MLTreeTaxonValidationReport:
    """Ensure an inferred ML tree contains exactly the taxa present in the input alignment."""
    manifest = load_engine_manifest(manifest_path)
    input_paths = [Path(path) for path in manifest["input_paths"]]
    output_paths = {key: Path(value) for key, value in dict(manifest["output_paths"]).items()}
    issues: list[str] = []
    if manifest.get("workflow") != "maximum-likelihood-tree":
        issues.append("manifest does not describe a maximum-likelihood-tree workflow")
    if not input_paths:
        issues.append("manifest does not include an input alignment path")
        expected_taxa: list[str] = []
    else:
        expected_taxa = sorted(record.identifier for record in load_fasta_alignment(input_paths[0]))
    tree_path = output_paths.get("tree")
    if tree_path is None or not tree_path.exists():
        issues.append("manifest tree output is missing")
        observed_taxa: list[str] = []
    else:
        observed_taxa = sorted(load_tree(tree_path).tip_names)
    missing_taxa = sorted(set(expected_taxa) - set(observed_taxa))
    unexpected_taxa = sorted(set(observed_taxa) - set(expected_taxa))
    if missing_taxa:
        issues.append("inferred tree is missing one or more expected taxa from the alignment")
    if unexpected_taxa:
        issues.append("inferred tree contains taxa not present in the input alignment")
    return MLTreeTaxonValidationReport(
        manifest_path=manifest_path,
        expected_taxa=expected_taxa,
        observed_taxa=observed_taxa,
        missing_taxa=missing_taxa,
        unexpected_taxa=unexpected_taxa,
        valid=not issues,
        issues=issues,
    )
