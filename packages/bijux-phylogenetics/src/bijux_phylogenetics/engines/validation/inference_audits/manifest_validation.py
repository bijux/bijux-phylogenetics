from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines.artifacts.iqtree import (
    parse_best_model_file,
    parse_iqtree_model_selection_summary,
    parse_log_likelihood_file,
    resolve_iqtree_model_sidecar,
)
from bijux_phylogenetics.engines.common import load_engine_manifest
from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.trees import load_tree_set

from .contracts import (
    BootstrapTreeSetValidationReport,
    MLTreeTaxonValidationReport,
    ModelSelectionValidationReport,
)


def validate_model_selection_against_engine_outputs(
    manifest_path: Path,
) -> ModelSelectionValidationReport:
    """Verify that the selected model exposed in workflow artifacts matches the engine outputs exactly."""
    manifest = load_engine_manifest(manifest_path)
    output_paths = {
        key: Path(value) for key, value in dict(manifest["output_paths"]).items()
    }
    iqtree_report = output_paths.get("iqtree_report")
    iqtree_log = output_paths.get("iqtree_log")
    selected_model_file = output_paths.get("selected_model")
    report_summary = (
        None
        if iqtree_report is None
        else parse_iqtree_model_selection_summary(
            iqtree_report_path=iqtree_report,
            model_sidecar_path=resolve_iqtree_model_sidecar(
                iqtree_report.with_suffix("")
            ),
        )
    )
    report_selected_model = (
        None if report_summary is None else report_summary.selected_model
    )
    report_selected_criterion = (
        None if report_summary is None else report_summary.selected_criterion
    )
    report_log_likelihood = (
        None if iqtree_report is None else parse_log_likelihood_file(iqtree_report)
    )
    if report_log_likelihood is None and iqtree_log is not None:
        report_log_likelihood = parse_log_likelihood_file(iqtree_log)
    artifact_selected_model = None
    if selected_model_file is not None and selected_model_file.exists():
        artifact_selected_model = (
            selected_model_file.read_text(encoding="utf-8").strip() or None
        )
    if report_selected_model is None and iqtree_report is not None:
        model_sidecar = resolve_iqtree_model_sidecar(iqtree_report.with_suffix(""))
        if model_sidecar is not None:
            report_selected_model = parse_best_model_file(model_sidecar)
    manifest_selected_model = manifest.get("selected_model")
    manifest_log_likelihood = manifest.get("log_likelihood")
    manifest_summary = manifest.get("model_selection_summary")
    manifest_selected_criterion = (
        None
        if not isinstance(manifest_summary, dict)
        or manifest_summary.get("selected_criterion") is None
        else str(manifest_summary["selected_criterion"])
    )
    issues: list[str] = []
    if manifest.get("workflow") != "model-selection":
        issues.append("manifest does not describe a model-selection workflow")
    if iqtree_log is None or not iqtree_log.exists():
        issues.append("iqtree log artifact is missing")
    if manifest_selected_model is None:
        issues.append("manifest selected_model field is missing")
    if manifest_log_likelihood is None:
        issues.append("manifest log_likelihood field is missing")
    if not isinstance(manifest_summary, dict):
        issues.append("manifest model_selection_summary field is missing")
    else:
        if int(manifest_summary.get("candidate_count", 0)) < 1:
            issues.append("manifest does not record any candidate substitution models")
        if manifest_summary.get("selected_criterion") is None:
            issues.append(
                "manifest model-selection summary does not record the selected criterion"
            )
    if report_selected_model is None:
        issues.append("engine report does not expose a parsable best-fit model")
    if report_summary is None or report_summary.candidate_count < 1:
        issues.append("engine report does not expose a parsable candidate-model table")
    elif (
        report_summary.best_model_aic is None
        or report_summary.best_model_aicc is None
        or report_summary.best_model_bic is None
    ):
        issues.append("engine report does not expose AIC, AICc, and BIC winners")
    if report_log_likelihood is None:
        issues.append("engine artifacts do not expose a parsable log-likelihood")
    if artifact_selected_model is None:
        issues.append("selected-model artifact is missing or blank")
    comparable = [
        value
        for value in (
            manifest_selected_model,
            report_selected_model,
            artifact_selected_model,
        )
        if value is not None
    ]
    if comparable and len(set(comparable)) > 1:
        issues.append(
            "selected model disagrees across manifest, report, and exported artifact"
        )
    if (
        manifest_log_likelihood is not None
        and report_log_likelihood is not None
        and float(manifest_log_likelihood) != report_log_likelihood
    ):
        issues.append(
            "log-likelihood disagrees between manifest and iqtree inference artifacts"
        )
    if (
        manifest_selected_criterion is not None
        and report_selected_criterion is not None
        and manifest_selected_criterion != report_selected_criterion
    ):
        issues.append(
            "selected criterion disagrees between manifest and iqtree inference artifacts"
        )
    return ModelSelectionValidationReport(
        manifest_path=manifest_path,
        manifest_selected_model=None
        if manifest_selected_model is None
        else str(manifest_selected_model),
        manifest_selected_criterion=manifest_selected_criterion,
        report_selected_model=report_selected_model,
        report_selected_criterion=report_selected_criterion,
        artifact_selected_model=artifact_selected_model,
        candidate_model_count=0
        if report_summary is None
        else report_summary.candidate_count,
        best_model_aic=None
        if report_summary is None
        else report_summary.best_model_aic,
        best_model_aicc=None
        if report_summary is None
        else report_summary.best_model_aicc,
        best_model_bic=None
        if report_summary is None
        else report_summary.best_model_bic,
        valid=not issues,
        issues=issues,
    )


def validate_ml_tree_contains_expected_taxa(
    manifest_path: Path,
) -> MLTreeTaxonValidationReport:
    """Ensure an inferred ML tree contains exactly the taxa present in the input alignment."""
    manifest = load_engine_manifest(manifest_path)
    input_paths = [Path(path) for path in manifest["input_paths"]]
    output_paths = {
        key: Path(value) for key, value in dict(manifest["output_paths"]).items()
    }
    issues: list[str] = []
    if manifest.get("workflow") != "maximum-likelihood-tree":
        issues.append("manifest does not describe a maximum-likelihood-tree workflow")
    if not input_paths:
        issues.append("manifest does not include an input alignment path")
        expected_taxa: list[str] = []
    else:
        expected_taxa = sorted(
            record.identifier for record in load_fasta_alignment(input_paths[0])
        )
    tree_path = output_paths.get("tree")
    if tree_path is None or not tree_path.exists():
        issues.append("manifest tree output is missing")
        observed_taxa: list[str] = []
    else:
        observed_taxa = sorted(load_tree(tree_path).tip_names)
    missing_taxa = sorted(set(expected_taxa) - set(observed_taxa))
    unexpected_taxa = sorted(set(observed_taxa) - set(expected_taxa))
    if missing_taxa:
        issues.append(
            "inferred tree is missing one or more expected taxa from the alignment"
        )
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


def validate_bootstrap_tree_set(path: Path) -> BootstrapTreeSetValidationReport:
    """Validate that every bootstrap tree parses and shares the same taxon set."""
    issues: list[str] = []
    try:
        report = load_tree_set(path)
    except Exception as error:
        return BootstrapTreeSetValidationReport(
            tree_set_path=path,
            tree_count=0,
            expected_taxa=[],
            valid=False,
            issues=[f"bootstrap tree set could not be parsed: {error}"],
        )
    taxa_sets = {tuple(record.taxa) for record in report.records}
    if len(taxa_sets) != 1:
        issues.append("bootstrap trees do not all contain the exact same taxon set")
    return BootstrapTreeSetValidationReport(
        tree_set_path=path,
        tree_count=report.tree_count,
        expected_taxa=report.shared_taxa,
        valid=not issues,
        issues=issues,
    )
