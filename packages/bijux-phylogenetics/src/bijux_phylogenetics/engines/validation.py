from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_descendant_taxa
from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_support_values,
    compare_tree_paths,
)
from bijux_phylogenetics.core.metadata import load_taxon_table
from bijux_phylogenetics.engines.common import (
    build_file_checksums,
    load_engine_manifest,
)
from bijux_phylogenetics.engines.iqtree_artifacts import (
    parse_best_model_file,
    parse_iqtree_model_selection_summary,
    parse_log_likelihood_file,
    resolve_iqtree_model_sidecar,
)
from bijux_phylogenetics.errors import InvalidAlignmentError
from bijux_phylogenetics.io.fasta import (
    load_fasta_alignment,
    summarize_alignment_readiness,
)
from bijux_phylogenetics.io.iqtree_support import (
    IqtreeBranchSupportLabel,
    parse_iqtree_branch_support_label,
    support_fraction,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.tree_set import load_tree_set


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
    manifest_selected_criterion: str | None
    report_selected_model: str | None
    report_selected_criterion: str | None
    artifact_selected_model: str | None
    candidate_model_count: int
    best_model_aic: str | None
    best_model_aicc: str | None
    best_model_bic: str | None
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


@dataclass(slots=True)
class MetadataClusterObservation:
    group: str
    tree_taxa: list[str]
    monophyletic: bool | None
    status: str
    note: str


@dataclass(slots=True)
class MetadataClusteringReport:
    tree_path: Path
    metadata_path: Path
    taxon_column: str
    group_column: str
    group_count: int
    monophyletic_group_count: int
    split_group_count: int
    observations: list[MetadataClusterObservation]


@dataclass(slots=True)
class InferenceFailureTaxonomyReport:
    workflow: str
    failure_category: str
    valid: bool
    issues: list[str]


@dataclass(slots=True)
class BootstrapTreeSetValidationReport:
    tree_set_path: Path
    tree_count: int
    expected_taxa: list[str]
    valid: bool
    issues: list[str]


@dataclass(slots=True)
class BootstrapSupportNode:
    node: str
    descendant_taxa: list[str]
    support: float
    support_fraction: float
    is_backbone: bool


@dataclass(slots=True)
class BootstrapSupportSummaryReport:
    tree_path: Path
    internal_node_count: int
    supported_node_count: int
    minimum_support: float | None
    maximum_support: float | None
    median_support: float | None
    weakly_supported_clade_count: int
    support_histogram: dict[str, int]
    nodes: list[BootstrapSupportNode]
    warnings: list[str]


@dataclass(slots=True)
class WeakBackboneReport:
    tree_path: Path
    threshold: float
    evaluated_backbone_node_count: int
    weak_backbone_node_count: int
    weak_nodes: list[BootstrapSupportNode]
    warnings: list[str]


@dataclass(slots=True)
class ShAlrtSupportNode:
    node: str
    descendant_taxa: list[str]
    sh_alrt_support: float | None
    sh_alrt_support_fraction: float | None
    ufboot_support: float | None
    ufboot_support_fraction: float | None
    is_backbone: bool
    sh_alrt_strong: bool
    ufboot_strong: bool
    conflicting_support_signal: bool
    support_agreement: str


@dataclass(slots=True)
class ShAlrtSupportSummaryReport:
    tree_path: Path
    internal_node_count: int
    annotated_node_count: int
    fully_scored_node_count: int
    minimum_sh_alrt_support: float | None
    maximum_sh_alrt_support: float | None
    minimum_ufboot_support: float | None
    maximum_ufboot_support: float | None
    weak_sh_alrt_clade_count: int
    weak_ufboot_clade_count: int
    conflicting_support_signal_count: int
    nodes: list[ShAlrtSupportNode]
    warnings: list[str]


@dataclass(slots=True)
class InferenceTreeComparisonReport:
    comparison_kind: str
    left_manifest_path: Path
    right_manifest_path: Path
    left_label: str
    right_label: str
    left_tree_path: Path
    right_tree_path: Path
    left_engine_name: str
    right_engine_name: str
    left_selected_model: str | None
    right_selected_model: str | None
    topology: object
    support: object
    branch_lengths: object
    warnings: list[str]


@dataclass(slots=True)
class InferenceOutputConsistencyReport:
    manifest_path: Path
    workflow: str
    failure_category: str
    current_output_checksum_match: bool
    valid: bool
    issues: list[str]


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
            model_sidecar_path=resolve_iqtree_model_sidecar(iqtree_report.with_suffix("")),
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
        candidate_model_count=0 if report_summary is None else report_summary.candidate_count,
        best_model_aic=None if report_summary is None else report_summary.best_model_aic,
        best_model_aicc=None if report_summary is None else report_summary.best_model_aicc,
        best_model_bic=None if report_summary is None else report_summary.best_model_bic,
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


def compare_inferred_tree_to_taxon_metadata(
    tree_path: Path,
    metadata_path: Path,
    *,
    group_column: str,
    taxon_column: str | None = None,
) -> MetadataClusteringReport:
    """Report whether metadata-defined biological groups cluster monophyletically in one inferred tree."""
    tree = load_tree(tree_path)
    table = load_taxon_table(metadata_path, taxon_column=taxon_column)
    if group_column not in table.columns:
        raise ValueError(f"metadata table does not contain column '{group_column}'")
    tree_taxa = set(tree.tip_names)
    group_taxa: dict[str, set[str]] = {}
    for row in table.rows:
        taxon = row[table.taxon_column]
        group = row[group_column].strip()
        if group and taxon in tree_taxa:
            group_taxa.setdefault(group, set()).add(taxon)
    node_taxa_sets = [set(node_descendant_taxa(node)) for node in tree.iter_nodes()]
    observations: list[MetadataClusterObservation] = []
    for group, taxa in sorted(group_taxa.items()):
        ordered_taxa = sorted(taxa)
        if len(ordered_taxa) < 2:
            observations.append(
                MetadataClusterObservation(
                    group=group,
                    tree_taxa=ordered_taxa,
                    monophyletic=None,
                    status="not_evaluable",
                    note="group has fewer than two taxa in the inferred tree",
                )
            )
            continue
        monophyletic = any(node_taxa == taxa for node_taxa in node_taxa_sets)
        observations.append(
            MetadataClusterObservation(
                group=group,
                tree_taxa=ordered_taxa,
                monophyletic=monophyletic,
                status="clusters_as_expected" if monophyletic else "split_unexpectedly",
                note=(
                    "all observed group members collapse to one internal clade"
                    if monophyletic
                    else "group members are distributed across multiple clades in the inferred tree"
                ),
            )
        )
    return MetadataClusteringReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        taxon_column=table.taxon_column,
        group_column=group_column,
        group_count=len(observations),
        monophyletic_group_count=sum(
            1 for row in observations if row.monophyletic is True
        ),
        split_group_count=sum(1 for row in observations if row.monophyletic is False),
        observations=observations,
    )


def classify_inference_workflow_failure(
    *,
    workflow: str,
    input_paths: list[Path],
    output_paths: dict[str, Path],
    run_exit_code: int | None = None,
) -> InferenceFailureTaxonomyReport:
    """Classify one workflow state into a stable inference-failure taxonomy."""
    issues: list[str] = []
    failure_category = "no_failure"
    if any(not path.exists() for path in input_paths):
        failure_category = "input_failure"
        issues.append("one or more declared input files are missing")
    else:
        for input_path in input_paths:
            if input_path.suffix.lower() in {".fasta", ".fa", ".fas", ".faa", ".fna"}:
                try:
                    load_fasta_alignment(input_path)
                except InvalidAlignmentError:
                    failure_category = "input_failure"
                    issues.append("one or more alignment inputs are invalid")
                    break
    if run_exit_code is not None and run_exit_code != 0:
        failure_category = "timeout" if run_exit_code == 124 else "engine_failure"
        issues.append(f"engine exited with code {run_exit_code}")
    missing_outputs = sorted(
        str(path) for path in output_paths.values() if not path.exists()
    )
    if missing_outputs:
        failure_category = "missing_output"
        issues.append("one or more expected outputs are missing")
    parse_failures: list[str] = []
    invalid_outputs: list[str] = []
    for key, path in output_paths.items():
        if not path.exists():
            continue
        if path.is_file() and not path.read_text(encoding="utf-8").strip():
            invalid_outputs.append(key)
            continue
        if path.suffix.lower() in {
            ".treefile",
            ".contree",
            ".nwk",
            ".tre",
            ".tree",
            ".ufboot",
        }:
            try:
                if key == "bootstrap_trees":
                    from bijux_phylogenetics.tree_set import load_tree_set

                    load_tree_set(path)
                else:
                    load_tree(path)
            except Exception:
                parse_failures.append(key)
    if invalid_outputs:
        failure_category = "invalid_output"
        issues.append("one or more outputs exist but are blank or unusable")
    if parse_failures:
        failure_category = "parse_failure"
        issues.append("one or more tree-like outputs could not be parsed")
    return InferenceFailureTaxonomyReport(
        workflow=workflow,
        failure_category=failure_category,
        valid=failure_category == "no_failure",
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


def summarize_bootstrap_support_distribution(
    tree_path: Path,
    *,
    weak_support_threshold: float = 70.0,
) -> BootstrapSupportSummaryReport:
    """Summarize internal-node support values and their distribution across one tree."""
    tree = load_tree(tree_path)
    nodes: list[BootstrapSupportNode] = []
    warnings: list[str] = []
    total_tip_count = tree.tip_count
    for node in tree.iter_nodes():
        if node.is_leaf():
            continue
        descendant_taxa = node_descendant_taxa(node)
        support_label = _parse_iqtree_support_label(node.name)
        support = None if support_label is None else support_label.ufboot_support
        if support is None:
            continue
        normalized_support_fraction = support_fraction(support)
        if normalized_support_fraction is None:
            continue
        nodes.append(
            BootstrapSupportNode(
                node="|".join(descendant_taxa)
                if descendant_taxa
                else (node.name or "<unnamed>"),
                descendant_taxa=descendant_taxa,
                support=support,
                support_fraction=normalized_support_fraction,
                is_backbone=len(descendant_taxa) >= max(2, total_tip_count // 2),
            )
        )
    histogram = {
        "lt50": sum(1 for node in nodes if node.support < 50.0),
        "50to69": sum(1 for node in nodes if 50.0 <= node.support < 70.0),
        "70to89": sum(1 for node in nodes if 70.0 <= node.support < 90.0),
        "ge90": sum(1 for node in nodes if node.support >= 90.0),
    }
    supports = sorted(node.support for node in nodes)
    if len(nodes) < sum(1 for node in tree.iter_nodes() if not node.is_leaf()):
        warnings.append(
            "one or more internal nodes did not expose numeric support labels"
        )
    if any(node.support < weak_support_threshold for node in nodes):
        warnings.append("one or more internal clades remain weakly supported")
    return BootstrapSupportSummaryReport(
        tree_path=tree_path,
        internal_node_count=sum(1 for node in tree.iter_nodes() if not node.is_leaf()),
        supported_node_count=len(nodes),
        minimum_support=None if not supports else supports[0],
        maximum_support=None if not supports else supports[-1],
        median_support=_median_support(supports),
        weakly_supported_clade_count=sum(
            1 for node in nodes if node.support < weak_support_threshold
        ),
        support_histogram=histogram,
        nodes=nodes,
        warnings=warnings,
    )


def detect_weakly_supported_backbone(
    tree_path: Path,
    *,
    threshold: float = 70.0,
) -> WeakBackboneReport:
    """Flag broad internal clades whose support falls below a declared backbone threshold."""
    summary = summarize_bootstrap_support_distribution(
        tree_path, weak_support_threshold=threshold
    )
    weak_nodes = [
        node for node in summary.nodes if node.is_backbone and node.support < threshold
    ]
    warnings = list(summary.warnings)
    if weak_nodes:
        warnings.append(
            "major internal branches remain weakly supported along the backbone"
        )
    return WeakBackboneReport(
        tree_path=tree_path,
        threshold=threshold,
        evaluated_backbone_node_count=sum(
            1 for node in summary.nodes if node.is_backbone
        ),
        weak_backbone_node_count=len(weak_nodes),
        weak_nodes=weak_nodes,
        warnings=warnings,
    )


def summarize_sh_alrt_support_distribution(
    tree_path: Path,
    *,
    sh_alrt_strong_threshold: float = 80.0,
    ufboot_strong_threshold: float = 95.0,
) -> ShAlrtSupportSummaryReport:
    """Summarize compound SH-aLRT/UFBoot branch labels across one tree."""
    tree = load_tree(tree_path)
    nodes: list[ShAlrtSupportNode] = []
    warnings: list[str] = []
    total_tip_count = tree.tip_count
    for node in tree.iter_nodes():
        if node.is_leaf():
            continue
        descendant_taxa = node_descendant_taxa(node)
        support_label = _parse_iqtree_support_label(node.name)
        if support_label is None:
            continue
        sh_alrt_value = support_label.sh_alrt_support
        ufboot_value = support_label.ufboot_support
        sh_alrt_strong = (
            sh_alrt_value is not None and sh_alrt_value >= sh_alrt_strong_threshold
        )
        ufboot_strong = (
            ufboot_value is not None and ufboot_value >= ufboot_strong_threshold
        )
        support_agreement = _support_agreement(
            sh_alrt_value=sh_alrt_value,
            ufboot_value=ufboot_value,
            sh_alrt_strong=sh_alrt_strong,
            ufboot_strong=ufboot_strong,
        )
        nodes.append(
            ShAlrtSupportNode(
                node="|".join(descendant_taxa)
                if descendant_taxa
                else (node.name or "<unnamed>"),
                descendant_taxa=descendant_taxa,
                sh_alrt_support=sh_alrt_value,
                sh_alrt_support_fraction=support_fraction(sh_alrt_value),
                ufboot_support=ufboot_value,
                ufboot_support_fraction=support_fraction(ufboot_value),
                is_backbone=len(descendant_taxa) >= max(2, total_tip_count // 2),
                sh_alrt_strong=sh_alrt_strong,
                ufboot_strong=ufboot_strong,
                conflicting_support_signal=support_agreement
                in {"sh_alrt_only", "ufboot_only"},
                support_agreement=support_agreement,
            )
        )
    if len(nodes) < sum(1 for node in tree.iter_nodes() if not node.is_leaf()):
        warnings.append(
            "one or more internal nodes did not expose parsable sh-alrt or ufboot labels"
        )
    if any(node.conflicting_support_signal for node in nodes):
        warnings.append(
            "one or more internal clades show conflicting sh-alrt and ufboot support signals"
        )
    if any(
        not node.sh_alrt_strong
        for node in nodes
        if node.sh_alrt_support is not None
    ):
        warnings.append("one or more internal clades remain weakly supported by sh-alrt")
    if any(
        not node.ufboot_strong
        for node in nodes
        if node.ufboot_support is not None
    ):
        warnings.append("one or more internal clades remain weakly supported by ufboot")
    sh_alrt_values = sorted(
        node.sh_alrt_support
        for node in nodes
        if node.sh_alrt_support is not None
    )
    ufboot_values = sorted(
        node.ufboot_support for node in nodes if node.ufboot_support is not None
    )
    return ShAlrtSupportSummaryReport(
        tree_path=tree_path,
        internal_node_count=sum(1 for node in tree.iter_nodes() if not node.is_leaf()),
        annotated_node_count=len(nodes),
        fully_scored_node_count=sum(
            1
            for node in nodes
            if node.sh_alrt_support is not None and node.ufboot_support is not None
        ),
        minimum_sh_alrt_support=None if not sh_alrt_values else sh_alrt_values[0],
        maximum_sh_alrt_support=None if not sh_alrt_values else sh_alrt_values[-1],
        minimum_ufboot_support=None if not ufboot_values else ufboot_values[0],
        maximum_ufboot_support=None if not ufboot_values else ufboot_values[-1],
        weak_sh_alrt_clade_count=sum(
            1
            for node in nodes
            if node.sh_alrt_support is not None and not node.sh_alrt_strong
        ),
        weak_ufboot_clade_count=sum(
            1
            for node in nodes
            if node.ufboot_support is not None and not node.ufboot_strong
        ),
        conflicting_support_signal_count=sum(
            1 for node in nodes if node.conflicting_support_signal
        ),
        nodes=nodes,
        warnings=warnings,
    )


def _parse_iqtree_support_label(
    value: str | None,
) -> IqtreeBranchSupportLabel | None:
    return parse_iqtree_branch_support_label(value)


def _support_agreement(
    *,
    sh_alrt_value: float | None,
    ufboot_value: float | None,
    sh_alrt_strong: bool,
    ufboot_strong: bool,
) -> str:
    if sh_alrt_value is None or ufboot_value is None:
        return "incomplete"
    if sh_alrt_strong and ufboot_strong:
        return "both_strong"
    if not sh_alrt_strong and not ufboot_strong:
        return "both_weak"
    if sh_alrt_strong:
        return "sh_alrt_only"
    return "ufboot_only"


def _median_support(values: list[float]) -> float | None:
    if not values:
        return None
    midpoint = len(values) // 2
    if len(values) % 2 == 1:
        return values[midpoint]
    return (values[midpoint - 1] + values[midpoint]) / 2.0


def compare_ml_trees_across_models(
    left_manifest_path: Path,
    right_manifest_path: Path,
) -> InferenceTreeComparisonReport:
    """Compare maximum-likelihood trees produced under different model choices."""
    return _compare_inference_trees(
        left_manifest_path,
        right_manifest_path,
        comparison_kind="model",
        left_label=_manifest_comparison_label(
            left_manifest_path, fallback="left-model"
        ),
        right_label=_manifest_comparison_label(
            right_manifest_path, fallback="right-model"
        ),
    )


def compare_inferred_trees_across_engines(
    left_manifest_path: Path,
    right_manifest_path: Path,
) -> InferenceTreeComparisonReport:
    """Compare inferred trees across two engine workflows."""
    left_manifest = load_engine_manifest(left_manifest_path)
    right_manifest = load_engine_manifest(right_manifest_path)
    return _compare_inference_trees(
        left_manifest_path,
        right_manifest_path,
        comparison_kind="engine",
        left_label=_display_engine_name(str(left_manifest["engine_name"])),
        right_label=_display_engine_name(str(right_manifest["engine_name"])),
    )


def _compare_inference_trees(
    left_manifest_path: Path,
    right_manifest_path: Path,
    *,
    comparison_kind: str,
    left_label: str,
    right_label: str,
) -> InferenceTreeComparisonReport:
    left_manifest = load_engine_manifest(left_manifest_path)
    right_manifest = load_engine_manifest(right_manifest_path)
    left_tree_path = _manifest_tree_output_path(left_manifest)
    right_tree_path = _manifest_tree_output_path(right_manifest)
    topology = compare_tree_paths(left_tree_path, right_tree_path)
    support = compare_support_values(left_tree_path, right_tree_path)
    branch_lengths = compare_branch_lengths(left_tree_path, right_tree_path)
    warnings: list[str] = []
    if not topology.topology_equal:
        warnings.append("inferred topologies differ across compared workflows")
    if topology.same_unrooted_topology and not topology.topology_equal:
        warnings.append(
            "compared workflows agree on unrooted splits but differ in rooting"
        )
    if topology.same_topology_different_branch_lengths:
        warnings.append(
            "compared workflows preserve topology but change branch-length interpretation"
        )
    return InferenceTreeComparisonReport(
        comparison_kind=comparison_kind,
        left_manifest_path=left_manifest_path,
        right_manifest_path=right_manifest_path,
        left_label=left_label,
        right_label=right_label,
        left_tree_path=left_tree_path,
        right_tree_path=right_tree_path,
        left_engine_name=str(left_manifest["engine_name"]),
        right_engine_name=str(right_manifest["engine_name"]),
        left_selected_model=_manifest_selected_model(left_manifest),
        right_selected_model=_manifest_selected_model(right_manifest),
        topology=topology,
        support=support,
        branch_lengths=branch_lengths,
        warnings=warnings,
    )


def _manifest_tree_output_path(manifest: dict[str, object]) -> Path:
    output_paths = {
        key: Path(value) for key, value in dict(manifest["output_paths"]).items()
    }
    tree_path = output_paths.get("tree")
    if tree_path is None:
        raise ValueError("manifest does not expose a tree output")
    return tree_path


def _manifest_selected_model(manifest: dict[str, object]) -> str | None:
    selected_model = manifest.get("selected_model")
    return None if selected_model is None else str(selected_model)


def _manifest_comparison_label(manifest_path: Path, *, fallback: str) -> str:
    manifest = load_engine_manifest(manifest_path)
    selected_model = _manifest_selected_model(manifest)
    if selected_model is not None:
        return selected_model
    return fallback


def _display_engine_name(raw: str) -> str:
    mapping = {
        "iqtree": "IQ-TREE",
        "iqtree2": "IQ-TREE",
        "fasttree": "FastTree",
        "mafft": "MAFFT",
        "trimal": "trimAl",
    }
    return mapping.get(raw.lower(), raw)


def validate_inference_engine_outputs(
    manifest_path: Path,
) -> InferenceOutputConsistencyReport:
    """Detect whether one engine workflow manifest and its current outputs still agree."""
    manifest = load_engine_manifest(manifest_path)
    workflow = str(manifest["workflow"])
    input_paths = [Path(path) for path in manifest["input_paths"]]
    output_paths = {
        key: Path(value) for key, value in dict(manifest["output_paths"]).items()
    }
    run_payload = dict(manifest["run"])
    failure = classify_inference_workflow_failure(
        workflow=workflow,
        input_paths=input_paths,
        output_paths=output_paths,
        run_exit_code=int(run_payload.get("exit_code", 0)),
    )
    issues = list(failure.issues)
    current_checksums = build_file_checksums(list(output_paths.values()))
    manifest_checksums = {
        str(key): str(value)
        for key, value in dict(manifest.get("output_checksums", {})).items()
    }
    current_output_checksum_match = current_checksums == manifest_checksums
    if not current_output_checksum_match:
        issues.append(
            "current output checksums do not match the recorded manifest outputs"
        )
    if workflow == "model-selection":
        model_validation = validate_model_selection_against_engine_outputs(
            manifest_path
        )
        issues.extend(model_validation.issues)
    elif workflow == "maximum-likelihood-tree":
        tree_validation = validate_ml_tree_contains_expected_taxa(manifest_path)
        issues.extend(tree_validation.issues)
        if output_paths.get("iqtree_log") is None:
            issues.append("maximum-likelihood manifest is missing the iqtree_log output")
        if manifest.get("selected_model") is None:
            issues.append(
                "maximum-likelihood manifest is missing the selected model"
            )
        if manifest.get("log_likelihood") is None:
            issues.append(
                "maximum-likelihood manifest is missing the log_likelihood field"
            )
    elif workflow == "bootstrap-support":
        bootstrap_path = output_paths.get("bootstrap_trees")
        if bootstrap_path is None:
            issues.append(
                "bootstrap-support manifest is missing the bootstrap_trees output"
            )
        else:
            bootstrap_validation = validate_bootstrap_tree_set(bootstrap_path)
            issues.extend(bootstrap_validation.issues)
        if output_paths.get("support_tree") is None:
            issues.append("bootstrap-support manifest is missing the support_tree output")
        if output_paths.get("support_table") is None:
            issues.append("bootstrap-support manifest is missing the support_table output")
        if output_paths.get("low_support_branches") is None:
            issues.append(
                "bootstrap-support manifest is missing the low_support_branches output"
            )
        if output_paths.get("support_histogram") is None:
            issues.append(
                "bootstrap-support manifest is missing the support_histogram output"
            )
        if output_paths.get("iqtree_log") is None:
            issues.append("bootstrap-support manifest is missing the iqtree_log output")
        if manifest.get("selected_model") is None:
            issues.append("bootstrap-support manifest is missing the selected model")
        if manifest.get("log_likelihood") is None:
            issues.append("bootstrap-support manifest is missing the log_likelihood field")
        iqtree_summary = manifest.get("iqtree_summary")
        if not isinstance(iqtree_summary, dict):
            issues.append("bootstrap-support manifest is missing the iqtree_summary")
        elif int(iqtree_summary.get("support_value_count", 0)) < 1:
            issues.append(
                "bootstrap-support manifest does not record parsed support values"
            )
        bootstrap_support_summary = manifest.get("bootstrap_support_summary")
        if not isinstance(bootstrap_support_summary, dict):
            issues.append(
                "bootstrap-support manifest is missing the bootstrap_support_summary"
            )
        else:
            if int(bootstrap_support_summary.get("supported_node_count", 0)) < 1:
                issues.append(
                    "bootstrap-support manifest does not record any supported internal nodes"
                )
            histogram = bootstrap_support_summary.get("support_histogram")
            if not isinstance(histogram, dict):
                issues.append(
                    "bootstrap-support manifest does not record the support histogram"
                )
            elif sorted(histogram) != ["50to69", "70to89", "ge90", "lt50"]:
                issues.append(
                    "bootstrap-support manifest records an incomplete support histogram"
                )
        weak_backbone_report = manifest.get("weak_backbone_report")
        if not isinstance(weak_backbone_report, dict):
            issues.append(
                "bootstrap-support manifest is missing the weak_backbone_report"
            )
        elif float(weak_backbone_report.get("threshold", 0.0)) <= 0.0:
            issues.append(
                "bootstrap-support manifest records an invalid weak-backbone threshold"
            )
    elif workflow == "bootstrap-consensus":
        if output_paths.get("iqtree_log") is None:
            issues.append(
                "bootstrap-consensus manifest is missing the iqtree_log output"
            )
        iqtree_summary = manifest.get("iqtree_summary")
        if not isinstance(iqtree_summary, dict):
            issues.append("bootstrap-consensus manifest is missing the iqtree_summary")
        elif int(iqtree_summary.get("support_value_count", 0)) < 1:
            issues.append(
                "bootstrap-consensus manifest does not record parsed support values"
            )
    elif workflow == "sh-alrt-support":
        bootstrap_path = output_paths.get("bootstrap_trees")
        if bootstrap_path is None:
            issues.append(
                "sh-alrt-support manifest is missing the bootstrap_trees output"
            )
        else:
            bootstrap_validation = validate_bootstrap_tree_set(bootstrap_path)
            issues.extend(bootstrap_validation.issues)
        if output_paths.get("support_tree") is None:
            issues.append("sh-alrt-support manifest is missing the support_tree output")
        if output_paths.get("support_table") is None:
            issues.append("sh-alrt-support manifest is missing the support_table output")
        if output_paths.get("conflicting_support_branches") is None:
            issues.append(
                "sh-alrt-support manifest is missing the conflicting_support_branches output"
            )
        if output_paths.get("iqtree_log") is None:
            issues.append("sh-alrt-support manifest is missing the iqtree_log output")
        if manifest.get("selected_model") is None:
            issues.append("sh-alrt-support manifest is missing the selected model")
        if manifest.get("log_likelihood") is None:
            issues.append("sh-alrt-support manifest is missing the log_likelihood field")
        iqtree_summary = manifest.get("iqtree_summary")
        if not isinstance(iqtree_summary, dict):
            issues.append("sh-alrt-support manifest is missing the iqtree_summary")
        elif int(iqtree_summary.get("support_value_count", 0)) < 1:
            issues.append(
                "sh-alrt-support manifest does not record parsed ufboot support values"
            )
        bootstrap_support_summary = manifest.get("bootstrap_support_summary")
        if not isinstance(bootstrap_support_summary, dict):
            issues.append(
                "sh-alrt-support manifest is missing the bootstrap_support_summary"
            )
        elif int(bootstrap_support_summary.get("supported_node_count", 0)) < 1:
            issues.append(
                "sh-alrt-support manifest does not record any parsed ufboot-supported nodes"
            )
        sh_alrt_support_summary = manifest.get("sh_alrt_support_summary")
        if not isinstance(sh_alrt_support_summary, dict):
            issues.append(
                "sh-alrt-support manifest is missing the sh_alrt_support_summary"
            )
        else:
            if int(sh_alrt_support_summary.get("annotated_node_count", 0)) < 1:
                issues.append(
                    "sh-alrt-support manifest does not record any parsed sh-alrt labels"
                )
            if int(sh_alrt_support_summary.get("fully_scored_node_count", 0)) < 1:
                issues.append(
                    "sh-alrt-support manifest does not record any jointly scored sh-alrt and ufboot branches"
                )
    return InferenceOutputConsistencyReport(
        manifest_path=manifest_path,
        workflow=workflow,
        failure_category=failure.failure_category,
        current_output_checksum_match=current_output_checksum_match,
        valid=not issues,
        issues=sorted(dict.fromkeys(issues)),
    )
