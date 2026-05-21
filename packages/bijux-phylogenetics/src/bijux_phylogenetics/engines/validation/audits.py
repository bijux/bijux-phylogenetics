from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_support_values,
    compare_tree_paths,
)
from bijux_phylogenetics.engines.artifacts.support import (
    BootstrapSupportNode,
    BootstrapSupportSummaryReport,
    FastTreeSupportNode,
    FastTreeSupportSummaryReport,
    ShAlrtSupportNode,
    ShAlrtSupportSummaryReport,
    WeakBackboneReport,
)
from bijux_phylogenetics.io.fasttree_support import (
    FastTreeBranchSupportLabel,
    parse_fasttree_branch_support_label,
)
from bijux_phylogenetics.io.iqtree_support import (
    IqtreeBranchSupportLabel,
    parse_iqtree_branch_support_label,
    support_fraction,
)
from bijux_phylogenetics.io.trees import load_tree

from .inference_audits import (
    BootstrapTreeSetValidationReport,
    InferenceFailureTaxonomyReport,
    InferenceOutputConsistencyReport,
    InferenceReadinessAuditReport,
    InferenceReadinessDecision,
    InferenceTreeComparisonReport,
    MetadataClusteringReport,
    MetadataClusterObservation,
    MLTreeTaxonValidationReport,
    ModelSelectionValidationReport,
    audit_alignment_inference_readiness,
    classify_inference_workflow_failure,
    compare_inferred_tree_to_taxon_metadata,
    validate_bootstrap_tree_set,
    validate_ml_tree_contains_expected_taxa,
    validate_model_selection_against_engine_outputs,
)


def summarize_bootstrap_support_distribution(
    tree_path: Path,
    *,
    weak_support_threshold: float = 70.0,
) -> BootstrapSupportSummaryReport:
    """Summarize internal-node support values and their distribution across one tree."""
    from bijux_phylogenetics.ancestral.common import node_descendant_taxa

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


def summarize_fasttree_support_distribution(
    tree_path: Path,
    *,
    weak_support_threshold: float = 0.7,
) -> FastTreeSupportSummaryReport:
    """Summarize FastTree SH-like local support labels across one tree."""
    from bijux_phylogenetics.ancestral.common import node_descendant_taxa

    tree = load_tree(tree_path)
    nodes: list[FastTreeSupportNode] = []
    warnings: list[str] = []
    total_tip_count = tree.tip_count
    internal_node_count = sum(1 for node in tree.iter_nodes() if not node.is_leaf())
    for node in tree.iter_nodes():
        if node.is_leaf():
            continue
        descendant_taxa = node_descendant_taxa(node)
        support_label = _parse_fasttree_support_label(node.name)
        if support_label is None:
            continue
        local_support = support_label.local_support
        nodes.append(
            FastTreeSupportNode(
                node="|".join(descendant_taxa)
                if descendant_taxa
                else (node.name or "<unnamed>"),
                descendant_taxa=descendant_taxa,
                local_support=local_support,
                support_fraction=local_support,
                is_backbone=len(descendant_taxa) >= max(2, total_tip_count // 2),
            )
        )
    histogram = {
        "lt0p5": sum(1 for node in nodes if node.local_support < 0.5),
        "0p5to0p69": sum(1 for node in nodes if 0.5 <= node.local_support < 0.7),
        "0p7to0p89": sum(1 for node in nodes if 0.7 <= node.local_support < 0.9),
        "ge0p9": sum(1 for node in nodes if node.local_support >= 0.9),
    }
    supports = sorted(node.local_support for node in nodes)
    if len(nodes) < internal_node_count:
        warnings.append(
            "one or more internal nodes did not expose parsable FastTree local support labels"
        )
    if any(node.local_support < weak_support_threshold for node in nodes):
        warnings.append(
            "one or more internal clades remain weakly supported by FastTree local support"
        )
    warnings.append(
        "FastTree is an approximately maximum-likelihood workflow and its SH-like local support values should be reviewed as approximate evidence"
    )
    return FastTreeSupportSummaryReport(
        tree_path=tree_path,
        internal_node_count=internal_node_count,
        annotated_node_count=len(nodes),
        minimum_local_support=None if not supports else supports[0],
        maximum_local_support=None if not supports else supports[-1],
        median_local_support=_median_support(supports),
        weakly_supported_clade_count=sum(
            1 for node in nodes if node.local_support < weak_support_threshold
        ),
        support_histogram=histogram,
        approximate_method=True,
        support_label_kind="sh-like-local-support",
        support_scale="proportion-0-to-1",
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
    from bijux_phylogenetics.ancestral.common import node_descendant_taxa

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
        not node.sh_alrt_strong for node in nodes if node.sh_alrt_support is not None
    ):
        warnings.append(
            "one or more internal clades remain weakly supported by sh-alrt"
        )
    if any(not node.ufboot_strong for node in nodes if node.ufboot_support is not None):
        warnings.append("one or more internal clades remain weakly supported by ufboot")
    sh_alrt_values = sorted(
        node.sh_alrt_support for node in nodes if node.sh_alrt_support is not None
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


def _parse_fasttree_support_label(
    value: str | None,
) -> FastTreeBranchSupportLabel | None:
    return parse_fasttree_branch_support_label(value)


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
    tree_path = (
        output_paths.get("tree")
        or output_paths.get("support_tree")
        or output_paths.get("consensus_tree")
    )
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
            issues.append(
                "maximum-likelihood manifest is missing the iqtree_log output"
            )
        if manifest.get("selected_model") is None:
            issues.append("maximum-likelihood manifest is missing the selected model")
        if manifest.get("log_likelihood") is None:
            issues.append(
                "maximum-likelihood manifest is missing the log_likelihood field"
            )
    elif workflow == "fast-approximate-tree":
        if output_paths.get("support_table") is None:
            issues.append(
                "fast-approximate-tree manifest is missing the support_table output"
            )
        if output_paths.get("low_support_branches") is None:
            issues.append(
                "fast-approximate-tree manifest is missing the low_support_branches output"
            )
        if output_paths.get("support_histogram") is None:
            issues.append(
                "fast-approximate-tree manifest is missing the support_histogram output"
            )
        fasttree_support_summary = manifest.get("fasttree_support_summary")
        if not isinstance(fasttree_support_summary, dict):
            issues.append(
                "fast-approximate-tree manifest is missing the fasttree_support_summary"
            )
        else:
            if int(fasttree_support_summary.get("annotated_node_count", 0)) < 1:
                issues.append(
                    "fast-approximate-tree manifest does not record any parsed FastTree local support labels"
                )
            histogram = fasttree_support_summary.get("support_histogram")
            if not isinstance(histogram, dict):
                issues.append(
                    "fast-approximate-tree manifest does not record the FastTree support histogram"
                )
            elif sorted(histogram) != ["0p5to0p69", "0p7to0p89", "ge0p9", "lt0p5"]:
                issues.append(
                    "fast-approximate-tree manifest records an incomplete FastTree support histogram"
                )
            if fasttree_support_summary.get("approximate_method") is not True:
                issues.append(
                    "fast-approximate-tree manifest does not record the FastTree approximation contract"
                )
            if (
                fasttree_support_summary.get("support_label_kind")
                != "sh-like-local-support"
            ):
                issues.append(
                    "fast-approximate-tree manifest does not record the FastTree support label kind"
                )
            if fasttree_support_summary.get("support_scale") != "proportion-0-to-1":
                issues.append(
                    "fast-approximate-tree manifest does not record the FastTree support scale"
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
            issues.append(
                "bootstrap-support manifest is missing the support_tree output"
            )
        if output_paths.get("support_table") is None:
            issues.append(
                "bootstrap-support manifest is missing the support_table output"
            )
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
            issues.append(
                "bootstrap-support manifest is missing the log_likelihood field"
            )
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
            issues.append(
                "sh-alrt-support manifest is missing the support_table output"
            )
        if output_paths.get("conflicting_support_branches") is None:
            issues.append(
                "sh-alrt-support manifest is missing the conflicting_support_branches output"
            )
        if output_paths.get("iqtree_log") is None:
            issues.append("sh-alrt-support manifest is missing the iqtree_log output")
        if manifest.get("selected_model") is None:
            issues.append("sh-alrt-support manifest is missing the selected model")
        if manifest.get("log_likelihood") is None:
            issues.append(
                "sh-alrt-support manifest is missing the log_likelihood field"
            )
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
