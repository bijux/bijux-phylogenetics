from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines.common import (
    build_file_checksums,
    load_engine_manifest,
)

from .contracts import InferenceOutputConsistencyReport
from .failure_taxonomy import classify_inference_workflow_failure
from .manifest_validation import (
    validate_bootstrap_tree_set,
    validate_ml_tree_contains_expected_taxa,
    validate_model_selection_against_engine_outputs,
)


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
