from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from bijux_phylogenetics.engines.common import load_engine_manifest


@dataclass(slots=True)
class TreeInferenceMethodsSummaryTextResult:
    output_path: Path
    title: str
    warning_count: int
    selected_model: str
    bootstrap_replicates: int
    trimmed_alignment_length: int
    supported_node_count: int
    text: str
    workflow_manifest_path: Path
    workflow_manifest: dict[str, object]


def _manifest_path_text(value: object | None) -> str:
    if value is None:
        return "not recorded"
    return str(value)


def _manifest_dict(payload: dict[str, object], key: str) -> dict[str, object]:
    value = payload.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"tree-inference workflow manifest is missing `{key}`")
    return value


def _manifest_nested_dict(
    payload: dict[str, object], *keys: str
) -> dict[str, object] | None:
    current: object = payload
    for key in keys:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
    return current if isinstance(current, dict) else None


def _manifest_list(payload: dict[str, object], key: str) -> list[object]:
    value = payload.get(key)
    return value if isinstance(value, list) else []


def _tree_inference_payload(
    *,
    workflow_report: Any | None,
    workflow_manifest_path: Path | None,
) -> tuple[dict[str, object], Path]:
    if (workflow_report is None) == (workflow_manifest_path is None):
        raise ValueError(
            "provide exactly one of workflow_report or workflow_manifest_path"
        )
    if workflow_report is not None:
        payload = asdict(workflow_report)
        manifest_path_value = payload.get("manifest_path")
        if manifest_path_value is None:
            raise ValueError(
                "tree-inference workflow report did not expose manifest_path"
            )
        return payload, Path(str(manifest_path_value))
    if workflow_manifest_path is None:
        raise ValueError(
            "workflow_manifest_path is required when workflow_report is absent"
        )
    return load_engine_manifest(workflow_manifest_path), workflow_manifest_path


def _tree_inference_warning_lines(payload: dict[str, object]) -> list[str]:
    support_summary = _manifest_nested_dict(payload, "support_summary")
    warning_texts = [
        str(item)
        for item in [
            *_manifest_list(payload, "warnings"),
            *(
                []
                if support_summary is None
                else _manifest_list(support_summary, "warnings")
            ),
        ]
        if str(item).strip()
    ]
    deduplicated: list[str] = []
    for warning in warning_texts:
        if warning not in deduplicated:
            deduplicated.append(warning)
    return deduplicated


def _output_paths_dict(workflow_payload: dict[str, object]) -> dict[str, object]:
    return _manifest_nested_dict(workflow_payload, "output_paths") or {}


def write_tree_inference_methods_summary_text(
    path: Path,
    *,
    workflow_report: Any | None = None,
    workflow_manifest_path: Path | None = None,
) -> TreeInferenceMethodsSummaryTextResult:
    """Write reviewer-facing methods text for one governed fasta-to-tree workflow."""
    payload, manifest_path = _tree_inference_payload(
        workflow_report=workflow_report,
        workflow_manifest_path=workflow_manifest_path,
    )
    workflow_name = str(payload.get("workflow", "fasta-to-tree"))
    if workflow_name != "fasta-to-tree":
        raise ValueError(
            f"tree-inference methods summary expects a fasta-to-tree workflow manifest, got `{workflow_name}`"
        )
    selected_model = payload.get("selected_model")
    if not isinstance(selected_model, str) or not selected_model.strip():
        raise ValueError(
            "tree-inference workflow manifest did not expose a selected model"
        )
    input_validation = _manifest_dict(payload, "input_validation")
    input_summary = _manifest_dict(input_validation, "summary")
    trimming_workflow = _manifest_dict(payload, "trimming_workflow")
    trimming_summary = _manifest_dict(trimming_workflow, "trimming_summary")
    model_selection_workflow = _manifest_dict(payload, "model_selection_workflow")
    model_selection_summary = _manifest_dict(
        model_selection_workflow, "model_selection_summary"
    )
    maximum_likelihood_workflow = _manifest_dict(payload, "maximum_likelihood_workflow")
    bootstrap_workflow = _manifest_dict(payload, "bootstrap_workflow")
    support_summary = _manifest_dict(payload, "support_summary")
    output_paths = _manifest_dict(payload, "output_paths")
    step_manifests = _manifest_dict(payload, "step_manifests")
    warnings = _tree_inference_warning_lines(payload)
    repaired_input_path = payload.get("prepared_input_path")
    raw_input_path = payload.get("input_path")
    raw_path_text = _manifest_path_text(raw_input_path)
    prepared_path_text = _manifest_path_text(repaired_input_path)
    input_repair = _manifest_nested_dict(payload, "input_repair")
    input_repair_lines = (
        [
            "- input FASTA was repaired before alignment",
            f"- prepared input path: `{prepared_path_text}`",
            "- normalized identifiers: "
            f"`{len(_manifest_list(input_repair, 'normalized_identifiers'))}`",
            "- removed invalid records: "
            f"`{len(_manifest_list(input_repair, 'removed_records'))}`",
        ]
        if input_repair is not None
        else [
            "- input FASTA was used directly without identifier or record repair",
            f"- prepared input path: `{prepared_path_text}`",
        ]
    )
    text = (
        "# Tree Inference Methods Summary\n\n"
        f"The workflow manifest `{manifest_path.name}` records one Bijux `{workflow_name}` run from raw FASTA input "
        f"`{Path(str(raw_input_path)).name}` through alignment, trimming, model selection, maximum-likelihood inference, "
        f"and bootstrap-supported tree finalization. The selected substitution model was `{selected_model}`, "
        f"the sequence type used for inference was `{payload.get('sequence_type')}`, and the final delivered tree was "
        f"`{Path(str(output_paths.get('tree'))).name}`.\n\n"
        "## Input And Alignment Preparation\n\n"
        f"- raw input path: `{raw_path_text}`\n"
        + "\n".join(input_repair_lines)
        + "\n"
        + f"- validated sequence count: `{input_summary.get('sequence_count')}`\n"
        + f"- total raw residue count: `{input_summary.get('total_residue_count')}`\n"
        + f"- inferred raw sequence alphabet: `{input_summary.get('inferred_alphabet')}`\n"
        + f"- alignment engine: `{_manifest_path_text(_manifest_dict(payload, 'alignment_workflow').get('engine_name'))}`\n"
        + f"- alignment mode: `{payload.get('alignment_mode')}`\n"
        + f"- aligned output path: `{_manifest_path_text(output_paths.get('alignment'))}`\n"
        + f"- trimming engine: `{_manifest_path_text(trimming_workflow.get('engine_name'))}`\n"
        + f"- trimming mode: `{payload.get('trimming_mode')}`\n"
        + f"- trimming gap threshold: `{format(float(payload.get('trim_gap_threshold', 0.0)), '.15g')}`\n"
        + f"- retained alignment length: `{trimming_summary.get('trimmed_alignment_length')}` of `{trimming_summary.get('input_alignment_length')}`\n"
        + f"- removed alignment sites: `{trimming_summary.get('removed_site_count')}`\n"
        + f"- trimmed alignment path: `{_manifest_path_text(output_paths.get('trimmed_alignment'))}`\n\n"
        "## Model Selection\n\n"
        + f"- model-selection engine: `{_manifest_path_text(model_selection_workflow.get('engine_name'))}`\n"
        + f"- candidate substitution models reviewed: `{model_selection_summary.get('candidate_count')}`\n"
        + f"- governing information criterion: `{_manifest_path_text(model_selection_summary.get('selected_criterion'))}`\n"
        + f"- selected substitution model: `{selected_model}`\n"
        + f"- iqtree random seed: `{payload.get('iqtree_seed')}`\n"
        + f"- iqtree threads: `{payload.get('iqtree_threads')}`\n"
        + f"- model-selection manifest: `{_manifest_path_text(step_manifests.get('model_selection'))}`\n\n"
        "## Maximum-Likelihood Inference\n\n"
        + f"- inference engine: `{_manifest_path_text(maximum_likelihood_workflow.get('engine_name'))}`\n"
        + f"- inference model: `{selected_model}`\n"
        + f"- maximum-likelihood log-likelihood: `{_manifest_path_text(maximum_likelihood_workflow.get('log_likelihood'))}`\n"
        + f"- unannotated maximum-likelihood tree artifact: `{_manifest_path_text(_output_paths_dict(maximum_likelihood_workflow).get('tree'))}`\n"
        + f"- inference manifest: `{_manifest_path_text(step_manifests.get('maximum_likelihood'))}`\n\n"
        "## Branch Support And Final Tree\n\n"
        + f"- support engine: `{_manifest_path_text(bootstrap_workflow.get('engine_name'))}`\n"
        + f"- support workflow: ultrafast bootstrap support on the same trimmed alignment under `{selected_model}`\n"
        + f"- bootstrap replicates: `{payload.get('bootstrap_replicates')}`\n"
        + f"- supported internal nodes: `{support_summary.get('supported_node_count')}` of `{support_summary.get('internal_node_count')}`\n"
        + f"- minimum/median/maximum support: `{_manifest_path_text(support_summary.get('minimum_support'))}` / `{_manifest_path_text(support_summary.get('median_support'))}` / `{_manifest_path_text(support_summary.get('maximum_support'))}`\n"
        + f"- weakly supported clade count: `{support_summary.get('weakly_supported_clade_count')}`\n"
        + f"- bootstrap-supported tree artifact: `{_manifest_path_text(_output_paths_dict(bootstrap_workflow).get('support_tree'))}`\n"
        + f"- bootstrap tree-set artifact: `{_manifest_path_text(_output_paths_dict(bootstrap_workflow).get('bootstrap_trees'))}`\n"
        + f"- support manifest: `{_manifest_path_text(step_manifests.get('bootstrap_support'))}`\n\n"
        "## Tree Processing And Traceability\n\n"
        + "- the final delivered tree is copied from the bootstrap-supported inference artifact so branch support remains attached to the reviewer-facing tree\n"
        + "- the workflow records a separate unannotated maximum-likelihood tree under engine-artifacts for audit and comparison\n"
        + "- no outgroup rooting, midpoint rerooting, or posterior summarization step is recorded in this fasta-to-tree workflow manifest\n"
        + f"- final tree path: `{_manifest_path_text(output_paths.get('tree'))}`\n"
        + f"- reviewer-facing model table: `{_manifest_path_text(output_paths.get('model_table'))}`\n"
        + f"- reviewer-facing support table: `{_manifest_path_text(output_paths.get('support_table'))}`\n"
        + f"- workflow manifest path: `{manifest_path}`\n\n"
        "## Workflow Warnings\n\n"
        + (
            "\n".join(f"- {warning}" for warning in warnings)
            if warnings
            else "- no workflow-level warning was recorded in the current manifest"
        )
        + "\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return TreeInferenceMethodsSummaryTextResult(
        output_path=path,
        title="Tree Inference Methods Summary",
        warning_count=len(warnings),
        selected_model=selected_model,
        bootstrap_replicates=int(payload.get("bootstrap_replicates", 0)),
        trimmed_alignment_length=int(
            trimming_summary.get("trimmed_alignment_length", 0)
        ),
        supported_node_count=int(support_summary.get("supported_node_count", 0)),
        text=text,
        workflow_manifest_path=manifest_path,
        workflow_manifest=payload,
    )
