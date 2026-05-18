from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from bijux_phylogenetics.core.alignment import AlignmentCleaningReport
from bijux_phylogenetics.core.ultrametric import APE_ULTRAMETRIC_TOLERANCE
from bijux_phylogenetics.diagnostics.validation import (
    LONG_BRANCH_OUTLIER_FACTOR,
    SHORT_BRANCH_OUTLIER_FACTOR,
    STAR_LIKE_FRACTION_THRESHOLD,
    TREE_IMBALANCE_WARNING_THRESHOLD,
    TreeForensicReport,
    TreeInspectionReport,
    TreeValidationReport,
    forensic_tree_path,
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.engines.common import load_engine_manifest
from bijux_phylogenetics.io.fasta import clean_alignment_with_profile


@dataclass(slots=True)
class AlignmentFilteringMethodsSummaryTextResult:
    output_path: Path
    title: str
    warning_count: int
    removed_site_count: int
    removed_sequence_count: int
    retained_sequence_count: int
    retained_alignment_length: int
    text: str
    cleaning: AlignmentCleaningReport


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


@dataclass(slots=True)
class TreeValidationMethodsSummaryTextResult:
    output_path: Path
    title: str
    warning_count: int
    blocked_context_count: int
    repair_item_count: int
    text: str
    validation: TreeValidationReport
    inspection: TreeInspectionReport
    forensic: TreeForensicReport


def _format_context_name(raw: str) -> str:
    return raw.replace("_", " ")


def _count_reason_values(values: list[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    return counts


def _filtering_policy_lines(cleaning: AlignmentCleaningReport) -> list[str]:
    profile = cleaning.profile
    lines = [
        f"- named profile: `{profile.name}`",
        f"- profile note: {profile.note}",
        (
            "- all-gap columns are removed"
            if profile.remove_all_gap_sites
            else "- all-gap columns are retained"
        ),
        (
            "- all-missing columns are removed"
            if profile.remove_all_missing_sites
            else "- all-missing columns are retained"
        ),
        (
            "- codon phase is preserved by expanding removed positions to full codons"
            if profile.preserve_codon_structure
            else "- codon phase preservation is disabled"
        ),
    ]
    if profile.site_missingness_threshold is not None:
        lines.append(
            "- site missingness threshold: "
            + f"`{format(profile.site_missingness_threshold, '.15g')}`"
        )
    else:
        lines.append("- site missingness threshold: not applied")
    if profile.sequence_missingness_threshold is not None:
        lines.append(
            "- sequence missingness threshold: "
            + f"`{format(profile.sequence_missingness_threshold, '.15g')}`"
        )
    else:
        lines.append("- sequence missingness threshold: not applied")
    return lines


def _filtering_removal_lines(cleaning: AlignmentCleaningReport) -> list[str]:
    removed_columns = cleaning.trim.removed_columns
    removed_sequences = cleaning.trim.removed_sequences
    column_reason_counts = _count_reason_values(
        [row.reason for row in removed_columns]
    )
    sequence_reason_counts = _count_reason_values(
        [row.reason for row in removed_sequences]
    )
    lines: list[str] = []
    if removed_columns:
        lines.append(
            "- removed sites: "
            + f"`{len(removed_columns)}`"
            + " ("
            + ", ".join(
                f"{reason}={count}"
                for reason, count in sorted(column_reason_counts.items())
            )
            + ")"
        )
    else:
        lines.append("- removed sites: `0`")
    if removed_sequences:
        sequence_details = ", ".join(
            f"`{row.identifier}` ({row.reason})" for row in removed_sequences
        )
        lines.append(
            "- removed sequences: "
            + f"`{len(removed_sequences)}`"
            + " ("
            + ", ".join(
                f"{reason}={count}"
                for reason, count in sorted(sequence_reason_counts.items())
            )
            + ")"
        )
        lines.append(f"- removed sequence identities: {sequence_details}")
    else:
        lines.append("- removed sequences: `0`")
    return lines


def _group_retention_lines(cleaning: AlignmentCleaningReport) -> list[str]:
    if not cleaning.group_retention:
        return ["- no metadata or trait group-retention audit was requested"]
    return [
        "- group retention "
        + f"`{row.column}={row.value}`: original `{row.original_count}`, retained `{row.retained_count}`, removed `{row.removed_count}`, removed fraction `{format(row.removed_fraction, '.15g')}`"
        for row in cleaning.group_retention
    ]


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
    assert workflow_manifest_path is not None
    return load_engine_manifest(workflow_manifest_path), workflow_manifest_path


def _tree_inference_warning_lines(payload: dict[str, object]) -> list[str]:
    support_summary = _manifest_nested_dict(payload, "support_summary")
    warning_texts = [
        str(item)
        for item in [
            *_manifest_list(payload, "warnings"),
            *([] if support_summary is None else _manifest_list(support_summary, "warnings")),
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
    maximum_likelihood_workflow = _manifest_dict(
        payload, "maximum_likelihood_workflow"
    )
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
            + f"`{len(_manifest_list(input_repair, 'normalized_identifiers'))}`",
            "- removed invalid records: "
            + f"`{len(_manifest_list(input_repair, 'removed_records'))}`",
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


def write_alignment_filtering_methods_summary_text(
    path: Path,
    *,
    alignment_path: Path,
    profile_name: str,
    group_table_path: Path | None = None,
    group_columns: list[str] | None = None,
) -> AlignmentFilteringMethodsSummaryTextResult:
    """Write reviewer-facing methods text for one profile-driven alignment cleaning pass."""
    _, cleaning = clean_alignment_with_profile(
        alignment_path,
        profile_name=profile_name,
        group_table_path=group_table_path,
        group_columns=group_columns,
    )
    trim = cleaning.trim
    comparison = cleaning.comparison
    warning_count = (
        len(cleaning.signal_warnings)
        + len(cleaning.warnings)
        + len(comparison.warnings)
    )
    text = (
        "# Alignment Filtering Methods Summary\n\n"
        f"The alignment `{alignment_path.name}` was cleaned with the Bijux profile-driven filtering surface using profile "
        f"`{cleaning.profile.name}`. The cleaning pass starts from the original aligned matrix and records all removed sites, "
        f"removed sequences, before-versus-after composition shifts, phylogenetic signal loss warnings, and optional group-retention bias checks.\n\n"
        "## Filtering Policy\n\n"
        + "\n".join(_filtering_policy_lines(cleaning))
        + "\n\n## Removed Content\n\n"
        + "\n".join(_filtering_removal_lines(cleaning))
        + "\n\n## Retained Dimensions\n\n"
        + f"- retained sequence count: `{trim.trimmed_sequence_count}` of `{trim.original_sequence_count}`\n"
        + f"- retained alignment length: `{trim.trimmed_alignment_length}` of `{trim.original_alignment_length}`\n"
        + f"- retained shared taxa: `{len(comparison.shared_taxa)}`\n"
        + f"- variable-site count before/after: `{comparison.left_variable_site_count}` -> `{comparison.right_variable_site_count}`\n"
        + f"- parsimony-informative-site count before/after: `{comparison.left_parsimony_informative_site_count}` -> `{comparison.right_parsimony_informative_site_count}`\n"
        + f"- missing-data fraction before/after: `{format(comparison.left_missing_data_fraction, '.15g')}` -> `{format(comparison.right_missing_data_fraction, '.15g')}`\n"
        + f"- gap fraction before/after: `{format(comparison.left_gap_fraction, '.15g')}` -> `{format(comparison.right_gap_fraction, '.15g')}`\n\n"
        "## Signal And Bias Checks\n\n"
        + (
            "\n".join(f"- {warning.message}" for warning in cleaning.signal_warnings)
            if cleaning.signal_warnings
            else "- no explicit phylogenetic signal-collapse warning was raised by the current cleaning pass"
        )
        + "\n"
        + (
            "\n".join(f"- {warning}" for warning in cleaning.warnings)
            if cleaning.warnings
            else "- no additional cleaning workflow warning was raised"
        )
        + "\n"
        + (
            "\n".join(f"- {warning}" for warning in comparison.warnings)
            if comparison.warnings
            else "- no cross-version alignment comparison warning was raised"
        )
        + "\n\n## Group Retention\n\n"
        + "\n".join(_group_retention_lines(cleaning))
        + "\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return AlignmentFilteringMethodsSummaryTextResult(
        output_path=path,
        title="Alignment Filtering Methods Summary",
        warning_count=warning_count,
        removed_site_count=len(trim.removed_columns),
        removed_sequence_count=len(trim.removed_sequences),
        retained_sequence_count=trim.trimmed_sequence_count,
        retained_alignment_length=trim.trimmed_alignment_length,
        text=text,
        cleaning=cleaning,
    )


def _repair_items(
    *,
    validation: TreeValidationReport,
    inspection: TreeInspectionReport,
) -> list[str]:
    items: list[str] = []
    if validation.missing_taxa:
        items.append(
            f"{validation.missing_taxa} unnamed tip label(s) require repair before deterministic downstream joins"
        )
    if validation.duplicate_taxa:
        items.append(
            "duplicate tip labels detected: "
            + ", ".join(f"`{label}`" for label in validation.duplicate_taxa)
        )
    if validation.unsafe_external_labels:
        items.append(
            "downstream-unsafe tip labels detected: "
            + ", ".join(
                f"`{label.raw_label}`" for label in validation.unsafe_external_labels
            )
        )
    if inspection.internal_label_conflicts:
        items.append(
            f"{len(inspection.internal_label_conflicts)} internal label conflict(s) require label interpretation review before publication-facing support summaries"
        )
    return items


def _downstream_consequences(forensic: TreeForensicReport) -> list[str]:
    context_lines: list[str] = []
    for context in forensic.branch_length_contexts:
        if context.allowed:
            detail = "allowed"
            if context.warnings:
                detail += " with warnings: " + "; ".join(context.warnings)
        else:
            detail = "blocked by " + "; ".join(context.blocked_by)
            if context.warnings:
                detail += " (warnings: " + "; ".join(context.warnings) + ")"
        context_lines.append(f"- `{_format_context_name(context.context)}`: {detail}")
    context_lines.append(
        "- `topology comparison`: "
        + (
            "allowed"
            if forensic.safe_for_topology_comparison
            else "blocked by syntax, duplicate-label, or unnamed-tip failures"
        )
    )
    context_lines.append(
        "- `visualization`: "
        + (
            "allowed"
            if forensic.safe_for_visualization
            else "blocked by syntax failures"
        )
    )
    context_lines.append(
        "- `publication`: "
        + (
            "allowed"
            if forensic.safe_for_publication
            else "blocked by biological safety or internal-label conflicts"
        )
    )
    return context_lines


def write_tree_validation_methods_summary_text(
    path: Path,
    *,
    tree_path: Path,
    source_format: str | None = None,
) -> TreeValidationMethodsSummaryTextResult:
    """Write reviewer-facing methods text for the tree-validation surface."""
    validation = validate_tree_path(tree_path, source_format=source_format)
    inspection = inspect_tree_path(tree_path, source_format=source_format)
    forensic = forensic_tree_path(tree_path, source_format=source_format)
    blocked_context_count = sum(
        1 for context in forensic.branch_length_contexts if not context.allowed
    )
    repair_items = _repair_items(validation=validation, inspection=inspection)
    warning_count = (
        len(validation.warnings)
        + len(inspection.warnings)
        + len(forensic.warnings)
        + len(repair_items)
    )
    consequences = _downstream_consequences(forensic)
    text = (
        "# Tree Validation Methods Summary\n\n"
        f"The tree `{tree_path.name}` was reviewed with Bijux tree validation, tree inspection, "
        f"and forensic downstream-safety surfaces. The parsed source format was `{validation.source_format}`, "
        f"the validation decision was `{validation.validity_decision}`, and the tree covered `{validation.tip_count}` tip(s) "
        f"and `{validation.internal_node_count}` internal node(s). Validation does not silently prune or repair taxa; "
        f"instead it records the exact blockers, warning details, and downstream analyses that remain unsafe.\n\n"
        "## Checks Performed\n\n"
        f"- structural integrity checks: syntax validity, cycle and duplicate-parentage detection, empty or degenerate roots, "
        f"duplicate tip labels, unnamed tips, and singleton internal nodes\n"
        f"- branch-length review: branch-length status `{validation.branch_length_status}`, zero-length branch count "
        f"`{validation.zero_length_branches}`, negative branch count `{validation.negative_branch_lengths}`, and missing internal or terminal lengths\n"
        f"- rooting and time-scale review: rooted classification `{validation.rooted}`, ultrametric status `{validation.ultrametric}`, "
        f"and root-state confidence `{validation.root_state_confidence.classification}`\n"
        f"- topology-resolution review: polytomy count `{validation.polytomy_count}` and biologically safe flag `{validation.biologically_safe}`\n"
        f"- internal-label audit: `{len(inspection.likely_support_labels)}` support-like internal label(s), "
        f"`{len(inspection.likely_named_internal_labels)}` named internal label(s), and "
        f"`{len(inspection.internal_label_conflicts)}` conflict(s)\n"
        f"- taxon-safety audit: `{len(validation.unsafe_external_labels)}` downstream-unsafe label(s) and "
        f"`{len(validation.taxon_identity_audit.suspicious_near_duplicates)}` suspicious near-duplicate label pair(s)\n"
        f"- shape heuristics: tree quality score `{inspection.tree_quality_score}`, long-branch outlier count "
        f"`{len(inspection.long_branch_outliers)}`, short-branch outlier count `{len(inspection.short_branch_outliers)}`, "
        f"star-like `{str(inspection.star_like).lower()}`, and comb-like `{str(inspection.comb_like).lower()}`\n\n"
        "## Thresholds\n\n"
        f"- ultrametricity is evaluated with the APE-compatible tolerance `{format(APE_ULTRAMETRIC_TOLERANCE, '.15g')}`\n"
        f"- unusually imbalanced trees are flagged when normalized Colless imbalance is at least "
        f"`{format(TREE_IMBALANCE_WARNING_THRESHOLD, '.15g')}`\n"
        f"- long terminal branches are flagged when they exceed `{format(LONG_BRANCH_OUTLIER_FACTOR, '.15g')}x` the median positive terminal branch length\n"
        f"- short nonzero branches are flagged when they fall below `{format(SHORT_BRANCH_OUTLIER_FACTOR, '.15g')}x` the median positive branch length\n"
        f"- star-like topologies are flagged when one node directly subtends at least `max(4, ceil({format(STAR_LIKE_FRACTION_THRESHOLD, '.15g')} * tip_count))` leaf children\n"
        "- support-like internal labels are interpreted against standard probability (`0-1`) or percentage (`1-100`) scales, and mixed or out-of-range scales are flagged\n\n"
        "## Exclusions And Repairs\n\n"
        + (
            "\n".join(f"- {item}" for item in repair_items)
            if repair_items
            else "- no taxa were excluded or flagged for repair by the current validation pass"
        )
        + "\n\n## Downstream Consequences\n\n"
        + "\n".join(consequences)
        + "\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return TreeValidationMethodsSummaryTextResult(
        output_path=path,
        title="Tree Validation Methods Summary",
        warning_count=warning_count,
        blocked_context_count=blocked_context_count,
        repair_item_count=len(repair_items),
        text=text,
        validation=validation,
        inspection=inspection,
        forensic=forensic,
    )
