from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

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
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    APE_ULTRAMETRIC_TOLERANCE,
)


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
