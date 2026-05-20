from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

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
from bijux_phylogenetics.io.fasta.cleaning import clean_alignment_with_profile


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
