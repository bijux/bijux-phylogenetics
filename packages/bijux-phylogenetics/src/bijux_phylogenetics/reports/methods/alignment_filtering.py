from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.io.fasta.cleaning import clean_alignment_with_profile
from bijux_phylogenetics.phylo.alignment import AlignmentCleaningReport


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
            f"`{format(profile.site_missingness_threshold, '.15g')}`"
        )
    else:
        lines.append("- site missingness threshold: not applied")
    if profile.sequence_missingness_threshold is not None:
        lines.append(
            "- sequence missingness threshold: "
            f"`{format(profile.sequence_missingness_threshold, '.15g')}`"
        )
    else:
        lines.append("- sequence missingness threshold: not applied")
    return lines


def _filtering_removal_lines(cleaning: AlignmentCleaningReport) -> list[str]:
    removed_columns = cleaning.trim.removed_columns
    removed_sequences = cleaning.trim.removed_sequences
    column_reason_counts = _count_reason_values([row.reason for row in removed_columns])
    sequence_reason_counts = _count_reason_values(
        [row.reason for row in removed_sequences]
    )
    lines: list[str] = []
    if removed_columns:
        lines.append(
            "- removed sites: "
            f"`{len(removed_columns)}`"
            " ("
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
            f"`{len(removed_sequences)}`"
            " ("
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
        f"`{row.column}={row.value}`: original `{row.original_count}`, retained `{row.retained_count}`, removed `{row.removed_count}`, removed fraction `{format(row.removed_fraction, '.15g')}`"
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
