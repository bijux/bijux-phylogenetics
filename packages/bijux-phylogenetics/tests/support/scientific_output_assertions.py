from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.validation import (
    ScientificOutputEquivalenceReport,
    compare_scientific_output,
)


def assert_scientific_outputs_equivalent(
    expected_root: Path,
    observed_root: Path,
) -> ScientificOutputEquivalenceReport:
    report = compare_scientific_output(expected_root, observed_root)
    assert report.equivalent, _format_equivalence_report(report)
    return report


def assert_selected_scientific_outputs_equivalent(
    expected_root: Path,
    observed_files: dict[str | Path, Path],
) -> list[ScientificOutputEquivalenceReport]:
    reports: list[ScientificOutputEquivalenceReport] = []
    for relative_name, observed_path in observed_files.items():
        relative_path = Path(relative_name)
        report = compare_scientific_output(expected_root / relative_path, observed_path)
        assert report.equivalent, _format_equivalence_report(report)
        reports.append(report)
    return reports


def _format_equivalence_report(report: ScientificOutputEquivalenceReport) -> str:
    lines = [
        "scientific output equivalence failed:",
        f"expected: {report.expected_root}",
        f"observed: {report.observed_root}",
    ]
    for issue in report.issues:
        lines.append(f"- {issue.relative_path}: {issue.kind}: {issue.detail}")
    return "\n".join(lines)
