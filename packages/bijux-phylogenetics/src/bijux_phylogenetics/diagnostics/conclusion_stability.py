from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import json
from pathlib import Path
import statistics

from bijux_phylogenetics.ancestral.discrete import DiscreteAncestralReport
from bijux_phylogenetics.ancestral.tree_set import DiscreteAncestralTreeSetReport
from bijux_phylogenetics.comparative.pgls import PGLSResult
from bijux_phylogenetics.comparative.pgls.posterior_tree import (
    PosteriorTreePGLSCoefficientSummaryRow,
    PosteriorTreePGLSReport,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.trees import (
    CladeFrequencyReport,
    CladeTableReport,
    CladeTableRow,
)

_STABLE_SCORE_THRESHOLD = 0.85
_WEAK_SCORE_THRESHOLD = 0.6


@dataclass(frozen=True, slots=True)
class KeyCladeStabilityRow:
    """One rooted clade scored across bootstrap trees and method variants."""

    clade_id: str
    descendant_taxa: tuple[str, ...]
    bootstrap_frequency: float
    method_presence_fraction: float
    combined_score: float
    stability_class: str
    evidence: str


@dataclass(frozen=True, slots=True)
class SupportValueStabilityRow:
    """One rooted clade support surface scored across bootstrap and method variants."""

    clade_id: str
    descendant_taxa: tuple[str, ...]
    baseline_support_fraction: float | None
    bootstrap_frequency: float
    method_presence_fraction: float
    mean_method_support_fraction: float | None
    maximum_support_delta: float | None
    method_support_consistency: float
    combined_score: float
    stability_class: str
    evidence: str


@dataclass(frozen=True, slots=True)
class AncestralStateStabilityRow:
    """One ancestral-state claim scored across bootstrap trees and method variants."""

    trait: str
    clade_id: str
    descendant_taxa: tuple[str, ...]
    baseline_state: str
    bootstrap_dominant_state: str | None
    bootstrap_state_consistency: float
    method_dominant_state: str | None
    method_state_consistency: float
    combined_score: float
    stability_class: str
    evidence: str


@dataclass(frozen=True, slots=True)
class ComparativeCoefficientStabilityRow:
    """One comparative coefficient claim scored across tree sets and method variants."""

    term: str
    baseline_estimate: float
    baseline_direction: str
    baseline_significant: bool
    bootstrap_direction_consistency: float
    bootstrap_significance_consistency: float
    method_direction_consistency: float
    method_significance_consistency: float
    combined_score: float
    stability_class: str
    evidence: str


@dataclass(frozen=True, slots=True)
class ConclusionStabilityConclusionRow:
    """One reviewer-facing conclusion scored and classified for stability."""

    category: str
    conclusion_id: str
    label: str
    combined_score: float
    stability_class: str
    evidence: str


@dataclass(frozen=True, slots=True)
class ConclusionStabilitySummary:
    """Compact summary across all scored conclusion surfaces."""

    stable_count: int
    weak_count: int
    unstable_count: int
    key_clade_count: int
    support_value_count: int
    ancestral_state_count: int
    comparative_coefficient_count: int


@dataclass(frozen=True, slots=True)
class ConclusionStabilityReport:
    """Stability scoring report across phylogenetic conclusion surfaces."""

    key_clade_rows: tuple[KeyCladeStabilityRow, ...]
    support_value_rows: tuple[SupportValueStabilityRow, ...]
    ancestral_state_rows: tuple[AncestralStateStabilityRow, ...]
    comparative_coefficient_rows: tuple[ComparativeCoefficientStabilityRow, ...]
    conclusion_rows: tuple[ConclusionStabilityConclusionRow, ...]
    summary: ConclusionStabilitySummary


def build_key_clade_stability_rows(
    *,
    baseline_clades: CladeTableReport,
    bootstrap_frequencies: CladeFrequencyReport,
    method_clade_reports: list[CladeTableReport],
) -> list[KeyCladeStabilityRow]:
    """Score rooted key clades across bootstrap trees and method variants."""
    baseline_rows = _informative_clade_rows(baseline_clades)
    bootstrap_frequency_by_clade = {
        row.clade: row.frequency for row in bootstrap_frequencies.clade_frequencies
    }
    method_clade_sets = [
        {row.clade_id for row in _informative_clade_rows(report)}
        for report in method_clade_reports
    ]
    method_tree_count = max(len(method_clade_sets), 1)
    rows: list[KeyCladeStabilityRow] = []
    for row in baseline_rows:
        method_presence_fraction = (
            sum(row.clade_id in clade_set for clade_set in method_clade_sets)
            / method_tree_count
        )
        bootstrap_frequency = bootstrap_frequency_by_clade.get(row.clade_id, 0.0)
        combined_score = _mean_score([bootstrap_frequency, method_presence_fraction])
        rows.append(
            KeyCladeStabilityRow(
                clade_id=row.clade_id,
                descendant_taxa=tuple(row.taxa),
                bootstrap_frequency=bootstrap_frequency,
                method_presence_fraction=method_presence_fraction,
                combined_score=combined_score,
                stability_class=_classify_stability(combined_score),
                evidence=(
                    f"bootstrap_frequency={_format_score(bootstrap_frequency)}; "
                    f"method_presence_fraction={_format_score(method_presence_fraction)}"
                ),
            )
        )
    rows.sort(key=lambda item: (item.combined_score, item.clade_id), reverse=True)
    return rows


def build_support_value_stability_rows(
    *,
    baseline_clades: CladeTableReport,
    bootstrap_frequencies: CladeFrequencyReport,
    method_clade_reports: list[CladeTableReport],
) -> list[SupportValueStabilityRow]:
    """Score clade support-value stability across bootstrap trees and method variants."""
    baseline_rows = _informative_clade_rows(baseline_clades)
    bootstrap_frequency_by_clade = {
        row.clade: row.frequency for row in bootstrap_frequencies.clade_frequencies
    }
    method_clade_maps = [_support_by_clade(report) for report in method_clade_reports]
    method_tree_count = max(len(method_clade_maps), 1)
    rows: list[SupportValueStabilityRow] = []
    for row in baseline_rows:
        baseline_support = row.support_fraction
        comparable_method_supports = [
            clade_map[row.clade_id]
            for clade_map in method_clade_maps
            if row.clade_id in clade_map and clade_map[row.clade_id] is not None
        ]
        method_presence_fraction = (
            sum(row.clade_id in clade_map for clade_map in method_clade_maps)
            / method_tree_count
        )
        if comparable_method_supports and baseline_support is not None:
            mean_abs_delta = statistics.fmean(
                abs(value - baseline_support) for value in comparable_method_supports
            )
            maximum_support_delta = max(
                abs(value - baseline_support) for value in comparable_method_supports
            )
        else:
            mean_abs_delta = 0.0
            maximum_support_delta = None
        method_support_consistency = max(
            0.0, method_presence_fraction * (1.0 - mean_abs_delta)
        )
        bootstrap_frequency = bootstrap_frequency_by_clade.get(row.clade_id, 0.0)
        combined_score = _mean_score([bootstrap_frequency, method_support_consistency])
        rows.append(
            SupportValueStabilityRow(
                clade_id=row.clade_id,
                descendant_taxa=tuple(row.taxa),
                baseline_support_fraction=baseline_support,
                bootstrap_frequency=bootstrap_frequency,
                method_presence_fraction=method_presence_fraction,
                mean_method_support_fraction=(
                    None
                    if not comparable_method_supports
                    else statistics.fmean(comparable_method_supports)
                ),
                maximum_support_delta=maximum_support_delta,
                method_support_consistency=method_support_consistency,
                combined_score=combined_score,
                stability_class=_classify_stability(combined_score),
                evidence=(
                    f"baseline_support={_format_optional_score(baseline_support)}; "
                    f"bootstrap_frequency={_format_score(bootstrap_frequency)}; "
                    f"method_support_consistency={_format_score(method_support_consistency)}"
                ),
            )
        )
    rows.sort(key=lambda item: (item.combined_score, item.clade_id), reverse=True)
    return rows


def build_ancestral_state_stability_rows(
    *,
    baseline_report: DiscreteAncestralReport,
    bootstrap_report: DiscreteAncestralTreeSetReport,
    method_reports: list[DiscreteAncestralReport],
) -> list[AncestralStateStabilityRow]:
    """Score discrete ancestral-state conclusions across bootstrap trees and method variants."""
    baseline_by_clade = _ancestral_state_by_clade(baseline_report)
    bootstrap_by_clade = {row.clade_id: row for row in bootstrap_report.clade_summaries}
    method_maps = [_ancestral_state_by_clade(report) for report in method_reports]
    method_tree_count = max(len(method_maps), 1)
    rows: list[AncestralStateStabilityRow] = []
    for clade_id, baseline in baseline_by_clade.items():
        bootstrap_row = bootstrap_by_clade.get(clade_id)
        if bootstrap_row is None or bootstrap_row.tree_presence_count == 0:
            bootstrap_consistency = 0.0
            bootstrap_dominant_state = None
        else:
            bootstrap_consistency = (
                bootstrap_row.state_distribution.get(baseline["state"], 0)
                / bootstrap_row.tree_presence_count
            )
            bootstrap_dominant_state = bootstrap_row.dominant_state
        method_states = [
            method_map[clade_id]["state"]
            for method_map in method_maps
            if clade_id in method_map
        ]
        method_state_consistency = (
            0.0
            if not method_states
            else sum(state == baseline["state"] for state in method_states)
            / method_tree_count
        )
        method_dominant_state = (
            None if not method_states else Counter(method_states).most_common(1)[0][0]
        )
        combined_score = _mean_score([bootstrap_consistency, method_state_consistency])
        rows.append(
            AncestralStateStabilityRow(
                trait=baseline_report.trait,
                clade_id=clade_id,
                descendant_taxa=tuple(baseline["taxa"]),
                baseline_state=baseline["state"],
                bootstrap_dominant_state=bootstrap_dominant_state,
                bootstrap_state_consistency=bootstrap_consistency,
                method_dominant_state=method_dominant_state,
                method_state_consistency=method_state_consistency,
                combined_score=combined_score,
                stability_class=_classify_stability(combined_score),
                evidence=(
                    f"baseline_state={baseline['state']}; "
                    f"bootstrap_state_consistency={_format_score(bootstrap_consistency)}; "
                    f"method_state_consistency={_format_score(method_state_consistency)}"
                ),
            )
        )
    rows.sort(
        key=lambda item: (item.combined_score, item.trait, item.clade_id),
        reverse=True,
    )
    return rows


def build_comparative_coefficient_stability_rows(
    *,
    baseline_result: PGLSResult,
    bootstrap_report: PosteriorTreePGLSReport,
    method_results: list[PGLSResult],
) -> list[ComparativeCoefficientStabilityRow]:
    """Score comparative coefficient conclusions across bootstrap trees and method variants."""
    baseline_by_term = {
        coefficient.name: coefficient
        for coefficient in baseline_result.coefficients
        if coefficient.name != "intercept"
    }
    bootstrap_by_term = {
        row.term: row for row in bootstrap_report.coefficient_summaries
    }
    method_by_term: dict[str, list[tuple[str, bool]]] = {}
    for result in method_results:
        for coefficient in result.coefficients:
            if coefficient.name == "intercept":
                continue
            method_by_term.setdefault(coefficient.name, []).append(
                (
                    _estimate_direction(coefficient.estimate),
                    coefficient.p_value <= 0.05,
                )
            )
    method_tree_count = max(len(method_results), 1)
    rows: list[ComparativeCoefficientStabilityRow] = []
    for term, coefficient in baseline_by_term.items():
        baseline_direction = _estimate_direction(coefficient.estimate)
        baseline_significant = coefficient.p_value <= 0.05
        bootstrap_summary = bootstrap_by_term.get(term)
        if bootstrap_summary is None:
            bootstrap_direction_consistency = 0.0
            bootstrap_significance_consistency = 0.0
        else:
            bootstrap_direction_consistency = _bootstrap_direction_fraction(
                bootstrap_summary, baseline_direction=baseline_direction
            )
            bootstrap_significance_consistency = _bootstrap_significance_fraction(
                bootstrap_summary, baseline_significant=baseline_significant
            )
        method_observations = method_by_term.get(term, [])
        if not method_observations:
            method_direction_consistency = 0.0
            method_significance_consistency = 0.0
        else:
            method_direction_consistency = (
                sum(
                    direction == baseline_direction
                    for direction, _significant in method_observations
                )
                / method_tree_count
            )
            method_significance_consistency = (
                sum(
                    significant == baseline_significant
                    for _direction, significant in method_observations
                )
                / method_tree_count
            )
        combined_score = _mean_score(
            [
                bootstrap_direction_consistency,
                bootstrap_significance_consistency,
                method_direction_consistency,
                method_significance_consistency,
            ]
        )
        rows.append(
            ComparativeCoefficientStabilityRow(
                term=term,
                baseline_estimate=coefficient.estimate,
                baseline_direction=baseline_direction,
                baseline_significant=baseline_significant,
                bootstrap_direction_consistency=bootstrap_direction_consistency,
                bootstrap_significance_consistency=bootstrap_significance_consistency,
                method_direction_consistency=method_direction_consistency,
                method_significance_consistency=method_significance_consistency,
                combined_score=combined_score,
                stability_class=_classify_stability(combined_score),
                evidence=(
                    f"baseline_direction={baseline_direction}; "
                    f"bootstrap_direction_consistency={_format_score(bootstrap_direction_consistency)}; "
                    f"method_direction_consistency={_format_score(method_direction_consistency)}"
                ),
            )
        )
    rows.sort(key=lambda item: (item.combined_score, item.term), reverse=True)
    return rows


def build_conclusion_stability_report(
    *,
    key_clade_rows: list[KeyCladeStabilityRow],
    support_value_rows: list[SupportValueStabilityRow],
    ancestral_state_rows: list[AncestralStateStabilityRow],
    comparative_coefficient_rows: list[ComparativeCoefficientStabilityRow],
) -> ConclusionStabilityReport:
    """Build one combined conclusion-stability ledger across all surfaces."""
    conclusion_rows: list[ConclusionStabilityConclusionRow] = []
    for row in key_clade_rows:
        conclusion_rows.append(
            ConclusionStabilityConclusionRow(
                category="key_clade",
                conclusion_id=row.clade_id,
                label=row.clade_id,
                combined_score=row.combined_score,
                stability_class=row.stability_class,
                evidence=row.evidence,
            )
        )
    for row in support_value_rows:
        conclusion_rows.append(
            ConclusionStabilityConclusionRow(
                category="support_value",
                conclusion_id=row.clade_id,
                label=row.clade_id,
                combined_score=row.combined_score,
                stability_class=row.stability_class,
                evidence=row.evidence,
            )
        )
    for row in ancestral_state_rows:
        conclusion_rows.append(
            ConclusionStabilityConclusionRow(
                category="ancestral_state",
                conclusion_id=f"{row.trait}:{row.clade_id}",
                label=f"{row.trait}::{row.baseline_state}",
                combined_score=row.combined_score,
                stability_class=row.stability_class,
                evidence=row.evidence,
            )
        )
    for row in comparative_coefficient_rows:
        conclusion_rows.append(
            ConclusionStabilityConclusionRow(
                category="comparative_coefficient",
                conclusion_id=row.term,
                label=row.term,
                combined_score=row.combined_score,
                stability_class=row.stability_class,
                evidence=row.evidence,
            )
        )
    summary = ConclusionStabilitySummary(
        stable_count=sum(row.stability_class == "stable" for row in conclusion_rows),
        weak_count=sum(row.stability_class == "weak" for row in conclusion_rows),
        unstable_count=sum(
            row.stability_class == "unstable" for row in conclusion_rows
        ),
        key_clade_count=len(key_clade_rows),
        support_value_count=len(support_value_rows),
        ancestral_state_count=len(ancestral_state_rows),
        comparative_coefficient_count=len(comparative_coefficient_rows),
    )
    return ConclusionStabilityReport(
        key_clade_rows=tuple(key_clade_rows),
        support_value_rows=tuple(support_value_rows),
        ancestral_state_rows=tuple(ancestral_state_rows),
        comparative_coefficient_rows=tuple(comparative_coefficient_rows),
        conclusion_rows=tuple(
            sorted(
                conclusion_rows,
                key=lambda row: (row.stability_class, row.combined_score, row.category),
                reverse=True,
            )
        ),
        summary=summary,
    )


def write_key_clade_stability_table(
    path: Path, rows: tuple[KeyCladeStabilityRow, ...]
) -> Path:
    """Write key-clade conclusion stability rows as TSV."""
    return write_taxon_rows(
        path,
        columns=[
            "clade_id",
            "descendant_taxa",
            "bootstrap_frequency",
            "method_presence_fraction",
            "combined_score",
            "stability_class",
            "evidence",
        ],
        rows=[
            {
                "clade_id": row.clade_id,
                "descendant_taxa": ",".join(row.descendant_taxa),
                "bootstrap_frequency": _format_score(row.bootstrap_frequency),
                "method_presence_fraction": _format_score(row.method_presence_fraction),
                "combined_score": _format_score(row.combined_score),
                "stability_class": row.stability_class,
                "evidence": row.evidence,
            }
            for row in rows
        ],
    )


def write_support_value_stability_table(
    path: Path, rows: tuple[SupportValueStabilityRow, ...]
) -> Path:
    """Write support-value stability rows as TSV."""
    return write_taxon_rows(
        path,
        columns=[
            "clade_id",
            "descendant_taxa",
            "baseline_support_fraction",
            "bootstrap_frequency",
            "method_presence_fraction",
            "mean_method_support_fraction",
            "maximum_support_delta",
            "method_support_consistency",
            "combined_score",
            "stability_class",
            "evidence",
        ],
        rows=[
            {
                "clade_id": row.clade_id,
                "descendant_taxa": ",".join(row.descendant_taxa),
                "baseline_support_fraction": _format_optional_score(
                    row.baseline_support_fraction
                ),
                "bootstrap_frequency": _format_score(row.bootstrap_frequency),
                "method_presence_fraction": _format_score(row.method_presence_fraction),
                "mean_method_support_fraction": _format_optional_score(
                    row.mean_method_support_fraction
                ),
                "maximum_support_delta": _format_optional_score(
                    row.maximum_support_delta
                ),
                "method_support_consistency": _format_score(
                    row.method_support_consistency
                ),
                "combined_score": _format_score(row.combined_score),
                "stability_class": row.stability_class,
                "evidence": row.evidence,
            }
            for row in rows
        ],
    )


def write_ancestral_state_stability_table(
    path: Path, rows: tuple[AncestralStateStabilityRow, ...]
) -> Path:
    """Write ancestral-state stability rows as TSV."""
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "clade_id",
            "descendant_taxa",
            "baseline_state",
            "bootstrap_dominant_state",
            "bootstrap_state_consistency",
            "method_dominant_state",
            "method_state_consistency",
            "combined_score",
            "stability_class",
            "evidence",
        ],
        rows=[
            {
                "trait": row.trait,
                "clade_id": row.clade_id,
                "descendant_taxa": ",".join(row.descendant_taxa),
                "baseline_state": row.baseline_state,
                "bootstrap_dominant_state": row.bootstrap_dominant_state or "",
                "bootstrap_state_consistency": _format_score(
                    row.bootstrap_state_consistency
                ),
                "method_dominant_state": row.method_dominant_state or "",
                "method_state_consistency": _format_score(row.method_state_consistency),
                "combined_score": _format_score(row.combined_score),
                "stability_class": row.stability_class,
                "evidence": row.evidence,
            }
            for row in rows
        ],
    )


def write_comparative_coefficient_stability_table(
    path: Path, rows: tuple[ComparativeCoefficientStabilityRow, ...]
) -> Path:
    """Write comparative coefficient stability rows as TSV."""
    return write_taxon_rows(
        path,
        columns=[
            "term",
            "baseline_estimate",
            "baseline_direction",
            "baseline_significant",
            "bootstrap_direction_consistency",
            "bootstrap_significance_consistency",
            "method_direction_consistency",
            "method_significance_consistency",
            "combined_score",
            "stability_class",
            "evidence",
        ],
        rows=[
            {
                "term": row.term,
                "baseline_estimate": _format_score(row.baseline_estimate),
                "baseline_direction": row.baseline_direction,
                "baseline_significant": str(row.baseline_significant).lower(),
                "bootstrap_direction_consistency": _format_score(
                    row.bootstrap_direction_consistency
                ),
                "bootstrap_significance_consistency": _format_score(
                    row.bootstrap_significance_consistency
                ),
                "method_direction_consistency": _format_score(
                    row.method_direction_consistency
                ),
                "method_significance_consistency": _format_score(
                    row.method_significance_consistency
                ),
                "combined_score": _format_score(row.combined_score),
                "stability_class": row.stability_class,
                "evidence": row.evidence,
            }
            for row in rows
        ],
    )


def write_conclusion_stability_summary_table(
    path: Path, report: ConclusionStabilityReport
) -> Path:
    """Write one compact summary of stability-class counts."""
    row = {
        "stable_count": str(report.summary.stable_count),
        "weak_count": str(report.summary.weak_count),
        "unstable_count": str(report.summary.unstable_count),
        "key_clade_count": str(report.summary.key_clade_count),
        "support_value_count": str(report.summary.support_value_count),
        "ancestral_state_count": str(report.summary.ancestral_state_count),
        "comparative_coefficient_count": str(
            report.summary.comparative_coefficient_count
        ),
    }
    return write_taxon_rows(path, columns=list(row.keys()), rows=[row])


def write_conclusion_stability_report_html(
    path: Path, report: ConclusionStabilityReport
) -> Path:
    """Render one reviewer-facing HTML stability report."""
    stable_rows = [
        row for row in report.conclusion_rows if row.stability_class == "stable"
    ]
    weak_rows = [row for row in report.conclusion_rows if row.stability_class == "weak"]
    unstable_rows = [
        row for row in report.conclusion_rows if row.stability_class == "unstable"
    ]
    return write_html_report(
        title="Bijux Conclusion Stability Report",
        sections=[
            ("stable-conclusions", _serialize_html_section(stable_rows)),
            ("weak-conclusions", _serialize_html_section(weak_rows)),
            (
                "unstable-conclusions",
                _serialize_html_section(unstable_rows),
            ),
        ],
        out_path=path,
        embedded_json={
            "summary": {
                "stable_count": report.summary.stable_count,
                "weak_count": report.summary.weak_count,
                "unstable_count": report.summary.unstable_count,
            },
            "stable_conclusions": [_asdict_conclusion(row) for row in stable_rows],
            "weak_conclusions": [_asdict_conclusion(row) for row in weak_rows],
            "unstable_conclusions": [_asdict_conclusion(row) for row in unstable_rows],
        },
        summary_metrics=[
            ("stable_count", report.summary.stable_count),
            ("weak_count", report.summary.weak_count),
            ("unstable_count", report.summary.unstable_count),
            ("key_clade_count", report.summary.key_clade_count),
            ("support_value_count", report.summary.support_value_count),
            ("ancestral_state_count", report.summary.ancestral_state_count),
            (
                "comparative_coefficient_count",
                report.summary.comparative_coefficient_count,
            ),
        ],
    )


def _informative_clade_rows(report: CladeTableReport) -> list[CladeTableRow]:
    total_taxa = max((row.taxon_count for row in report.rows), default=0)
    return [
        row
        for row in report.rows
        if row.node_kind != "tip" and 1 < row.taxon_count < total_taxa
    ]


def _support_by_clade(report: CladeTableReport) -> dict[str, float | None]:
    return {
        row.clade_id: row.support_fraction for row in _informative_clade_rows(report)
    }


def _ancestral_state_by_clade(
    report: DiscreteAncestralReport,
) -> dict[str, dict[str, object]]:
    rows: dict[str, dict[str, object]] = {}
    for estimate in report.estimates:
        if estimate.is_tip:
            continue
        clade_id = _clade_id_from_taxa(estimate.descendant_taxa)
        rows[clade_id] = {
            "state": estimate.most_likely_state,
            "taxa": list(estimate.descendant_taxa),
        }
    return rows


def _bootstrap_direction_fraction(
    row: PosteriorTreePGLSCoefficientSummaryRow, *, baseline_direction: str
) -> float:
    if baseline_direction == "positive":
        return row.positive_tree_count / row.tree_fit_count
    if baseline_direction == "negative":
        return row.negative_tree_count / row.tree_fit_count
    return row.zero_tree_count / row.tree_fit_count


def _bootstrap_significance_fraction(
    row: PosteriorTreePGLSCoefficientSummaryRow, *, baseline_significant: bool
) -> float:
    significant_fraction = row.significant_tree_count / row.tree_fit_count
    return significant_fraction if baseline_significant else 1.0 - significant_fraction


def _estimate_direction(estimate: float) -> str:
    if abs(estimate) <= 1e-12:
        return "zero"
    return "positive" if estimate > 0.0 else "negative"


def _classify_stability(score: float) -> str:
    if score >= _STABLE_SCORE_THRESHOLD:
        return "stable"
    if score >= _WEAK_SCORE_THRESHOLD:
        return "weak"
    return "unstable"


def _clade_id_from_taxa(taxa: list[str]) -> str:
    return "|".join(sorted(taxa))


def _mean_score(values: list[float]) -> float:
    return statistics.fmean(values) if values else 0.0


def _format_score(value: float) -> str:
    return format(value, ".15g")


def _format_optional_score(value: float | None) -> str:
    return "" if value is None else _format_score(value)


def _asdict_conclusion(row: ConclusionStabilityConclusionRow) -> dict[str, object]:
    return {
        "category": row.category,
        "conclusion_id": row.conclusion_id,
        "label": row.label,
        "combined_score": row.combined_score,
        "stability_class": row.stability_class,
        "evidence": row.evidence,
    }


def _serialize_html_section(rows: list[ConclusionStabilityConclusionRow]) -> str:
    return json.dumps(
        [_asdict_conclusion(row) for row in rows],
        indent=2,
        sort_keys=True,
    )
