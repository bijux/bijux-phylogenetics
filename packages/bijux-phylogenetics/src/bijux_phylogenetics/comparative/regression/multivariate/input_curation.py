from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.pgls import (
    PGLSInputReport,
    PGLSTaxonExclusion,
    inspect_pgls_inputs,
)
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .contracts import MultivariateTaxonExclusion


def build_response_formula(response: str, predictors: list[str]) -> str:
    return f"{response} ~ {' + '.join(predictors)}"


def inspect_response_model(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str,
    predictors: list[str],
    taxon_column: str,
) -> PGLSInputReport:
    return inspect_pgls_inputs(
        tree_path,
        traits_path,
        formula=build_response_formula(response, predictors),
        taxon_column=taxon_column,
    )


def raise_for_input_blockers(reports: list[PGLSInputReport]) -> None:
    blockers: list[str] = []
    for report in reports:
        for blocker in report.blockers:
            blockers.append(f"{report.formula.response}: {blocker}")
    if blockers:
        raise ComparativeMethodError("; ".join(blockers))


def shared_analysis_taxa(reports: list[PGLSInputReport]) -> list[str]:
    if not reports:
        return []
    shared = set(reports[0].analysis_taxa)
    for report in reports[1:]:
        shared &= set(report.analysis_taxa)
    return sorted(shared)


def build_shared_taxon_exclusions(
    *,
    overlap_taxa: list[str],
    responses: list[str],
    overlap_reports: list[PGLSInputReport],
) -> list[MultivariateTaxonExclusion]:
    exclusions_by_response: dict[str, dict[str, PGLSTaxonExclusion]] = {
        report.formula.response: {
            exclusion.taxon: exclusion
            for exclusion in report.formula_audit.excluded_taxa
        }
        for report in overlap_reports
    }
    excluded_rows: list[MultivariateTaxonExclusion] = []
    for taxon in overlap_taxa:
        missing_columns: set[str] = set()
        blocking_responses: list[str] = []
        invalid_details: list[str] = []
        missing_details: list[str] = []
        other_details: list[str] = []
        for response in responses:
            exclusion = exclusions_by_response.get(response, {}).get(taxon)
            if exclusion is None:
                continue
            blocking_responses.append(response)
            if exclusion.reason == "missing_value":
                parsed_missing = parse_missing_columns(exclusion.details)
                missing_columns.update(parsed_missing)
                missing_details.append(exclusion.details)
            elif exclusion.reason == "non_numeric_or_invalid_value":
                invalid_details.append(exclusion.details)
            else:
                other_details.append(exclusion.details)
        if not blocking_responses:
            continue
        if invalid_details:
            reason = "invalid_required_values"
            details = "; ".join(sorted(set(invalid_details)))
        elif missing_columns:
            reason = "missing_required_values"
            details = "; ".join(sorted(set(missing_details)))
        else:
            reason = "excluded_from_shared_complete_case"
            details = "; ".join(sorted(set(other_details)))
        excluded_rows.append(
            MultivariateTaxonExclusion(
                taxon=taxon,
                reason=reason,
                missing_columns=sorted(missing_columns),
                blocking_responses=blocking_responses,
                details=details,
            )
        )
    return excluded_rows


def parse_missing_columns(details: str) -> list[str]:
    prefix = "taxon is missing required value(s): "
    if not details.startswith(prefix):
        return []
    missing = details.removeprefix(prefix)
    return [column.strip() for column in missing.split(",") if column.strip()]


__all__ = [
    "build_response_formula",
    "build_shared_taxon_exclusions",
    "inspect_response_model",
    "parse_missing_columns",
    "raise_for_input_blockers",
    "shared_analysis_taxa",
]
