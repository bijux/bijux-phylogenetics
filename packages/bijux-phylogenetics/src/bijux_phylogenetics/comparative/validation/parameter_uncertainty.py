from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.comparative.continuous import (
    fit_brownian_motion_model,
    fit_ornstein_uhlenbeck_model,
)
from bijux_phylogenetics.comparative.validation.reference_examples import (
    validate_comparative_reference_examples,
)


@dataclass(slots=True)
class ComparativeParameterIntervalAuditRow:
    """Audit row showing whether an external reference estimate is covered by one interval."""

    model: str
    parameter: str
    estimate: float
    lower_95: float
    upper_95: float
    reference_estimate: float
    interval_method: str
    contains_reference: bool
    reaches_search_boundary: bool
    boundary_note: str | None


@dataclass(slots=True)
class ComparativeParameterUncertaintyAudit:
    """Audit of BM/OU interval behavior against trusted reference estimates."""

    rows: list[ComparativeParameterIntervalAuditRow]
    warnings: list[str]
    all_reference_estimates_covered: bool


def audit_comparative_parameter_uncertainty() -> ComparativeParameterUncertaintyAudit:
    """Audit BM/OU parameter intervals against external reference estimates."""
    reference = {
        observation.case: observation
        for observation in validate_comparative_reference_examples().observations
    }
    root = Path(__file__).resolve().parents[4] / "tests/fixtures"
    tree = root / "trees/example_tree.nwk"
    traits = root / "metadata/example_traits_comparative.tsv"
    brownian = fit_brownian_motion_model(tree, traits, trait="response")
    ou = fit_ornstein_uhlenbeck_model(tree, traits, trait="response")

    rows: list[ComparativeParameterIntervalAuditRow] = []
    for interval in brownian.confidence_intervals:
        reference_estimate = reference["brownian-example-tree"].expected_parameters[
            interval.name
        ]
        rows.append(
            ComparativeParameterIntervalAuditRow(
                model="brownian",
                parameter=interval.name,
                estimate=interval.estimate,
                lower_95=interval.lower_95,
                upper_95=interval.upper_95,
                reference_estimate=reference_estimate,
                interval_method=interval.method,
                contains_reference=interval.lower_95
                <= reference_estimate
                <= interval.upper_95,
                reaches_search_boundary=False,
                boundary_note=None,
            )
        )
    for interval in ou.confidence_intervals:
        if interval.name not in {"alpha", "theta"}:
            continue
        reference_estimate = reference["ou-example-tree-grid"].expected_parameters[
            interval.name
        ]
        reaches_boundary = interval.name == "alpha" and math.isclose(
            interval.upper_95, interval.estimate, abs_tol=1e-9
        )
        boundary_note = (
            "supported alpha interval reaches the upper grid boundary and should be interpreted with the boundary warning"
            if reaches_boundary
            else None
        )
        rows.append(
            ComparativeParameterIntervalAuditRow(
                model="ou",
                parameter=interval.name,
                estimate=interval.estimate,
                lower_95=interval.lower_95,
                upper_95=interval.upper_95,
                reference_estimate=reference_estimate,
                interval_method=interval.method,
                contains_reference=interval.lower_95
                <= reference_estimate
                <= interval.upper_95,
                reaches_search_boundary=reaches_boundary,
                boundary_note=boundary_note,
            )
        )
    warnings = [row.boundary_note for row in rows if row.boundary_note]
    return ComparativeParameterUncertaintyAudit(
        rows=rows,
        warnings=warnings,
        all_reference_estimates_covered=all(row.contains_reference for row in rows),
    )
