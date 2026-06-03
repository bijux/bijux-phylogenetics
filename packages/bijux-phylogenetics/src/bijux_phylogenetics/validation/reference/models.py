from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class ReferenceFixtureCheck:
    """Observed-versus-expected validation for one checked-in reference fixture."""

    goal_id: int
    suite: str
    name: str
    fixture_paths: list[Path]
    passed: bool
    expected: dict[str, object]
    observed: dict[str, object]
    notes: list[str]


@dataclass(slots=True)
class ReferenceValidationSuiteReport:
    """One reviewer-facing suite of reference validations."""

    goal_id: int
    suite: str
    reviewer_goal: str
    passed: bool
    fixture_count: int
    passed_fixture_count: int
    failed_fixture_count: int
    fixtures: list[ReferenceFixtureCheck]
    coverage_notes: list[str]
    limitations: list[str]


@dataclass(slots=True)
class CoreWorkflowValidationRow:
    """One Level 1 workflow with fixture coverage and trust notes."""

    workflow: str
    fixture_suite_names: list[str]
    fixture_count: int
    expected_outputs: list[str]
    limitations: list[str]
    passed: bool
    notes: list[str]


@dataclass(frozen=True, slots=True)
class CoreWorkflowFailureCase:
    """Known workflow failure or warning case with expected behavior."""

    workflow: str
    fixture_name: str
    outcome_kind: str
    observed_code: str
    observed_summary: str
    passed: bool


@dataclass(frozen=True, slots=True)
class WorkflowMaturityClassification:
    """Reviewer-facing maturity label for one core workflow."""

    workflow: str
    maturity: str
    rationale: list[str]
    outstanding_risks: list[str]


@dataclass(slots=True)
class CoreWorkflowValidationReport:
    """Aggregate validation report for the Level 1 trust surface."""

    suites: list[ReferenceValidationSuiteReport]
    workflows: list[CoreWorkflowValidationRow]
    failure_gallery: list[CoreWorkflowFailureCase]
    maturity_classifications: list[WorkflowMaturityClassification]
    total_fixture_count: int
    passed_fixture_count: int
    failed_fixture_count: int
    limitations: list[str]


@dataclass(slots=True)
class LevelOneReleaseGateDecision:
    """Gate decision for whether the example Level 1 workflow is review-ready."""

    decision: str
    rationale: list[str]
    retained_taxa: list[str]
    excluded_taxa: list[str]
    blocked_analyses: list[str]
    allowed_analyses: list[str]
    reviewer_visible_warnings: list[str]


@dataclass(slots=True)
class LevelOneReleaseGateReport:
    """Integrated release gate built around the checked-in workflow fixtures."""

    fixtures_root: Path
    validation: CoreWorkflowValidationReport
    dataset_readiness_decision: str
    dataset_blockers: list[str]
    dataset_warnings: list[str]
    exclusion_causes: dict[str, list[str]]
    taxon_first_loss_stage: dict[str, str | None]
    gate: LevelOneReleaseGateDecision
