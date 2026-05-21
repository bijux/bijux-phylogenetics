from __future__ import annotations

from pathlib import Path
import tempfile

from .models import ReferenceFixtureCheck, ReferenceValidationSuiteReport


def default_fixtures_root() -> Path:
    return Path(__file__).resolve().parents[4] / "tests" / "fixtures"


def fixture(root: Path, *parts: str) -> Path:
    return root.joinpath(*parts)


def normalize(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): normalize(item) for key, item in value.items()}
    if isinstance(value, list):
        return [normalize(item) for item in value]
    return value


def check(
    *,
    goal_id: int,
    suite: str,
    name: str,
    fixture_paths: list[Path],
    expected: dict[str, object],
    observed: dict[str, object],
    notes: list[str] | None = None,
) -> ReferenceFixtureCheck:
    normalized_expected = normalize(expected)
    normalized_observed = normalize(observed)
    return ReferenceFixtureCheck(
        goal_id=goal_id,
        suite=suite,
        name=name,
        fixture_paths=fixture_paths,
        passed=normalized_expected == normalized_observed,
        expected=normalized_expected,
        observed=normalized_observed,
        notes=[] if notes is None else list(notes),
    )


def suite_report(
    *,
    goal_id: int,
    suite: str,
    reviewer_goal: str,
    fixtures: list[ReferenceFixtureCheck],
    coverage_notes: list[str],
    limitations: list[str],
) -> ReferenceValidationSuiteReport:
    passed_fixture_count = sum(
        1 for checked_fixture in fixtures if checked_fixture.passed
    )
    failed_fixture_count = len(fixtures) - passed_fixture_count
    return ReferenceValidationSuiteReport(
        goal_id=goal_id,
        suite=suite,
        reviewer_goal=reviewer_goal,
        passed=failed_fixture_count == 0,
        fixture_count=len(fixtures),
        passed_fixture_count=passed_fixture_count,
        failed_fixture_count=failed_fixture_count,
        fixtures=fixtures,
        coverage_notes=coverage_notes,
        limitations=limitations,
    )


def error_observation(error: Exception) -> dict[str, object]:
    return {
        "error_type": type(error).__name__,
        "error_code": getattr(error, "code", "unknown"),
        "message": str(error),
    }


def temp_reference_dir(name: str) -> Path:
    """Return a deterministic temporary directory for reference checks."""
    return Path(tempfile.gettempdir()) / name
