from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re


PORTABILITY_ALLOWED_REPO_PREFIXES = (
    ".github/",
    "apis/",
    "artifacts/",
    "block-payloads/",
    "configs/",
    "docs/",
    "evidence-book/",
    "index/",
    "makes/",
    "mkdocs.yml",
    "packages/",
    "README.md",
    "reference/",
    "studies/",
    "provenance/",
)
PORTABILITY_ALLOWED_EXTERNAL_PREFIXES = ("external:",)
PORTABILITY_TRACKED_KEYS = {
    "build_script_path",
    "checked_fixture_path",
    "code_path",
    "governed_code_paths",
    "governed_reload_inputs",
    "locator",
    "path",
    "python_reference_script",
    "r_reference_script",
    "relative_path",
    "source_basis_locators",
}
LINE_ANCHOR_PATTERN = re.compile(r"#L\d+(?:-L\d+)?$")


@dataclass(frozen=True, slots=True)
class EvidencePathIssue:
    relative_file_path: str
    json_pointer: str
    value: str
    issue_kind: str
    message: str


@dataclass(frozen=True, slots=True)
class EvidencePathValue:
    json_pointer: str
    key: str
    value: str
    locator_kind: str


def classify_locator_kind(value: str) -> str:
    if value.startswith("/Users/"):
        return "workstation_absolute"
    if Path(value).is_absolute():
        return "absolute_path"
    if value.startswith("../") or "/../" in value or value.startswith("..\\") or "\\..\\" in value:
        return "parent_traversal"
    if any(value.startswith(prefix) for prefix in PORTABILITY_ALLOWED_EXTERNAL_PREFIXES):
        return "external_locator"
    normalized = _strip_line_anchor(value)
    if any(normalized == prefix.rstrip("/") or normalized.startswith(prefix) for prefix in PORTABILITY_ALLOWED_REPO_PREFIXES):
        return "repo_relative"
    if "/" in normalized or normalized.endswith((".json", ".md", ".py", ".R", ".csv", ".nwk", ".tsv")):
        return "suspicious_path_like"
    return "non_path_text"


def is_portable_locator(value: str) -> bool:
    return classify_locator_kind(value) in {"external_locator", "repo_relative"}


def collect_payload_path_values(payload: object) -> list[EvidencePathValue]:
    collected: list[EvidencePathValue] = []
    _collect_payload_path_values(payload, json_pointer="$", parent_key=None, collected=collected)
    return collected


def audit_payload_path_values(payload: object, *, relative_file_path: str) -> list[EvidencePathIssue]:
    issues: list[EvidencePathIssue] = []
    for entry in collect_payload_path_values(payload):
        if entry.locator_kind == "non_path_text":
            continue
        if entry.locator_kind == "workstation_absolute":
            issues.append(
                EvidencePathIssue(
                    relative_file_path=relative_file_path,
                    json_pointer=entry.json_pointer,
                    value=entry.value,
                    issue_kind="workstation_absolute",
                    message="durable evidence path values must not embed workstation-local absolute paths",
                )
            )
        elif entry.locator_kind == "absolute_path":
            issues.append(
                EvidencePathIssue(
                    relative_file_path=relative_file_path,
                    json_pointer=entry.json_pointer,
                    value=entry.value,
                    issue_kind="absolute_path",
                    message="durable evidence path values must be repository-relative or explicitly externalized",
                )
            )
        elif entry.locator_kind == "parent_traversal":
            issues.append(
                EvidencePathIssue(
                    relative_file_path=relative_file_path,
                    json_pointer=entry.json_pointer,
                    value=entry.value,
                    issue_kind="parent_traversal",
                    message="durable evidence path values must not rely on parent-directory traversal or sibling-repository guesses",
                )
            )
        elif entry.locator_kind == "suspicious_path_like":
            issues.append(
                EvidencePathIssue(
                    relative_file_path=relative_file_path,
                    json_pointer=entry.json_pointer,
                    value=entry.value,
                    issue_kind="suspicious_path_like",
                    message="durable evidence path values must use repository-relative prefixes or explicit external locators",
                )
            )
    return issues


def load_json_payload(path: Path) -> dict[str, object] | list[object]:
    return json.loads(path.read_text(encoding="utf-8"))


def render_portability_rules_markdown() -> str:
    lines = [
        "# Portability Rules",
        "",
        "Checked-in evidence paths must stay portable and reviewer-readable.",
        "",
        "Allowed locator forms:",
        "",
        "- repository-relative paths rooted at governed prefixes such as `evidence-book/`, `packages/`, `docs/`, or `artifacts/`",
        "- explicit external locators such as `external:lund/...`",
        "",
        "Forbidden locator forms:",
        "",
        "- workstation-local absolute paths tied to one machine layout",
        "- parent traversal such as `../...` that guesses sibling repositories or ambient workspace shape",
        "- ambiguous path-like strings that are neither explicit external locators nor rooted repository-relative paths",
        "",
        "Portable checked-in reports and plots must reference their governed sources through those same locator forms.",
        "",
    ]
    return "\n".join(lines)


def _collect_payload_path_values(
    payload: object,
    *,
    json_pointer: str,
    parent_key: str | None,
    collected: list[EvidencePathValue],
) -> None:
    if isinstance(payload, dict):
        for key, value in payload.items():
            child_pointer = f"{json_pointer}/{key}"
            if isinstance(value, str) and _should_track_string(key, value):
                collected.append(
                    EvidencePathValue(
                        json_pointer=child_pointer,
                        key=key,
                        value=value,
                        locator_kind=classify_locator_kind(value),
                    )
                )
            else:
                _collect_payload_path_values(
                    value,
                    json_pointer=child_pointer,
                    parent_key=key,
                    collected=collected,
                )
        return
    if isinstance(payload, list):
        for index, value in enumerate(payload):
            child_pointer = f"{json_pointer}/{index}"
            if isinstance(value, str) and _should_track_string(parent_key, value):
                collected.append(
                    EvidencePathValue(
                        json_pointer=child_pointer,
                        key=parent_key or "list-item",
                        value=value,
                        locator_kind=classify_locator_kind(value),
                    )
                )
            else:
                _collect_payload_path_values(
                    value,
                    json_pointer=child_pointer,
                    parent_key=parent_key,
                    collected=collected,
                )


def _should_track_string(key: str | None, value: str) -> bool:
    if key in PORTABILITY_TRACKED_KEYS:
        return True
    if value.startswith(PORTABILITY_ALLOWED_EXTERNAL_PREFIXES):
        return True
    if value.startswith(PORTABILITY_ALLOWED_REPO_PREFIXES):
        return True
    if value.startswith("/Users/") or Path(value).is_absolute():
        return True
    if value.startswith("../") or "/../" in value or value.startswith("..\\") or "\\..\\" in value:
        return True
    return False


def _strip_line_anchor(value: str) -> str:
    match = LINE_ANCHOR_PATTERN.search(value)
    if match is None:
        return value
    return value[: match.start()]
