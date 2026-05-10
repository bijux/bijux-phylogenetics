from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
import subprocess  # nosec B404


FRESHNESS_REPORT_JSON = "freshness-report.json"
FRESHNESS_REPORT_MARKDOWN = "freshness-report.md"


def _load_json(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object at {path}")
    return payload


def _git_last_change_date(repo_root: Path, path: str) -> str | None:
    candidate = repo_root / path
    if not candidate.exists():
        return None
    result = subprocess.run(  # nosec B603
        ["git", "-C", str(repo_root), "log", "-1", "--format=%cs", "--", path],
        check=False,
        capture_output=True,
        text=True,
    )
    value = result.stdout.strip()
    if result.returncode != 0 or not value:
        return None
    return value


def _tracked_change_dates(repo_root: Path, manifest: dict[str, object]) -> list[str]:
    freshness = manifest.get("freshness", {})
    if not isinstance(freshness, dict):
        return []
    locators = [
        str(path)
        for path in freshness.get("governed_code_paths", [])
        if isinstance(path, str)
    ]
    locators.extend(
        str(locator)
        for locator in freshness.get("source_basis_locators", [])
        if isinstance(locator, str) and ":" not in locator
    )
    change_dates = []
    for locator in locators:
        change_date = _git_last_change_date(repo_root, locator)
        if change_date is not None:
            change_dates.append(change_date)
    return sorted(set(change_dates))


def build_evidence_freshness_report(repo_root: Path) -> dict[str, object]:
    repo_root = Path(repo_root)
    book_root = repo_root / "evidence-book"
    entries = []
    counts: Counter[str] = Counter()
    for manifest_path in sorted(book_root.glob("studies/*/evidence-*/manifest.json")):
        manifest = _load_json(manifest_path)
        change_dates = _tracked_change_dates(repo_root, manifest)
        last_generated_on = str(manifest["freshness"]["last_generated_on"])
        latest_change_on = change_dates[-1] if change_dates else None
        if latest_change_on is None:
            status = "source_unresolved"
        elif last_generated_on >= latest_change_on:
            status = "current"
        else:
            status = "stale"
        counts.update([status])
        entries.append(
            {
                "study_id": manifest["study_id"],
                "evidence_id": manifest["evidence_id"],
                "relative_path": manifest_path.relative_to(book_root).parent.as_posix(),
                "last_generated_on": last_generated_on,
                "latest_governed_change_on": latest_change_on,
                "freshness_status": status,
                "tracked_input_count": len(change_dates),
            }
        )
    return {
        "schema_version": 1,
        "bundle_count": len(entries),
        "freshness_status_counts": dict(sorted(counts.items())),
        "entries": entries,
    }


def render_evidence_freshness_report(payload: dict[str, object]) -> str:
    lines = [
        "# Evidence Freshness Report",
        "",
        f"- bundles: `{payload['bundle_count']}`",
        "",
        "## Status Counts",
        "",
    ]
    for status, count in payload["freshness_status_counts"].items():
        lines.append(f"- `{status}`: `{count}`")
    lines.extend(
        [
            "",
            "## Bundles",
            "",
        ]
    )
    for entry in payload["entries"]:
        lines.append(
            f"- `{entry['study_id']}/{entry['evidence_id']}` — `{entry['freshness_status']}`"
        )
        lines.append(
            f"  Generated: `{entry['last_generated_on']}`; latest governed change: `{entry['latest_governed_change_on']}`"
        )
    lines.append("")
    return "\n".join(lines)


def encode_evidence_freshness_report(payload: dict[str, object]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"
