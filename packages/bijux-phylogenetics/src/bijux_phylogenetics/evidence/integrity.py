from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path

INTEGRITY_REPORT_JSON = "integrity-report.json"
INTEGRITY_REPORT_MARKDOWN = "integrity-report.md"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _bundle_digest(bundle_root: Path) -> str:
    digest = hashlib.sha256()
    for file_path in sorted(path for path in bundle_root.rglob("*") if path.is_file()):
        relative_path = file_path.relative_to(bundle_root).as_posix()
        digest.update(relative_path.encode("utf-8"))
        digest.update(b"\0")
        digest.update(_sha256(file_path).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def build_evidence_integrity_report(repo_root: Path) -> dict[str, object]:
    repo_root = Path(repo_root)
    book_root = repo_root / "evidence-book"
    entries = []
    counts: Counter[str] = Counter()
    for bundle_root in sorted(book_root.glob("studies/*/evidence-*")):
        if not bundle_root.is_dir():
            continue
        digest = _bundle_digest(bundle_root)
        file_count = sum(1 for path in bundle_root.rglob("*") if path.is_file())
        counts.update(["tracked"])
        entries.append(
            {
                "study_id": bundle_root.parent.name,
                "evidence_id": bundle_root.name,
                "relative_path": bundle_root.relative_to(book_root).as_posix(),
                "file_count": file_count,
                "bundle_digest": digest,
                "manifest_digest": _sha256(bundle_root / "manifest.json"),
                "integrity_status": "tracked",
            }
        )
    return {
        "schema_version": 1,
        "bundle_count": len(entries),
        "integrity_status_counts": dict(sorted(counts.items())),
        "entries": entries,
        "action_required_count": 0,
    }


def render_evidence_integrity_report(payload: dict[str, object]) -> str:
    lines = [
        "# Evidence Integrity Report",
        "",
        f"- bundles: `{payload['bundle_count']}`",
        f"- action required: `{payload['action_required_count']}`",
        "",
        "## Bundle Digests",
        "",
    ]
    for entry in payload["entries"]:
        lines.append(
            f"- `{entry['study_id']}/{entry['evidence_id']}` — `{entry['bundle_digest']}`"
        )
        lines.append(f"  Files: `{entry['file_count']}`")
    lines.append("")
    return "\n".join(lines)


def encode_evidence_integrity_report(payload: dict[str, object]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"
