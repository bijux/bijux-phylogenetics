from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import json
from pathlib import Path
import re


EVIDENCE_BOOK_DIRNAME = "evidence-book"
EVIDENCE_STUDIES_DIRNAME = "studies"
EVIDENCE_INDEX_DIRNAME = "index"
EVIDENCE_STUDY_MANIFEST = "study.json"
EVIDENCE_BUNDLE_MANIFEST = "manifest.json"
EVIDENCE_CATALOG = "catalog.md"
EVIDENCE_INDEX = "evidence-index.json"
EVIDENCE_ID_PATTERN = re.compile(r"^evidence-\d{3}$")
STUDY_ID_PATTERN = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")

STUDY_REQUIRED_KEYS = {
    "study_id",
    "study_title",
    "summary",
    "owner_package",
}
EVIDENCE_REQUIRED_KEYS = {
    "schema_version",
    "study_id",
    "evidence_id",
    "evidence_title",
    "summary",
    "owner_package",
    "source_basis",
    "claim_tags",
    "verdict",
    "limitations",
}
VERDICT_REQUIRED_KEYS = {"status", "summary"}


@dataclass(frozen=True, slots=True)
class EvidenceBookValidationIssue:
    path: Path
    message: str


@dataclass(slots=True)
class EvidenceBookValidationReport:
    root: Path
    valid: bool
    issues: list[EvidenceBookValidationIssue]
    bundle_paths: list[Path]


def evidence_book_root(repo_root: Path) -> Path:
    return Path(repo_root) / EVIDENCE_BOOK_DIRNAME


def _relative_to(root: Path, path: Path) -> Path:
    try:
        return path.relative_to(root)
    except ValueError:
        return path


def _load_json(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("expected a JSON object")
    return payload


def _study_paths(book_root: Path) -> list[Path]:
    studies_root = book_root / EVIDENCE_STUDIES_DIRNAME
    if not studies_root.exists():
        return []
    return sorted(path for path in studies_root.iterdir() if path.is_dir())


def _bundle_paths(book_root: Path) -> list[Path]:
    bundle_paths: list[Path] = []
    for study_root in _study_paths(book_root):
        bundle_paths.extend(
            sorted(
                path
                for path in study_root.iterdir()
                if path.is_dir() and EVIDENCE_ID_PATTERN.fullmatch(path.name)
            )
        )
    return bundle_paths


def _validate_study_manifest(
    book_root: Path, study_root: Path, issues: list[EvidenceBookValidationIssue]
) -> dict[str, object] | None:
    if not STUDY_ID_PATTERN.fullmatch(study_root.name):
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, study_root),
                "study directory name must be a durable kebab-case identifier",
            )
        )
    readme_path = study_root / "README.md"
    if not readme_path.exists():
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, study_root),
                "study directory must include README.md",
            )
        )
    manifest_path = study_root / EVIDENCE_STUDY_MANIFEST
    if not manifest_path.exists():
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, study_root),
                f"study directory must include {EVIDENCE_STUDY_MANIFEST}",
            )
        )
        return None

    try:
        manifest = _load_json(manifest_path)
    except (ValueError, json.JSONDecodeError) as error:
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, manifest_path),
                f"invalid study manifest: {error}",
            )
        )
        return None

    missing_keys = sorted(STUDY_REQUIRED_KEYS - manifest.keys())
    if missing_keys:
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, manifest_path),
                f"study manifest missing keys: {', '.join(missing_keys)}",
            )
        )
    if manifest.get("study_id") != study_root.name:
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, manifest_path),
                "study manifest study_id must match directory name",
            )
        )
    return manifest


def _validate_bundle_manifest(
    book_root: Path,
    study_root: Path,
    bundle_root: Path,
    issues: list[EvidenceBookValidationIssue],
) -> dict[str, object] | None:
    manifest_path = bundle_root / EVIDENCE_BUNDLE_MANIFEST
    readme_path = bundle_root / "README.md"
    if not readme_path.exists():
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, bundle_root),
                "evidence bundle must include README.md",
            )
        )
    if not manifest_path.exists():
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, bundle_root),
                f"evidence bundle must include {EVIDENCE_BUNDLE_MANIFEST}",
            )
        )
        return None

    try:
        manifest = _load_json(manifest_path)
    except (ValueError, json.JSONDecodeError) as error:
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, manifest_path),
                f"invalid evidence manifest: {error}",
            )
        )
        return None

    missing_keys = sorted(EVIDENCE_REQUIRED_KEYS - manifest.keys())
    if missing_keys:
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, manifest_path),
                f"evidence manifest missing keys: {', '.join(missing_keys)}",
            )
        )
    if manifest.get("study_id") != study_root.name:
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, manifest_path),
                "evidence manifest study_id must match parent study directory",
            )
        )
    if manifest.get("evidence_id") != bundle_root.name:
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, manifest_path),
                "evidence manifest evidence_id must match bundle directory name",
            )
        )

    verdict = manifest.get("verdict")
    if not isinstance(verdict, dict):
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, manifest_path),
                "evidence manifest verdict must be a JSON object",
            )
        )
    else:
        missing_verdict_keys = sorted(VERDICT_REQUIRED_KEYS - verdict.keys())
        if missing_verdict_keys:
            issues.append(
                EvidenceBookValidationIssue(
                    _relative_to(book_root, manifest_path),
                    "evidence manifest verdict missing keys: "
                    + ", ".join(missing_verdict_keys),
                )
            )

    source_basis = manifest.get("source_basis")
    if not isinstance(source_basis, list) or not source_basis:
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, manifest_path),
                "evidence manifest source_basis must be a non-empty list",
            )
        )
    claim_tags = manifest.get("claim_tags")
    if not isinstance(claim_tags, list):
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, manifest_path),
                "evidence manifest claim_tags must be a list",
            )
        )
    limitations = manifest.get("limitations")
    if not isinstance(limitations, list):
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(book_root, manifest_path),
                "evidence manifest limitations must be a list",
            )
        )
    return manifest


def validate_evidence_book(repo_root: Path) -> EvidenceBookValidationReport:
    root = evidence_book_root(repo_root)
    issues: list[EvidenceBookValidationIssue] = []
    if not root.exists():
        issues.append(
            EvidenceBookValidationIssue(Path(EVIDENCE_BOOK_DIRNAME), "directory not found")
        )
        return EvidenceBookValidationReport(
            root=root, valid=False, issues=issues, bundle_paths=[]
        )

    studies_root = root / EVIDENCE_STUDIES_DIRNAME
    index_root = root / EVIDENCE_INDEX_DIRNAME
    if not studies_root.exists():
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(root, studies_root),
                f"missing {EVIDENCE_STUDIES_DIRNAME}/ directory",
            )
        )
    if not index_root.exists():
        issues.append(
            EvidenceBookValidationIssue(
                _relative_to(root, index_root),
                f"missing {EVIDENCE_INDEX_DIRNAME}/ directory",
            )
        )

    bundle_paths: list[Path] = []
    for study_root in _study_paths(root):
        _validate_study_manifest(root, study_root, issues)
        for child in sorted(path for path in study_root.iterdir() if path.is_dir()):
            if child.name in {"reference", "data", "figures", "provenance"}:
                continue
            if EVIDENCE_ID_PATTERN.fullmatch(child.name):
                bundle_paths.append(child)
                _validate_bundle_manifest(root, study_root, child, issues)
            elif child.name != "__pycache__":
                issues.append(
                    EvidenceBookValidationIssue(
                        _relative_to(root, child),
                        "unexpected directory inside study; evidence bundles must use evidence-### names",
                    )
                )

    return EvidenceBookValidationReport(
        root=root,
        valid=not issues,
        issues=issues,
        bundle_paths=sorted(bundle_paths),
    )


def build_evidence_book_index(repo_root: Path) -> dict[str, object]:
    root = evidence_book_root(repo_root)
    report = validate_evidence_book(repo_root)
    if not report.valid:
        messages = "; ".join(
            f"{issue.path.as_posix()}: {issue.message}" for issue in report.issues
        )
        raise ValueError(f"evidence-book is invalid: {messages}")

    evidence_entries: list[dict[str, object]] = []
    studies: list[dict[str, object]] = []
    for study_root in _study_paths(root):
        study_manifest = _load_json(study_root / EVIDENCE_STUDY_MANIFEST)
        study_entries: list[dict[str, object]] = []
        for bundle_root in sorted(
            path
            for path in study_root.iterdir()
            if path.is_dir() and EVIDENCE_ID_PATTERN.fullmatch(path.name)
        ):
            manifest = _load_json(bundle_root / EVIDENCE_BUNDLE_MANIFEST)
            verdict = manifest["verdict"]
            entry = {
                "study_id": study_manifest["study_id"],
                "study_title": study_manifest["study_title"],
                "evidence_id": manifest["evidence_id"],
                "evidence_title": manifest["evidence_title"],
                "owner_package": manifest["owner_package"],
                "relative_path": _relative_to(root, bundle_root).as_posix(),
                "claim_tags": manifest["claim_tags"],
                "verdict_status": verdict["status"],
                "verdict_summary": verdict["summary"],
            }
            study_entries.append(entry)
            evidence_entries.append(entry)
        studies.append(
            {
                "study_id": study_manifest["study_id"],
                "study_title": study_manifest["study_title"],
                "summary": study_manifest["summary"],
                "owner_package": study_manifest["owner_package"],
                "bundle_count": len(study_entries),
                "evidence": study_entries,
            }
        )

    verdict_counts = Counter(
        str(entry["verdict_status"]) for entry in evidence_entries
    )
    return {
        "schema_version": 1,
        "root": EVIDENCE_BOOK_DIRNAME,
        "study_count": len(studies),
        "evidence_count": len(evidence_entries),
        "verdict_counts": dict(sorted(verdict_counts.items())),
        "studies": studies,
        "evidence": evidence_entries,
    }


def render_evidence_catalog(index_payload: dict[str, object]) -> str:
    lines = [
        "# Evidence Catalog",
        "",
        "This catalog is generated from `index/evidence-index.json` and lists every",
        "governed evidence bundle in the repository evidence-book.",
        "",
        f"- studies: `{index_payload['study_count']}`",
        f"- evidence bundles: `{index_payload['evidence_count']}`",
        "",
    ]
    verdict_counts = index_payload.get("verdict_counts", {})
    if isinstance(verdict_counts, dict) and verdict_counts:
        lines.append("## Verdict Counts")
        lines.append("")
        for verdict, count in verdict_counts.items():
            lines.append(f"- `{verdict}`: `{count}`")
        lines.append("")

    lines.append("## Studies")
    lines.append("")
    for study in index_payload["studies"]:
        lines.append(f"### {study['study_title']}")
        lines.append("")
        lines.append(f"- study id: `{study['study_id']}`")
        lines.append(f"- owner package: `{study['owner_package']}`")
        lines.append(f"- bundle count: `{study['bundle_count']}`")
        lines.append(f"- summary: {study['summary']}")
        lines.append("")
        for entry in study["evidence"]:
            lines.append(
                f"- `{entry['evidence_id']}` — {entry['evidence_title']} "
                f"(`{entry['verdict_status']}`)"
            )
            lines.append(f"  Path: `{entry['relative_path']}`")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def write_evidence_book_index(repo_root: Path) -> tuple[Path, Path]:
    root = evidence_book_root(repo_root)
    index_root = root / EVIDENCE_INDEX_DIRNAME
    index_root.mkdir(parents=True, exist_ok=True)
    payload = build_evidence_book_index(repo_root)
    index_path = index_root / EVIDENCE_INDEX
    catalog_path = index_root / EVIDENCE_CATALOG
    index_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    catalog_path.write_text(
        render_evidence_catalog(payload),
        encoding="utf-8",
    )
    return index_path, catalog_path
