from __future__ import annotations

import json
from pathlib import Path
import re

EVIDENCE_BOOK_DIRNAME = "evidence-book"
EVIDENCE_STUDIES_DIRNAME = "studies"
EVIDENCE_INDEX_DIRNAME = "index"
EVIDENCE_BUNDLE_MANIFEST = "manifest.json"
EVIDENCE_CATALOG = "catalog.md"
EVIDENCE_INDEX = "evidence-index.json"
EVIDENCE_CLAIM_MAP = "claim-map.json"
EVIDENCE_PARITY_DASHBOARD = "parity-dashboard.json"
EVIDENCE_PARITY_SUMMARY = "parity-dashboard.md"
EVIDENCE_MISMATCH_ARCHIVE = "mismatch-archive.json"
EVIDENCE_MISMATCH_SUMMARY = "mismatch-archive.md"
EVIDENCE_VERDICT_WORKFLOWS = "verdict-workflows.json"
EVIDENCE_VERDICT_WORKFLOWS_SUMMARY = "verdict-workflows.md"
EVIDENCE_FALSE_CONFIDENCE_AUDIT = "false-confidence-audit.json"
EVIDENCE_FALSE_CONFIDENCE_SUMMARY = "false-confidence-audit.md"
EVIDENCE_SCIENTIFIC_DEBT_REGISTER = "scientific-debt-register.json"
EVIDENCE_SCIENTIFIC_DEBT_SUMMARY = "scientific-debt-register.md"
EVIDENCE_PORTABILITY_AUDIT = "portability-audit.json"
EVIDENCE_PORTABILITY_SUMMARY = "portability-audit.md"
EVIDENCE_FRAGILE_EXAMPLE_AUDIT = "fragile-example-audit.json"
EVIDENCE_FRAGILE_EXAMPLE_SUMMARY = "fragile-example-audit.md"
EVIDENCE_REGENERATION_CONTRACT = "regeneration-contract.json"
EVIDENCE_REGENERATION_SUMMARY = "regeneration-contract.md"
EVIDENCE_CLAIM_REAUDIT = "claim-reaudit.json"
EVIDENCE_CLAIM_REAUDIT_SUMMARY = "claim-reaudit.md"
EVIDENCE_ANALYTICAL_SURFACE_COVERAGE = "analytical-surface-coverage.json"
EVIDENCE_ANALYTICAL_SURFACE_COVERAGE_SUMMARY = "analytical-surface-coverage.md"
EVIDENCE_CLOSURE_CRITERIA = "closure-criteria.json"
EVIDENCE_CLOSURE_CRITERIA_SUMMARY = "closure-criteria.md"
EVIDENCE_MATURITY_SCORECARD = "evidence-maturity-scorecard.json"
EVIDENCE_MATURITY_SCORECARD_SUMMARY = "evidence-maturity-scorecard.md"
EVIDENCE_REVIEW_RITUAL = "evidence-review-ritual.json"
EVIDENCE_REVIEW_RITUAL_SUMMARY = "evidence-review-ritual.md"
EVIDENCE_COMPLETION_GATES = "completion-gates.json"
EVIDENCE_COMPLETION_GATES_SUMMARY = "completion-gates.md"
EVIDENCE_ID_PATTERN = re.compile(r"^evidence-\d{3}$")
STUDY_ID_PATTERN = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")


def evidence_book_root(repo_root: Path) -> Path:
    return Path(repo_root) / EVIDENCE_BOOK_DIRNAME


def relative_to(root: Path, path: Path) -> Path:
    try:
        return path.relative_to(root)
    except ValueError:
        return path


def load_json(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("expected a JSON object")
    return payload


def load_json_list(path: Path) -> list[object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, list):
        raise ValueError("expected a JSON array")
    return payload


def load_bundle_claim_rows(bundle_root: Path) -> list[dict[str, object]]:
    claims_path = bundle_root / "claims.json"
    legacy_claim_verdicts_path = bundle_root / "claim_verdicts.json"
    if claims_path.exists():
        payload = load_json(claims_path)
        claims = payload.get("claims")
        if not isinstance(claims, list):
            raise ValueError("claims.json must contain a claims list")
        return [row for row in claims if isinstance(row, dict)]
    if legacy_claim_verdicts_path.exists():
        rows = load_json_list(legacy_claim_verdicts_path)
        return [row for row in rows if isinstance(row, dict)]
    return []


def study_paths(book_root: Path) -> list[Path]:
    studies_root = book_root / EVIDENCE_STUDIES_DIRNAME
    if not studies_root.exists():
        return []
    return sorted(path for path in studies_root.iterdir() if path.is_dir())


def bundle_paths(book_root: Path) -> list[Path]:
    discovered_bundle_paths: list[Path] = []
    for study_root in study_paths(book_root):
        discovered_bundle_paths.extend(
            sorted(
                path
                for path in study_root.iterdir()
                if path.is_dir() and EVIDENCE_ID_PATTERN.fullmatch(path.name)
            )
        )
    return discovered_bundle_paths
