from __future__ import annotations

import json
from pathlib import Path
import re

from ..bundle_contracts import (
    REQUIRED_BUNDLE_LOCAL_ARTIFACTS,
    REQUIRED_BUNDLE_RESULT_ARTIFACTS,
)
from ..coverage import COVERAGE_GAPS_JSON, COVERAGE_GAPS_MARKDOWN
from ..freshness import FRESHNESS_REPORT_JSON, FRESHNESS_REPORT_MARKDOWN
from ..integrity import INTEGRITY_REPORT_JSON, INTEGRITY_REPORT_MARKDOWN
from ..portability import audit_payload_path_values
from ..reviewer import REVIEWER_SUMMARY_JSON, REVIEWER_SUMMARY_MARKDOWN
from ..study_contracts import (
    ALLOWED_STUDY_ROOT_DIRS,
    ALLOWED_STUDY_ROOT_FILES,
    load_study_contract,
)
from ..teaching import (
    ALLOWED_COMPARISON_MODES,
    ALLOWED_STUDY_CATEGORIES,
    TEACHING_AND_MIGRATION_INDEX_FILENAME,
    TEACHING_AND_MIGRATION_SUMMARY_FILENAME,
    study_metadata,
    teaching_study_ids,
)
from .layout import (
    EVIDENCE_ANALYTICAL_SURFACE_COVERAGE,
    EVIDENCE_ANALYTICAL_SURFACE_COVERAGE_SUMMARY,
    EVIDENCE_BOOK_DIRNAME,
    EVIDENCE_BUNDLE_MANIFEST,
    EVIDENCE_CATALOG,
    EVIDENCE_CLAIM_MAP,
    EVIDENCE_CLAIM_REAUDIT,
    EVIDENCE_CLAIM_REAUDIT_SUMMARY,
    EVIDENCE_CLOSURE_CRITERIA,
    EVIDENCE_CLOSURE_CRITERIA_SUMMARY,
    EVIDENCE_COMPLETION_GATES,
    EVIDENCE_COMPLETION_GATES_SUMMARY,
    EVIDENCE_FALSE_CONFIDENCE_AUDIT,
    EVIDENCE_FALSE_CONFIDENCE_SUMMARY,
    EVIDENCE_FRAGILE_EXAMPLE_AUDIT,
    EVIDENCE_FRAGILE_EXAMPLE_SUMMARY,
    EVIDENCE_ID_PATTERN,
    EVIDENCE_INDEX,
    EVIDENCE_INDEX_DIRNAME,
    EVIDENCE_MATURITY_SCORECARD,
    EVIDENCE_MATURITY_SCORECARD_SUMMARY,
    EVIDENCE_MISMATCH_ARCHIVE,
    EVIDENCE_MISMATCH_SUMMARY,
    EVIDENCE_PARITY_DASHBOARD,
    EVIDENCE_PARITY_SUMMARY,
    EVIDENCE_PORTABILITY_AUDIT,
    EVIDENCE_PORTABILITY_SUMMARY,
    EVIDENCE_REGENERATION_CONTRACT,
    EVIDENCE_REGENERATION_SUMMARY,
    EVIDENCE_REVIEW_RITUAL,
    EVIDENCE_REVIEW_RITUAL_SUMMARY,
    EVIDENCE_SCIENTIFIC_DEBT_REGISTER,
    EVIDENCE_SCIENTIFIC_DEBT_SUMMARY,
    EVIDENCE_STUDIES_DIRNAME,
    EVIDENCE_VERDICT_WORKFLOWS,
    EVIDENCE_VERDICT_WORKFLOWS_SUMMARY,
    STUDY_ID_PATTERN,
    evidence_book_root,
    load_bundle_claim_rows,
    load_json,
    relative_to,
    study_paths,
)
from .models import EvidenceBookValidationIssue, EvidenceBookValidationReport

EVIDENCE_REQUIRED_KEYS = {
    "schema_version",
    "study_id",
    "evidence_id",
    "evidence_title",
    "summary",
    "owner_package",
    "claim_ids",
    "source_basis",
    "freshness",
    "ownership",
    "claim_tags",
    "comparison_mode",
    "verdict",
    "limitations",
}
VERDICT_REQUIRED_KEYS = {"status", "summary"}
FRESHNESS_REQUIRED_KEYS = {
    "last_generated_on",
    "governed_code_paths",
    "source_basis_locators",
}
OWNERSHIP_REQUIRED_KEYS = {"owner_package", "analytical_surfaces"}
ALLOWED_EVIDENCE_VERDICT_STATUSES = {
    "matched",
    "matched_with_tolerance",
    "mismatch_explained",
    "mismatch_unexplained",
    "not_comparable",
}


def scan_for_local_path_leaks(
    book_root: Path, issues: list[EvidenceBookValidationIssue]
) -> None:
    for path in sorted(
        candidate for candidate in book_root.rglob("*") if candidate.is_file()
    ):
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        if "/Users/" in text:
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, path),
                    "durable evidence content must not contain workstation-local /Users/ paths",
                )
            )


def scan_for_portability_issues(
    book_root: Path, issues: list[EvidenceBookValidationIssue]
) -> None:
    for path in sorted(
        candidate for candidate in book_root.rglob("*.json") if candidate.is_file()
    ):
        relative_path = relative_to(book_root, path).as_posix()
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        if not isinstance(payload, (dict, list)):
            continue
        for issue in audit_payload_path_values(
            payload, relative_file_path=relative_path
        ):
            issues.append(
                EvidenceBookValidationIssue(
                    Path(issue.relative_file_path),
                    f"{issue.json_pointer}: {issue.message}",
                )
            )


def validate_study_contract(
    book_root: Path, study_root: Path, issues: list[EvidenceBookValidationIssue]
) -> dict[str, object] | None:
    if not STUDY_ID_PATTERN.fullmatch(study_root.name):
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, study_root),
                "study directory name must be a durable kebab-case identifier",
            )
        )
    readme_path = study_root / "README.md"
    if not readme_path.exists():
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, study_root),
                "study directory must include README.md",
            )
        )
    allowed_entries = {
        *ALLOWED_STUDY_ROOT_FILES,
        *ALLOWED_STUDY_ROOT_DIRS,
        *[
            child.name
            for child in study_root.iterdir()
            if child.is_dir() and EVIDENCE_ID_PATTERN.fullmatch(child.name)
        ],
    }
    for child in sorted(study_root.iterdir()):
        if child.name == "__pycache__":
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, child),
                    "study root must not contain __pycache__",
                )
            )
            continue
        if child.name not in allowed_entries:
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, child),
                    "study root may only contain README.md, datasets/, reference/, provenance/, and evidence-### directories",
                )
            )

    provenance_root = study_root / "provenance"
    if not provenance_root.is_dir():
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, provenance_root),
                "study directory must include provenance/",
            )
        )
    dataset_registry_path = study_root / "datasets" / "registry.json"
    if not dataset_registry_path.is_file():
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, dataset_registry_path),
                "study directory must include datasets/registry.json",
            )
        )
    provenance_files = sorted(
        path for path in provenance_root.glob("*.json") if path.is_file()
    )
    if len(provenance_files) != 1:
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, provenance_root),
                "study provenance/ must contain exactly one descriptor json file",
            )
        )

    contract = load_study_contract(study_root)
    study_categories = contract.get("study_categories")
    if (
        not isinstance(study_categories, list)
        or not study_categories
        or not all(isinstance(value, str) and value for value in study_categories)
    ):
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, study_root),
                "study contract study_categories must be a non-empty list of strings",
            )
        )
    elif not set(study_categories) <= ALLOWED_STUDY_CATEGORIES:
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, study_root),
                "study contract study_categories must use governed category names",
            )
        )
    return contract


def validate_bundle_manifest(
    book_root: Path,
    study_root: Path,
    bundle_root: Path,
    issues: list[EvidenceBookValidationIssue],
    *,
    require_generated_outputs: bool,
) -> dict[str, object] | None:
    manifest_path = bundle_root / EVIDENCE_BUNDLE_MANIFEST
    readme_path = bundle_root / "README.md"
    allowed_root_entries = {
        "README.md",
        "analysis.py",
        "checks.json",
        "claims.json",
        "inputs.manifest.json",
        "manifest.json",
        "provenance.json",
        "reference.R",
        "report.md",
        "results",
    }
    if not readme_path.exists():
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, bundle_root),
                "evidence bundle must include README.md",
            )
        )
    if not manifest_path.exists():
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, bundle_root),
                f"evidence bundle must include {EVIDENCE_BUNDLE_MANIFEST}",
            )
        )
        return None

    try:
        manifest = load_json(manifest_path)
    except (ValueError, json.JSONDecodeError) as error:
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, manifest_path),
                f"invalid evidence manifest: {error}",
            )
        )
        return None

    for child in sorted(bundle_root.iterdir()):
        if child.name == "__pycache__":
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, child),
                    "evidence bundle root must not contain __pycache__",
                )
            )
            continue
        if child.name not in allowed_root_entries:
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, child),
                    "evidence bundle root must use the governed symmetric layout and move machine outputs under results/",
                )
            )

    missing_keys = sorted(EVIDENCE_REQUIRED_KEYS - manifest.keys())
    if missing_keys:
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, manifest_path),
                f"evidence manifest missing keys: {', '.join(missing_keys)}",
            )
        )
    if manifest.get("study_id") != study_root.name:
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, manifest_path),
                "evidence manifest study_id must match parent study directory",
            )
        )
    if manifest.get("evidence_id") != bundle_root.name:
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, manifest_path),
                "evidence manifest evidence_id must match bundle directory name",
            )
        )
    for filename in REQUIRED_BUNDLE_LOCAL_ARTIFACTS:
        artifact_path = bundle_root / filename
        if not artifact_path.exists():
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, artifact_path),
                    f"evidence bundle must include governed local artifact {filename}",
                )
            )
    for relative_path in REQUIRED_BUNDLE_RESULT_ARTIFACTS:
        result_path = bundle_root / relative_path
        if not result_path.exists():
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, result_path),
                    f"evidence bundle must include governed result surface {relative_path}",
                )
            )
    if require_generated_outputs:
        reviewer_summary_json_path = bundle_root / REVIEWER_SUMMARY_JSON
        reviewer_summary_markdown_path = bundle_root / REVIEWER_SUMMARY_MARKDOWN
        if not reviewer_summary_json_path.exists():
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, reviewer_summary_json_path),
                    "evidence bundle must include a governed reviewer-summary.json output",
                )
            )
        if not reviewer_summary_markdown_path.exists():
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, reviewer_summary_markdown_path),
                    "evidence bundle must include a governed reviewer-summary.md output",
                )
            )

    verdict = manifest.get("verdict")
    if not isinstance(verdict, dict):
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, manifest_path),
                "evidence manifest verdict must be a JSON object",
            )
        )
    else:
        missing_verdict_keys = sorted(VERDICT_REQUIRED_KEYS - verdict.keys())
        if missing_verdict_keys:
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, manifest_path),
                    "evidence manifest verdict missing keys: "
                    + ", ".join(missing_verdict_keys),
                )
            )
        status = verdict.get("status")
        if status not in ALLOWED_EVIDENCE_VERDICT_STATUSES:
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, manifest_path),
                    "evidence manifest verdict status must be one of: "
                    + ", ".join(sorted(ALLOWED_EVIDENCE_VERDICT_STATUSES)),
                )
            )

    source_basis = manifest.get("source_basis")
    if not isinstance(source_basis, list) or not source_basis:
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, manifest_path),
                "evidence manifest source_basis must be a non-empty list",
            )
        )
    claim_tags = manifest.get("claim_tags")
    if not isinstance(claim_tags, list):
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, manifest_path),
                "evidence manifest claim_tags must be a list",
            )
        )
    claim_ids = manifest.get("claim_ids")
    if (
        not isinstance(claim_ids, list)
        or not claim_ids
        or not all(isinstance(value, str) and value for value in claim_ids)
    ):
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, manifest_path),
                "evidence manifest claim_ids must be a non-empty list of strings",
            )
        )
    else:
        claim_rows = load_bundle_claim_rows(bundle_root)
        if not claim_rows:
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, bundle_root),
                    "evidence bundle must include claims.json or claim_verdicts.json for declared claim_ids",
                )
            )
        else:
            row_claim_ids = {
                row.get("claim_id")
                for row in claim_rows
                if isinstance(row.get("claim_id"), str)
            }
            if set(claim_ids) != row_claim_ids:
                issues.append(
                    EvidenceBookValidationIssue(
                        relative_to(book_root, bundle_root),
                        "evidence bundle claim rows must match manifest claim_ids exactly",
                    )
                )

    freshness = manifest.get("freshness")
    if not isinstance(freshness, dict):
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, manifest_path),
                "evidence manifest freshness must be a JSON object",
            )
        )
    else:
        missing_freshness_keys = sorted(FRESHNESS_REQUIRED_KEYS - freshness.keys())
        if missing_freshness_keys:
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, manifest_path),
                    "evidence manifest freshness missing keys: "
                    + ", ".join(missing_freshness_keys),
                )
            )
        last_generated_on = freshness.get("last_generated_on")
        if not isinstance(last_generated_on, str) or not re.fullmatch(
            r"\d{4}-\d{2}-\d{2}", last_generated_on
        ):
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, manifest_path),
                    "evidence manifest freshness last_generated_on must use YYYY-MM-DD",
                )
            )
        governed_code_paths = freshness.get("governed_code_paths")
        if (
            not isinstance(governed_code_paths, list)
            or not governed_code_paths
            or not all(
                isinstance(value, str) and value for value in governed_code_paths
            )
        ):
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, manifest_path),
                    "evidence manifest freshness governed_code_paths must be a non-empty list of strings",
                )
            )
        source_basis_locators = freshness.get("source_basis_locators")
        if (
            not isinstance(source_basis_locators, list)
            or not source_basis_locators
            or not all(
                isinstance(value, str) and value for value in source_basis_locators
            )
        ):
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, manifest_path),
                    "evidence manifest freshness source_basis_locators must be a non-empty list of strings",
                )
            )
        elif isinstance(source_basis, list):
            bundle_locators = {
                entry.get("locator")
                for entry in source_basis
                if isinstance(entry, dict) and isinstance(entry.get("locator"), str)
            }
            if not set(source_basis_locators) <= bundle_locators:
                issues.append(
                    EvidenceBookValidationIssue(
                        relative_to(book_root, manifest_path),
                        "evidence manifest freshness source_basis_locators must refer to source_basis locators",
                    )
                )

    ownership = manifest.get("ownership")
    if not isinstance(ownership, dict):
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, manifest_path),
                "evidence manifest ownership must be a JSON object",
            )
        )
    else:
        missing_ownership_keys = sorted(OWNERSHIP_REQUIRED_KEYS - ownership.keys())
        if missing_ownership_keys:
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, manifest_path),
                    "evidence manifest ownership missing keys: "
                    + ", ".join(missing_ownership_keys),
                )
            )
        if ownership.get("owner_package") != manifest.get("owner_package"):
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, manifest_path),
                    "evidence manifest ownership owner_package must match owner_package",
                )
            )
        analytical_surfaces = ownership.get("analytical_surfaces")
        if (
            not isinstance(analytical_surfaces, list)
            or not analytical_surfaces
            or not all(
                isinstance(value, str) and value for value in analytical_surfaces
            )
        ):
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(book_root, manifest_path),
                    "evidence manifest ownership analytical_surfaces must be a non-empty list of strings",
                )
            )
    limitations = manifest.get("limitations")
    if not isinstance(limitations, list):
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, manifest_path),
                "evidence manifest limitations must be a list",
            )
        )
    comparison_mode = manifest.get("comparison_mode")
    if comparison_mode not in ALLOWED_COMPARISON_MODES:
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(book_root, manifest_path),
                "evidence manifest comparison_mode must be one of: "
                + ", ".join(sorted(ALLOWED_COMPARISON_MODES)),
            )
        )
    return manifest


def validate_evidence_book(
    repo_root: Path,
    *,
    require_index_outputs: bool = True,
    require_generated_bundle_outputs: bool = True,
) -> EvidenceBookValidationReport:
    root = evidence_book_root(repo_root)
    issues: list[EvidenceBookValidationIssue] = []
    if not root.exists():
        issues.append(
            EvidenceBookValidationIssue(
                Path(EVIDENCE_BOOK_DIRNAME), "directory not found"
            )
        )
        return EvidenceBookValidationReport(
            root=root, valid=False, issues=issues, bundle_paths=[]
        )

    studies_root = root / EVIDENCE_STUDIES_DIRNAME
    index_root = root / EVIDENCE_INDEX_DIRNAME
    if not studies_root.exists():
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(root, studies_root),
                f"missing {EVIDENCE_STUDIES_DIRNAME}/ directory",
            )
        )
    if not index_root.exists():
        issues.append(
            EvidenceBookValidationIssue(
                relative_to(root, index_root),
                f"missing {EVIDENCE_INDEX_DIRNAME}/ directory",
            )
        )

    bundle_paths: list[Path] = []
    for study_root in study_paths(root):
        study_manifest = validate_study_contract(root, study_root, issues)
        if (
            isinstance(study_manifest, dict)
            and str(study_manifest["study_id"]) in teaching_study_ids()
        ):
            expected_categories = set(
                study_metadata(str(study_manifest["study_id"]))["study_categories"]
            )
            actual_categories = set(study_manifest.get("study_categories", []))
            if actual_categories != expected_categories:
                issues.append(
                    EvidenceBookValidationIssue(
                        relative_to(root, study_root / "README.md"),
                        "teaching study categories must match governed teaching metadata",
                    )
                )
        for child in sorted(path for path in study_root.iterdir() if path.is_dir()):
            if child.name in ALLOWED_STUDY_ROOT_DIRS:
                continue
            if EVIDENCE_ID_PATTERN.fullmatch(child.name):
                bundle_paths.append(child)
                validate_bundle_manifest(
                    root,
                    study_root,
                    child,
                    issues,
                    require_generated_outputs=require_generated_bundle_outputs,
                )
            elif child.name != "__pycache__":
                issues.append(
                    EvidenceBookValidationIssue(
                        relative_to(root, child),
                        "unexpected directory inside study; evidence bundles must use evidence-### names",
                    )
                )

    index_path = index_root / EVIDENCE_INDEX
    catalog_path = index_root / EVIDENCE_CATALOG
    if index_root.exists() and require_index_outputs:
        if not index_path.exists():
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(root, index_path),
                    f"missing {EVIDENCE_INDEX_DIRNAME}/{EVIDENCE_INDEX}",
                )
            )
        if not catalog_path.exists():
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(root, catalog_path),
                    f"missing {EVIDENCE_INDEX_DIRNAME}/{EVIDENCE_CATALOG}",
                )
            )
        claim_map_path = index_root / EVIDENCE_CLAIM_MAP
        if not claim_map_path.exists():
            issues.append(
                EvidenceBookValidationIssue(
                    relative_to(root, claim_map_path),
                    f"missing {EVIDENCE_INDEX_DIRNAME}/{EVIDENCE_CLAIM_MAP}",
                )
            )
        expected_index_outputs = (
            EVIDENCE_PARITY_DASHBOARD,
            EVIDENCE_PARITY_SUMMARY,
            EVIDENCE_MISMATCH_ARCHIVE,
            EVIDENCE_MISMATCH_SUMMARY,
            EVIDENCE_VERDICT_WORKFLOWS,
            EVIDENCE_VERDICT_WORKFLOWS_SUMMARY,
            EVIDENCE_FALSE_CONFIDENCE_AUDIT,
            EVIDENCE_FALSE_CONFIDENCE_SUMMARY,
            EVIDENCE_SCIENTIFIC_DEBT_REGISTER,
            EVIDENCE_SCIENTIFIC_DEBT_SUMMARY,
            EVIDENCE_PORTABILITY_AUDIT,
            EVIDENCE_PORTABILITY_SUMMARY,
            EVIDENCE_FRAGILE_EXAMPLE_AUDIT,
            EVIDENCE_FRAGILE_EXAMPLE_SUMMARY,
            EVIDENCE_REGENERATION_CONTRACT,
            EVIDENCE_REGENERATION_SUMMARY,
            TEACHING_AND_MIGRATION_INDEX_FILENAME,
            TEACHING_AND_MIGRATION_SUMMARY_FILENAME,
            EVIDENCE_CLAIM_REAUDIT,
            EVIDENCE_CLAIM_REAUDIT_SUMMARY,
            EVIDENCE_ANALYTICAL_SURFACE_COVERAGE,
            EVIDENCE_ANALYTICAL_SURFACE_COVERAGE_SUMMARY,
            EVIDENCE_CLOSURE_CRITERIA,
            EVIDENCE_CLOSURE_CRITERIA_SUMMARY,
            EVIDENCE_MATURITY_SCORECARD,
            EVIDENCE_MATURITY_SCORECARD_SUMMARY,
            EVIDENCE_REVIEW_RITUAL,
            EVIDENCE_REVIEW_RITUAL_SUMMARY,
            EVIDENCE_COMPLETION_GATES,
            EVIDENCE_COMPLETION_GATES_SUMMARY,
            COVERAGE_GAPS_JSON,
            COVERAGE_GAPS_MARKDOWN,
            FRESHNESS_REPORT_JSON,
            FRESHNESS_REPORT_MARKDOWN,
            INTEGRITY_REPORT_JSON,
            INTEGRITY_REPORT_MARKDOWN,
        )
        for relative_name in expected_index_outputs:
            output_path = index_root / relative_name
            if not output_path.exists():
                issues.append(
                    EvidenceBookValidationIssue(
                        relative_to(root, output_path),
                        f"missing {EVIDENCE_INDEX_DIRNAME}/{relative_name}",
                    )
                )
        if index_path.exists():
            try:
                index_payload = load_json(index_path)
            except (ValueError, json.JSONDecodeError) as error:
                issues.append(
                    EvidenceBookValidationIssue(
                        relative_to(root, index_path),
                        f"invalid evidence index: {error}",
                    )
                )
            else:
                discovered_paths = {
                    relative_to(root, bundle_root).as_posix()
                    for bundle_root in bundle_paths
                }
                indexed_paths = {
                    str(entry.get("relative_path"))
                    for entry in index_payload.get("evidence", [])
                    if isinstance(entry, dict)
                }
                if discovered_paths != indexed_paths:
                    issues.append(
                        EvidenceBookValidationIssue(
                            relative_to(root, index_path),
                            "evidence index must cover every discovered evidence bundle exactly once",
                        )
                    )
        false_confidence_audit_path = index_root / EVIDENCE_FALSE_CONFIDENCE_AUDIT
        if false_confidence_audit_path.exists():
            try:
                false_confidence_payload = load_json(false_confidence_audit_path)
            except (ValueError, json.JSONDecodeError) as error:
                issues.append(
                    EvidenceBookValidationIssue(
                        relative_to(root, false_confidence_audit_path),
                        f"invalid false-confidence audit: {error}",
                    )
                )
            else:
                if false_confidence_payload.get("action_required_count") != 0:
                    issues.append(
                        EvidenceBookValidationIssue(
                            relative_to(root, false_confidence_audit_path),
                            "false-confidence audit must not leave action_required entries unresolved",
                        )
                    )
        portability_audit_path = index_root / EVIDENCE_PORTABILITY_AUDIT
        if portability_audit_path.exists():
            try:
                portability_payload = load_json(portability_audit_path)
            except (ValueError, json.JSONDecodeError) as error:
                issues.append(
                    EvidenceBookValidationIssue(
                        relative_to(root, portability_audit_path),
                        f"invalid portability audit: {error}",
                    )
                )
            else:
                if portability_payload.get("action_required_count") != 0:
                    issues.append(
                        EvidenceBookValidationIssue(
                            relative_to(root, portability_audit_path),
                            "portability audit must not leave action_required entries unresolved",
                        )
                    )

    scan_for_local_path_leaks(root, issues)
    scan_for_portability_issues(root, issues)
    return EvidenceBookValidationReport(
        root=root,
        valid=not issues,
        issues=issues,
        bundle_paths=bundle_paths,
    )
