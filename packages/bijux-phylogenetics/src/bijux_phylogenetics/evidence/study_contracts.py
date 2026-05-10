from __future__ import annotations

import json
from pathlib import Path

from .teaching import study_metadata, teaching_study_ids


ALLOWED_STUDY_ROOT_DIRS = {"datasets", "reference", "provenance"}
ALLOWED_STUDY_ROOT_FILES = {"README.md"}
DEFAULT_OWNER_PACKAGE = "bijux-phylogenetics"


def _load_json(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object at {path}")
    return payload


def _repo_root_from_study(study_root: Path) -> Path:
    return study_root.parents[2]


def _read_markdown_title(readme_path: Path, *, fallback: str) -> str:
    for raw_line in readme_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line.startswith("# "):
            return line.removeprefix("# ").strip() or fallback
    return fallback


def _read_markdown_summary(readme_path: Path, *, fallback: str) -> str:
    lines = readme_path.read_text(encoding="utf-8").splitlines()
    after_heading = False
    paragraph_lines: list[str] = []
    for raw_line in lines:
        line = raw_line.strip()
        if not after_heading:
            if line.startswith("# "):
                after_heading = True
            continue
        if not line:
            if paragraph_lines:
                break
            continue
        if line.startswith("#"):
            break
        paragraph_lines.append(line)
    return " ".join(paragraph_lines) if paragraph_lines else fallback


def _study_source_basis(
    provenance_payload: dict[str, object] | None,
) -> list[dict[str, object]]:
    if provenance_payload is None:
        return []
    sources = provenance_payload.get("sources")
    if not isinstance(sources, list):
        return []
    rows: list[dict[str, object]] = []
    for source in sources:
        if not isinstance(source, dict):
            continue
        locator = source.get("locator")
        if not isinstance(locator, str) or not locator:
            continue
        row = {
            "kind": source.get("kind", "source"),
            "label": source.get("label", locator),
            "locator": locator,
        }
        rows.append(row)
    return rows


def _study_scope(
    provenance_payload: dict[str, object] | None,
    dataset_payload: dict[str, object] | None,
) -> dict[str, object]:
    coverage_focus: set[str] = set()
    untouched_source_locators: list[str] = []
    if provenance_payload is not None:
        for source in provenance_payload.get("sources", []):
            if not isinstance(source, dict):
                continue
            locator = source.get("locator")
            if source.get("read_only") is True and isinstance(locator, str) and locator:
                untouched_source_locators.append(locator)
            for concept_tag in source.get("concept_tags", []):
                if isinstance(concept_tag, str) and concept_tag:
                    coverage_focus.add(concept_tag)
    if dataset_payload is not None:
        for dataset in dataset_payload.get("datasets", []):
            if not isinstance(dataset, dict):
                continue
            dataset_id = dataset.get("dataset_id")
            if isinstance(dataset_id, str) and dataset_id:
                coverage_focus.add(dataset_id)
    return {
        "coverage_focus": sorted(coverage_focus),
        "untouched_source_locators": untouched_source_locators,
    }


def load_study_contract(study_root: Path) -> dict[str, object]:
    study_id = study_root.name
    repo_root = _repo_root_from_study(study_root)
    readme_path = study_root / "README.md"
    fallback_title = study_id.replace("-", " ").title()
    fallback_summary = f"Governed evidence study for {study_id}."
    title = (
        _read_markdown_title(readme_path, fallback=fallback_title)
        if readme_path.exists()
        else fallback_title
    )
    summary = (
        _read_markdown_summary(readme_path, fallback=fallback_summary)
        if readme_path.exists()
        else fallback_summary
    )
    provenance_paths = sorted(
        path for path in (study_root / "provenance").glob("*.json") if path.is_file()
    )
    provenance_payload = (
        _load_json(provenance_paths[0]) if len(provenance_paths) == 1 else None
    )
    dataset_registry_path = study_root / "datasets" / "registry.json"
    dataset_payload = _load_json(dataset_registry_path) if dataset_registry_path.is_file() else None
    if study_id in teaching_study_ids():
        metadata = study_metadata(study_id)
        study_categories = list(metadata["study_categories"])
    else:
        study_categories = ["scientific-validation"]
    return {
        "study_id": study_id,
        "study_title": title,
        "summary": summary,
        "owner_package": DEFAULT_OWNER_PACKAGE,
        "study_categories": study_categories,
        "source_intake_policy": ""
        if provenance_payload is None
        else str(provenance_payload.get("intake_policy", "")),
        "provenance_descriptor_locator": ""
        if len(provenance_paths) != 1
        else provenance_paths[0].relative_to(repo_root).as_posix(),
        "dataset_registry_locator": ""
        if not dataset_registry_path.is_file()
        else dataset_registry_path.relative_to(repo_root).as_posix(),
        "source_basis": _study_source_basis(provenance_payload),
        "study_scope": _study_scope(provenance_payload, dataset_payload),
    }
