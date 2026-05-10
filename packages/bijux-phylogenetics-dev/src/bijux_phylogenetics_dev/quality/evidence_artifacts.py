"""Sync and validate governed local artifact surfaces for each evidence bundle."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from bijux_phylogenetics.evidence.bundle_artifacts import (
    ARTIFACT_JSON_FILENAMES,
    REQUIRED_BUNDLE_LOCAL_ARTIFACTS,
    build_bundle_local_artifacts,
    iter_bundle_roots,
)


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object at {path}")
    return payload


def _selected_bundle_roots(
    repo_root: Path,
    *,
    study_id: str | None = None,
    evidence_id: str | None = None,
) -> list[Path]:
    roots = iter_bundle_roots(repo_root)
    if study_id is not None:
        roots = [root for root in roots if root.parent.name == study_id]
    if evidence_id is not None:
        roots = [root for root in roots if root.name == evidence_id]
    return roots


def sync_evidence_artifacts(
    repo_root: Path,
    *,
    study_id: str | None = None,
    evidence_id: str | None = None,
) -> list[Path]:
    repo_root = repo_root.resolve()
    written: list[Path] = []
    for bundle_root in _selected_bundle_roots(
        repo_root, study_id=study_id, evidence_id=evidence_id
    ):
        artifacts = build_bundle_local_artifacts(repo_root, bundle_root)
        for filename, payload in artifacts.items():
            target = bundle_root / filename
            if filename in ARTIFACT_JSON_FILENAMES:
                _write_json(target, payload)
            else:
                if not isinstance(payload, str):
                    raise TypeError(f"expected text payload for {filename}")
                _write_text(target, payload)
            written.append(target)
    return written


def check_evidence_artifacts(
    repo_root: Path,
    *,
    study_id: str | None = None,
    evidence_id: str | None = None,
) -> list[str]:
    repo_root = repo_root.resolve()
    mismatches: list[str] = []
    for bundle_root in _selected_bundle_roots(
        repo_root, study_id=study_id, evidence_id=evidence_id
    ):
        expected = build_bundle_local_artifacts(repo_root, bundle_root)
        for filename in REQUIRED_BUNDLE_LOCAL_ARTIFACTS:
            target = bundle_root / filename
            if not target.is_file():
                mismatches.append(
                    f"{target.relative_to(repo_root)}: missing governed local artifact"
                )
                continue
            if filename in ARTIFACT_JSON_FILENAMES:
                actual = _load_json(target)
                if actual != expected[filename]:
                    mismatches.append(
                        f"{target.relative_to(repo_root)}: stale governed local artifact"
                    )
            else:
                actual_text = target.read_text(encoding="utf-8")
                expected_text = expected[filename]
                if not isinstance(expected_text, str):
                    raise TypeError(f"expected text payload for {filename}")
                if actual_text != expected_text:
                    mismatches.append(
                        f"{target.relative_to(repo_root)}: stale governed local artifact"
                    )
    return mismatches


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sync or validate bundle-local evidence artifact surfaces."
    )
    parser.add_argument("command", choices=("sync", "check"))
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--study-id", default="")
    parser.add_argument("--evidence-id", default="")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    study_id = args.study_id or None
    evidence_id = args.evidence_id or None
    if args.command == "sync":
        written = sync_evidence_artifacts(
            repo_root,
            study_id=study_id,
            evidence_id=evidence_id,
        )
        print(
            json.dumps(
                {
                    "written": [
                        path.relative_to(repo_root).as_posix() for path in written
                    ]
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    mismatches = check_evidence_artifacts(
        repo_root,
        study_id=study_id,
        evidence_id=evidence_id,
    )
    if mismatches:
        raise SystemExit(
            "evidence artifact check failed:\n" + "\n".join(mismatches)
        )
    print(
        json.dumps(
            {
                "status": "ok",
                "bundle_count": len(
                    _selected_bundle_roots(
                        repo_root,
                        study_id=study_id,
                        evidence_id=evidence_id,
                    )
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
