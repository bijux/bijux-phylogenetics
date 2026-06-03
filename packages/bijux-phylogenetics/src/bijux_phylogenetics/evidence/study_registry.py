from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from bijux_phylogenetics.runtime.errors import EvidenceContractError

from .studies.primate_longevity_signal import (
    EVIDENCE_ID as PRIMATE_PCM1_SUMMARY_EVIDENCE_ID,
)
from .studies.primate_longevity_signal import (
    STUDY_ID as PRIMATE_PCM1_STUDY_ID,
)
from .studies.primate_longevity_signal import (
    refresh_primate_summary_bundle,
)
from .studies.primate_pcm1_component_bundles import (
    COMPONENT_BUNDLE_DEFINITIONS,
    build_primate_pcm1_component_bundles,
)
from .studies.primate_pgls_and_signal import (
    BUNDLE_DEFINITIONS,
    SUMMARY_EVIDENCE_ID,
    build_primate_pgls_signal_bundle,
)
from .studies.primate_pgls_and_signal import (
    STUDY_ID as PRIMATE_PCM2_STUDY_ID,
)
from .study_contracts import load_study_contract


@dataclass(frozen=True, slots=True)
class EvidenceStudyRegistration:
    study_id: str
    study_title: str
    build_script_path: str | None
    supported_evidence_ids: tuple[str, ...]
    supports_partial_rerun: bool


@dataclass(slots=True)
class EvidenceStudyBuildReport:
    study_id: str
    build_script_path: str | None
    updated_paths: list[str]


@dataclass(slots=True)
class EvidenceStudyRerunReport:
    study_id: str
    selected_evidence_ids: list[str]
    updated_paths: list[str]


def _repo_root(path: Path | str) -> Path:
    return Path(path)


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _primate_pcm1_supported_evidence_ids() -> tuple[str, ...]:
    return tuple(
        sorted(
            str(definition["evidence_id"])
            for definition in COMPONENT_BUNDLE_DEFINITIONS
        )
    )


def _primate_pcm2_supported_evidence_ids() -> tuple[str, ...]:
    return tuple(
        sorted(
            {
                SUMMARY_EVIDENCE_ID,
                *(str(definition["evidence_id"]) for definition in BUNDLE_DEFINITIONS),
            }
        )
    )


def study_registrations(repo_root: Path | str) -> tuple[EvidenceStudyRegistration, ...]:
    repo_root = _repo_root(repo_root)
    return (
        EvidenceStudyRegistration(
            study_id=PRIMATE_PCM1_STUDY_ID,
            study_title=str(
                load_study_contract(
                    repo_root / "evidence-book" / "studies" / PRIMATE_PCM1_STUDY_ID
                )["study_title"]
            ),
            build_script_path=None,
            supported_evidence_ids=_primate_pcm1_supported_evidence_ids(),
            supports_partial_rerun=True,
        ),
        EvidenceStudyRegistration(
            study_id=PRIMATE_PCM2_STUDY_ID,
            study_title=str(
                load_study_contract(
                    repo_root / "evidence-book" / "studies" / PRIMATE_PCM2_STUDY_ID
                )["study_title"]
            ),
            build_script_path=None,
            supported_evidence_ids=_primate_pcm2_supported_evidence_ids(),
            supports_partial_rerun=True,
        ),
    )


def get_study_registration(
    repo_root: Path | str, study_id: str
) -> EvidenceStudyRegistration:
    for registration in study_registrations(repo_root):
        if registration.study_id == study_id:
            return registration
    raise EvidenceContractError(f"unknown evidence study: {study_id}")


def build_registered_study(
    repo_root: Path | str, study_id: str
) -> EvidenceStudyBuildReport:
    repo_root = _repo_root(repo_root)
    registration = get_study_registration(repo_root, study_id)
    rerun_report = rerun_selected_evidence(
        repo_root, study_id, list(registration.supported_evidence_ids)
    )
    return EvidenceStudyBuildReport(
        study_id=registration.study_id,
        build_script_path=registration.build_script_path,
        updated_paths=rerun_report.updated_paths,
    )


def _rerun_primate_pcm1_selection(
    repo_root: Path, evidence_ids: list[str]
) -> list[str]:
    study_root = repo_root / "evidence-book" / "studies" / PRIMATE_PCM1_STUDY_ID
    updated_paths: list[str] = []
    bundles = build_primate_pcm1_component_bundles(repo_root)
    for evidence_id in evidence_ids:
        bundle = bundles[evidence_id]
        bundle_root = study_root / evidence_id
        _write_json(bundle_root / "manifest.json", bundle["manifest"])
        _write_json(bundle_root / "claims.json", bundle["claims"])
        _write_json(
            bundle_root / "results" / bundle["report_filename"],
            bundle["report_payload"],
        )
        (bundle_root / "README.md").write_text(bundle["readme"], encoding="utf-8")
        updated_paths.extend(
            [
                str((bundle_root / "manifest.json").relative_to(repo_root)),
                str((bundle_root / "claims.json").relative_to(repo_root)),
                str(
                    (bundle_root / "results" / bundle["report_filename"]).relative_to(
                        repo_root
                    )
                ),
                str((bundle_root / "README.md").relative_to(repo_root)),
            ]
        )
    summary_bundle = refresh_primate_summary_bundle(repo_root)
    summary_bundle_root = study_root / PRIMATE_PCM1_SUMMARY_EVIDENCE_ID
    _write_json(summary_bundle_root / "manifest.json", summary_bundle["manifest"])
    _write_json(summary_bundle_root / "claims.json", summary_bundle["claims"])
    updated_paths.extend(
        [
            str((summary_bundle_root / "manifest.json").relative_to(repo_root)),
            str((summary_bundle_root / "claims.json").relative_to(repo_root)),
        ]
    )
    return sorted(set(updated_paths))


def _rerun_primate_pcm2_selection(
    repo_root: Path, evidence_ids: list[str]
) -> list[str]:
    study_root = repo_root / "evidence-book" / "studies" / PRIMATE_PCM2_STUDY_ID
    updated_paths: list[str] = []
    for evidence_id in evidence_ids:
        bundle = build_primate_pgls_signal_bundle(repo_root, evidence_id)
        bundle_root = study_root / evidence_id
        bundle_root.mkdir(parents=True, exist_ok=True)
        _write_json(bundle_root / "manifest.json", bundle["manifest"])
        _write_json(bundle_root / "claims.json", bundle["claims"])
        _write_json(
            bundle_root / "results" / bundle["report_filename"],
            bundle["report_payload"],
        )
        (bundle_root / "README.md").write_text(bundle["readme"], encoding="utf-8")
        updated_paths.extend(
            [
                str((bundle_root / "manifest.json").relative_to(repo_root)),
                str((bundle_root / "claims.json").relative_to(repo_root)),
                str(
                    (bundle_root / "results" / bundle["report_filename"]).relative_to(
                        repo_root
                    )
                ),
                str((bundle_root / "README.md").relative_to(repo_root)),
            ]
        )
        if evidence_id == "evidence-001":
            _write_json(
                bundle_root / "results" / "scalar-parity-table.json",
                bundle["scalar_parity_table"],
            )
            (bundle_root / "results" / "scalar-parity-table.md").write_text(
                bundle["scalar_parity_markdown"],
                encoding="utf-8",
            )
            updated_paths.extend(
                [
                    str(
                        (
                            bundle_root / "results" / "scalar-parity-table.json"
                        ).relative_to(repo_root)
                    ),
                    str(
                        (
                            bundle_root / "results" / "scalar-parity-table.md"
                        ).relative_to(repo_root)
                    ),
                ]
            )
    return sorted(set(updated_paths))


def rerun_selected_evidence(
    repo_root: Path | str, study_id: str, evidence_ids: list[str]
) -> EvidenceStudyRerunReport:
    repo_root = _repo_root(repo_root)
    registration = get_study_registration(repo_root, study_id)
    if not registration.supports_partial_rerun:
        raise EvidenceContractError(
            f"study does not support partial evidence rerun: {study_id}"
        )
    selected = sorted(set(evidence_ids))
    unsupported = sorted(set(selected) - set(registration.supported_evidence_ids))
    if unsupported:
        raise EvidenceContractError(
            f"unsupported evidence id(s) for {study_id}: {', '.join(unsupported)}"
        )
    if study_id == PRIMATE_PCM1_STUDY_ID:
        updated_paths = _rerun_primate_pcm1_selection(repo_root, selected)
    elif study_id == PRIMATE_PCM2_STUDY_ID:
        updated_paths = _rerun_primate_pcm2_selection(repo_root, selected)
    else:
        raise EvidenceContractError(
            f"partial evidence rerun not implemented for {study_id}"
        )
    return EvidenceStudyRerunReport(
        study_id=study_id,
        selected_evidence_ids=selected,
        updated_paths=updated_paths,
    )
