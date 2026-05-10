from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import subprocess  # nosec B404
import sys

from bijux_phylogenetics.errors import EvidenceContractError

from .studies.comparative_trust_boundaries import (
    STUDY_ID as COMPARATIVE_STUDY_ID,
    build_comparative_trust_boundaries_bundles,
    build_comparative_trust_boundaries_claim_registry,
    build_comparative_trust_boundaries_family_index,
    build_comparative_trust_boundaries_provenance,
    build_comparative_trust_boundaries_source_fragment_map,
    render_comparative_trust_boundaries_study_manifest,
    render_comparative_trust_boundaries_study_readme,
    write_json as write_comparative_json,
    write_weak_signal_traits_table,
)
from .studies.primate_longevity_signal import (
    EVIDENCE_ID as PRIMATE_PCM1_SUMMARY_EVIDENCE_ID,
    STUDY_ID as PRIMATE_PCM1_STUDY_ID,
    build_primate_claim_registry,
    build_primate_family_index,
    build_primate_parity_policy,
    build_primate_source_fragment_map,
    build_primate_summary_bundle_claims,
)
from .studies.primate_pcm1_component_bundles import (
    build_primate_data_preparation_bundle_index,
    build_primate_pcm1_component_bundles,
    build_primate_structural_parity_table,
    render_primate_data_preparation_bundle_index_markdown,
    render_primate_structural_parity_table_markdown,
)
from .studies.primate_pgls_and_signal import (
    STUDY_ID as PRIMATE_PCM2_STUDY_ID,
    build_primate_pgls_signal_bundles,
    build_primate_pgls_signal_claim_registry,
    build_primate_pgls_signal_evidence_registry,
    build_primate_pgls_signal_external_sources,
    build_primate_pgls_signal_family_index,
    build_primate_pgls_signal_parity_policy,
    build_primate_pgls_signal_source_fragment_map,
    render_primate_pgls_signal_study_manifest,
    render_primate_pgls_signal_study_readme,
)


@dataclass(frozen=True, slots=True)
class EvidenceStudyRegistration:
    study_id: str
    study_title: str
    build_script_path: str
    supported_evidence_ids: tuple[str, ...]
    supports_partial_rerun: bool


@dataclass(slots=True)
class EvidenceStudyBuildReport:
    study_id: str
    build_script_path: str
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


def study_registrations(repo_root: Path | str) -> tuple[EvidenceStudyRegistration, ...]:
    repo_root = _repo_root(repo_root)
    return (
        EvidenceStudyRegistration(
            study_id=COMPARATIVE_STUDY_ID,
            study_title="Comparative trust boundary evidence study",
            build_script_path="evidence-book/studies/comparative-trust-boundaries/build_evidence.py",
            supported_evidence_ids=tuple(
                sorted(build_comparative_trust_boundaries_bundles(repo_root))
            ),
            supports_partial_rerun=True,
        ),
        EvidenceStudyRegistration(
            study_id=PRIMATE_PCM1_STUDY_ID,
            study_title="Primate longevity signal evidence study",
            build_script_path="evidence-book/studies/primate-longevity-signal/build_evidence.py",
            supported_evidence_ids=tuple(
                sorted(
                    evidence_id
                    for evidence_id in build_primate_pcm1_component_bundles(repo_root)
                    if evidence_id != PRIMATE_PCM1_SUMMARY_EVIDENCE_ID
                )
            ),
            supports_partial_rerun=True,
        ),
        EvidenceStudyRegistration(
            study_id=PRIMATE_PCM2_STUDY_ID,
            study_title="Primate PGLS and signal evidence study",
            build_script_path="evidence-book/studies/primate-pgls-and-signal/build_evidence.py",
            supported_evidence_ids=tuple(
                sorted(build_primate_pgls_signal_bundles(repo_root))
            ),
            supports_partial_rerun=True,
        ),
        EvidenceStudyRegistration(
            study_id="taxon-trust",
            study_title="Taxon trust evidence study",
            build_script_path="evidence-book/studies/taxon-trust/build_evidence.py",
            supported_evidence_ids=(),
            supports_partial_rerun=False,
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
    script_path = repo_root / registration.build_script_path
    subprocess.run(  # nosec B603
        [sys.executable, str(script_path)],
        check=True,
        cwd=str(repo_root),
    )
    return EvidenceStudyBuildReport(
        study_id=registration.study_id,
        build_script_path=registration.build_script_path,
        updated_paths=[registration.build_script_path],
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
        _write_json(bundle_root / bundle["report_filename"], bundle["report_payload"])
        (bundle_root / "README.md").write_text(bundle["readme"], encoding="utf-8")
        updated_paths.extend(
            [
                str((bundle_root / "manifest.json").relative_to(repo_root)),
                str((bundle_root / "claims.json").relative_to(repo_root)),
                str((bundle_root / bundle["report_filename"]).relative_to(repo_root)),
                str((bundle_root / "README.md").relative_to(repo_root)),
            ]
        )
    _write_json(
        study_root / "source-fragment-map.json",
        build_primate_source_fragment_map(repo_root),
    )
    _write_json(study_root / "family-index.json", build_primate_family_index(repo_root))
    _write_json(
        study_root / "claim-registry.json", build_primate_claim_registry(repo_root)
    )
    _write_json(
        study_root / "parity-policy.json", build_primate_parity_policy(repo_root)
    )
    _write_json(
        study_root / PRIMATE_PCM1_SUMMARY_EVIDENCE_ID / "claims.json",
        build_primate_summary_bundle_claims(repo_root),
    )
    structural_table = build_primate_structural_parity_table(repo_root)
    _write_json(study_root / "structural-parity-table.json", structural_table)
    (study_root / "structural-parity-table.md").write_text(
        render_primate_structural_parity_table_markdown(structural_table),
        encoding="utf-8",
    )
    data_preparation = build_primate_data_preparation_bundle_index(repo_root)
    _write_json(study_root / "data-preparation-bundles.json", data_preparation)
    (study_root / "data-preparation-bundles.md").write_text(
        render_primate_data_preparation_bundle_index_markdown(data_preparation),
        encoding="utf-8",
    )
    updated_paths.extend(
        [
            "evidence-book/studies/primate-longevity-signal/source-fragment-map.json",
            "evidence-book/studies/primate-longevity-signal/family-index.json",
            "evidence-book/studies/primate-longevity-signal/claim-registry.json",
            "evidence-book/studies/primate-longevity-signal/parity-policy.json",
            "evidence-book/studies/primate-longevity-signal/evidence-001/claims.json",
            "evidence-book/studies/primate-longevity-signal/structural-parity-table.json",
            "evidence-book/studies/primate-longevity-signal/structural-parity-table.md",
            "evidence-book/studies/primate-longevity-signal/data-preparation-bundles.json",
            "evidence-book/studies/primate-longevity-signal/data-preparation-bundles.md",
        ]
    )
    return sorted(set(updated_paths))


def _rerun_primate_pcm2_selection(
    repo_root: Path, evidence_ids: list[str]
) -> list[str]:
    study_root = repo_root / "evidence-book" / "studies" / PRIMATE_PCM2_STUDY_ID
    updated_paths: list[str] = []
    (study_root / "README.md").write_text(
        render_primate_pgls_signal_study_readme(repo_root),
        encoding="utf-8",
    )
    _write_json(
        study_root / "study.json",
        render_primate_pgls_signal_study_manifest(repo_root),
    )
    _write_json(
        study_root / "evidence-registry.json",
        build_primate_pgls_signal_evidence_registry(repo_root),
    )
    bundles = build_primate_pgls_signal_bundles(repo_root)
    for evidence_id in evidence_ids:
        bundle = bundles[evidence_id]
        bundle_root = study_root / evidence_id
        bundle_root.mkdir(parents=True, exist_ok=True)
        _write_json(bundle_root / "manifest.json", bundle["manifest"])
        _write_json(bundle_root / "claims.json", bundle["claims"])
        _write_json(bundle_root / bundle["report_filename"], bundle["report_payload"])
        (bundle_root / "README.md").write_text(bundle["readme"], encoding="utf-8")
        updated_paths.extend(
            [
                str((bundle_root / "manifest.json").relative_to(repo_root)),
                str((bundle_root / "claims.json").relative_to(repo_root)),
                str((bundle_root / bundle["report_filename"]).relative_to(repo_root)),
                str((bundle_root / "README.md").relative_to(repo_root)),
            ]
        )
        if evidence_id == "evidence-001":
            _write_json(
                bundle_root / "scalar-parity-table.json",
                bundle["scalar_parity_table"],
            )
            (bundle_root / "scalar-parity-table.md").write_text(
                bundle["scalar_parity_markdown"],
                encoding="utf-8",
            )
            updated_paths.extend(
                [
                    str(
                        (bundle_root / "scalar-parity-table.json").relative_to(
                            repo_root
                        )
                    ),
                    str(
                        (bundle_root / "scalar-parity-table.md").relative_to(repo_root)
                    ),
                ]
            )
    provenance_root = study_root / "provenance"
    provenance_root.mkdir(exist_ok=True)
    _write_json(
        provenance_root / "lund-course-sources.json",
        build_primate_pgls_signal_external_sources(),
    )
    _write_json(
        study_root / "source-fragment-map.json",
        build_primate_pgls_signal_source_fragment_map(),
    )
    _write_json(
        study_root / "family-index.json",
        build_primate_pgls_signal_family_index(repo_root),
    )
    _write_json(
        study_root / "parity-policy.json", build_primate_pgls_signal_parity_policy()
    )
    _write_json(
        study_root / "claim-registry.json",
        build_primate_pgls_signal_claim_registry(repo_root),
    )
    updated_paths.extend(
        [
            "evidence-book/studies/primate-pgls-and-signal/README.md",
            "evidence-book/studies/primate-pgls-and-signal/study.json",
            "evidence-book/studies/primate-pgls-and-signal/evidence-registry.json",
            "evidence-book/studies/primate-pgls-and-signal/provenance/lund-course-sources.json",
            "evidence-book/studies/primate-pgls-and-signal/source-fragment-map.json",
            "evidence-book/studies/primate-pgls-and-signal/family-index.json",
            "evidence-book/studies/primate-pgls-and-signal/parity-policy.json",
            "evidence-book/studies/primate-pgls-and-signal/claim-registry.json",
        ]
    )
    return sorted(set(updated_paths))


def _rerun_comparative_selection(repo_root: Path, evidence_ids: list[str]) -> list[str]:
    study_root = repo_root / "evidence-book" / "studies" / COMPARATIVE_STUDY_ID
    updated_paths: list[str] = []
    (study_root / "README.md").write_text(
        render_comparative_trust_boundaries_study_readme(),
        encoding="utf-8",
    )
    write_comparative_json(
        study_root / "study.json",
        render_comparative_trust_boundaries_study_manifest(),
    )
    write_weak_signal_traits_table(repo_root)
    provenance_root = study_root / "provenance"
    provenance_root.mkdir(exist_ok=True)
    write_comparative_json(
        provenance_root / "runtime-sources.json",
        build_comparative_trust_boundaries_provenance(),
    )
    write_comparative_json(
        study_root / "source-fragment-map.json",
        build_comparative_trust_boundaries_source_fragment_map(),
    )
    write_comparative_json(
        study_root / "family-index.json",
        build_comparative_trust_boundaries_family_index(repo_root),
    )
    write_comparative_json(
        study_root / "claim-registry.json",
        build_comparative_trust_boundaries_claim_registry(repo_root),
    )
    bundles = build_comparative_trust_boundaries_bundles(repo_root)
    for evidence_id in evidence_ids:
        bundle = bundles[evidence_id]
        bundle_root = study_root / evidence_id
        write_comparative_json(bundle_root / "manifest.json", bundle["manifest"])
        write_comparative_json(bundle_root / "claims.json", bundle["claims"])
        write_comparative_json(
            bundle_root / bundle["report_filename"], bundle["report_payload"]
        )
        (bundle_root / "README.md").write_text(bundle["readme"], encoding="utf-8")
        updated_paths.extend(
            [
                str((bundle_root / "manifest.json").relative_to(repo_root)),
                str((bundle_root / "claims.json").relative_to(repo_root)),
                str((bundle_root / bundle["report_filename"]).relative_to(repo_root)),
                str((bundle_root / "README.md").relative_to(repo_root)),
            ]
        )
    updated_paths.extend(
        [
            "evidence-book/studies/comparative-trust-boundaries/README.md",
            "evidence-book/studies/comparative-trust-boundaries/study.json",
            "evidence-book/studies/comparative-trust-boundaries/provenance/runtime-sources.json",
            "evidence-book/studies/comparative-trust-boundaries/source-fragment-map.json",
            "evidence-book/studies/comparative-trust-boundaries/family-index.json",
            "evidence-book/studies/comparative-trust-boundaries/claim-registry.json",
            "evidence-book/studies/comparative-trust-boundaries/evidence-002/weak_signal_traits.tsv",
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
    elif study_id == COMPARATIVE_STUDY_ID:
        updated_paths = _rerun_comparative_selection(repo_root, selected)
    else:
        raise EvidenceContractError(
            f"partial evidence rerun not implemented for {study_id}"
        )
    return EvidenceStudyRerunReport(
        study_id=study_id,
        selected_evidence_ids=selected,
        updated_paths=updated_paths,
    )
