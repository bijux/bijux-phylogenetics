from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess  # nosec B404
import sys

from bijux_phylogenetics.evidence.book import write_evidence_book_index
from bijux_phylogenetics.evidence.studies.primate_pgls_and_signal import (
    STUDY_ID,
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
from bijux_phylogenetics_dev.quality.evidence_artifacts import sync_evidence_artifacts
from bijux_phylogenetics_dev.quality.evidence_inputs import sync_inputs_manifests


REPO_ROOT = Path(__file__).resolve().parents[3]
STUDY_ROOT = Path(__file__).resolve().parent
REFERENCE_ROOT = STUDY_ROOT / "reference"
PROVENANCE_ROOT = STUDY_ROOT / "provenance"
R_SCRIPT = REFERENCE_ROOT / "primate_pgls_and_signal_reference_r.R"
R_OUTPUT = REFERENCE_ROOT / "reference_results.json"


def _run_r_reference_script() -> None:
    env = dict(os.environ)
    default_r_lib = REPO_ROOT / "artifacts" / "root" / "r-library"
    if "R_LIBS_USER" not in env and default_r_lib.exists():
        env["R_LIBS_USER"] = str(default_r_lib)
    subprocess.run(  # nosec B603
        ["Rscript", str(R_SCRIPT), str(REPO_ROOT), str(R_OUTPUT)],
        check=True,
        cwd=str(REPO_ROOT),
        env=env,
    )

def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    REFERENCE_ROOT.mkdir(parents=True, exist_ok=True)
    PROVENANCE_ROOT.mkdir(parents=True, exist_ok=True)
    _run_r_reference_script()

    (STUDY_ROOT / "README.md").write_text(
        render_primate_pgls_signal_study_readme(REPO_ROOT),
        encoding="utf-8",
    )
    _write_json(
        STUDY_ROOT / "study.json",
        render_primate_pgls_signal_study_manifest(REPO_ROOT),
    )
    _write_json(
        STUDY_ROOT / "evidence-registry.json",
        build_primate_pgls_signal_evidence_registry(REPO_ROOT),
    )
    _write_json(
        PROVENANCE_ROOT / "lund-course-sources.json",
        build_primate_pgls_signal_external_sources(),
    )
    _write_json(
        STUDY_ROOT / "source-fragment-map.json",
        build_primate_pgls_signal_source_fragment_map(),
    )
    _write_json(
        STUDY_ROOT / "family-index.json",
        build_primate_pgls_signal_family_index(REPO_ROOT),
    )
    _write_json(
        STUDY_ROOT / "parity-policy.json",
        build_primate_pgls_signal_parity_policy(),
    )
    _write_json(
        STUDY_ROOT / "claim-registry.json",
        build_primate_pgls_signal_claim_registry(REPO_ROOT),
    )

    for evidence_id, bundle in build_primate_pgls_signal_bundles(REPO_ROOT).items():
        bundle_root = STUDY_ROOT / evidence_id
        bundle_root.mkdir(parents=True, exist_ok=True)
        _write_json(bundle_root / "manifest.json", bundle["manifest"])
        _write_json(bundle_root / "claims.json", bundle["claims"])
        _write_json(bundle_root / bundle["report_filename"], bundle["report_payload"])
        (bundle_root / "README.md").write_text(bundle["readme"], encoding="utf-8")
        if evidence_id == "evidence-001":
            _write_json(
                bundle_root / "scalar-parity-table.json", bundle["scalar_parity_table"]
            )
            (bundle_root / "scalar-parity-table.md").write_text(
                bundle["scalar_parity_markdown"], encoding="utf-8"
            )

    sync_evidence_artifacts(REPO_ROOT, study_id=STUDY_ID)
    sync_inputs_manifests(REPO_ROOT, study_id=STUDY_ID)
    write_evidence_book_index(REPO_ROOT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
