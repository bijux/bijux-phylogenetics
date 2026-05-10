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
    build_primate_pgls_signal_external_sources,
    build_primate_pgls_signal_family_index,
    build_primate_pgls_signal_parity_policy,
    build_primate_pgls_signal_source_fragment_map,
)


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


def _study_manifest() -> dict[str, object]:
    return {
        "study_id": STUDY_ID,
        "study_title": "Primate PGLS and signal evidence study",
        "summary": "Governed parity study for the regression and signal sections of the Lund primate comparative lecture, with explicit coverage boundaries for the remaining EB and ancestral fragments.",
        "owner_package": "bijux-phylogenetics",
        "study_categories": ["teaching-study", "migration-study"],
        "confidence_posture": "governed-parity-in-progress",
        "coverage_boundary_evidence_ids": ["evidence-006"],
    }


def _study_readme() -> str:
    lines = [
        "# Primate PGLS And Signal",
        "",
        "This study turns the regression and signal sections of the Lund primate",
        "comparative lecture into governed Evidence IDs backed by checked-in R",
        "reference outputs and canonical `bijux-phylogenetics` reproductions.",
        "",
        "It is intentionally strict about confidence posture:",
        "",
        "- baseline GLS, Pagel-lambda PGLS, signal testing, and scalar diagnostics",
        "  are backed by governed parity bundles",
        "- transformed-tree, EB-model, and ancestral-mode fragments are kept visible",
        "  as explicit coverage boundaries instead of being implied as validated",
        "",
        "Current bundles:",
        "",
        "- `evidence-001` reload semantics and study summary",
        "- `evidence-002` baseline GLS parity",
        "- `evidence-003` Pagel-lambda regression parity",
        "- `evidence-004` phylogenetic signal parity",
        "- `evidence-005` residual diagnostics parity",
        "- `evidence-006` coverage boundaries for the still-open lecture fragments",
        "",
    ]
    return "\n".join(lines)


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    REFERENCE_ROOT.mkdir(parents=True, exist_ok=True)
    PROVENANCE_ROOT.mkdir(parents=True, exist_ok=True)
    _run_r_reference_script()

    (STUDY_ROOT / "README.md").write_text(_study_readme(), encoding="utf-8")
    _write_json(STUDY_ROOT / "study.json", _study_manifest())
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

    write_evidence_book_index(REPO_ROOT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
