from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.evidence.book import write_evidence_book_index
from bijux_phylogenetics.evidence.studies.comparative_trust_boundaries import (
    STUDY_ID,
    build_comparative_trust_boundaries_bundles,
    build_comparative_trust_boundaries_claim_registry,
    build_comparative_trust_boundaries_family_index,
    build_comparative_trust_boundaries_provenance,
    build_comparative_trust_boundaries_source_fragment_map,
    render_comparative_trust_boundaries_study_manifest,
    render_comparative_trust_boundaries_study_readme,
    write_json,
    write_weak_signal_traits_table,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
STUDY_ROOT = Path(__file__).resolve().parent
PROVENANCE_ROOT = STUDY_ROOT / "provenance"


def main() -> int:
    PROVENANCE_ROOT.mkdir(parents=True, exist_ok=True)
    write_weak_signal_traits_table(REPO_ROOT)
    (STUDY_ROOT / "README.md").write_text(
        render_comparative_trust_boundaries_study_readme(),
        encoding="utf-8",
    )
    write_json(STUDY_ROOT / "study.json", render_comparative_trust_boundaries_study_manifest())
    write_json(
        PROVENANCE_ROOT / "runtime-sources.json",
        build_comparative_trust_boundaries_provenance(),
    )
    write_json(
        STUDY_ROOT / "source-fragment-map.json",
        build_comparative_trust_boundaries_source_fragment_map(),
    )
    write_json(
        STUDY_ROOT / "family-index.json",
        build_comparative_trust_boundaries_family_index(REPO_ROOT),
    )
    write_json(
        STUDY_ROOT / "claim-registry.json",
        build_comparative_trust_boundaries_claim_registry(REPO_ROOT),
    )

    for evidence_id, bundle in build_comparative_trust_boundaries_bundles(REPO_ROOT).items():
        bundle_root = STUDY_ROOT / evidence_id
        bundle_root.mkdir(parents=True, exist_ok=True)
        write_json(bundle_root / "manifest.json", bundle["manifest"])
        write_json(bundle_root / "claims.json", bundle["claims"])
        write_json(bundle_root / bundle["report_filename"], bundle["report_payload"])
        (bundle_root / "README.md").write_text(bundle["readme"], encoding="utf-8")

    write_evidence_book_index(REPO_ROOT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
