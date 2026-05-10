from __future__ import annotations

from bijux_phylogenetics.evidence.teaching import (
    build_migration_guide,
    build_student_safe_reproducibility_contract,
    build_teaching_and_migration_index,
    build_teaching_guide,
    render_migration_guide_markdown,
    render_student_safe_reproducibility_markdown,
    render_teaching_and_migration_index_markdown,
    render_teaching_guide_markdown,
)


def _study_manifest() -> dict[str, object]:
    return {
        "study_id": "primate-pgls-and-signal",
        "study_title": "Primate PGLS and signal evidence study",
        "study_categories": ["teaching-study", "migration-study"],
    }


def _family_index() -> dict[str, object]:
    return {
        "families": [
            {
                "family_id": "baseline-regression",
                "family_title": "Baseline regression",
                "family_verdict": "matched",
                "coverage_status": "covered",
                "evidence_ids": ["evidence-002"],
                "fragment_ids": ["baseline-gls-fit"],
            },
            {
                "family_id": "diagnostics",
                "family_title": "Diagnostics",
                "family_verdict": "matched_with_tolerance",
                "coverage_status": "covered",
                "evidence_ids": ["evidence-005"],
                "fragment_ids": ["estimated-lambda-diagnostics"],
            },
        ]
    }


def _source_fragment_map() -> dict[str, object]:
    return {
        "fragments": [
            {
                "fragment_id": "baseline-gls-fit",
                "fragment_title": "Non-phylogenetic GLS fit",
                "script_locators": ["external:lund/pcm2-modes-pgls/script#L122-L136"],
            },
            {
                "fragment_id": "pagel-lambda-regression",
                "fragment_title": "Fixed-lambda and estimated-lambda PGLS fits",
                "script_locators": ["external:lund/pcm2-modes-pgls/script#L138-L179"],
            },
            {
                "fragment_id": "phylogenetic-signal-test",
                "fragment_title": "Intercept-only PGLS and lambda-zero likelihood-ratio testing",
                "script_locators": ["external:lund/pcm2-modes-pgls/script#L181-L192"],
            },
            {
                "fragment_id": "estimated-lambda-diagnostics",
                "fragment_title": "Estimated-lambda residual and fitted diagnostics",
                "script_locators": ["external:lund/pcm2-modes-pgls/script#L168-L179"],
            },
        ]
    }


def _bundle_manifests() -> list[dict[str, object]]:
    return [
        {"comparison_mode": "direct_parity"},
        {"comparison_mode": "direct_parity"},
        {"comparison_mode": "direct_parity"},
    ]


def test_teaching_guide_builds_family_narratives_and_topic_tags() -> None:
    payload = build_teaching_guide(_study_manifest(), _family_index(), _source_fragment_map())

    assert payload["study_id"] == "primate-pgls-and-signal"
    assert payload["family_count"] == 2
    assert "gls" in payload["concept_tags"]
    assert payload["families"][0]["teaching_narrative"]
    assert "external:lund/pcm2-modes-pgls/script#L122-L136" in payload["families"][0]["source_locators"]
    assert "Baseline regression" in render_teaching_guide_markdown(payload)


def test_migration_guide_builds_side_by_side_examples() -> None:
    payload = build_migration_guide(
        _study_manifest(),
        _source_fragment_map(),
        _bundle_manifests(),
    )

    assert payload["comparison_mode_counts"] == {"direct_parity": 3}
    assert payload["example_count"] == 4
    assert payload["examples"][0]["r_source_locators"]
    assert payload["examples"][0]["bijux_locators"]
    assert "Why migrate" in render_migration_guide_markdown(payload)


def test_student_safe_reproducibility_contract_surfaces_portable_scope() -> None:
    payload = build_student_safe_reproducibility_contract(_study_manifest())

    assert payload["entrypoint_command"].startswith(
        "UV_PROJECT_ENVIRONMENT=artifacts/root/venv"
    )
    assert "no workstation-local /Users paths" in payload["forbidden_assumptions"]
    assert "Student-Safe Reproducibility" in render_student_safe_reproducibility_markdown(
        payload
    )


def test_teaching_and_migration_index_aggregates_guides() -> None:
    teaching_payload = build_teaching_guide(
        _study_manifest(),
        _family_index(),
        _source_fragment_map(),
    )
    migration_payload = build_migration_guide(
        _study_manifest(),
        _source_fragment_map(),
        _bundle_manifests(),
    )
    contract_payload = build_student_safe_reproducibility_contract(_study_manifest())

    index_payload = build_teaching_and_migration_index(
        [teaching_payload],
        [migration_payload],
        [contract_payload],
    )

    assert index_payload["teaching_study_count"] == 1
    assert index_payload["migration_study_count"] == 1
    assert index_payload["course_topic_counts"]["gls"] == 1
    assert "Teaching And Migration Evidence" in render_teaching_and_migration_index_markdown(
        index_payload
    )
