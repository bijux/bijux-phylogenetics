from __future__ import annotations

from collections import Counter
from copy import deepcopy

TEACHING_GUIDE_FILENAME = "teaching-guide.json"
TEACHING_GUIDE_MARKDOWN_FILENAME = "teaching-guide.md"
MIGRATION_GUIDE_FILENAME = "migration-guide.json"
MIGRATION_GUIDE_MARKDOWN_FILENAME = "migration-guide.md"
STUDENT_SAFE_REPRODUCIBILITY_FILENAME = "student-safe-reproducibility.json"
STUDENT_SAFE_REPRODUCIBILITY_MARKDOWN_FILENAME = "student-safe-reproducibility.md"
TEACHING_AND_MIGRATION_INDEX_FILENAME = "teaching-and-migration.json"
TEACHING_AND_MIGRATION_SUMMARY_FILENAME = "teaching-and-migration.md"

ALLOWED_STUDY_CATEGORIES = {
    "migration-study",
    "scientific-validation",
    "teaching-study",
}
ALLOWED_COMPARISON_MODES = {
    "bijux_native_reinterpretation",
    "direct_parity",
}

TEACHING_STUDY_METADATA = {
    "primate-longevity-signal": {
        "course_material_label": "Lund PCM1 plots and signal lecture",
        "course_material_locator": "external:lund/pcm1-plots-signal/script",
        "study_categories": ["teaching-study", "migration-study"],
        "student_safe_reproducibility": {
            "supported_scope": (
                "regenerate the governed teaching and migration summaries from the "
                "checked-in evidence bundle and study indexes"
            ),
            "entrypoint_command": (
                "UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 "
                'python -c "from pathlib import Path; '
                "from bijux_phylogenetics.evidence.book import write_evidence_book_index; "
                "write_evidence_book_index(Path('.'))\""
            ),
            "portable_prerequisites": [
                "Python 3.11",
                "UV project environment at artifacts/root/venv",
            ],
            "forbidden_assumptions": [
                "no workstation-local /Users paths",
                "no sibling repository checkout assumptions",
                "no hidden lecture data paths outside the repository",
            ],
            "expected_paths": [
                "evidence-book/index/teaching-and-migration.json",
            ],
            "not_claimed": [
                "raw Lund workbook reconstruction is not re-executed from inside this repository",
                "rendered-figure equivalence is still outside the current teaching trust claim",
            ],
        },
        "families": {
            "workflow-contracts": {
                "concept_tags": [
                    "reproducibility",
                    "package-loading",
                    "course-context",
                ],
                "teaching_narrative": (
                    "Students see which package and environment assumptions belong to the "
                    "lecture setup before any numerical claim is made."
                ),
            },
            "data-preparation": {
                "concept_tags": [
                    "data-import",
                    "type-repair",
                    "missing-data",
                    "species-aggregation",
                ],
                "teaching_narrative": (
                    "The lecture data-cleaning path is broken into explicit steps so "
                    "students can inspect what must match before comparative inference "
                    "is trusted."
                ),
                "migration_examples": [
                    {
                        "example_id": "pcm1-data-preparation",
                        "family_id": "data-preparation",
                        "fragment_id": "primate-data-preprocessing",
                        "bijux_locators": [
                            "bijux_phylogenetics.evidence.studies.primate_pcm1_component_bundles:build_primate_pcm1_component_bundles",
                            "bijux_phylogenetics.evidence.studies.primate_longevity_signal:build_primate_source_fragment_map",
                        ],
                        "bijux_summary": (
                            "Bijux records the cleaned primate table as a governed "
                            "evidence surface instead of burying preprocessing behind "
                            "later model output."
                        ),
                        "why_migrate_note": (
                            "The Python-side evidence makes type repair, missing-data "
                            "accounting, and grouped species decisions machine-checkable."
                        ),
                    }
                ],
            },
            "tree-operations": {
                "concept_tags": [
                    "tree-import",
                    "pruning",
                    "topology",
                    "tree-data-alignment",
                ],
                "teaching_narrative": (
                    "Tree loading, pruning, and topology operations are separated from "
                    "plotting so students can review structure before interpretation."
                ),
                "migration_examples": [
                    {
                        "example_id": "pcm1-tree-operations",
                        "family_id": "tree-operations",
                        "fragment_id": "tree-import-and-pruning",
                        "bijux_locators": [
                            "bijux_phylogenetics.phylo.pruning:prune_tree_to_requested_taxa",
                            "bijux_phylogenetics.evidence.studies.primate_pcm1_component_bundles:build_primate_pcm1_component_bundles",
                        ],
                        "bijux_summary": (
                            "The tree-trimming and alignment surface is governed by "
                            "structural parity outputs rather than visual inspection."
                        ),
                        "why_migrate_note": (
                            "Students can move from interactive R tree surgery to a "
                            "scriptable tree-review workflow with explicit structural checks."
                        ),
                    }
                ],
            },
            "visual-surfaces": {
                "concept_tags": ["plotting", "figure-boundaries", "teaching-visuals"],
                "teaching_narrative": (
                    "Plotting-oriented blocks remain visible as course material while the "
                    "repo stays honest that rendered-figure equivalence is not yet claimed."
                ),
            },
            "simulation-inputs": {
                "concept_tags": ["simulation", "random-seeds", "signal-inputs"],
                "teaching_narrative": (
                    "The lecture's seeded random inputs are frozen so comparisons focus on "
                    "signal-fitting behavior rather than hidden simulation drift."
                ),
            },
            "comparative-signal": {
                "concept_tags": ["pagel-lambda", "likelihood-ratio", "signal-testing"],
                "teaching_narrative": (
                    "Students can follow how lambda estimation and lambda-zero testing are "
                    "checked numerically instead of being accepted from screenshots or prose."
                ),
                "migration_examples": [
                    {
                        "example_id": "pcm1-signal-fitting",
                        "family_id": "comparative-signal",
                        "fragment_id": "primate-lambda-fit",
                        "bijux_locators": [
                            "bijux_phylogenetics.comparative.signal:estimate_pagels_lambda",
                            "bijux_phylogenetics.validation:build_scientific_validation_report",
                        ],
                        "bijux_summary": (
                            "Bijux turns the lecture's lambda-fitting path into a tracked "
                            "parity bundle with explicit tolerance rules."
                        ),
                        "why_migrate_note": (
                            "The migration value is not shorter syntax; it is a governed "
                            "record of where the fitted values align and where boundaries remain."
                        ),
                    }
                ],
            },
            "ancestral-reconstruction": {
                "concept_tags": [
                    "ancestral-states",
                    "node-estimates",
                    "confidence-intervals",
                ],
                "teaching_narrative": (
                    "Ancestral-state outputs stay tied to explicit node estimates and "
                    "interval summaries so reviewers can inspect the exact numerical claim."
                ),
                "migration_examples": [
                    {
                        "example_id": "pcm1-ancestral-reconstruction",
                        "family_id": "ancestral-reconstruction",
                        "fragment_id": "continuous-ancestral-point-estimates",
                        "bijux_locators": [
                            "bijux_phylogenetics.validation:build_scientific_validation_report",
                            "bijux_phylogenetics.evidence.studies.primate_longevity_signal:build_primate_scalar_parity_table",
                        ],
                        "bijux_summary": (
                            "Ancestral-state checks are rendered as parity rows and "
                            "reviewable evidence rather than isolated notebook output."
                        ),
                        "why_migrate_note": (
                            "The Python-side evidence helps students compare node-level "
                            "claims without manually reading long console output."
                        ),
                    }
                ],
            },
            "artifact-provenance": {
                "concept_tags": ["processed-exports", "workspace-artifacts", "handoff"],
                "teaching_narrative": (
                    "Saved tables and trees are treated as governed handoff artifacts so "
                    "downstream lessons do not depend on undocumented files."
                ),
            },
        },
    },
    "primate-pgls-and-signal": {
        "course_material_label": "Lund PCM2 modes and PGLS lecture",
        "course_material_locator": "external:lund/pcm2-modes-pgls/script",
        "study_categories": ["teaching-study", "migration-study"],
        "student_safe_reproducibility": {
            "supported_scope": (
                "regenerate the governed PCM2 study outputs from the checked-in study "
                "sources and repository-managed reference script"
            ),
            "entrypoint_command": (
                "UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 "
                "python -m bijux_phylogenetics.command_line evidence book build "
                "--study-id primate-pgls-and-signal"
            ),
            "portable_prerequisites": [
                "Python 3.11",
                "UV project environment at artifacts/root/venv",
                "Rscript on PATH",
                "optional R_LIBS_USER rooted at artifacts/root/r-library",
            ],
            "forbidden_assumptions": [
                "no workstation-local /Users paths",
                "no sibling repository checkout assumptions",
                "no hidden manual patching of reference outputs",
            ],
            "expected_paths": [
                "evidence-book/studies/primate-pgls-and-signal/evidence-001/results/scalar-parity-table.json",
                "evidence-book/index/teaching-and-migration.json",
            ],
            "not_claimed": [
                "EB and ancestral-mode parity remains an explicit coverage boundary",
                "plot rendering is summarized as diagnostics rather than figure equivalence",
            ],
        },
        "families": {
            "workflow-contracts": {
                "concept_tags": [
                    "workspace-reload",
                    "course-context",
                    "reproducibility",
                ],
                "teaching_narrative": (
                    "The lecture's one-line reload assumption is made explicit so students "
                    "know which objects and paths the rest of the workflow depends on."
                ),
            },
            "baseline-regression": {
                "concept_tags": ["gls", "baseline-model", "regression"],
                "teaching_narrative": (
                    "Baseline GLS is isolated before phylogenetic covariance is introduced, "
                    "which makes the teaching sequence and the trust sequence line up."
                ),
                "migration_examples": [
                    {
                        "example_id": "pcm2-baseline-gls",
                        "family_id": "baseline-regression",
                        "fragment_id": "baseline-gls-fit",
                        "bijux_locators": [
                            "bijux_phylogenetics.comparative.pgls:run_pgls",
                            "bijux_phylogenetics.evidence.studies.primate_pgls_and_signal:build_primate_pgls_signal_bundles",
                        ],
                        "bijux_summary": (
                            "The baseline regression surface stays visible as a governed "
                            "comparison instead of disappearing once lambda-based models start."
                        ),
                        "why_migrate_note": (
                            "R users can see the exact non-phylogenetic baseline they already "
                            "know, then compare how Bijux layers explicit evidence on top."
                        ),
                    }
                ],
            },
            "phylogenetic-regression": {
                "concept_tags": [
                    "pgls",
                    "pagel-lambda",
                    "fixed-lambda",
                    "estimated-lambda",
                ],
                "teaching_narrative": (
                    "Fixed-lambda and estimated-lambda regression are kept in one evidence "
                    "family so the lecture transition from GLS to PGLS remains reviewable."
                ),
                "migration_examples": [
                    {
                        "example_id": "pcm2-pgls",
                        "family_id": "phylogenetic-regression",
                        "fragment_id": "pagel-lambda-regression",
                        "bijux_locators": [
                            "bijux_phylogenetics.comparative.pgls:run_pgls",
                            "bijux_phylogenetics.comparative.signal:estimate_pagels_lambda",
                        ],
                        "bijux_summary": (
                            "Bijux records the regression and lambda-fitting outputs with "
                            "explicit tolerance rules and verdicts."
                        ),
                        "why_migrate_note": (
                            "The migration value is that regression parity and tolerance logic "
                            "become explicit review artifacts instead of silent assumptions."
                        ),
                    }
                ],
            },
            "phylogenetic-signal": {
                "concept_tags": [
                    "phylogenetic-signal",
                    "lambda-zero",
                    "likelihood-ratio",
                ],
                "teaching_narrative": (
                    "Signal testing is treated as its own teaching family so the lecture's "
                    "intercept-only workflow does not get mixed into broader model-comparison claims."
                ),
                "migration_examples": [
                    {
                        "example_id": "pcm2-signal-test",
                        "family_id": "phylogenetic-signal",
                        "fragment_id": "phylogenetic-signal-test",
                        "bijux_locators": [
                            "bijux_phylogenetics.comparative.signal:compute_phylogenetic_signal_test",
                            "bijux_phylogenetics.comparative.signal:estimate_pagels_lambda",
                        ],
                        "bijux_summary": (
                            "The intercept-only signal workflow is carried into a governed "
                            "Python surface with explicit verdicts and tolerances."
                        ),
                        "why_migrate_note": (
                            "Students can compare the same signal decision path while seeing "
                            "which numerical differences remain bounded and reviewable."
                        ),
                    }
                ],
            },
            "diagnostics": {
                "concept_tags": [
                    "residuals",
                    "qq",
                    "heteroscedasticity",
                    "diagnostics",
                ],
                "teaching_narrative": (
                    "Residual and QQ diagnostics are converted into machine-recorded summary "
                    "values so the teaching surface is auditable even without figure matching."
                ),
                "migration_examples": [
                    {
                        "example_id": "pcm2-diagnostics",
                        "family_id": "diagnostics",
                        "fragment_id": "estimated-lambda-diagnostics",
                        "bijux_locators": [
                            "bijux_phylogenetics.evidence.studies.primate_pgls_and_signal:build_primate_pgls_signal_scalar_parity_table",
                            "bijux_phylogenetics.comparative.pgls:run_pgls",
                        ],
                        "bijux_summary": (
                            "Bijux exposes diagnostics as recorded scalar checks rather than "
                            "requiring manual plot comparison to understand parity."
                        ),
                        "why_migrate_note": (
                            "The migration benefit is reviewer-visible diagnostics that can be "
                            "versioned and compared, not just re-plotted."
                        ),
                    }
                ],
            },
            "coverage-boundaries": {
                "concept_tags": ["eb", "ou", "ancestral-gaps", "coverage-boundary"],
                "teaching_narrative": (
                    "Open EB and ancestral-mode gaps stay explicit so students are not taught "
                    "a false sense of completeness."
                ),
            },
        },
    },
}


def teaching_study_ids() -> set[str]:
    return set(TEACHING_STUDY_METADATA)


def study_metadata(study_id: str) -> dict[str, object]:
    return deepcopy(TEACHING_STUDY_METADATA[study_id])


def _fallback_family_metadata(
    family: dict[str, object],
    fragments_by_id: dict[str, dict[str, object]],
) -> dict[str, object]:
    family_id = str(family["family_id"])
    concept_tags = sorted(
        {
            *[
                str(fragments_by_id[fragment_id].get("concept_family"))
                for fragment_id in family.get("fragment_ids", [])
                if (
                    isinstance(fragment_id, str)
                    and isinstance(fragments_by_id.get(fragment_id), dict)
                    and isinstance(
                        fragments_by_id[fragment_id].get("concept_family"), str
                    )
                )
            ],
            *[token for token in family_id.split("-") if token],
        }
    )
    if not concept_tags:
        concept_tags = ["teaching-review-required"]
    return {
        "concept_tags": concept_tags,
        "teaching_narrative": (
            "This evidence family is indexed for teaching review, but its dedicated "
            "teaching narrative has not been curated yet. Review the bundle-level "
            "reports and source fragments directly until the family guide is strengthened."
        ),
    }


def build_teaching_guide(
    study_manifest: dict[str, object],
    family_index: dict[str, object],
    source_fragment_map: dict[str, object],
) -> dict[str, object]:
    metadata = study_metadata(str(study_manifest["study_id"]))
    fragments_by_id = {
        str(fragment["fragment_id"]): fragment
        for fragment in source_fragment_map.get("fragments", [])
        if isinstance(fragment, dict)
    }
    families = []
    concept_counter: Counter[str] = Counter()
    for family in family_index.get("families", []):
        if not isinstance(family, dict):
            continue
        family_id = str(family["family_id"])
        family_metadata = metadata["families"].get(family_id)
        if not isinstance(family_metadata, dict):
            family_metadata = _fallback_family_metadata(family, fragments_by_id)
        concept_tags = list(family_metadata["concept_tags"])
        concept_counter.update(concept_tags)
        family_fragments = [
            fragments_by_id[fragment_id]
            for fragment_id in family.get("fragment_ids", [])
            if fragment_id in fragments_by_id
        ]
        families.append(
            {
                "family_id": family_id,
                "family_title": family["family_title"],
                "family_verdict": family["family_verdict"],
                "coverage_status": family.get("coverage_status", "covered"),
                "concept_tags": concept_tags,
                "teaching_narrative": family_metadata["teaching_narrative"],
                "evidence_ids": family.get("evidence_ids", []),
                "fragment_ids": family.get("fragment_ids", []),
                "source_locators": [
                    locator
                    for fragment in family_fragments
                    for locator in fragment.get("script_locators", [])
                ],
            }
        )
    return {
        "schema_version": 1,
        "study_id": study_manifest["study_id"],
        "study_title": study_manifest["study_title"],
        "study_categories": study_manifest["study_categories"],
        "course_material_label": metadata["course_material_label"],
        "course_material_locator": metadata["course_material_locator"],
        "concept_tags": sorted(concept_counter),
        "family_count": len(families),
        "families": families,
    }


def render_teaching_guide_markdown(payload: dict[str, object]) -> str:
    lines = [
        f"# {payload['study_title']} Teaching Guide",
        "",
        f"- study id: `{payload['study_id']}`",
        f"- categories: `{', '.join(payload['study_categories'])}`",
        f"- course source: `{payload['course_material_locator']}`",
        f"- concept tags: `{', '.join(payload['concept_tags'])}`",
        "",
    ]
    for family in payload["families"]:
        lines.append(f"## {family['family_title']}")
        lines.append("")
        lines.append(f"- family id: `{family['family_id']}`")
        lines.append(f"- verdict: `{family['family_verdict']}`")
        lines.append(f"- coverage: `{family['coverage_status']}`")
        lines.append(f"- concept tags: `{', '.join(family['concept_tags'])}`")
        lines.append(f"- evidence ids: `{', '.join(family['evidence_ids'])}`")
        lines.append(f"- fragment ids: `{', '.join(family['fragment_ids'])}`")
        lines.append(f"- teaching narrative: {family['teaching_narrative']}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_migration_guide(
    study_manifest: dict[str, object],
    source_fragment_map: dict[str, object],
    bundle_manifests: list[dict[str, object]],
) -> dict[str, object]:
    metadata = study_metadata(str(study_manifest["study_id"]))
    comparison_mode_counts = Counter(
        str(manifest["comparison_mode"]) for manifest in bundle_manifests
    )
    fragment_map = {
        str(fragment["fragment_id"]): fragment
        for fragment in source_fragment_map.get("fragments", [])
        if isinstance(fragment, dict)
    }
    examples = []
    for family_id, family_metadata in metadata["families"].items():
        for example in family_metadata.get("migration_examples", []):
            fragment = fragment_map[str(example["fragment_id"])]
            examples.append(
                {
                    "example_id": example["example_id"],
                    "family_id": family_id,
                    "fragment_id": example["fragment_id"],
                    "r_fragment_title": fragment["fragment_title"],
                    "r_source_locators": fragment.get("script_locators", []),
                    "bijux_locators": example["bijux_locators"],
                    "bijux_summary": example["bijux_summary"],
                    "comparison_mode": "direct_parity",
                    "why_migrate_note": example["why_migrate_note"],
                }
            )
    return {
        "schema_version": 1,
        "study_id": study_manifest["study_id"],
        "study_title": study_manifest["study_title"],
        "course_material_locator": metadata["course_material_locator"],
        "comparison_mode_counts": dict(sorted(comparison_mode_counts.items())),
        "example_count": len(examples),
        "examples": examples,
    }


def render_migration_guide_markdown(payload: dict[str, object]) -> str:
    lines = [
        f"# {payload['study_title']} Migration Guide",
        "",
        f"- study id: `{payload['study_id']}`",
        f"- course source: `{payload['course_material_locator']}`",
        "",
        "## Comparison Modes",
        "",
    ]
    for mode, count in payload["comparison_mode_counts"].items():
        lines.append(f"- `{mode}`: `{count}`")
    lines.append("")
    lines.append("## Side-By-Side Examples")
    lines.append("")
    for example in payload["examples"]:
        lines.append(f"### {example['r_fragment_title']}")
        lines.append("")
        lines.append(f"- example id: `{example['example_id']}`")
        lines.append(f"- family id: `{example['family_id']}`")
        lines.append(f"- fragment id: `{example['fragment_id']}`")
        lines.append(f"- comparison mode: `{example['comparison_mode']}`")
        lines.append(f"- R locators: `{', '.join(example['r_source_locators'])}`")
        lines.append(f"- Bijux locators: `{', '.join(example['bijux_locators'])}`")
        lines.append(f"- Bijux summary: {example['bijux_summary']}")
        lines.append(f"- Why migrate: {example['why_migrate_note']}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_student_safe_reproducibility_contract(
    study_manifest: dict[str, object],
) -> dict[str, object]:
    metadata = study_metadata(str(study_manifest["study_id"]))
    contract = metadata["student_safe_reproducibility"]
    return {
        "schema_version": 1,
        "study_id": study_manifest["study_id"],
        "study_title": study_manifest["study_title"],
        "study_categories": study_manifest["study_categories"],
        "supported_scope": contract["supported_scope"],
        "entrypoint_command": contract["entrypoint_command"],
        "portable_prerequisites": contract["portable_prerequisites"],
        "forbidden_assumptions": contract["forbidden_assumptions"],
        "expected_paths": contract["expected_paths"],
        "not_claimed": contract["not_claimed"],
    }


def render_student_safe_reproducibility_markdown(payload: dict[str, object]) -> str:
    lines = [
        f"# {payload['study_title']} Student-Safe Reproducibility",
        "",
        f"- study id: `{payload['study_id']}`",
        f"- categories: `{', '.join(payload['study_categories'])}`",
        f"- supported scope: {payload['supported_scope']}",
        f"- entrypoint: `{payload['entrypoint_command']}`",
        "",
        "## Portable Prerequisites",
        "",
    ]
    for item in payload["portable_prerequisites"]:
        lines.append(f"- {item}")
    lines.append("")
    lines.append("## Forbidden Assumptions")
    lines.append("")
    for item in payload["forbidden_assumptions"]:
        lines.append(f"- {item}")
    lines.append("")
    lines.append("## Expected Paths")
    lines.append("")
    for item in payload["expected_paths"]:
        lines.append(f"- `{item}`")
    lines.append("")
    lines.append("## Not Claimed")
    lines.append("")
    for item in payload["not_claimed"]:
        lines.append(f"- {item}")
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_teaching_and_migration_index(
    teaching_guides: list[dict[str, object]],
    migration_guides: list[dict[str, object]],
    reproducibility_contracts: list[dict[str, object]],
) -> dict[str, object]:
    reproducibility_by_study = {
        str(contract["study_id"]): contract for contract in reproducibility_contracts
    }
    teaching_studies = []
    migration_studies = []
    concept_counter: Counter[str] = Counter()
    for guide in teaching_guides:
        concept_counter.update(str(tag) for tag in guide.get("concept_tags", []))
        teaching_studies.append(
            {
                "study_id": guide["study_id"],
                "study_title": guide["study_title"],
                "study_categories": guide["study_categories"],
                "course_material_locator": guide["course_material_locator"],
                "family_count": guide["family_count"],
                "concept_tags": guide["concept_tags"],
            }
        )
    for guide in migration_guides:
        migration_studies.append(
            {
                "study_id": guide["study_id"],
                "study_title": guide["study_title"],
                "comparison_mode_counts": guide["comparison_mode_counts"],
                "example_count": guide["example_count"],
                "supported_scope": reproducibility_by_study[str(guide["study_id"])][
                    "supported_scope"
                ],
            }
        )
    return {
        "schema_version": 1,
        "teaching_study_count": len(teaching_studies),
        "migration_study_count": len(migration_studies),
        "course_topic_counts": dict(sorted(concept_counter.items())),
        "teaching_studies": teaching_studies,
        "migration_studies": migration_studies,
    }


def render_teaching_and_migration_index_markdown(payload: dict[str, object]) -> str:
    lines = [
        "# Teaching And Migration Evidence",
        "",
        "This landing page separates Lund-derived teaching and migration studies",
        "from the broader scientific-validation surfaces in the evidence-book.",
        "",
        f"- teaching studies: `{payload['teaching_study_count']}`",
        f"- migration studies: `{payload['migration_study_count']}`",
        "",
    ]
    if payload["course_topic_counts"]:
        lines.append("## Course Topics")
        lines.append("")
        for topic, count in payload["course_topic_counts"].items():
            lines.append(f"- `{topic}`: `{count}`")
        lines.append("")
    lines.append("## Teaching Studies")
    lines.append("")
    for study in payload["teaching_studies"]:
        lines.append(f"- `{study['study_id']}` — {study['study_title']}")
        lines.append(f"  Course source: `{study['course_material_locator']}`")
        lines.append(f"  Families: `{study['family_count']}`")
        lines.append(f"  Concept tags: `{', '.join(study['concept_tags'])}`")
    lines.append("")
    lines.append("## Migration Studies")
    lines.append("")
    for study in payload["migration_studies"]:
        counts = ", ".join(
            f"{mode}={count}" for mode, count in study["comparison_mode_counts"].items()
        )
        lines.append(f"- `{study['study_id']}` — {study['study_title']}")
        lines.append(f"  Comparison modes: {counts}")
        lines.append(f"  Side-by-side examples: `{study['example_count']}`")
        lines.append(f"  Supported scope: {study['supported_scope']}")
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"
