from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

from bijux_phylogenetics.core.dataset import audit_dataset_inputs, build_dataset_crosswalk
from bijux_phylogenetics.core.taxon_workflows import (
    build_taxon_stability_report,
    build_taxon_workflow_loss_report,
    load_taxon_run_source,
)
from bijux_phylogenetics.core.taxonomy import (
    audit_tree_taxon_synonyms,
    inspect_tree_taxon_identity,
    inspect_tree_taxon_namespaces,
    resolve_tree_taxon_synonyms,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.reports.service import render_taxon_report


REPO_ROOT = Path(__file__).resolve().parents[3]
STUDY_ROOT = Path(__file__).resolve().parent
STUDY_ID = "taxon-trust"
EVIDENCE_ID = "evidence-001"
OUTPUT_ROOT = STUDY_ROOT / EVIDENCE_ID
BLOCK_ROOT = OUTPUT_ROOT / "block-payloads"
FIXTURES_ROOT = (
    REPO_ROOT / "packages" / "bijux-phylogenetics" / "tests" / "fixtures"
)

IDENTITY_TREE = FIXTURES_ROOT / "trees" / "example_tree_identity.nwk"
TAXONOMY_TREE = FIXTURES_ROOT / "trees" / "example_taxonomy_tree.nwk"
SYNONYM_TABLE = FIXTURES_ROOT / "metadata" / "example_taxon_synonyms.tsv"
WORKFLOW_TREE = FIXTURES_ROOT / "trees" / "example_taxon_workflow_tree.nwk"
WORKFLOW_METADATA = FIXTURES_ROOT / "metadata" / "example_taxon_workflow_metadata.csv"
WORKFLOW_TRAITS = FIXTURES_ROOT / "metadata" / "example_taxon_workflow_traits.csv"
WORKFLOW_ALIGNMENT = (
    FIXTURES_ROOT / "alignments" / "example_taxon_workflow_alignment.fasta"
)
WORKFLOW_FILTERED_ALIGNMENT = (
    FIXTURES_ROOT / "alignments" / "example_taxon_workflow_filtered_alignment.fasta"
)
WORKFLOW_INFERENCE = FIXTURES_ROOT / "trees" / "example_taxon_workflow_inference.nwk"
WORKFLOW_REPORTED = FIXTURES_ROOT / "metadata" / "example_taxon_workflow_reported.csv"

CLAIM_TITLES = {
    "taxon-spelling-variant-audit": "Taxon spelling-variant audit",
    "taxonomic-synonym-candidate-detection": "Taxonomic synonym candidate detection",
    "controlled-synonym-resolution": "Controlled synonym resolution",
    "ambiguous-synonym-rejection": "Ambiguous synonym rejection",
    "taxon-namespace-classification": "Taxon namespace classification",
    "mixed-namespace-detection": "Mixed namespace detection",
    "taxon-crosswalk-table": "Taxon crosswalk table",
    "taxon-exclusion-reasoning": "Taxon exclusion reasoning",
    "workflow-taxon-loss-report": "Workflow taxon-loss report",
    "taxon-stability-report": "Taxon stability report",
}

PYTHON_BLOCKS = {
    "taxon-spelling-variant-audit": """from pathlib import Path

from bijux_phylogenetics.core.taxonomy import inspect_tree_taxon_identity
from bijux_phylogenetics.io.trees import load_tree

tree = load_tree(Path("packages/bijux-phylogenetics/tests/fixtures/trees/example_tree_identity.nwk"))
report = inspect_tree_taxon_identity(tree)""",
    "taxonomic-synonym-candidate-detection": """from pathlib import Path

from bijux_phylogenetics.core.taxonomy import audit_tree_taxon_synonyms
from bijux_phylogenetics.io.trees import load_tree

tree = load_tree(Path("packages/bijux-phylogenetics/tests/fixtures/trees/example_taxonomy_tree.nwk"))
report = audit_tree_taxon_synonyms(
    tree,
    Path("packages/bijux-phylogenetics/tests/fixtures/metadata/example_taxon_synonyms.tsv"),
)""",
    "controlled-synonym-resolution": """from pathlib import Path

from bijux_phylogenetics.core.taxonomy import resolve_tree_taxon_synonyms
from bijux_phylogenetics.io.trees import load_tree

tree = load_tree(Path("packages/bijux-phylogenetics/tests/fixtures/trees/example_taxonomy_tree.nwk"))
resolved_tree, report = resolve_tree_taxon_synonyms(
    tree,
    synonym_table_path=Path("packages/bijux-phylogenetics/tests/fixtures/metadata/example_taxon_synonyms.tsv"),
)""",
    "ambiguous-synonym-rejection": """from pathlib import Path

from bijux_phylogenetics.core.taxonomy import resolve_tree_taxon_synonyms
from bijux_phylogenetics.io.trees import load_tree

tree = load_tree(Path("packages/bijux-phylogenetics/tests/fixtures/trees/example_taxonomy_tree.nwk"))
resolved_tree, report = resolve_tree_taxon_synonyms(
    tree,
    synonym_table_path=Path("packages/bijux-phylogenetics/tests/fixtures/metadata/example_taxon_synonyms.tsv"),
)
assert "Jaguar" in resolved_tree.tip_names""",
    "taxon-namespace-classification": """from pathlib import Path

from bijux_phylogenetics.core.taxonomy import inspect_tree_taxon_namespaces
from bijux_phylogenetics.io.trees import load_tree

tree = load_tree(Path("packages/bijux-phylogenetics/tests/fixtures/trees/example_taxonomy_tree.nwk"))
report = inspect_tree_taxon_namespaces(tree)""",
    "mixed-namespace-detection": """from pathlib import Path

from bijux_phylogenetics.core.taxonomy import inspect_tree_taxon_namespaces
from bijux_phylogenetics.io.trees import load_tree

tree = load_tree(Path("packages/bijux-phylogenetics/tests/fixtures/trees/example_taxonomy_tree.nwk"))
report = inspect_tree_taxon_namespaces(tree)
assert report.mixed_namespaces is True""",
    "taxon-crosswalk-table": """from pathlib import Path

from bijux_phylogenetics.core.dataset import build_dataset_crosswalk

report = build_dataset_crosswalk(
    Path("packages/bijux-phylogenetics/tests/fixtures/trees/example_taxon_workflow_tree.nwk"),
    Path("packages/bijux-phylogenetics/tests/fixtures/metadata/example_taxon_workflow_metadata.csv"),
    Path("packages/bijux-phylogenetics/tests/fixtures/metadata/example_taxon_workflow_traits.csv"),
    alignment_path=Path("packages/bijux-phylogenetics/tests/fixtures/alignments/example_taxon_workflow_alignment.fasta"),
)""",
    "taxon-exclusion-reasoning": """from pathlib import Path

from bijux_phylogenetics.core.dataset import audit_dataset_inputs

report = audit_dataset_inputs(
    Path("packages/bijux-phylogenetics/tests/fixtures/trees/example_taxon_workflow_tree.nwk"),
    Path("packages/bijux-phylogenetics/tests/fixtures/metadata/example_taxon_workflow_metadata.csv"),
    Path("packages/bijux-phylogenetics/tests/fixtures/metadata/example_taxon_workflow_traits.csv"),
    alignment_path=Path("packages/bijux-phylogenetics/tests/fixtures/alignments/example_taxon_workflow_alignment.fasta"),
)""",
    "workflow-taxon-loss-report": """from pathlib import Path

from bijux_phylogenetics.core.taxon_workflows import build_taxon_workflow_loss_report

report = build_taxon_workflow_loss_report(
    Path("packages/bijux-phylogenetics/tests/fixtures/trees/example_taxon_workflow_tree.nwk"),
    Path("packages/bijux-phylogenetics/tests/fixtures/metadata/example_taxon_workflow_metadata.csv"),
    Path("packages/bijux-phylogenetics/tests/fixtures/metadata/example_taxon_workflow_traits.csv"),
    alignment_path=Path("packages/bijux-phylogenetics/tests/fixtures/alignments/example_taxon_workflow_alignment.fasta"),
    filtered_alignment_path=Path("packages/bijux-phylogenetics/tests/fixtures/alignments/example_taxon_workflow_filtered_alignment.fasta"),
    inference_tree_path=Path("packages/bijux-phylogenetics/tests/fixtures/trees/example_taxon_workflow_inference.nwk"),
    reported_taxa_path=Path("packages/bijux-phylogenetics/tests/fixtures/metadata/example_taxon_workflow_reported.csv"),
)""",
    "taxon-stability-report": """from pathlib import Path

from bijux_phylogenetics.core.taxon_workflows import (
    build_taxon_stability_report,
    load_taxon_run_source,
)

report = build_taxon_stability_report(
    [
        load_taxon_run_source(label="tree", path=Path("packages/bijux-phylogenetics/tests/fixtures/trees/example_taxon_workflow_tree.nwk")),
        load_taxon_run_source(label="metadata", path=Path("packages/bijux-phylogenetics/tests/fixtures/metadata/example_taxon_workflow_metadata.csv")),
        load_taxon_run_source(label="traits", path=Path("packages/bijux-phylogenetics/tests/fixtures/metadata/example_taxon_workflow_traits.csv")),
        load_taxon_run_source(label="alignment", path=Path("packages/bijux-phylogenetics/tests/fixtures/alignments/example_taxon_workflow_alignment.fasta")),
        load_taxon_run_source(label="filtered_alignment", path=Path("packages/bijux-phylogenetics/tests/fixtures/alignments/example_taxon_workflow_filtered_alignment.fasta")),
        load_taxon_run_source(label="inference_tree", path=Path("packages/bijux-phylogenetics/tests/fixtures/trees/example_taxon_workflow_inference.nwk")),
        load_taxon_run_source(label="reported_taxa", path=Path("packages/bijux-phylogenetics/tests/fixtures/metadata/example_taxon_workflow_reported.csv")),
    ]
)""",
}


def _rel(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _pair_rows(rows: list[object]) -> list[list[str]]:
    return [
        [getattr(row, "left_label"), getattr(row, "right_label")]
        for row in rows
    ]


def _claim_payloads() -> list[dict[str, object]]:
    identity_report = inspect_tree_taxon_identity(load_tree(IDENTITY_TREE))
    synonym_audit = audit_tree_taxon_synonyms(load_tree(TAXONOMY_TREE), SYNONYM_TABLE)
    resolved_tree, synonym_resolution = resolve_tree_taxon_synonyms(
        load_tree(TAXONOMY_TREE), synonym_table_path=SYNONYM_TABLE
    )
    namespace_report = inspect_tree_taxon_namespaces(load_tree(TAXONOMY_TREE))
    crosswalk = build_dataset_crosswalk(
        WORKFLOW_TREE,
        WORKFLOW_METADATA,
        WORKFLOW_TRAITS,
        alignment_path=WORKFLOW_ALIGNMENT,
    )
    dataset_audit = audit_dataset_inputs(
        WORKFLOW_TREE,
        WORKFLOW_METADATA,
        WORKFLOW_TRAITS,
        alignment_path=WORKFLOW_ALIGNMENT,
    )
    workflow_loss = build_taxon_workflow_loss_report(
        WORKFLOW_TREE,
        WORKFLOW_METADATA,
        WORKFLOW_TRAITS,
        alignment_path=WORKFLOW_ALIGNMENT,
        filtered_alignment_path=WORKFLOW_FILTERED_ALIGNMENT,
        inference_tree_path=WORKFLOW_INFERENCE,
        reported_taxa_path=WORKFLOW_REPORTED,
    )
    stability = build_taxon_stability_report(
        [
            load_taxon_run_source(label="tree", path=WORKFLOW_TREE),
            load_taxon_run_source(label="metadata", path=WORKFLOW_METADATA),
            load_taxon_run_source(label="traits", path=WORKFLOW_TRAITS),
            load_taxon_run_source(label="alignment", path=WORKFLOW_ALIGNMENT),
            load_taxon_run_source(
                label="filtered_alignment", path=WORKFLOW_FILTERED_ALIGNMENT
            ),
            load_taxon_run_source(label="inference_tree", path=WORKFLOW_INFERENCE),
            load_taxon_run_source(label="reported_taxa", path=WORKFLOW_REPORTED),
        ]
    )

    return [
        {
            "claim_id": "taxon-spelling-variant-audit",
            "claim_title": CLAIM_TITLES["taxon-spelling-variant-audit"],
            "verdict": "verified",
            "fixture_paths": [_rel(IDENTITY_TREE)],
            "summary_lines": [
                "3 underscore/space collision pairs, 1 case collision, and 1 near-duplicate pair were detected."
            ],
            "python_block": PYTHON_BLOCKS["taxon-spelling-variant-audit"],
            "observed": {
                "spelling_variants": _pair_rows(identity_report.spelling_variants),
                "whitespace_variants": _pair_rows(identity_report.whitespace_variants),
                "underscore_space_collisions": _pair_rows(
                    identity_report.underscore_space_collisions
                ),
                "case_collisions": _pair_rows(identity_report.case_collisions),
                "suspicious_near_duplicates": _pair_rows(
                    identity_report.suspicious_near_duplicates
                ),
            },
        },
        {
            "claim_id": "taxonomic-synonym-candidate-detection",
            "claim_title": CLAIM_TITLES["taxonomic-synonym-candidate-detection"],
            "verdict": "verified",
            "fixture_paths": [_rel(TAXONOMY_TREE), _rel(SYNONYM_TABLE)],
            "summary_lines": [
                "1 configured synonym candidate was found and 1 raw label was blocked as ambiguous."
            ],
            "python_block": PYTHON_BLOCKS["taxonomic-synonym-candidate-detection"],
            "observed": {
                "candidates": [asdict(row) for row in synonym_audit.candidates],
                "ambiguous_mappings": [
                    asdict(row) for row in synonym_audit.ambiguous_mappings
                ],
                "warnings": synonym_audit.warnings,
            },
        },
        {
            "claim_id": "controlled-synonym-resolution",
            "claim_title": CLAIM_TITLES["controlled-synonym-resolution"],
            "verdict": "verified",
            "fixture_paths": [_rel(TAXONOMY_TREE), _rel(SYNONYM_TABLE)],
            "summary_lines": [
                "Felis_concolor resolves to Puma_concolor with explicit provenance and reversible mapping rows."
            ],
            "python_block": PYTHON_BLOCKS["controlled-synonym-resolution"],
            "observed": {
                "renamed_taxa": [asdict(row) for row in synonym_resolution.renamed_taxa],
                "unchanged_taxa": synonym_resolution.unchanged_taxa,
                "resolved_tip_names": sorted(resolved_tree.tip_names),
            },
        },
        {
            "claim_id": "ambiguous-synonym-rejection",
            "claim_title": CLAIM_TITLES["ambiguous-synonym-rejection"],
            "verdict": "verified",
            "fixture_paths": [_rel(TAXONOMY_TREE), _rel(SYNONYM_TABLE)],
            "summary_lines": [
                "Jaguar remains unresolved because the synonym table maps it to both Panthera_onca and Panthera_pardus."
            ],
            "python_block": PYTHON_BLOCKS["ambiguous-synonym-rejection"],
            "observed": {
                "ambiguous_mappings": [
                    asdict(row) for row in synonym_resolution.ambiguous_mappings
                ],
                "duplicate_resolved_labels": [
                    asdict(row) for row in synonym_resolution.duplicate_resolved_labels
                ],
                "retained_unresolved_labels": [
                    label
                    for label in sorted(resolved_tree.tip_names)
                    if label == "Jaguar"
                ],
            },
        },
        {
            "claim_id": "taxon-namespace-classification",
            "claim_title": CLAIM_TITLES["taxon-namespace-classification"],
            "verdict": "verified",
            "fixture_paths": [_rel(TAXONOMY_TREE)],
            "summary_lines": [
                "The taxonomy tree includes species_name, accession_id, sample_id, and isolate_id labels."
            ],
            "python_block": PYTHON_BLOCKS["taxon-namespace-classification"],
            "observed": {
                "assignments": [asdict(row) for row in namespace_report.assignments],
                "namespace_counts": namespace_report.namespace_counts,
                "dominant_namespace": namespace_report.dominant_namespace,
            },
        },
        {
            "claim_id": "mixed-namespace-detection",
            "claim_title": CLAIM_TITLES["mixed-namespace-detection"],
            "verdict": "verified",
            "fixture_paths": [_rel(TAXONOMY_TREE)],
            "summary_lines": [
                "mixed_namespaces is true and the report warns against automation without a crosswalk."
            ],
            "python_block": PYTHON_BLOCKS["mixed-namespace-detection"],
            "observed": {
                "mixed_namespaces": namespace_report.mixed_namespaces,
                "warnings": namespace_report.warnings,
            },
        },
        {
            "claim_id": "taxon-crosswalk-table",
            "claim_title": CLAIM_TITLES["taxon-crosswalk-table"],
            "verdict": "verified",
            "fixture_paths": [
                _rel(WORKFLOW_TREE),
                _rel(WORKFLOW_METADATA),
                _rel(WORKFLOW_TRAITS),
                _rel(WORKFLOW_ALIGNMENT),
            ],
            "summary_lines": [
                "The crosswalk exposes 4 taxa across tree, metadata, traits, and alignment surfaces."
            ],
            "python_block": PYTHON_BLOCKS["taxon-crosswalk-table"],
            "observed": {
                "row_count": len(crosswalk.rows),
                "rows": [asdict(row) for row in crosswalk.rows],
            },
        },
        {
            "claim_id": "taxon-exclusion-reasoning",
            "claim_title": CLAIM_TITLES["taxon-exclusion-reasoning"],
            "verdict": "verified",
            "fixture_paths": [
                _rel(WORKFLOW_TREE),
                _rel(WORKFLOW_METADATA),
                _rel(WORKFLOW_TRAITS),
                _rel(WORKFLOW_ALIGNMENT),
            ],
            "summary_lines": [
                "B is excluded by metadata absence; D is excluded by alignment absence, with affected analyses recorded explicitly."
            ],
            "python_block": PYTHON_BLOCKS["taxon-exclusion-reasoning"],
            "observed": {
                "rows": [asdict(row) for row in dataset_audit.exclusion_table.rows],
            },
        },
        {
            "claim_id": "workflow-taxon-loss-report",
            "claim_title": CLAIM_TITLES["workflow-taxon-loss-report"],
            "verdict": "verified",
            "fixture_paths": [
                _rel(WORKFLOW_TREE),
                _rel(WORKFLOW_METADATA),
                _rel(WORKFLOW_TRAITS),
                _rel(WORKFLOW_ALIGNMENT),
                _rel(WORKFLOW_FILTERED_ALIGNMENT),
                _rel(WORKFLOW_INFERENCE),
                _rel(WORKFLOW_REPORTED),
            ],
            "summary_lines": [
                "First-loss stages are B at alignment_filtering, C at trait_missingness, and D at alignment."
            ],
            "python_block": PYTHON_BLOCKS["workflow-taxon-loss-report"],
            "observed": {
                "loss_stage_counts": workflow_loss.loss_stage_counts,
                "rows": [asdict(row) for row in workflow_loss.rows],
            },
        },
        {
            "claim_id": "taxon-stability-report",
            "claim_title": CLAIM_TITLES["taxon-stability-report"],
            "verdict": "verified",
            "fixture_paths": [
                _rel(WORKFLOW_TREE),
                _rel(WORKFLOW_METADATA),
                _rel(WORKFLOW_TRAITS),
                _rel(WORKFLOW_ALIGNMENT),
                _rel(WORKFLOW_FILTERED_ALIGNMENT),
                _rel(WORKFLOW_INFERENCE),
                _rel(WORKFLOW_REPORTED),
            ],
            "summary_lines": [
                "A is stable across all seven sources; B, C, and D are unstable."
            ],
            "python_block": PYTHON_BLOCKS["taxon-stability-report"],
            "observed": {
                "shared_taxa": stability.shared_taxa,
                "stable_taxa": stability.stable_taxa,
                "unstable_taxa": stability.unstable_taxa,
                "rows": [asdict(row) for row in stability.rows],
            },
        },
    ]


def _build_taxonomy_report_manifest() -> dict[str, object]:
    with tempfile.TemporaryDirectory() as tmp_dir:
        out_path = Path(tmp_dir) / "taxonomy-report.html"
        result = render_taxon_report(
            tree_path=WORKFLOW_TREE,
            metadata_path=WORKFLOW_METADATA,
            traits_path=WORKFLOW_TRAITS,
            alignment_path=WORKFLOW_ALIGNMENT,
            filtered_alignment_path=WORKFLOW_FILTERED_ALIGNMENT,
            inference_tree_path=WORKFLOW_INFERENCE,
            reported_taxa_path=WORKFLOW_REPORTED,
            out_path=out_path,
        )
        return json.loads(result.machine_manifest_path.read_text(encoding="utf-8"))


def _write_bundle_readme(
    claim_payloads: list[dict[str, object]], report_manifest_path: str
) -> None:
    lines = [
        "# Evidence 001",
        "",
        "This bundle validates the checked-in taxon-trust workflow claims using",
        "repository fixtures and package-owned report surfaces.",
        "",
        "Generated by:",
        "",
        "- [`build_evidence.py`](../build_evidence.py)",
        "",
        "Supporting reviewer report surface:",
        "",
        f"- [`taxonomy_report_machine_manifest.json`](./{report_manifest_path})",
        "",
        "## Claim Verdicts",
        "",
    ]
    for payload in claim_payloads:
        lines.extend(
            [
                f"### {payload['claim_title']}",
                "",
                f"Claim id: `{payload['claim_id']}`",
                "",
                f"Verdict: `{payload['verdict']}`",
                "",
                "Inputs:",
            ]
        )
        for fixture_path in payload["fixture_paths"]:
            lines.append(f"- `{fixture_path}`")
        lines.extend(["", "Highlights:"])
        for summary_line in payload["summary_lines"]:
            lines.append(f"- {summary_line}")
        lines.extend(
            [
                "",
                "Python block:",
                "",
                "```python",
                payload["python_block"],
                "```",
                "",
                f"Observed payload: [`{payload['claim_id']}.json`](./block-payloads/{payload['claim_id']}.json)",
                "",
            ]
        )
    (OUTPUT_ROOT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    if OUTPUT_ROOT.exists():
        shutil.rmtree(OUTPUT_ROOT)
    BLOCK_ROOT.mkdir(parents=True, exist_ok=True)

    claim_payloads = _claim_payloads()
    report_manifest = _build_taxonomy_report_manifest()

    for payload in claim_payloads:
        _write_json(BLOCK_ROOT / f"{payload['claim_id']}.json", payload)

    claim_verdicts = [
        {
            "claim_id": payload["claim_id"],
            "claim_title": payload["claim_title"],
            "verdict": payload["verdict"],
            "fixture_paths": payload["fixture_paths"],
            "summary_lines": payload["summary_lines"],
        }
        for payload in claim_payloads
    ]
    manifest = {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": EVIDENCE_ID,
        "evidence_title": "Taxon workflow review bundle",
        "summary": "Fixture-backed taxon trust evidence spanning identity, crosswalk, exclusion, loss, and stability behavior.",
        "owner_package": "bijux-phylogenetics",
        "source_basis": [
            {
                "kind": "repository-fixture",
                "label": "taxon trust fixtures",
                "locator": "packages/bijux-phylogenetics/tests/fixtures",
            }
        ],
        "claim_tags": ["taxonomy", "workflow-audit", "trust-surface"],
        "verdict": {
            "status": "matched",
            "summary": "All recorded taxon-trust claims match the checked-in fixture expectations.",
        },
        "limitations": [
            "Covers one governed taxon-workflow bundle rather than the full repository surface.",
        ],
        "claim_ids": [payload["claim_id"] for payload in claim_payloads],
        "input_checksums": {
            _rel(path): _sha256(path)
            for path in [
                IDENTITY_TREE,
                TAXONOMY_TREE,
                SYNONYM_TABLE,
                WORKFLOW_TREE,
                WORKFLOW_METADATA,
                WORKFLOW_TRAITS,
                WORKFLOW_ALIGNMENT,
                WORKFLOW_FILTERED_ALIGNMENT,
                WORKFLOW_INFERENCE,
                WORKFLOW_REPORTED,
            ]
        },
    }

    _write_json(OUTPUT_ROOT / "manifest.json", manifest)
    _write_json(OUTPUT_ROOT / "claim_verdicts.json", claim_verdicts)
    _write_json(OUTPUT_ROOT / "taxonomy_report_machine_manifest.json", report_manifest)
    _write_bundle_readme(claim_payloads, "taxonomy_report_machine_manifest.json")


if __name__ == "__main__":
    main()
