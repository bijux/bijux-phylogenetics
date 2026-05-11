from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from .primate_longevity_signal import STUDY_ID

REFERENCE_SCRIPT_PATH = (
    "evidence-book/studies/primate-longevity-signal/reference/"
    "primate_lifespan_signal_reference_r.R"
)
REFERENCE_BUNDLE_ID = "evidence-001"
REFERENCE_BUNDLE_ROOT = (
    Path("evidence-book") / "studies" / "primate-longevity-signal" / REFERENCE_BUNDLE_ID
)
COMPONENT_BUNDLE_DEFINITIONS = [
    {
        "evidence_id": "evidence-002",
        "report_filename": "workbook-loading-parity.json",
        "title": "Primate workbook loading parity bundle",
        "summary": "Governed evidence for the raw workbook intake boundary and the resulting processed table contract.",
        "claim_id": "pcm1-workbook-loading-parity",
        "claim_title": "PCM1 workbook loading contract is governed explicitly",
        "claim_summary": "The lecture workbook boundary is represented explicitly and the resulting processed table contract matches the checked-in governed reference artifact.",
        "verdict": "matched",
        "claim_tags": ["teaching", "parity", "data-preparation", "workbook-loading"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["data-preparation"],
        "limitations": [
            "The raw workbook remains an external teaching source; the Python-side review starts from the governed processed reference artifact."
        ],
        "source_fragments": ["primate-data-preprocessing"],
        "reference_line_specs": ["33-34", "77"],
        "source_basis": [
            {
                "kind": "external-course-material",
                "label": "lund pcm1 raw primate workbook",
                "locator": "external:lund/pcm1-plots-signal/data/primate_raw.xlsx",
            },
            {
                "kind": "repository-reference",
                "label": "governed processed primate table",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv",
            },
        ],
    },
    {
        "evidence_id": "evidence-003",
        "report_filename": "type-repair-parity.json",
        "title": "Primate type repair parity bundle",
        "summary": "Governed evidence for the factor and numeric coercion contract before grouped downstream analysis.",
        "claim_id": "pcm1-type-repair-parity",
        "claim_title": "PCM1 type repair contract is explicit and reviewable",
        "claim_summary": "The reference factor and numeric coercion rules are made explicit and the governed processed artifact satisfies the repaired column contract.",
        "verdict": "matched",
        "claim_tags": ["teaching", "parity", "data-preparation", "type-repair"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["data-preparation"],
        "limitations": [
            "The Python-side review validates the repaired processed artifact rather than replaying the external workbook import directly."
        ],
        "source_fragments": ["primate-data-preprocessing"],
        "reference_line_specs": ["39-40", "59-71"],
        "source_basis": [
            {
                "kind": "repository-reference",
                "label": "governed processed primate table",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv",
            }
        ],
    },
    {
        "evidence_id": "evidence-004",
        "report_filename": "missing-data-accounting-parity.json",
        "title": "Primate missing-data accounting parity bundle",
        "summary": "Governed evidence for how missing values are handled before the processed primate table is used downstream.",
        "claim_id": "pcm1-missing-data-accounting-parity",
        "claim_title": "PCM1 missing-data accounting is visible before downstream inference",
        "claim_summary": "Missing-value handling is surfaced explicitly at the preprocessing boundary instead of being hidden beneath later comparative analyses.",
        "verdict": "matched",
        "claim_tags": ["teaching", "parity", "data-preparation", "missing-data"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["data-preparation"],
        "limitations": [
            "The raw workbook remains external, so the bundle records the governed missing-value contract at the processed artifact boundary."
        ],
        "source_fragments": ["primate-data-preprocessing"],
        "reference_line_specs": ["41-47"],
        "source_basis": [
            {
                "kind": "repository-reference",
                "label": "governed processed primate table",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv",
            }
        ],
    },
    {
        "evidence_id": "evidence-005",
        "report_filename": "duplicate-species-aggregation-parity.json",
        "title": "Primate duplicate-species aggregation parity bundle",
        "summary": "Governed evidence for duplicate-species consolidation before processed trait analysis.",
        "claim_id": "pcm1-duplicate-species-aggregation-parity",
        "claim_title": "PCM1 duplicate-species aggregation is explicit and quantified",
        "claim_summary": "The raw-to-processed row contraction and the resulting unique-species contract are broken out as standalone evidence instead of being buried inside later model checks.",
        "verdict": "matched",
        "claim_tags": [
            "teaching",
            "parity",
            "data-preparation",
            "duplicate-aggregation",
        ],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["data-preparation"],
        "limitations": [
            "The raw workbook is external, so duplicate consolidation is demonstrated through the governed row-count and species-count contract."
        ],
        "source_fragments": ["primate-data-preprocessing"],
        "reference_line_specs": ["41-47"],
        "source_basis": [
            {
                "kind": "repository-reference",
                "label": "governed processed primate table",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv",
            }
        ],
    },
    {
        "evidence_id": "evidence-006",
        "report_filename": "tree-import-parity.json",
        "title": "Primate tree import parity bundle",
        "summary": "Governed evidence for reading the primate tree and preserving the expected original and trimmed tree contract.",
        "claim_id": "pcm1-tree-import-parity",
        "claim_title": "PCM1 tree import behavior is governed explicitly",
        "claim_summary": "Original and trimmed tree imports are surfaced as their own evidence boundary with explicit tip-count and branch-length expectations.",
        "verdict": "matched",
        "claim_tags": ["teaching", "parity", "tree-operations", "tree-import"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["tree-preprocessing"],
        "limitations": [
            "This bundle focuses on import and preserved tree shape; downstream correspondence and diagnostics are broken out separately."
        ],
        "source_fragments": ["tree-import-and-pruning"],
        "reference_line_specs": ["79-85"],
        "source_basis": [
            {
                "kind": "external-course-material",
                "label": "lund pcm1 original primate tree",
                "locator": "external:lund/pcm1-plots-signal/data/primatetree.nex",
            },
            {
                "kind": "repository-reference",
                "label": "governed trimmed tree export",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_trimmed_primatetree.nwk",
            },
        ],
    },
    {
        "evidence_id": "evidence-007",
        "report_filename": "tree-diagnostics-parity.json",
        "title": "Primate tree diagnostics parity bundle",
        "summary": "Governed evidence for rootedness, binary-tree status, ultrametricity, and node-label interpretation.",
        "claim_id": "pcm1-tree-diagnostics-parity",
        "claim_title": "PCM1 tree diagnostics are separated from later comparative results",
        "claim_summary": "Tree rootedness, binary structure, ultrametricity, and node-label interpretation are recorded as standalone trust surfaces.",
        "verdict": "matched",
        "claim_tags": ["teaching", "parity", "tree-operations", "tree-diagnostics"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["tree-validation"],
        "limitations": [
            "This bundle records the checked-in diagnostic contract; visual tree surfaces remain intentionally separate."
        ],
        "source_fragments": ["tree-import-and-pruning", "extract-clade-node-77"],
        "reference_line_specs": ["80-83", "87-95", "109"],
        "source_basis": [
            {
                "kind": "repository-reference",
                "label": "governed tree-processing comparison payload",
                "locator": "evidence-book/studies/primate-longevity-signal/evidence-001/results/block-payloads/tree-import-and-pruning.json",
            }
        ],
    },
    {
        "evidence_id": "evidence-008",
        "report_filename": "tree-data-correspondence-parity.json",
        "title": "Primate tree-data correspondence parity bundle",
        "summary": "Governed evidence for missing-tip trimming, tip-order alignment, and node-aware tree-data joining.",
        "claim_id": "pcm1-tree-data-correspondence-parity",
        "claim_title": "PCM1 tree-data correspondence is explicit before model fitting",
        "claim_summary": "Missing-tip handling, aligned species ordering, and joined node identity are isolated as prerequisite trust surfaces before downstream inference.",
        "verdict": "matched",
        "claim_tags": [
            "teaching",
            "parity",
            "data-preparation",
            "tree-data-correspondence",
        ],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["tree-data-alignment"],
        "limitations": [
            "This bundle covers correspondence and join readiness only; later lambda and ancestral-state results are separate evidence."
        ],
        "source_fragments": [
            "tree-import-and-pruning",
            "tip-order-alignment",
            "treeio-node-mapping-and-join",
        ],
        "reference_line_specs": ["81-83", "99-107"],
        "source_basis": [
            {
                "kind": "repository-reference",
                "label": "governed processed primate table",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv",
            },
            {
                "kind": "repository-reference",
                "label": "governed trimmed tree export",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_trimmed_primatetree.nwk",
            },
        ],
    },
    {
        "evidence_id": "evidence-009",
        "report_filename": "processed-export-parity.json",
        "title": "Primate processed export parity bundle",
        "summary": "Governed evidence for the exported processed CSV and trimmed tree artifacts consumed by later teaching steps.",
        "claim_id": "pcm1-processed-export-parity",
        "claim_title": "PCM1 processed exports are governed as explicit downstream inputs",
        "claim_summary": "The processed CSV and trimmed tree exports are treated as first-class governed artifacts rather than hidden scaffolding for later parity claims.",
        "verdict": "matched",
        "claim_tags": ["teaching", "parity", "artifact-provenance", "processed-export"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["artifact-provenance"],
        "limitations": [
            "The governed export contract records the downstream artifacts and their parity expectations, not the rendered plot outputs that consume them."
        ],
        "source_fragments": ["processed-analysis-artifacts"],
        "reference_line_specs": ["77", "85"],
        "source_basis": [
            {
                "kind": "repository-reference",
                "label": "governed processed primate table",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv",
            },
            {
                "kind": "repository-reference",
                "label": "governed trimmed tree export",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_trimmed_primatetree.nwk",
            },
        ],
    },
]

NUMERIC_COLUMNS = [
    "body_mass",
    "gestation",
    "home_range",
    "longevity",
    "social_group_size",
]
CATEGORICAL_COLUMNS = ["sex_dimorphism", "mating_system"]
TEXT_COLUMNS = ["family", "species"]
DATA_PREPARATION_BUNDLE_IDS = [
    "evidence-002",
    "evidence-003",
    "evidence-004",
    "evidence-005",
    "evidence-008",
]


def _study_root(repo_root: Path) -> Path:
    return Path(repo_root) / "evidence-book" / "studies" / "primate-longevity-signal"


def _base_bundle_root(repo_root: Path) -> Path:
    return _study_root(repo_root) / REFERENCE_BUNDLE_ID


def _read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _reference_script_locators(line_specs: list[str]) -> list[str]:
    locators: list[str] = []
    for spec in line_specs:
        for part in spec.split(","):
            normalized = part.strip()
            if not normalized:
                continue
            if "-" in normalized:
                start_text, end_text = normalized.split("-", maxsplit=1)
                locators.append(
                    f"{REFERENCE_SCRIPT_PATH}#L{int(start_text)}-L{int(end_text)}"
                )
            else:
                locators.append(f"{REFERENCE_SCRIPT_PATH}#L{int(normalized)}")
    return locators


def _load_context(repo_root: Path) -> dict[str, object]:
    bundle_root = _base_bundle_root(repo_root)
    return {
        "r_results": _read_json(bundle_root / "results" / "r_reference_results.json"),
        "bijux_results": _read_json(
            bundle_root / "results" / "bijux_reference_results.json"
        ),
        "reference_rows": _read_csv_rows(
            _study_root(repo_root) / "datasets" / "reference_primate.csv"
        ),
        "block_payloads": {
            path.stem: _read_json(path)
            for path in sorted(
                (bundle_root / "results" / "block-payloads").glob("*.json")
            )
        },
    }


def _missing_counts(rows: list[dict[str, str]]) -> dict[str, int]:
    counts = dict.fromkeys(rows[0].keys(), 0)
    for row in rows:
        for column, value in row.items():
            if value == "":
                counts[column] += 1
    return counts


def _report_payload_for_spec(
    spec: dict[str, object], context: dict[str, object]
) -> dict[str, object]:
    r_results = context["r_results"]
    bijux_results = context["bijux_results"]
    reference_rows = context["reference_rows"]
    block_payloads = context["block_payloads"]
    evidence_id = spec["evidence_id"]

    if evidence_id == "evidence-002":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "raw_workbook_locator": "external:lund/pcm1-plots-signal/data/primate_raw.xlsx",
            "governed_processed_reference_locator": (
                "evidence-book/studies/primate-longevity-signal/evidence-001/"
                "reference_primate.csv"
            ),
            "raw_row_count": r_results["data_processing"]["raw_row_count"],
            "processed_row_count": r_results["data_processing"]["processed_row_count"],
            "processed_species_count": r_results["data_processing"][
                "processed_species_count"
            ],
            "checked_in_processed_matches_reference": r_results["data_processing"][
                "checked_in_processed_matches_reference"
            ],
            "reference_column_order": list(reference_rows[0].keys()),
            "boundary_note": (
                "The raw workbook remains external; the governed processed reference table is the reproducible repository-side handoff."
            ),
        }
    if evidence_id == "evidence-003":
        numeric_parse_failures = {
            column: [
                row["species"]
                for row in reference_rows
                if row[column] and _safe_float(row[column]) is None
            ]
            for column in NUMERIC_COLUMNS
        }
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "factor_columns": ["sex_dimorphism"],
            "numeric_columns": NUMERIC_COLUMNS,
            "text_columns": TEXT_COLUMNS + ["mating_system"],
            "numeric_parse_failures": numeric_parse_failures,
            "numeric_columns_parse_without_error": all(
                not failures for failures in numeric_parse_failures.values()
            ),
            "sex_dimorphism_values": sorted(
                {row["sex_dimorphism"] for row in reference_rows}
            ),
            "mating_system_values_sample": sorted(
                {row["mating_system"] for row in reference_rows}
            )[:6],
            "checked_in_processed_matches_reference": r_results["data_processing"][
                "checked_in_processed_matches_reference"
            ],
        }
    if evidence_id == "evidence-004":
        missing_counts = _missing_counts(reference_rows)
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "missing_counts_by_column": missing_counts,
            "columns_with_missing_values": [
                column for column, count in missing_counts.items() if count > 0
            ],
            "rows_with_any_missing_values": sum(
                1
                for row in reference_rows
                if any(value == "" for value in row.values())
            ),
            "na_rm_reference_rule": "across(body_mass:social_group_size, ~mean(.x, na.rm = TRUE))",
            "post_repair_row_count": len(reference_rows),
        }
    if evidence_id == "evidence-005":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "raw_row_count": r_results["data_processing"]["raw_row_count"],
            "processed_row_count": r_results["data_processing"]["processed_row_count"],
            "processed_species_count": r_results["data_processing"][
                "processed_species_count"
            ],
            "raw_to_processed_row_delta": r_results["data_processing"]["raw_row_count"]
            - r_results["data_processing"]["processed_row_count"],
            "duplicate_species_contract_count": r_results["data_processing"][
                "duplicate_species_after_grouping"
            ],
            "duplicates_remaining_after_grouping": 0,
        }
    if evidence_id == "evidence-006":
        tree_processing = block_payloads["tree-import-and-pruning"]["evidence"]
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "original_tree_locator": "external:lund/pcm1-plots-signal/data/primatetree.nex",
            "trimmed_tree_locator": "external:lund/pcm1-plots-signal/data/trimmed_primatetree.nex",
            "r_original_tip_count": tree_processing["r"]["original_tip_count"],
            "bijux_original_tip_count": tree_processing["bijux"]["original_tree"][
                "inspect"
            ]["tip_count"],
            "r_trimmed_tip_count": tree_processing["r"]["trimmed_tip_count"],
            "bijux_trimmed_tip_count": tree_processing["bijux"]["trimmed_tree"][
                "inspect"
            ]["tip_count"],
            "original_tree_has_branch_lengths": tree_processing["bijux"][
                "original_tree"
            ]["inspect"]["has_branch_lengths"],
            "trimmed_tree_has_branch_lengths": tree_processing["bijux"]["trimmed_tree"][
                "inspect"
            ]["has_branch_lengths"],
        }
    if evidence_id == "evidence-007":
        tree_processing = block_payloads["tree-import-and-pruning"]["evidence"]
        extract_clade = block_payloads["extract-clade-node-77"]["evidence"]
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "rooted": tree_processing["r"]["rooted"],
            "binary": tree_processing["r"]["binary"],
            "ultrametric": tree_processing["r"]["ultrametric"],
            "bijux_rooted": bijux_results["data_tree_alignment"]["readiness"]["rooted"],
            "bijux_binary": bijux_results["data_tree_alignment"]["readiness"]["binary"],
            "node_label_contract": {
                "r_source_node_numeric": extract_clade["r_source_node_numeric"],
                "r_source_node_label": extract_clade["r_source_node_label"],
                "bijux_matched_node_name": extract_clade["bijux_matched_node_name"],
                "same_taxa": extract_clade["same_taxa"],
            },
        }
    if evidence_id == "evidence-008":
        tree_processing = block_payloads["tree-import-and-pruning"]["evidence"]
        tip_order = block_payloads["tip-order-alignment"]["evidence"]
        tree_join = block_payloads["treeio-node-mapping-and-join"]["evidence"]
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "missing_tips": tree_processing["r"]["missing_tips"],
            "checked_in_trimmed_tip_set_matches_reference": tree_processing["r"][
                "checked_in_trimmed_tip_set_matches_reference"
            ],
            "checked_in_trimmed_tip_order_matches_reference": tree_processing["r"][
                "checked_in_trimmed_tip_order_matches_reference"
            ],
            "aligned_species_equals_tip_order": tip_order["r_aligned_equals_tip_order"],
            "first_six_species_match": tip_order["first_six_species_match"],
            "nodeid_examples_r": tree_join["nodeid_examples_r"],
            "nodeid_examples_bijux": tree_join["nodeid_examples_bijux"],
        }
    if evidence_id == "evidence-009":
        tree_processing = block_payloads["tree-import-and-pruning"]["evidence"]
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "processed_csv_locator": r_results["data_processing"][
                "checked_in_processed_path"
            ],
            "governed_processed_reference_locator": r_results["data_processing"][
                "reference_processed_path"
            ],
            "trimmed_tree_locator": "external:lund/pcm1-plots-signal/data/trimmed_primatetree.nex",
            "governed_trimmed_tree_reference_locator": (
                "evidence-book/studies/primate-longevity-signal/evidence-001/"
                "reference_trimmed_primatetree.nwk"
            ),
            "checked_in_processed_matches_reference": r_results["data_processing"][
                "checked_in_processed_matches_reference"
            ],
            "checked_in_trimmed_tip_set_matches_reference": tree_processing["r"][
                "checked_in_trimmed_tip_set_matches_reference"
            ],
            "checked_in_trimmed_tip_order_matches_reference": tree_processing["r"][
                "checked_in_trimmed_tip_order_matches_reference"
            ],
        }
    raise ValueError(f"unsupported primate component bundle: {evidence_id}")


def _safe_float(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _build_claims_payload(spec: dict[str, object]) -> dict[str, object]:
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": spec["evidence_id"],
        "claim_count": 1,
        "claims": [
            {
                "claim_id": spec["claim_id"],
                "claim_title": spec["claim_title"],
                "summary": spec["claim_summary"],
                "verdict": spec["verdict"],
                "evidence_ids": [spec["evidence_id"]],
                "source_fragments": spec["source_fragments"],
            }
        ],
    }


def _build_manifest(
    spec: dict[str, object], *, report_filename: str, report_payload: dict[str, object]
) -> dict[str, object]:
    source_basis = list(spec["source_basis"]) + [
        {
            "kind": "repository-reference",
            "label": f"{spec['title']} report payload",
            "locator": (
                f"evidence-book/studies/primate-longevity-signal/{spec['evidence_id']}/"
                f"{report_filename}"
            ),
        }
    ]
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": spec["evidence_id"],
        "evidence_title": spec["title"],
        "summary": spec["summary"],
        "owner_package": "bijux-phylogenetics",
        "claim_ids": [spec["claim_id"]],
        "source_basis": source_basis,
        "freshness": {
            "last_generated_on": "2026-05-10",
            "governed_code_paths": [
                "packages/bijux-phylogenetics/src/bijux_phylogenetics/evidence/studies/primate_pcm1_component_bundles.py"
            ],
            "source_basis_locators": [entry["locator"] for entry in source_basis],
        },
        "ownership": {
            "owner_package": "bijux-phylogenetics",
            "analytical_surfaces": spec["analytical_surfaces"],
        },
        "claim_tags": spec["claim_tags"],
        "comparison_mode": spec["comparison_mode"],
        "verdict": {
            "status": spec["verdict"],
            "summary": spec["claim_summary"],
        },
        "limitations": spec["limitations"],
        "source_fragments": spec["source_fragments"],
        "reference_script_locators": _reference_script_locators(
            spec["reference_line_specs"]
        ),
        "reference_bundle_locator": REFERENCE_BUNDLE_ROOT.as_posix(),
        "supporting_report_locator": (
            f"evidence-book/studies/primate-longevity-signal/{spec['evidence_id']}/"
            f"{report_filename}"
        ),
        "report_keys": sorted(report_payload.keys()),
    }


def _render_bundle_readme(
    spec: dict[str, object], *, report_filename: str, manifest: dict[str, object]
) -> str:
    lines = [
        f"# {spec['title']}",
        "",
        spec["summary"],
        "",
        f"- evidence id: `{spec['evidence_id']}`",
        f"- claim id: `{spec['claim_id']}`",
        f"- verdict: `{spec['verdict']}`",
        f"- source fragments: `{', '.join(spec['source_fragments'])}`",
        "",
        "Governed files:",
        "",
        f"- `{report_filename}`",
        "- `claims.json`",
        "- `manifest.json`",
        "",
        "Reference script locators:",
        "",
    ]
    lines.extend(f"- `{locator}`" for locator in manifest["reference_script_locators"])
    lines.extend(
        [
            "",
            "Limitations:",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in spec["limitations"])
    lines.append("")
    return "\n".join(lines)


def build_primate_pcm1_component_bundles(
    repo_root: Path,
) -> dict[str, dict[str, object]]:
    context = _load_context(repo_root)
    bundles: dict[str, dict[str, object]] = {}
    for spec in COMPONENT_BUNDLE_DEFINITIONS:
        report_payload = _report_payload_for_spec(spec, context)
        manifest = _build_manifest(
            spec,
            report_filename=spec["report_filename"],
            report_payload=report_payload,
        )
        bundles[spec["evidence_id"]] = {
            "manifest": manifest,
            "claims": _build_claims_payload(spec),
            "report_filename": spec["report_filename"],
            "report_payload": report_payload,
            "readme": _render_bundle_readme(
                spec,
                report_filename=spec["report_filename"],
                manifest=manifest,
            ),
        }
    return bundles


def build_primate_data_preparation_bundle_index(repo_root: Path) -> dict[str, object]:
    bundles = build_primate_pcm1_component_bundles(repo_root)
    entries: list[dict[str, object]] = []
    for evidence_id in DATA_PREPARATION_BUNDLE_IDS:
        manifest = bundles[evidence_id]["manifest"]
        entries.append(
            {
                "evidence_id": evidence_id,
                "title": manifest["evidence_title"],
                "claim_id": manifest["claim_ids"][0],
                "relative_path": (
                    Path("studies") / "primate-longevity-signal" / evidence_id
                ).as_posix(),
                "claim_tags": manifest["claim_tags"],
                "source_fragments": manifest["source_fragments"],
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "bundle_count": len(entries),
        "bundle_ids": DATA_PREPARATION_BUNDLE_IDS,
        "bundles": entries,
        "review_rule": (
            "Downstream lambda and ancestral-state evidence must rely on these preprocessing bundles instead of hiding data preparation inside later model comparisons."
        ),
    }


def render_primate_data_preparation_bundle_index_markdown(
    payload: dict[str, object],
) -> str:
    lines = [
        "# Primate Data-Preparation Parity Bundles",
        "",
        payload["review_rule"],
        "",
        f"Bundles: `{payload['bundle_count']}`",
        "",
    ]
    for entry in payload["bundles"]:
        lines.append(
            f"- `{entry['evidence_id']}` — {entry['title']} "
            f"(`{', '.join(entry['source_fragments'])}`)"
        )
    lines.append("")
    return "\n".join(lines)


def build_primate_structural_parity_table(repo_root: Path) -> dict[str, object]:
    context = _load_context(repo_root)
    block_payloads = context["block_payloads"]
    r_results = context["r_results"]
    bijux_results = context["bijux_results"]
    tree_processing = block_payloads["tree-import-and-pruning"]["evidence"]
    unrooted = block_payloads["unrooted-tree-demo"]["evidence"]
    extract_clade = block_payloads["extract-clade-node-77"]["evidence"]
    rotate_nodes = block_payloads["rotate-nodes-behavior"]["evidence"]
    tip_order = block_payloads["tip-order-alignment"]["evidence"]
    tree_join = block_payloads["treeio-node-mapping-and-join"]["evidence"]

    def row(
        *,
        row_id: str,
        evidence_id: str,
        metric_name: str,
        r_value: object,
        bijux_value: object,
    ) -> dict[str, object]:
        matched = r_value == bijux_value
        return {
            "row_id": row_id,
            "evidence_id": evidence_id,
            "metric_name": metric_name,
            "r_value": r_value,
            "bijux_value": bijux_value,
            "verdict": "matched" if matched else "mismatch_unexplained",
        }

    rows = [
        row(
            row_id="original-tip-count",
            evidence_id="evidence-006",
            metric_name="original_tip_count",
            r_value=tree_processing["r"]["original_tip_count"],
            bijux_value=tree_processing["bijux"]["original_tree"]["inspect"][
                "tip_count"
            ],
        ),
        row(
            row_id="trimmed-tip-count",
            evidence_id="evidence-006",
            metric_name="trimmed_tip_count",
            r_value=tree_processing["r"]["trimmed_tip_count"],
            bijux_value=tree_processing["bijux"]["trimmed_tree"]["inspect"][
                "tip_count"
            ],
        ),
        row(
            row_id="rooted",
            evidence_id="evidence-007",
            metric_name="rooted",
            r_value=tree_processing["r"]["rooted"],
            bijux_value=bijux_results["data_tree_alignment"]["readiness"]["rooted"],
        ),
        row(
            row_id="binary",
            evidence_id="evidence-007",
            metric_name="binary",
            r_value=tree_processing["r"]["binary"],
            bijux_value=bijux_results["data_tree_alignment"]["readiness"]["binary"],
        ),
        row(
            row_id="ultrametric",
            evidence_id="evidence-007",
            metric_name="ultrametric",
            r_value=tree_processing["r"]["ultrametric"],
            bijux_value=tree_processing["bijux"]["trimmed_tree"]["validate"][
                "ultrametric"
            ],
        ),
        row(
            row_id="node-label-match",
            evidence_id="evidence-007",
            metric_name="source_node_label",
            r_value=extract_clade["r_source_node_label"],
            bijux_value=extract_clade["bijux_matched_node_name"],
        ),
        row(
            row_id="missing-tip-set",
            evidence_id="evidence-008",
            metric_name="missing_tips",
            r_value="|".join(sorted(r_results["tree_processing"]["missing_tips"])),
            bijux_value="|".join(
                sorted(tree_processing["bijux"]["pruning"]["removed_taxa"])
            ),
        ),
        row(
            row_id="tip-order-aligned",
            evidence_id="evidence-008",
            metric_name="aligned_species_equals_tip_order",
            r_value=tip_order["r_aligned_equals_tip_order"],
            bijux_value=tip_order["bijux_aligned_equals_tip_order"],
        ),
        row(
            row_id="first-six-tip-order",
            evidence_id="evidence-008",
            metric_name="aligned_species_first_6",
            r_value="|".join(
                r_results["data_tree_alignment"]["aligned_species_first_6"]
            ),
            bijux_value="|".join(
                bijux_results["data_tree_alignment"]["aligned_species_first_6"]
            ),
        ),
        row(
            row_id="joined-tip-count",
            evidence_id="evidence-008",
            metric_name="joined_tip_count",
            r_value=r_results["data_tree_alignment"]["joined_tip_count"],
            bijux_value=tree_join["bijux_analysis_taxa_count"],
        ),
        row(
            row_id="nodeid-pan-paniscus",
            evidence_id="evidence-008",
            metric_name="nodeid:Pan_paniscus",
            r_value=tree_join["nodeid_examples_r"]["Pan_paniscus"],
            bijux_value=tree_join["nodeid_examples_bijux"]["Pan_paniscus"],
        ),
        row(
            row_id="nodeid-hylobates-lar",
            evidence_id="evidence-008",
            metric_name="nodeid:Hylobates_lar",
            r_value=tree_join["nodeid_examples_r"]["Hylobates_lar"],
            bijux_value=tree_join["nodeid_examples_bijux"]["Hylobates_lar"],
        ),
        row(
            row_id="nodeid-node32",
            evidence_id="evidence-008",
            metric_name="nodeid:Node32",
            r_value=tree_join["nodeid_examples_r"]["Node32"],
            bijux_value=tree_join["nodeid_examples_bijux"]["Node32"],
        ),
        row(
            row_id="extract-clade-tip-count",
            evidence_id="evidence-007",
            metric_name="extract_clade_tip_count",
            r_value=extract_clade["r_tip_count"],
            bijux_value=extract_clade["bijux_tip_count"],
        ),
        row(
            row_id="unrooted-tip-count",
            evidence_id="evidence-007",
            metric_name="unrooted_tip_count",
            r_value=unrooted["r"]["tip_count"],
            bijux_value=unrooted["bijux"]["tip_count"],
        ),
        row(
            row_id="rotate-single-tip-order-sha256",
            evidence_id="evidence-007",
            metric_name="rotate_single_tip_order_sha256",
            r_value=_sha256_text("|".join(rotate_nodes["r"]["tip_order"])),
            bijux_value=_sha256_text("|".join(rotate_nodes["bijux"]["tip_order"])),
        ),
        row(
            row_id="rotate-all-tip-order-sha256",
            evidence_id="evidence-007",
            metric_name="rotate_all_tip_order_sha256",
            r_value=_sha256_text("|".join(rotate_nodes["r_all"]["tip_order"])),
            bijux_value=_sha256_text("|".join(rotate_nodes["bijux_all"]["tip_order"])),
        ),
    ]
    verdict_counts: dict[str, int] = {}
    for entry in rows:
        verdict = entry["verdict"]
        verdict_counts[verdict] = verdict_counts.get(verdict, 0) + 1
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "row_count": len(rows),
        "verdict_counts": dict(sorted(verdict_counts.items())),
        "rows": rows,
    }


def render_primate_structural_parity_table_markdown(payload: dict[str, object]) -> str:
    lines = [
        "# Primate Structural Parity Table",
        "",
        "| Row | Evidence | Metric | R | Bijux | Verdict |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload["rows"]:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["row_id"]),
                    str(row["evidence_id"]),
                    str(row["metric_name"]),
                    str(row["r_value"]),
                    str(row["bijux_value"]),
                    str(row["verdict"]),
                ]
            )
            + " |"
        )
    lines.extend(["", f"Rows: `{payload['row_count']}`", ""])
    return "\n".join(lines)
