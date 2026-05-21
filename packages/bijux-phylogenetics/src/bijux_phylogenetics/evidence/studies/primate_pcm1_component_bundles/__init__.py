from __future__ import annotations

from datetime import date
import hashlib
from pathlib import Path

from .definitions import CATEGORICAL_COLUMNS
from .definitions import COMPONENT_BUNDLE_DEFINITIONS
from .definitions import DATA_PREPARATION_BUNDLE_IDS
from .definitions import NUMERIC_COLUMNS
from .definitions import REFERENCE_BUNDLE_ID
from .definitions import REFERENCE_BUNDLE_ROOT
from .definitions import REFERENCE_SCRIPT_PATH
from .definitions import STUDY_ID
from .definitions import TEXT_COLUMNS
from .report_payloads import build_component_report_payload
from .study_context import load_reference_context
from .study_context import reference_bundle_root
from .study_context import reference_script_locators
from .study_context import study_root


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
            "last_generated_on": date.today().isoformat(),
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
        "reference_script_locators": reference_script_locators(
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
    context = load_reference_context(repo_root)
    bundles: dict[str, dict[str, object]] = {}
    for spec in COMPONENT_BUNDLE_DEFINITIONS:
        report_payload = build_component_report_payload(spec, context)
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
    context = load_reference_context(repo_root)
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
