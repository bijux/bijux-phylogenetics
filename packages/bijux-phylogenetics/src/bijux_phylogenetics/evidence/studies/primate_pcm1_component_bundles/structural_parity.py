from __future__ import annotations

import hashlib
from pathlib import Path

from .definitions import STUDY_ID
from .study_context import load_reference_context


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
            r_value=sha256_text("|".join(rotate_nodes["r"]["tip_order"])),
            bijux_value=sha256_text("|".join(rotate_nodes["bijux"]["tip_order"])),
        ),
        row(
            row_id="rotate-all-tip-order-sha256",
            evidence_id="evidence-007",
            metric_name="rotate_all_tip_order_sha256",
            r_value=sha256_text("|".join(rotate_nodes["r_all"]["tip_order"])),
            bijux_value=sha256_text("|".join(rotate_nodes["bijux_all"]["tip_order"])),
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


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()
