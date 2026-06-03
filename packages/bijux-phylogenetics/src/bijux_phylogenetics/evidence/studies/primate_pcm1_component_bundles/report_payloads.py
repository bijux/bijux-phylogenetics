from __future__ import annotations

from .definitions import NUMERIC_COLUMNS, STUDY_ID, TEXT_COLUMNS
from .study_context import missing_counts


def build_component_report_payload(
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
                if row[column] and safe_float(row[column]) is None
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
        column_missing_counts = missing_counts(reference_rows)
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "missing_counts_by_column": column_missing_counts,
            "columns_with_missing_values": [
                column for column, count in column_missing_counts.items() if count > 0
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


def safe_float(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None
