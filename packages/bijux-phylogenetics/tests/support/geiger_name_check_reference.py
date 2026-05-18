from __future__ import annotations

GEIGER_NAME_CHECK_REFERENCE_PAYLOADS: dict[str, dict[str, object]] = {
    "example_tree_example_traits": {
        "reference_outcome": "mismatch",
        "tree_not_data": ["D"],
        "data_not_tree": ["E"],
    },
    "example_tree_reordered_traits": {
        "reference_outcome": "OK",
        "tree_not_data": [],
        "data_not_tree": [],
    },
    "case_sensitive_probe": {
        "reference_outcome": "mismatch",
        "tree_not_data": ["A"],
        "data_not_tree": ["a"],
    },
}
