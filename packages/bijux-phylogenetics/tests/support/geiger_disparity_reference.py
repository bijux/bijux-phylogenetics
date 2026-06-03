from __future__ import annotations

GEIGER_DISPARITY_REFERENCE_PAYLOADS: dict[str, dict[str, object]] = {
    "example_tree_response_univariate": {
        "trait_columns": ["response"],
        "ape_node_id_order": [5, 6, 7],
        "descendant_taxa": [
            ["A", "B", "C", "D"],
            ["A", "B"],
            ["C", "D"],
        ],
        "clade_disparity": [2.16666666666667, 2.25, 2.25],
    },
    "example_tree_response_predictor_one_multivariate": {
        "trait_columns": ["response", "predictor_one"],
        "ape_node_id_order": [5, 6, 7],
        "descendant_taxa": [
            ["A", "B", "C", "D"],
            ["A", "B"],
            ["C", "D"],
        ],
        "clade_disparity": [5.5, 3.25, 3.25],
    },
}
