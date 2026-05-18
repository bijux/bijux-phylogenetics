from __future__ import annotations

GEIGER_TREEDATA_REFERENCE_PAYLOADS: dict[str, dict[str, object]] = {
    "example_tree_example_traits": {
        "aligned_taxa": ["A", "B", "C"],
        "dropped_tree_taxa": ["D"],
        "dropped_trait_taxa": ["E"],
    },
    "example_tree_reordered_traits": {
        "aligned_taxa": ["A", "B", "C", "D"],
        "dropped_tree_taxa": [],
        "dropped_trait_taxa": [],
    },
    "example_tree_traits_validate": {
        "aligned_taxa": ["A", "B", "C", "D"],
        "dropped_tree_taxa": [],
        "dropped_trait_taxa": [],
    },
}
