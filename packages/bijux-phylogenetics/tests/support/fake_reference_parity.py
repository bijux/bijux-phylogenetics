from __future__ import annotations

from pathlib import Path

from tests.support.fake_external_engines import write_executable


def fake_ape_rscript(
    path: Path,
    *,
    ape_available: bool = True,
    summary_overrides: dict[str, object] | None = None,
    normalized_tree_overrides: dict[str, str] | None = None,
) -> Path:
    summary_payload = repr(summary_overrides or {})
    normalized_tree_payload = repr(normalized_tree_overrides or {})
    return write_executable(
        path,
        f"""#!/usr/bin/env python3
import csv
import json
import math
import sys
from pathlib import Path

from Bio import Phylo

TABULAR_CASES = {{
    "dna-base-frequency-lowercase": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "state_count": 17,
        }},
        "rows_name": "base-frequency.tsv",
        "rows": [
            {{"state": "a", "frequency": 6 / 18}},
            {{"state": "c", "frequency": 2 / 18}},
            {{"state": "g", "frequency": 3 / 18}},
            {{"state": "t", "frequency": 5 / 18}},
            {{"state": "r", "frequency": 0.0}},
            {{"state": "m", "frequency": 0.0}},
            {{"state": "w", "frequency": 0.0}},
            {{"state": "s", "frequency": 0.0}},
            {{"state": "k", "frequency": 0.0}},
            {{"state": "y", "frequency": 0.0}},
            {{"state": "v", "frequency": 0.0}},
            {{"state": "h", "frequency": 0.0}},
            {{"state": "d", "frequency": 0.0}},
            {{"state": "b", "frequency": 0.0}},
            {{"state": "n", "frequency": 2 / 18}},
            {{"state": "-", "frequency": 0.0}},
            {{"state": "?", "frequency": 0.0}},
        ],
    }},
    "dna-base-frequency-ambiguity": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "state_count": 17,
        }},
        "rows_name": "base-frequency.tsv",
        "rows": [
            {{"state": "a", "frequency": 3 / 18}},
            {{"state": "c", "frequency": 3 / 18}},
            {{"state": "g", "frequency": 3 / 18}},
            {{"state": "t", "frequency": 3 / 18}},
            {{"state": "r", "frequency": 1 / 18}},
            {{"state": "m", "frequency": 0.0}},
            {{"state": "w", "frequency": 0.0}},
            {{"state": "s", "frequency": 0.0}},
            {{"state": "k", "frequency": 0.0}},
            {{"state": "y", "frequency": 0.0}},
            {{"state": "v", "frequency": 0.0}},
            {{"state": "h", "frequency": 0.0}},
            {{"state": "d", "frequency": 0.0}},
            {{"state": "b", "frequency": 0.0}},
            {{"state": "n", "frequency": 1 / 18}},
            {{"state": "-", "frequency": 1 / 18}},
            {{"state": "?", "frequency": 3 / 18}},
        ],
    }},
    "dna-raw-distance-clean": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 8,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.125}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.5}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.625}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.125}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.625}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.5}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.5}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.625}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.125}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.625}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.5}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.125}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-raw-distance-gaps": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 1 / 6}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 1 / 2}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 1 / 6}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 1 / 2}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 1 / 2}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 1 / 2}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 1 / 4}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 1 / 4}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-raw-distance-identical": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 8,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 1 / 8}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 1 / 8}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 1 / 8}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 1 / 8}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 1 / 8}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 1 / 8}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 1 / 8}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 1 / 8}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 1 / 8}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 1 / 8}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-raw-distance-high-divergence": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 4,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 1.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.25}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 1.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.75}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.25}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.75}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-raw-distance-missing-data": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-translation-valid-frame": {{
        "summary": {{
            "sequence_count": 3,
            "translated_length": 3,
            "stop_codon_count": 0,
        }},
        "rows_name": "translation.tsv",
        "rows": [
            {{"identifier": "valid_a", "amino_acid_sequence": "MEL"}},
            {{"identifier": "valid_b", "amino_acid_sequence": "MKM"}},
            {{"identifier": "valid_c", "amino_acid_sequence": "MTW"}},
        ],
    }},
    "dna-translation-internal-stop": {{
        "summary": {{
            "sequence_count": 1,
            "translated_length": 3,
            "stop_codon_count": 1,
        }},
        "rows_name": "translation.tsv",
        "rows": [
            {{"identifier": "internal_stop", "amino_acid_sequence": "M*W"}},
        ],
    }},
    "dna-translation-terminal-stop": {{
        "summary": {{
            "sequence_count": 1,
            "translated_length": 3,
            "stop_codon_count": 1,
        }},
        "rows_name": "translation.tsv",
        "rows": [
            {{"identifier": "terminal_stop", "amino_acid_sequence": "ME*"}},
        ],
    }},
}}

SUMMARY_OVERRIDES = {summary_payload}
NORMALIZED_TREE_OVERRIDES = {normalized_tree_payload}
APE_AVAILABLE = {str(ape_available)}

def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\\n", encoding="utf-8")

def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\\t")
        writer.writeheader()
        writer.writerows(rows)

def normalize_label(value):
    return "" if value is None else str(value)

def parse_support_label(value):
    text = normalize_label(value)
    if not text:
        return ""
    if "/" in text:
        try:
            return float(text.split("/")[-1])
        except ValueError:
            return ""
    try:
        return float(text)
    except ValueError:
        return ""

def descendant_taxa(clade):
    return sorted(terminal.name for terminal in clade.get_terminals() if terminal.name)

def clade_rows(tree, tree_index_value):
    rows = []
    for clade in tree.find_clades(order="preorder"):
        taxa = descendant_taxa(clade)
        node_kind = "tip" if clade.is_terminal() else "internal"
        if clade == tree.root:
            node_kind = "root"
        if clade.is_terminal():
            node_label = normalize_label(clade.name)
        else:
            node_label = normalize_label(clade.name if clade.name is not None else getattr(clade, "confidence", None))
        rows.append(
            {{
                "tree_index": tree_index_value,
                "node_kind": node_kind,
                "clade_id": "|".join(taxa),
                "node_label": node_label,
                "taxon_count": len(taxa),
                "taxa": "|".join(taxa),
                "support": parse_support_label(node_label),
                "branch_length": "" if clade.branch_length is None else clade.branch_length,
            }}
        )
    node_order = {{"root": 0, "internal": 1, "tip": 2}}
    return sorted(
        rows,
        key=lambda row: (
            0 if row["tree_index"] == "" else int(row["tree_index"]),
            node_order.get(row["node_kind"], 9),
            row["clade_id"],
            row["node_label"],
        ),
    )

def is_rooted_tree(tree):
    return len(getattr(tree.root, "clades", [])) == 2

def iter_internal_clades_preorder(clade):
    yield clade
    for child in getattr(clade, "clades", []):
        if child.is_terminal():
            continue
        yield from iter_internal_clades_preorder(child)

def node_depth_lookup(tree):
    lookup = {{id(tree.root): 0.0}}
    def walk(clade):
        base_depth = lookup[id(clade)]
        for child in getattr(clade, "clades", []):
            lookup[id(child)] = base_depth + float(child.branch_length or 0.0)
            walk(child)
    walk(tree.root)
    return lookup

def node_depth_rows(tree):
    depths = node_depth_lookup(tree)
    rows = []
    tip_clades = list(tree.get_terminals())
    internal_clades = list(iter_internal_clades_preorder(tree.root))
    for node_id, clade in enumerate(tip_clades, start=1):
        rows.append(
            {{
                "node_id": node_id,
                "node_kind": "tip",
                "node_label": normalize_label(clade.name),
                "descendant_taxa": "|".join(descendant_taxa(clade)),
                "branch_length_depth": depths[id(clade)],
                "branch_length": ""
                if clade.branch_length is None
                else clade.branch_length,
            }}
        )
    for offset, clade in enumerate(internal_clades, start=1):
        node_id = len(tip_clades) + offset
        node_label = normalize_label(
            clade.name if clade.name is not None else getattr(clade, "confidence", None)
        )
        rows.append(
            {{
                "node_id": node_id,
                "node_kind": "root" if clade is tree.root else "internal",
                "node_label": node_label,
                "descendant_taxa": "|".join(descendant_taxa(clade)),
                "branch_length_depth": depths[id(clade)],
                "branch_length": ""
                if clade.branch_length is None
                else clade.branch_length,
            }}
        )
    return rows

def branching_time_rows(tree):
    node_depth_table = node_depth_rows(tree)
    tip_depths = [
        row["branch_length_depth"]
        for row in node_depth_table
        if row["node_kind"] == "tip"
    ]
    root_age = max(tip_depths)
    return [
        {{
            "node_id": row["node_id"],
            "node_kind": row["node_kind"],
            "node_label": row["node_label"],
            "descendant_taxa": row["descendant_taxa"],
            "node_depth": row["branch_length_depth"],
            "branching_time": root_age - row["branch_length_depth"],
        }}
        for row in node_depth_table
        if row["node_kind"] != "tip"
    ]

def ultrametric_rows(tree):
    node_depth_table = node_depth_rows(tree)
    return [
        {{
            "node_id": row["node_id"],
            "tip_label": row["node_label"],
            "root_to_tip_depth": row["branch_length_depth"],
            "deviation_from_mean_depth": None,
            "deviation_from_min_depth": None,
            "deviation_from_max_depth": None,
            "is_offending_taxon": False,
        }}
        for row in node_depth_table
        if row["node_kind"] == "tip"
    ]

def signature_id(taxa):
    return "|".join(sorted(set(taxa)))

def rooted_topology_signatures(tree):
    total_tip_count = len(list(tree.get_terminals()))
    signatures = {{}}
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root or clade.is_terminal():
            continue
        taxa = descendant_taxa(clade)
        if len(taxa) <= 1 or len(taxa) >= total_tip_count:
            continue
        signatures[signature_id(taxa)] = taxa
    return signatures

def canonical_unrooted_signature(taxa, all_taxa):
    selected = sorted(set(taxa))
    complement = sorted(taxon for taxon in all_taxa if taxon not in selected)
    if len(selected) < len(complement):
        return selected
    if len(complement) < len(selected):
        return complement
    if "|".join(selected) <= "|".join(complement):
        return selected
    return complement

def unrooted_topology_signatures(tree):
    all_taxa = descendant_taxa(tree.root)
    signatures = {{}}
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root or clade.is_terminal():
            continue
        selected = canonical_unrooted_signature(descendant_taxa(clade), all_taxa)
        if len(selected) <= 1 or len(selected) >= len(all_taxa):
            continue
        signatures[signature_id(selected)] = selected
    return signatures

def topology_distance_rows(left_tree, right_tree, rf_mode):
    if rf_mode == "rooted":
        left_signatures = rooted_topology_signatures(left_tree)
        right_signatures = rooted_topology_signatures(right_tree)
        split_kind = "clade"
    else:
        left_signatures = unrooted_topology_signatures(left_tree)
        right_signatures = unrooted_topology_signatures(right_tree)
        split_kind = "split"
    left_ids = set(left_signatures)
    right_ids = set(right_signatures)
    shared_ids = left_ids & right_ids
    left_only_ids = left_ids - right_ids
    right_only_ids = right_ids - left_ids
    all_ids = sorted(
        left_ids | right_ids,
        key=lambda split_id: (len(split_id.split("|")), split_id),
    )
    rows = [
        {{
            "split_id": split_id,
            "split_kind": split_kind,
            "comparison_status": (
                "shared"
                if split_id in shared_ids
                else "left_only"
                if split_id in left_only_ids
                else "right_only"
            ),
            "taxon_count": len((left_signatures.get(split_id) or right_signatures[split_id])),
            "descendant_taxa": "|".join(left_signatures.get(split_id) or right_signatures[split_id]),
            "left_present": split_id in left_ids,
            "right_present": split_id in right_ids,
        }}
        for split_id in all_ids
    ]
    return rows, len(left_ids), len(right_ids), len(shared_ids), len(left_only_ids), len(right_only_ids)

def informative_clades(tree, shared_taxa):
    total_tip_count = len(shared_taxa)
    clades = []
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root or clade.is_terminal():
            continue
        taxa = frozenset(descendant_taxa(clade))
        if 1 < len(taxa) < total_tip_count and taxa <= shared_taxa:
            clades.append(taxa)
    return clades

def split_counts(tree_set):
    counts = {{}}
    for tree in tree_set:
        for split_id in unrooted_topology_signatures(tree):
            counts[split_id] = counts.get(split_id, 0) + 1
    return counts

def canonical_bipartition(taxa, all_taxa):
    return "|".join(canonical_unrooted_signature(sorted(taxa), sorted(all_taxa)))

def clade_support_status(supporting_tree_count, tree_count, node_kind, unscored_reason=None):
    if node_kind == "root":
        return (
            "fixed",
            "the root spans the full compatible taxon set and is present in every comparison tree",
        )
    if supporting_tree_count is None:
        if unscored_reason == "absent-root-split":
            return (
                "not-counted",
                "ape::prop.clades leaves this root-adjacent split unscored when the comparison tree set never realizes the matching bipartition",
            )
        return (
            "not-counted",
            "ape::prop.clades leaves this root-adjacent clade unscored because its complement is a singleton tip",
        )
    if supporting_tree_count == 0:
        return (
            "absent",
            "the reference clade is absent from the comparison tree set",
        )
    if supporting_tree_count == tree_count:
        return (
            "fixed",
            "the reference clade is present in every comparison tree",
        )
    return (
        "partial-support",
        "the reference clade is present in only a subset of comparison trees",
    )

def prop_clades_rows(reference_tree, comparison_trees):
    reference_taxa = frozenset(terminal.name for terminal in reference_tree.get_terminals())
    comparison_taxa = frozenset(terminal.name for terminal in comparison_trees[0].get_terminals())
    if any(
        frozenset(terminal.name for terminal in tree.get_terminals()) != comparison_taxa
        for tree in comparison_trees[1:]
    ):
        raise ValueError(
            "reference tree support mapping requires all comparison trees to share the exact same taxon set"
        )
    if reference_taxa != comparison_taxa:
        raise ValueError(
            "reference tree and comparison tree set must share the exact same taxon set"
        )

    tree_count = len(comparison_trees)
    clade_counts = {{}}
    for tree in comparison_trees:
        for clade in informative_clades(tree, comparison_taxa):
            clade_id = "|".join(sorted(clade))
            clade_counts[clade_id] = clade_counts.get(clade_id, 0) + 1
    split_count_lookup = split_counts(comparison_trees)
    depths = node_depth_lookup(reference_tree)
    rows = []
    supported_clade_count = 0
    absent_clade_count = 0
    unscored_clade_count = 0
    root_children = list(getattr(reference_tree.root, "clades", []))
    tip_count = len(list(reference_tree.get_terminals()))
    internal_clades = list(iter_internal_clades_preorder(reference_tree.root))

    for offset, clade in enumerate(internal_clades, start=1):
        node_id = tip_count + offset
        taxa = descendant_taxa(clade)
        clade_id = "|".join(taxa)
        node_kind = "root" if clade is reference_tree.root else "internal"
        if clade is reference_tree.root:
            supporting_tree_count = tree_count
            clade_frequency = 1.0
            support_percent = 100.0
            unscored_reason = None
        elif len(taxa) == len(reference_taxa) - 1:
            supporting_tree_count = None
            clade_frequency = ""
            support_percent = ""
            unscored_reason = "singleton-complement"
            unscored_clade_count += 1
        elif clade in root_children:
            split_support = split_count_lookup.get(
                canonical_bipartition(taxa, reference_taxa),
                0,
            )
            if split_support == 0:
                supporting_tree_count = None
                clade_frequency = ""
                support_percent = ""
                unscored_reason = "absent-root-split"
                unscored_clade_count += 1
            else:
                supporting_tree_count = split_support
                clade_frequency = supporting_tree_count / tree_count
                support_percent = clade_frequency * 100.0
                unscored_reason = None
                supported_clade_count += 1
        else:
            supporting_tree_count = clade_counts.get(clade_id, 0)
            clade_frequency = supporting_tree_count / tree_count
            support_percent = clade_frequency * 100.0
            unscored_reason = None
            if supporting_tree_count == 0:
                absent_clade_count += 1
            else:
                supported_clade_count += 1

        support_status, explanation = clade_support_status(
            supporting_tree_count,
            tree_count,
            node_kind,
            unscored_reason,
        )
        node_label = normalize_label(
            clade.name if clade.name is not None else getattr(clade, "confidence", None)
        )
        rows.append(
            {{
                "node_id": node_id,
                "node_kind": node_kind,
                "node_label": node_label,
                "descendant_taxa": clade_id,
                "supporting_tree_count": ""
                if supporting_tree_count is None
                else supporting_tree_count,
                "clade_frequency": clade_frequency,
                "support_percent": support_percent,
                "support_status": support_status,
                "explanation": explanation,
                "reference_branch_length": ""
                if clade.branch_length is None
                else clade.branch_length,
                "reference_root_depth": depths[id(clade)],
            }}
        )

    summary = {{
        "tree_count": tree_count,
        "shared_taxa": sorted(reference_taxa),
        "shared_taxon_count": len(reference_taxa),
        "internal_node_count": len(internal_clades),
        "supported_clade_count": supported_clade_count,
        "absent_clade_count": absent_clade_count,
        "unscored_clade_count": unscored_clade_count,
    }}
    return summary, rows

def matrix_rank(matrix, tolerance=1e-12):
    working = [list(map(float, row)) for row in matrix]
    row_count = len(working)
    column_count = len(working[0]) if working else 0
    rank = 0
    pivot_row = 0
    for pivot_column in range(column_count):
        candidate_row = max(
            range(pivot_row, row_count),
            key=lambda index: abs(working[index][pivot_column]),
            default=None,
        )
        if candidate_row is None:
            break
        pivot_value = working[candidate_row][pivot_column]
        if abs(pivot_value) <= tolerance:
            continue
        working[pivot_row], working[candidate_row] = (
            working[candidate_row],
            working[pivot_row],
        )
        pivot = working[pivot_row][pivot_column]
        working[pivot_row] = [value / pivot for value in working[pivot_row]]
        for row_index in range(row_count):
            if row_index == pivot_row:
                continue
            factor = working[row_index][pivot_column]
            if abs(factor) <= tolerance:
                continue
            working[row_index] = [
                row_value - factor * pivot_value
                for row_value, pivot_value in zip(
                    working[row_index], working[pivot_row], strict=True
                )
            ]
        rank += 1
        pivot_row += 1
        if pivot_row == row_count:
            break
    return rank

def invert_matrix(matrix):
    size = len(matrix)
    augmented = [
        [float(value) for value in row] + [1.0 if row_index == column_index else 0.0 for column_index in range(size)]
        for row_index, row in enumerate(matrix)
    ]
    for pivot_index in range(size):
        pivot_row = max(
            range(pivot_index, size),
            key=lambda row_index: abs(augmented[row_index][pivot_index]),
        )
        pivot_value = augmented[pivot_row][pivot_index]
        if abs(pivot_value) <= 1e-12:
            raise ValueError("matrix is singular")
        if pivot_row != pivot_index:
            augmented[pivot_index], augmented[pivot_row] = augmented[pivot_row], augmented[pivot_index]
        pivot_value = augmented[pivot_index][pivot_index]
        augmented[pivot_index] = [value / pivot_value for value in augmented[pivot_index]]
        for row_index in range(size):
            if row_index == pivot_index:
                continue
            factor = augmented[row_index][pivot_index]
            if abs(factor) <= 1e-15:
                continue
            augmented[row_index] = [
                row_value - factor * pivot_value
                for row_value, pivot_value in zip(
                    augmented[row_index], augmented[pivot_index], strict=True
                )
            ]
    return [row[size:] for row in augmented]

def matrix_infinity_norm(matrix):
    return max((sum(abs(value) for value in row) for row in matrix), default=0.0)

def symmetric_matrix_eigenvalues(matrix, tolerance=1e-15, max_iterations=10000):
    size = len(matrix)
    if size == 0:
        return []
    if size == 1:
        return [float(matrix[0][0])]
    working = [list(map(float, row)) for row in matrix]
    for _ in range(max_iterations):
        pivot_row = 0
        pivot_column = 1
        pivot_value = 0.0
        for row_index in range(size):
            for column_index in range(row_index + 1, size):
                candidate = abs(working[row_index][column_index])
                if candidate > pivot_value:
                    pivot_row = row_index
                    pivot_column = column_index
                    pivot_value = candidate
        if pivot_value <= tolerance:
            return [working[index][index] for index in range(size)]
        app = working[pivot_row][pivot_row]
        aqq = working[pivot_column][pivot_column]
        apq = working[pivot_row][pivot_column]
        tau = (aqq - app) / (2.0 * apq)
        tangent = (
            math.copysign(1.0, tau) / (abs(tau) + math.sqrt(1.0 + tau * tau))
            if abs(tau) > tolerance
            else 1.0
        )
        cosine = 1.0 / math.sqrt(1.0 + tangent * tangent)
        sine = tangent * cosine
        for index in range(size):
            if index in (pivot_row, pivot_column):
                continue
            left = working[index][pivot_row]
            right = working[index][pivot_column]
            working[index][pivot_row] = working[pivot_row][index] = cosine * left - sine * right
            working[index][pivot_column] = working[pivot_column][index] = sine * left + cosine * right
        working[pivot_row][pivot_row] = (
            cosine * cosine * app
            - 2.0 * sine * cosine * apq
            + sine * sine * aqq
        )
        working[pivot_column][pivot_column] = (
            sine * sine * app
            + 2.0 * sine * cosine * apq
            + cosine * cosine * aqq
        )
        working[pivot_row][pivot_column] = 0.0
        working[pivot_column][pivot_row] = 0.0
    raise ValueError("symmetric eigenvalue iteration did not converge")

def symmetric_matrix_condition_number(matrix, tolerance=1e-12):
    singular_values = sorted(
        abs(value)
        for value in symmetric_matrix_eigenvalues(matrix, tolerance=tolerance)
    )
    if not singular_values:
        return 0.0
    if singular_values[0] <= tolerance:
        return math.inf
    return singular_values[-1] / singular_values[0]

def matrix_log_determinant(matrix):
    size = len(matrix)
    working = [list(map(float, row)) for row in matrix]
    sign = 1.0
    log_abs_det = 0.0
    for pivot_index in range(size):
        pivot_row = max(
            range(pivot_index, size),
            key=lambda row_index: abs(working[row_index][pivot_index]),
        )
        pivot_value = working[pivot_row][pivot_index]
        if abs(pivot_value) <= 1e-12:
            raise ValueError("matrix determinant is zero")
        if pivot_row != pivot_index:
            working[pivot_index], working[pivot_row] = working[pivot_row], working[pivot_index]
            sign *= -1.0
        pivot_value = working[pivot_index][pivot_index]
        if pivot_value < 0:
            sign *= -1.0
        log_abs_det += math.log(abs(pivot_value))
        for row_index in range(pivot_index + 1, size):
            factor = working[row_index][pivot_index] / pivot_value
            if abs(factor) <= 1e-15:
                continue
            for column_index in range(pivot_index, size):
                working[row_index][column_index] -= factor * working[pivot_index][column_index]
    if sign <= 0:
        raise ValueError("matrix determinant is not positive")
    return log_abs_det

def ancestor_depths(tree):
    lookup = {{}}
    def walk(clade, depth, path):
        current = dict(path)
        current[id(clade)] = depth
        if clade.is_terminal():
            lookup[clade.name] = current
            return
        for child in clade.clades:
            walk(child, depth + (child.branch_length or 0.0), current)
    walk(tree.root, 0.0, {{}})
    return lookup

runner_path = Path(sys.argv[1])
case_path = Path(sys.argv[2])
output_root = Path(sys.argv[3])
output_root.mkdir(parents=True, exist_ok=True)
case_payload = json.loads(case_path.read_text(encoding="utf-8"))
execution_path = output_root / "reference-execution.json"

if not APE_AVAILABLE:
    write_json(
        execution_path,
        {{
            "status": "unavailable",
            "mismatch_reason": "ape_package_unavailable",
            "message": "ape is not installed in the fake reference environment",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": None,
        }},
    )
    raise SystemExit(0)

case_id = case_payload["case_id"]
if case_payload["operation"] in {{
    "read-tree-structure",
    "write-tree-structure",
    "root-tree-outgroup",
    "unroot-tree",
    "drop-tree-taxa",
    "keep-tree-taxa",
    "extract-tree-clade",
    "read-tree-set-structure",
    "write-tree-set-structure",
}}:
    try:
        if case_payload["operation"] in {{"read-tree-structure", "write-tree-structure", "root-tree-outgroup", "unroot-tree", "drop-tree-taxa", "keep-tree-taxa", "extract-tree-clade"}}:
            newick_path = output_root / "normalized-tree.nwk"
            if case_payload["operation"] == "root-tree-outgroup":
                outgroup_taxa = tuple(case_payload.get("outgroup_taxa", []))
                if outgroup_taxa == ("Z",):
                    raise ValueError("specified outgroup not in labels of the tree")
                if outgroup_taxa == ("B", "D"):
                    raise ValueError("the specified outgroup is not monophyletic")
                if outgroup_taxa == ("D",):
                    newick_text = "(((A:0.2,B:0.2):0.7,C:0.1):0,D:0.1);\\n"
                elif outgroup_taxa == ("C", "D"):
                    newick_text = "((A:0.2,B:0.2):0.7,(C:0.1,D:0.1):0);\\n"
                else:
                    newick_text = Path(case_payload["input_fixture"]).read_text(encoding="utf-8")
                newick_path.write_text(newick_text, encoding="utf-8")
                tree = Phylo.read(newick_path, "newick")
            elif case_payload["operation"] == "unroot-tree":
                if case_id == "unroot-tree-balanced-rooted":
                    newick_text = "(A:0.1,B:0.1,(C:0.2,D:0.2):0.3);\\n"
                elif case_id == "unroot-tree-rootable":
                    newick_text = "(A:0.2,B:0.2,(C:0.1,D:0.1):0.7);\\n"
                elif case_id == "unroot-tree-after-outgroup-rooting":
                    newick_text = "((A:0.2,B:0.2):0.7,C:0.1,D:0.1);\\n"
                else:
                    newick_text = Path(case_payload["input_fixture"]).read_text(encoding="utf-8")
                newick_path.write_text(newick_text, encoding="utf-8")
                tree = Phylo.read(newick_path, "newick")
            elif case_payload["operation"] == "drop-tree-taxa":
                if case_id == "drop-tip-rooted-single":
                    newick_text = "((A:0.1,B:0.1):0.2,C:0.3);\\n"
                    dropped_taxa = ["D"]
                    absent_requested_taxa = []
                elif case_id == "drop-tip-rooted-multiple":
                    newick_text = "(A:0.3,C:0.3);\\n"
                    dropped_taxa = ["B", "D"]
                    absent_requested_taxa = []
                elif case_id == "drop-tip-root-change-after-outgroup-rooting":
                    newick_text = "((A:0.2,B:0.2):0.7,C:0.1);\\n"
                    dropped_taxa = ["D"]
                    absent_requested_taxa = []
                elif case_id == "drop-tip-unrooted-three-tip":
                    newick_text = "(A:0.1,B:0.2,C:0.3);\\n"
                    dropped_taxa = ["D"]
                    absent_requested_taxa = []
                elif case_id == "drop-tip-unrooted-two-tip":
                    newick_text = "(A:0.1,B:0.2);\\n"
                    dropped_taxa = ["C", "D"]
                    absent_requested_taxa = []
                else:
                    newick_text = Path(case_payload["input_fixture"]).read_text(encoding="utf-8")
                    dropped_taxa = []
                    absent_requested_taxa = ["Z"]
                newick_path.write_text(newick_text, encoding="utf-8")
                tree = Phylo.read(newick_path, "newick")
            elif case_payload["operation"] == "keep-tree-taxa":
                if case_id == "keep-tip-rooted-selected-two":
                    newick_text = "(A:0.3,C:0.3);\\n"
                    requested_taxa = ["A", "C"]
                    dropped_taxa = ["B", "D"]
                elif case_id == "keep-tip-rooted-order-insensitive":
                    newick_text = "(A:0.3,C:0.3);\\n"
                    requested_taxa = ["A", "C"]
                    dropped_taxa = ["B", "D"]
                elif case_id == "keep-tip-root-change-after-outgroup-rooting":
                    newick_text = "((A:0.2,B:0.2):0.7,C:0.1);\\n"
                    requested_taxa = ["A", "B", "C"]
                    dropped_taxa = ["D"]
                elif case_id == "keep-tip-unrooted-three-tip":
                    newick_text = "(A:0.1,B:0.2,C:0.3);\\n"
                    requested_taxa = ["A", "B", "C"]
                    dropped_taxa = ["D"]
                elif case_id == "keep-tip-unrooted-two-tip":
                    newick_text = "(A:0.1,B:0.2);\\n"
                    requested_taxa = ["A", "B"]
                    dropped_taxa = ["C", "D"]
                else:
                    newick_text = Path(case_payload["input_fixture"]).read_text(encoding="utf-8")
                    requested_taxa = sorted(set(case_payload.get("requested_taxa", [])))
                    dropped_taxa = []
                newick_path.write_text(newick_text, encoding="utf-8")
                tree = Phylo.read(newick_path, "newick")
            elif case_payload["operation"] == "extract-tree-clade":
                if case_id == "extract-clade-tip-node-invalid":
                    raise ValueError("node number must be greater than the number of tips")
                elif case_id == "extract-clade-node-out-of-bounds":
                    raise IndexError("subscript out of bounds")
                elif case_id == "extract-clade-root":
                    newick_text = "((A:0.1,B:0.1)Mammals:0.2,(C:0.2,D:0.2)Birds:0.1)Root;\\n"
                    requested_node_id = 5
                    matched_node_id = 5
                    matched_node_name = "Root"
                elif case_id == "extract-clade-mammals":
                    newick_text = "(A:0.1,B:0.1)Mammals;\\n"
                    requested_node_id = 6
                    matched_node_id = 6
                    matched_node_name = "Mammals"
                else:
                    newick_text = "(C:0.2,D:0.2)Birds;\\n"
                    requested_node_id = 7
                    matched_node_id = 7
                    matched_node_name = "Birds"
                newick_path.write_text(newick_text, encoding="utf-8")
                tree = Phylo.read(newick_path, "newick")
            else:
                tree = Phylo.read(case_payload["input_fixture"], "newick")
            summary = {{
                "tree_count": 1,
                "tip_count": len(tree.get_terminals()),
                "internal_node_count": len(tree.get_nonterminals()),
                "edge_count": len(tree.get_terminals()) + len(tree.get_nonterminals()) - 1,
                "rooted": is_rooted_tree(tree),
                "tip_labels": [terminal.name for terminal in tree.get_terminals()],
                "branch_length_count": sum(
                    1 for clade in tree.find_clades(order="preorder")
                    if clade is not tree.root and clade.branch_length is not None
                ),
            }}
            if case_payload["operation"] == "drop-tree-taxa":
                summary["dropped_taxa"] = dropped_taxa
                summary["absent_requested_taxa"] = absent_requested_taxa
            if case_payload["operation"] == "keep-tree-taxa":
                summary["requested_taxa"] = requested_taxa
                summary["dropped_taxa"] = dropped_taxa
            if case_payload["operation"] == "extract-tree-clade":
                summary["requested_node_id"] = requested_node_id
                summary["matched_node_id"] = matched_node_id
                summary["matched_node_name"] = matched_node_name
            rows = clade_rows(tree, "")
            summary.update(SUMMARY_OVERRIDES)
            summary_path = output_root / "summary.json"
            clades_path = output_root / "clades.tsv"
            write_json(summary_path, summary)
            write_tsv(clades_path, rows)
            if case_payload["operation"] in {{"root-tree-outgroup", "unroot-tree", "drop-tree-taxa", "keep-tree-taxa", "extract-tree-clade"}}:
                pass
            elif case_id == "read-tree-quoted-taxon-labels":
                newick_path.write_text(
                    "('Homo_sapiens':0.1,'Mus_musculus':0.2,'A.B-1':0.3);\\n",
                    encoding="utf-8",
                )
            elif case_id in NORMALIZED_TREE_OVERRIDES:
                newick_path.write_text(
                    NORMALIZED_TREE_OVERRIDES[case_id],
                    encoding="utf-8",
                )
            else:
                newick_path.write_text(
                    Path(case_payload["input_fixture"]).read_text(encoding="utf-8"),
                    encoding="utf-8",
                )
            outputs = {{
                "summary_json": str(summary_path),
                "clades": str(clades_path),
                "normalized_tree": str(newick_path),
            }}
        else:
            trees = list(Phylo.parse(case_payload["input_fixture"], "newick"))
            if not trees:
                raise ValueError("tree set contains no trees")
            rows = []
            for index, tree in enumerate(trees, start=1):
                rows.extend(clade_rows(tree, index))
            shared_tip_labels = sorted(set(t.name for t in trees[0].get_terminals()))
            summary = {{
                "tree_count": len(trees),
                "source_format": "newick",
                "tree_indices": list(range(1, len(trees) + 1)),
                "shared_tip_labels": shared_tip_labels,
                "unique_tip_label_count": len(shared_tip_labels),
            }}
            summary.update(SUMMARY_OVERRIDES)
            summary_path = output_root / "summary.json"
            clades_path = output_root / "clades.tsv"
            tree_set_path = output_root / "normalized-tree-set.nwk"
            write_json(summary_path, summary)
            write_tsv(clades_path, rows)
            tree_set_path.write_text(
                Path(case_payload["input_fixture"]).read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            outputs = {{
                "summary_json": str(summary_path),
                "clades": str(clades_path),
                "normalized_tree_set": str(tree_set_path),
            }}
    except Exception as error:
        write_json(
            execution_path,
            {{
                "status": "failed",
                "mismatch_reason": "reference_execution_failed",
                "error_type": "TreeRootingError" if case_payload["operation"] == "root-tree-outgroup" else "TreeParseError",
                "message": str(error),
                "case_id": case_payload["case_id"],
                "function_name": case_payload["function_name"],
                "input_fixture": case_payload["input_fixture"],
                "r_version": "4.6.0",
                "ape_version": "5.0.0",
            }},
        )
        raise SystemExit(0)

    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": outputs,
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "get-tree-mrca":
    try:
        if case_id == "get-mrca-missing-tip":
            raise ValueError("missing value where TRUE/FALSE needed")
        if case_id == "get-mrca-balanced-two-tip":
            summary = {{
                "requested_taxa": ["A", "B"],
                "unique_requested_taxa": ["A", "B"],
                "duplicate_requested_taxa": [],
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B"],
                "matched_extra_taxa": [],
                "matched_tip_count": 2,
                "is_root": False,
            }}
        elif case_id == "get-mrca-balanced-full-tip-set":
            summary = {{
                "requested_taxa": ["A", "B", "C", "D"],
                "unique_requested_taxa": ["A", "B", "C", "D"],
                "duplicate_requested_taxa": [],
                "matched_node_id": 5,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C", "D"],
                "matched_extra_taxa": [],
                "matched_tip_count": 4,
                "is_root": True,
            }}
        elif case_id == "get-mrca-balanced-duplicate-request":
            summary = {{
                "requested_taxa": ["A", "A", "B"],
                "unique_requested_taxa": ["A", "B"],
                "duplicate_requested_taxa": ["A"],
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B"],
                "matched_extra_taxa": [],
                "matched_tip_count": 2,
                "is_root": False,
            }}
        elif case_id == "get-mrca-pectinate-many-tip":
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C"],
                "matched_extra_taxa": [],
                "matched_tip_count": 3,
                "is_root": False,
            }}
        elif case_id == "get-mrca-rooted-polytomy":
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C"],
                "matched_extra_taxa": [],
                "matched_tip_count": 3,
                "is_root": False,
            }}
        else:
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C"],
                "matched_extra_taxa": [],
                "matched_tip_count": 3,
                "is_root": False,
            }}
    except Exception as error:
        write_json(
            execution_path,
            {{
                "status": "failed",
                "mismatch_reason": "reference_execution_failed",
                "error_type": "TreeMrcaError",
                "message": str(error),
                "case_id": case_payload["case_id"],
                "function_name": case_payload["function_name"],
                "input_fixture": case_payload["input_fixture"],
                "r_version": "4.6.0",
                "ape_version": "5.0.0",
            }},
        )
        raise SystemExit(0)

    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    write_json(summary_path, summary)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{"summary_json": str(summary_path)}},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "assess-tree-monophyly":
    try:
        if case_id == "is-monophyletic-all-missing-rerooted":
            raise ValueError("specified outgroup not in labels of the tree")
        if case_id == "is-monophyletic-rooted-two-tip":
            summary = {{
                "requested_taxa": ["A", "B"],
                "unique_requested_taxa": ["A", "B"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B"],
                "reroot": False,
                "rooted": True,
                "monophyletic": True,
                "complementary_clade_used": False,
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B"],
                "matched_extra_taxa": [],
                "matched_tip_count": 2,
                "is_root": False,
            }}
        elif case_id == "is-monophyletic-rooted-three-tip-reroot-false":
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B", "C"],
                "reroot": False,
                "rooted": True,
                "monophyletic": False,
                "complementary_clade_used": False,
                "matched_node_id": 5,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C", "D"],
                "matched_extra_taxa": ["D"],
                "matched_tip_count": 4,
                "is_root": True,
            }}
        elif case_id == "is-monophyletic-rooted-three-tip-reroot-true":
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B", "C"],
                "reroot": True,
                "rooted": True,
                "monophyletic": True,
                "complementary_clade_used": True,
                "matched_node_id": 5,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C", "D"],
                "matched_extra_taxa": ["D"],
                "matched_tip_count": 4,
                "is_root": True,
            }}
        elif case_id == "is-monophyletic-rooted-full-tip-set":
            summary = {{
                "requested_taxa": ["A", "B", "C", "D"],
                "unique_requested_taxa": ["A", "B", "C", "D"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B", "C", "D"],
                "reroot": False,
                "rooted": True,
                "monophyletic": True,
                "complementary_clade_used": False,
                "matched_node_id": 5,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C", "D"],
                "matched_extra_taxa": [],
                "matched_tip_count": 4,
                "is_root": True,
            }}
        elif case_id == "is-monophyletic-rooted-mixed-missing":
            summary = {{
                "requested_taxa": ["A", "Z"],
                "unique_requested_taxa": ["A", "Z"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": ["Z"],
                "present_requested_taxa": ["A"],
                "reroot": False,
                "rooted": True,
                "monophyletic": True,
                "complementary_clade_used": False,
                "matched_node_id": 1,
                "matched_node_name": "A",
                "matched_taxa": ["A"],
                "matched_extra_taxa": [],
                "matched_tip_count": 1,
                "is_root": False,
            }}
        elif case_id == "is-monophyletic-unrooted-two-tip":
            summary = {{
                "requested_taxa": ["A", "B"],
                "unique_requested_taxa": ["A", "B"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B"],
                "reroot": True,
                "rooted": False,
                "monophyletic": False,
                "complementary_clade_used": False,
                "matched_node_id": 5,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C", "D"],
                "matched_extra_taxa": ["C", "D"],
                "matched_tip_count": 4,
                "is_root": True,
            }}
        elif case_id == "is-monophyletic-unrooted-three-tip":
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B", "C"],
                "reroot": True,
                "rooted": False,
                "monophyletic": True,
                "complementary_clade_used": True,
                "matched_node_id": 5,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C", "D"],
                "matched_extra_taxa": ["D"],
                "matched_tip_count": 4,
                "is_root": True,
            }}
        elif case_id == "is-monophyletic-after-outgroup-rooting":
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B", "C"],
                "reroot": False,
                "rooted": True,
                "monophyletic": True,
                "complementary_clade_used": False,
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C"],
                "matched_extra_taxa": [],
                "matched_tip_count": 3,
                "is_root": False,
            }}
        else:
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B", "C"],
                "reroot": False,
                "rooted": True,
                "monophyletic": True,
                "complementary_clade_used": False,
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C"],
                "matched_extra_taxa": [],
                "matched_tip_count": 3,
                "is_root": False,
            }}
    except Exception as error:
        write_json(
            execution_path,
            {{
                "status": "failed",
                "mismatch_reason": "reference_execution_failed",
                "error_type": "TreeMonophylyError",
                "message": str(error),
                "case_id": case_payload["case_id"],
                "function_name": case_payload["function_name"],
                "input_fixture": case_payload["input_fixture"],
                "r_version": "4.6.0",
                "ape_version": "5.0.0",
            }},
        )
        raise SystemExit(0)

    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    write_json(summary_path, summary)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{"summary_json": str(summary_path)}},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-consensus":
    if case_id == "consensus-mismatched-taxon-set":
        write_json(
            execution_path,
            {{
                "status": "failed",
                "mismatch_reason": "reference_execution_failed",
                "error_type": "ConsensusTreeError",
                "message": "consensus requires all trees to share the exact same taxon set",
                "case_id": case_payload["case_id"],
                "function_name": case_payload["function_name"],
                "input_fixture": case_payload["input_fixture"],
                "r_version": "4.6.0",
                "ape_version": "5.0.0",
            }},
        )
        raise SystemExit(0)

    if case_id == "consensus-majority-conflicting-four-taxon":
        newick_text = "((A:0.1,B:0.1)66.6666666666667:0.2,(C:0.1,D:0.1)66.6666666666667:0.2);\\n"
        rows = [
            {{"clade": "A|B", "tree_count": 2, "frequency": 2 / 3}},
            {{"clade": "A|C", "tree_count": 1, "frequency": 1 / 3}},
            {{"clade": "B|D", "tree_count": 1, "frequency": 1 / 3}},
            {{"clade": "C|D", "tree_count": 2, "frequency": 2 / 3}},
        ]
        summary = {{
            "tree_count": 3,
            "shared_taxa": ["A", "B", "C", "D"],
            "shared_taxon_count": 4,
            "tip_count": 4,
            "rooted": False,
            "consensus_method": "majority-rule",
            "consensus_threshold": 0.5,
            "included_clade_count": 1,
            "clade_frequency_count": 4,
        }}
    elif case_id == "consensus-strict-conflicting-four-taxon":
        newick_text = "(A:0.1,B:0.1,C:0.1,D:0.1);\\n"
        rows = [
            {{"clade": "A|B", "tree_count": 2, "frequency": 2 / 3}},
            {{"clade": "A|C", "tree_count": 1, "frequency": 1 / 3}},
            {{"clade": "B|D", "tree_count": 1, "frequency": 1 / 3}},
            {{"clade": "C|D", "tree_count": 2, "frequency": 2 / 3}},
        ]
        summary = {{
            "tree_count": 3,
            "shared_taxa": ["A", "B", "C", "D"],
            "shared_taxon_count": 4,
            "tip_count": 4,
            "rooted": False,
            "consensus_method": "strict",
            "consensus_threshold": 1.0,
            "included_clade_count": 0,
            "clade_frequency_count": 4,
        }}
    elif case_id == "consensus-majority-posterior-six-taxon":
        newick_text = "(A:1,B:1,(C:1,D:1)60:1.66666666666667,(E:1,F:1)60:1.66666666666667);\\n"
        rows = [
            {{"clade": "A|B", "tree_count": 2, "frequency": 0.4}},
            {{"clade": "A|B|C|D", "tree_count": 2, "frequency": 0.4}},
            {{"clade": "A|B|D|E", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "A|B|E|F", "tree_count": 2, "frequency": 0.4}},
            {{"clade": "A|C", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "A|D", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "A|E", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "B|D", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "B|E", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "B|F", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "C|D", "tree_count": 3, "frequency": 0.6}},
            {{"clade": "C|F", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "E|F", "tree_count": 3, "frequency": 0.6}},
        ]
        summary = {{
            "tree_count": 5,
            "shared_taxa": ["A", "B", "C", "D", "E", "F"],
            "shared_taxon_count": 6,
            "tip_count": 6,
            "rooted": False,
            "consensus_method": "majority-rule",
            "consensus_threshold": 0.5,
            "included_clade_count": 2,
            "clade_frequency_count": 13,
        }}
    else:
        raise ValueError(f"unsupported fake consensus case: {{case_id}}")

    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / "clade-frequencies.tsv"
    newick_path = output_root / "normalized-tree.nwk"
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    newick_path.write_text(newick_text, encoding="utf-8")
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "clade_frequencies": str(rows_path),
                "normalized_tree": str(newick_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-clade-support":
    try:
        reference_tree = Phylo.read(case_payload["reference_tree_path"], "newick")
        comparison_trees = list(Phylo.parse(case_payload["input_fixture"], "newick"))
        if not comparison_trees:
            raise ValueError("tree set contains no trees")
        summary, rows = prop_clades_rows(reference_tree, comparison_trees)
    except Exception as error:
        write_json(
            execution_path,
            {{
                "status": "failed",
                "mismatch_reason": "reference_execution_failed",
                "error_type": "PropCladesError",
                "message": str(error),
                "case_id": case_payload["case_id"],
                "function_name": case_payload["function_name"],
                "input_fixture": case_payload["input_fixture"],
                "r_version": "4.6.0",
                "ape_version": "5.0.0",
            }},
        )
        raise SystemExit(0)

    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / "support-table.tsv"
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "support_table": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-tip-distance":
    tree = Phylo.read(case_payload["input_fixture"], "newick")
    tip_labels = [terminal.name for terminal in tree.get_terminals()]
    rows = []
    for left in tip_labels:
        for right in tip_labels:
            rows.append(
                {{
                    "left_identifier": left,
                    "right_identifier": right,
                    "distance": tree.distance(left, right),
                }}
            )
    summary = {{
        "tip_count": len(tip_labels),
        "rooted": is_rooted_tree(tree),
        "tip_labels": tip_labels,
        "pair_count": len(rows),
        "diagonal_zero": True,
        "symmetric": True,
        "complete_branch_lengths": True,
        "missing_branch_length_policy": "error",
    }}
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / "tip-distance-long.tsv"
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "tip_distance_long": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-topology-distance":
    trees = list(Phylo.parse(case_payload["input_fixture"], "newick"))
    if len(trees) != 2:
        raise ValueError("ape topology-distance parity fixtures must contain exactly two trees")
    left_tree, right_tree = trees
    rf_mode = case_payload.get("rf_mode", "rooted")
    (
        rows,
        left_split_count,
        right_split_count,
        shared_split_count,
        left_only_split_count,
        right_only_split_count,
    ) = topology_distance_rows(left_tree, right_tree, rf_mode)
    summary = {{
        "tip_count": len(descendant_taxa(left_tree.root)),
        "shared_taxa": descendant_taxa(left_tree.root),
        "left_only_taxa": [],
        "right_only_taxa": [],
        "taxon_overlap_policy": "require-identical",
        "rf_mode": rf_mode,
        "rooted_left": is_rooted_tree(left_tree),
        "rooted_right": is_rooted_tree(right_tree),
        "polytomy_present_left": any(len(clade.clades) > 2 for clade in left_tree.find_clades()),
        "polytomy_present_right": any(len(clade.clades) > 2 for clade in right_tree.find_clades()),
        "left_split_count": left_split_count,
        "right_split_count": right_split_count,
        "shared_split_count": shared_split_count,
        "left_only_split_count": left_only_split_count,
        "right_only_split_count": right_only_split_count,
        "robinson_foulds_distance": left_only_split_count + right_only_split_count,
        "normalized_robinson_foulds": (
            0.0
            if (left_split_count + right_split_count) == 0
            else (left_only_split_count + right_only_split_count)
            / (left_split_count + right_split_count)
        ),
        "topology_equal": (left_only_split_count + right_only_split_count) == 0,
    }}
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / "split-table.tsv"
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "split_table": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-brownian-covariance":
    tree = Phylo.read(case_payload["input_fixture"], "newick")
    tip_labels = [terminal.name for terminal in tree.get_terminals()]
    depth_lookup = ancestor_depths(tree)
    covariance = []
    for left in tip_labels:
        row = []
        left_path = depth_lookup[left]
        for right in tip_labels:
            right_path = depth_lookup[right]
            shared_depth = max(
                left_path[node_id]
                for node_id in set(left_path) & set(right_path)
            )
            row.append(shared_depth)
        covariance.append(row)
    root_depths = [covariance[index][index] for index in range(len(tip_labels))]
    branch_lengths = [
        clade.branch_length
        for clade in tree.find_clades(order="preorder")
        if clade is not tree.root and clade.branch_length is not None
    ]
    covariance_rank = matrix_rank(covariance)
    singular = covariance_rank < len(covariance)
    try:
        raw_log_determinant = matrix_log_determinant(covariance)
        positive_definite = True
    except ValueError:
        raw_log_determinant = None
        positive_definite = False
    condition_number = None if singular else symmetric_matrix_condition_number(covariance)
    near_singular = singular or (
        condition_number is not None and condition_number >= 1e12
    )
    matrix_rows = []
    long_rows = []
    for row_index, left in enumerate(tip_labels):
        matrix_row = {{"taxon": left}}
        for column_index, right in enumerate(tip_labels):
            value = covariance[row_index][column_index]
            matrix_row[right] = value
            long_rows.append(
                {{
                    "left_taxon": left,
                    "right_taxon": right,
                    "shared_ancestry_covariance": value,
                }}
            )
        matrix_rows.append(matrix_row)
    summary = {{
        "tip_count": len(tip_labels),
        "rooted": is_rooted_tree(tree),
        "tip_labels": tip_labels,
        "pair_count": len(long_rows),
        "tree_is_ultrametric": max(root_depths) - min(root_depths) <= 1e-12,
        "minimum_root_to_tip_depth": min(root_depths),
        "maximum_root_to_tip_depth": max(root_depths),
        "minimum_branch_length": min(branch_lengths),
        "maximum_branch_length": max(branch_lengths),
        "matrix_dimension": len(covariance),
        "matrix_rank": covariance_rank,
        "singular": singular,
        "near_singular": near_singular,
        "positive_definite": positive_definite,
        "condition_number": condition_number,
        "raw_log_determinant": raw_log_determinant,
    }}
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    matrix_path = output_root / "covariance-matrix.tsv"
    rows_path = output_root / "covariance-long.tsv"
    write_json(summary_path, summary)
    write_tsv(matrix_path, matrix_rows)
    write_tsv(rows_path, long_rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "covariance_matrix": str(matrix_path),
                "covariance_long": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-node-depth":
    tree = Phylo.read(case_payload["input_fixture"], "newick")
    rows = node_depth_rows(tree)
    tip_depths = [
        row["branch_length_depth"] for row in rows if row["node_kind"] == "tip"
    ]
    internal_depths = [
        row["branch_length_depth"] for row in rows if row["node_kind"] != "tip"
    ]
    summary = {{
        "node_count": len(rows),
        "tip_count": len(tip_depths),
        "internal_node_count": len(internal_depths),
        "rooted": is_rooted_tree(tree),
        "tip_labels": [terminal.name for terminal in tree.get_terminals()],
        "tree_is_ultrametric": (
            abs(max(tip_depths) - min(tip_depths)) <= 1e-12 if tip_depths else True
        ),
        "zero_branch_length_count": sum(
            1
            for clade in tree.find_clades(order="preorder")
            if clade is not tree.root and clade.branch_length == 0.0
        ),
        "minimum_tip_depth": min(tip_depths),
        "maximum_tip_depth": max(tip_depths),
        "minimum_internal_depth": min(internal_depths),
        "maximum_internal_depth": max(internal_depths),
    }}
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / "node-depths.tsv"
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "node_depths": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-branching-times":
    tree = Phylo.read(case_payload["input_fixture"], "newick")
    depth_lookup = node_depth_lookup(tree)
    tip_depths = [
        depth_lookup[id(terminal)]
        for terminal in tree.get_terminals()
    ]
    rows = branching_time_rows(tree)
    summary = {{
        "internal_node_count": len(rows),
        "rooted": is_rooted_tree(tree),
        "tip_labels": [terminal.name for terminal in tree.get_terminals()],
        "tree_is_ultrametric": (
            abs(max(tip_depths) - min(tip_depths)) <= 1e-12 if tip_depths else True
        ),
        "root_age": max(tip_depths),
        "zero_branch_length_count": sum(
            1
            for clade in tree.find_clades(order="preorder")
            if clade is not tree.root and clade.branch_length == 0.0
        ),
        "minimum_tip_depth": min(tip_depths),
        "maximum_tip_depth": max(tip_depths),
        "max_tip_depth_deviation": max(tip_depths) - min(tip_depths),
        "tolerance": 1e-12,
    }}
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / "branching-times.tsv"
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "branching_times": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-ultrametricity":
    tree = Phylo.read(case_payload["input_fixture"], "newick")
    rows = ultrametric_rows(tree)
    tip_depths = [row["root_to_tip_depth"] for row in rows]
    minimum_tip_depth = min(tip_depths)
    maximum_tip_depth = max(tip_depths)
    mean_tip_depth = sum(tip_depths) / len(tip_depths)
    max_tip_depth_deviation = maximum_tip_depth - minimum_tip_depth
    option = int(case_payload.get("ultrametric_option", 1))
    if option == 1:
        if math.isclose(maximum_tip_depth, 0.0, abs_tol=1e-15):
            criterion_value = 0.0 if math.isclose(max_tip_depth_deviation, 0.0, abs_tol=1e-15) else math.inf
        else:
            criterion_value = max_tip_depth_deviation / maximum_tip_depth
        criterion_name = "scaled-range"
    else:
        criterion_name = "variance"
        if len(tip_depths) <= 1:
            criterion_value = 0.0
        else:
            criterion_value = sum((depth - mean_tip_depth) ** 2 for depth in tip_depths) / (len(tip_depths) - 1)
    offending_taxa = sorted(
        {{
            row["tip_label"]
            for row in rows
            if math.isclose(row["root_to_tip_depth"], minimum_tip_depth, abs_tol=1e-12)
            or math.isclose(row["root_to_tip_depth"], maximum_tip_depth, abs_tol=1e-12)
        }}
    )
    if math.isclose(max_tip_depth_deviation, 0.0, abs_tol=1e-12):
        offending_taxa = []
    rows = [
        {{
            **row,
            "deviation_from_mean_depth": abs(row["root_to_tip_depth"] - mean_tip_depth),
            "deviation_from_min_depth": row["root_to_tip_depth"] - minimum_tip_depth,
            "deviation_from_max_depth": maximum_tip_depth - row["root_to_tip_depth"],
            "is_offending_taxon": row["tip_label"] in offending_taxa,
        }}
        for row in rows
    ]
    summary = {{
        "tip_count": len(rows),
        "rooted": is_rooted_tree(tree),
        "tip_labels": [terminal.name for terminal in tree.get_terminals()],
        "ultrametric": criterion_value <= case_payload["tolerance"],
        "criterion_name": criterion_name,
        "criterion_value": criterion_value,
        "tolerance": case_payload["tolerance"],
        "option": option,
        "minimum_tip_depth": minimum_tip_depth,
        "maximum_tip_depth": maximum_tip_depth,
        "mean_tip_depth": mean_tip_depth,
        "max_tip_depth_deviation": max_tip_depth_deviation,
        "root_age": maximum_tip_depth,
        "offending_taxa": offending_taxa,
    }}
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / "ultrametric-diagnostics.tsv"
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "ultrametric_diagnostics": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_id in TABULAR_CASES:
    payload = TABULAR_CASES[case_id]
    summary = dict(payload["summary"])
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / payload["rows_name"]
    write_json(summary_path, summary)
    write_tsv(rows_path, payload["rows"])
    outputs = {{"summary_json": str(summary_path)}}
    if payload["rows_name"] == "base-frequency.tsv":
        outputs["base_frequency"] = str(rows_path)
    elif payload["rows_name"] == "distance-matrix.tsv":
        outputs["distance_matrix"] = str(rows_path)
    else:
        outputs["translation"] = str(rows_path)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": outputs,
        }},
    )
    raise SystemExit(0)

write_json(
    execution_path,
    {{
        "status": "failed",
        "mismatch_reason": "unsupported_operation",
        "message": f"unsupported ape parity operation: {{case_payload['operation']}}",
        "case_id": case_payload["case_id"],
        "function_name": case_payload["function_name"],
        "input_fixture": case_payload["input_fixture"],
        "r_version": "4.6.0",
        "ape_version": "5.0.0",
    }},
)
""",
    )
