from __future__ import annotations

from pathlib import Path

from tests.support.fake_external_engines import write_executable


def fake_ape_rscript(
    path: Path,
    *,
    ape_available: bool = True,
    summary_overrides: dict[str, object] | None = None,
) -> Path:
    summary_payload = repr(summary_overrides or {})
    return write_executable(
        path,
        f"""#!/usr/bin/env python3
import csv
import json
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
    "read-tree-set-structure",
    "write-tree-set-structure",
}}:
    try:
        if case_payload["operation"] in {{"read-tree-structure", "write-tree-structure"}}:
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
            rows = clade_rows(tree, "")
            summary.update(SUMMARY_OVERRIDES)
            summary_path = output_root / "summary.json"
            clades_path = output_root / "clades.tsv"
            newick_path = output_root / "normalized-tree.nwk"
            write_json(summary_path, summary)
            write_tsv(clades_path, rows)
            if case_id == "read-tree-quoted-taxon-labels":
                newick_path.write_text(
                    "('Homo_sapiens':0.1,'Mus_musculus':0.2,'A.B-1':0.3);\\n",
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
            write_json(summary_path, summary)
            write_tsv(clades_path, rows)
            outputs = {{
                "summary_json": str(summary_path),
                "clades": str(clades_path),
            }}
    except Exception as error:
        write_json(
            execution_path,
            {{
                "status": "failed",
                "mismatch_reason": "reference_execution_failed",
                "error_type": "TreeParseError",
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
