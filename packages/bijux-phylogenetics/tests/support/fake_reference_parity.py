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

TREE_SUMMARIES = {{
    "read-tree-balanced-rooted-ultrametric": {{
        "tip_count": 4,
        "internal_node_count": 3,
        "edge_count": 6,
        "rooted": True,
        "tip_labels": ["A", "B", "C", "D"],
        "branch_length_count": 6,
    }},
    "read-tree-unrooted-branch-length": {{
        "tip_count": 4,
        "internal_node_count": 1,
        "edge_count": 4,
        "rooted": False,
        "tip_labels": ["A", "B", "C", "D"],
        "branch_length_count": 4,
    }},
    "read-tree-quoted-taxon-labels": {{
        "tip_count": 3,
        "internal_node_count": 1,
        "edge_count": 3,
        "rooted": False,
        "tip_labels": ["'Homo sapiens'", "'Mus musculus'", "'A.B-1'"],
        "branch_length_count": 3,
    }},
}}

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
if case_id in TREE_SUMMARIES:
    summary = dict(TREE_SUMMARIES[case_id])
    summary.update(SUMMARY_OVERRIDES)
    tip_rows = [{{"position": index, "label": label}} for index, label in enumerate(summary["tip_labels"], start=1)]
    summary_path = output_root / "summary.json"
    tips_path = output_root / "tips.tsv"
    newick_path = output_root / "normalized-tree.nwk"
    write_json(summary_path, summary)
    write_tsv(tips_path, tip_rows)
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
                "tip_table": str(tips_path),
                "normalized_tree": str(newick_path),
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
