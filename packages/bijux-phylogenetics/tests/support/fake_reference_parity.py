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
import json
import sys
from pathlib import Path

CASE_SUMMARIES = {{
    "read-tree-example-rooted": {{
        "tip_count": 4,
        "internal_node_count": 3,
        "edge_count": 6,
        "rooted": True,
        "tip_labels": ["A", "B", "C", "D"],
        "branch_length_count": 6,
    }},
    "read-tree-example-unrooted": {{
        "tip_count": 4,
        "internal_node_count": 1,
        "edge_count": 4,
        "rooted": False,
        "tip_labels": ["A", "B", "C", "D"],
        "branch_length_count": 4,
    }},
}}

SUMMARY_OVERRIDES = {summary_payload}
APE_AVAILABLE = {str(ape_available)}

def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\\n", encoding="utf-8")

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

summary = dict(CASE_SUMMARIES[case_payload["case_id"]])
summary.update(SUMMARY_OVERRIDES)
tip_rows = [{{"position": index, "label": label}} for index, label in enumerate(summary["tip_labels"], start=1)]
summary_path = output_root / "summary.json"
tips_path = output_root / "tips.tsv"
newick_path = output_root / "normalized-tree.nwk"
write_json(summary_path, summary)
tips_path.write_text(
    "position\\tlabel\\n"
    + "".join(f"{{row['position']}}\\t{{row['label']}}\\n" for row in tip_rows),
    encoding="utf-8",
)
newick_path.write_text(Path(case_payload["input_fixture"]).read_text(encoding="utf-8"), encoding="utf-8")
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
""",
    )
