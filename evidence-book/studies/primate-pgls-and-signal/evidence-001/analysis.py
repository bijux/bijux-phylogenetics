from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

REPO_ROOT = Path(__file__).resolve().parents[4]
BUNDLE_ROOT = Path(__file__).resolve().parent
RESULTS_ROOT = BUNDLE_ROOT / 'results'
STUDY_ID = "primate-pgls-and-signal"
EVIDENCE_ID = "evidence-001"
COMPARISON_MODE = "direct_parity"
PRIMARY_OUTPUTS = ['evidence-book/studies/primate-pgls-and-signal/evidence-001/rdata-reload-semantics.json', 'evidence-book/studies/primate-pgls-and-signal/evidence-001/scalar-parity-table.json', 'evidence-book/studies/primate-pgls-and-signal/evidence-001/scalar-parity-table.md']
BUILD_SCRIPT = 'evidence-book/studies/primate-pgls-and-signal/build_evidence.py'

def main() -> None:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    payload = {
        "study_id": STUDY_ID,
        "evidence_id": EVIDENCE_ID,
        "comparison_mode": COMPARISON_MODE,
        "build_script": BUILD_SCRIPT,
    }
    if BUILD_SCRIPT is not None:
        subprocess.run(
            [sys.executable, str(REPO_ROOT / BUILD_SCRIPT)],
            cwd=str(REPO_ROOT),
            check=True,
        )
        payload["execution_mode"] = "study_build_wrapper"
    else:
        payload["execution_mode"] = "bundle_contract_only"
    payload["primary_outputs"] = [
        path for path in PRIMARY_OUTPUTS if (REPO_ROOT / path).is_file()
    ]
    output_path = RESULTS_ROOT / 'analysis-run.json'
    output_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + '\n',
        encoding='utf-8',
    )
    print(json.dumps(payload, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
