from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

REPO_ROOT = Path(__file__).resolve().parents[4]
BUNDLE_ROOT = Path(__file__).resolve().parent
STUDY_ID = "primate-longevity-signal"
EVIDENCE_ID = "evidence-002"
COMPARISON_MODE = "direct_parity"
PRIMARY_OUTPUTS = ['evidence-book/studies/primate-longevity-signal/evidence-002/workbook-loading-parity.json']
BUILD_SCRIPT = 'evidence-book/studies/primate-longevity-signal/build_evidence.py'

def main() -> None:
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
    print(json.dumps(payload, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()
