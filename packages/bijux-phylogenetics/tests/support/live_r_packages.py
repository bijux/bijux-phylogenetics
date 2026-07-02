from __future__ import annotations

import os
from pathlib import Path
import subprocess


def r_package_available(rscript: str, package_name: str) -> bool:
    repository_root = Path(__file__).resolve().parents[3]
    environment = dict(os.environ)
    r_library = repository_root / "artifacts" / "r-lib"
    if r_library.is_dir():
        environment["R_LIBS_USER"] = str(r_library)
    result = subprocess.run(
        [
            rscript,
            "-e",
            f"cat(requireNamespace('{package_name}', quietly=TRUE), '\\n')",
        ],
        capture_output=True,
        check=False,
        cwd=repository_root,
        env=environment,
        text=True,
    )
    return result.returncode == 0 and result.stdout.strip() == "TRUE"
