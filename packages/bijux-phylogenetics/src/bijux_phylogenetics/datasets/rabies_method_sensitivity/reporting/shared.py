from __future__ import annotations

import os
from pathlib import Path


def _relative_bundle_path(base_path: Path, value: Path) -> str:
    return Path(os.path.relpath(value, start=base_path.parent)).as_posix()


def _format_float(value: float) -> str:
    return format(value, ".12g")
