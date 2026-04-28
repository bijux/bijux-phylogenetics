from __future__ import annotations

from pathlib import Path


def load_nexus(path: Path) -> None:
    """Reserved NEXUS loader entrypoint."""
    raise NotImplementedError(f"NEXUS parsing is not implemented yet for {path}")

