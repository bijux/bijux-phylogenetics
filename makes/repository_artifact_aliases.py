#!/usr/bin/env python3
"""Bootstrap repository artifact alias materialization from the dev package."""

from __future__ import annotations

from pathlib import Path
import sys


def _dev_src_path() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "packages"
        / "bijux-phylogenetics-dev"
        / "src"
    )


def main() -> int:
    dev_src_path = _dev_src_path()
    sys.path.insert(0, str(dev_src_path))
    from bijux_phylogenetics_dev.workspace.artifact_aliases import main as package_main

    return package_main(sys.argv[1:])


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
