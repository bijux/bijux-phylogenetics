from __future__ import annotations

from pathlib import Path

MEBIBYTE = 1024 * 1024
MINIMUM_MEMORY_MIB = 1024
MINIMUM_SCRATCH_MIB = 256
MINIMUM_WALLCLOCK_MINUTES = 20


def directory_bytes(path: Path) -> int:
    if not path.is_dir():
        return 0
    return sum(
        candidate.stat().st_size for candidate in path.rglob("*") if candidate.is_file()
    )


def format_float(value: float) -> str:
    text = f"{value:.2f}"
    return text.rstrip("0").rstrip(".") if "." in text else text


def format_slurm_time(minutes: int) -> str:
    hours, remaining_minutes = divmod(minutes, 60)
    return f"{hours:02d}:{remaining_minutes:02d}:00"


def round_up(value: int, *, quantum: int) -> int:
    return ((value + quantum - 1) // quantum) * quantum


def write_tsv(path: Path, rows: list[list[str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )
    return path
