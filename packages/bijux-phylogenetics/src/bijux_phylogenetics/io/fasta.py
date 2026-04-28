from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.core.alignment import AlignmentSummary


def summarise_fasta(path: Path) -> AlignmentSummary:
    """Summarise a FASTA alignment without loading a heavy dependency."""
    if not path.exists():
        raise FileNotFoundError(f"alignment file not found: {path}")

    lengths: list[int] = []
    current: list[str] = []
    count = 0

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    lengths.append(len("".join(current)))
                    current = []
                count += 1
                continue
            current.append(line)

    if current:
        lengths.append(len("".join(current)))

    if count == 0:
        raise ValueError(f"alignment contains no FASTA records: {path}")

    return AlignmentSummary(
        path=path,
        sequence_count=count,
        min_sequence_length=min(lengths),
        max_sequence_length=max(lengths),
    )

