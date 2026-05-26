from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.engines import run_alignment_trimming
from bijux_phylogenetics.io.fasta import load_fasta_alignment

pytestmark = pytest.mark.engine_contract

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def _write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _fake_trimal(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

if "--version" in sys.argv:
    print("trimAl v2.0")
    raise SystemExit(0)

args = sys.argv[1:]
input_path = Path(args[args.index("-in") + 1])
output_path = Path(args[args.index("-out") + 1])
records = []
identifier = None
sequence = []
for raw_line in input_path.read_text(encoding="utf-8").splitlines():
    line = raw_line.strip()
    if not line:
        continue
    if line.startswith(">"):
        if identifier is not None:
            records.append((identifier, "".join(sequence)))
        identifier = line[1:]
        sequence = []
    else:
        sequence.append(line)
if identifier is not None:
    records.append((identifier, "".join(sequence)))
output_path.parent.mkdir(parents=True, exist_ok=True)
with output_path.open("w", encoding="utf-8") as handle:
    for identifier, sequence in records:
        handle.write(f">{identifier}\\n{sequence[:-1]}\\n")
print("warning: trimal fixture gap-threshold trimmed one trailing site", file=sys.stderr)
""",
    )


def test_run_alignment_trimming_writes_retained_site_summary_output(
    tmp_path: Path,
) -> None:
    executable = _fake_trimal(tmp_path / "trimal-fixture")
    input_path = fixture("alignments/example_alignment_trim.fasta")
    output_path = tmp_path / "trimmed.fasta"

    report = run_alignment_trimming(
        input_path,
        output_path,
        executable=executable,
        gap_threshold=0.2,
    )

    input_alignment_length = len(load_fasta_alignment(input_path)[0].sequence)
    assert report.output_paths["trimming_summary"].exists()
    trimming_summary_rows = (
        report.output_paths["trimming_summary"].read_text(encoding="utf-8").splitlines()
    )
    assert trimming_summary_rows[0] == "metric\tvalue"
    assert "mode\tgap-threshold" in trimming_summary_rows
    assert "gap_threshold\t0.200000" in trimming_summary_rows
    assert f"input_alignment_length\t{input_alignment_length}" in trimming_summary_rows
    assert (
        f"trimmed_alignment_length\t{input_alignment_length - 1}"
        in trimming_summary_rows
    )
    assert f"retained_site_count\t{input_alignment_length - 1}" in trimming_summary_rows
    assert "removed_site_count\t1" in trimming_summary_rows
