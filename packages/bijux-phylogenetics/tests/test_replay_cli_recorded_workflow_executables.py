from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.engines import run_fasta_to_tree_workflow

pytestmark = pytest.mark.engine_contract

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def _write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _fake_mafft(path: Path, *, version: str = "7.999") -> Path:
    return _write_executable(
        path,
        f"""#!/usr/bin/env python3
import sys
from pathlib import Path

if "--version" in sys.argv:
    print("mafft v{version}", file=sys.stderr)
    raise SystemExit(0)

input_path = Path(sys.argv[-1])
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
width = max(len(row[1]) for row in records)
for identifier, sequence in records:
    print(f">{{identifier}}")
    print(sequence.ljust(width, "-"))
""",
    )


def _fake_trimal(path: Path, *, version: str = "2.0") -> Path:
    return _write_executable(
        path,
        f"""#!/usr/bin/env python3
import sys
from pathlib import Path

if "--version" in sys.argv:
    print("trimAl v{version}")
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
with output_path.open("w", encoding="utf-8") as handle:
    for identifier, sequence in records:
        handle.write(f">{{identifier}}\\n{{sequence[:-1]}}\\n")
""",
    )


def _fake_iqtree(path: Path, *, version: str = "2.9.9") -> Path:
    return _write_executable(
        path,
        f"""#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    print("IQ-TREE multicore version {version}")
    raise SystemExit(0)

prefix = Path(args[args.index("-pre") + 1]) if "-pre" in args else Path("iqtree")
prefix.parent.mkdir(parents=True, exist_ok=True)
if "-m" in args and args[args.index("-m") + 1] == "MF":
    prefix.with_suffix(".iqtree").write_text(
        " No. Model         -LnL         df  AIC          AICc         BIC\\n"
        "  1  GTR+G         123.456      12  270.912      330.912      272.912\\n"
        "  2  HKY+G         124.000      10  268.000      320.000      269.000\\n"
        "  3  JC            130.500      5   271.000      300.000      271.500\\n"
        "Akaike Information Criterion:           HKY+G\\n"
        "Corrected Akaike Information Criterion: JC\\n"
        "Bayesian Information Criterion:         GTR+G\\n"
        "Best-fit model according to BIC: GTR+G\\n"
        "Log-likelihood of the tree: -123.456\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".log").write_text(
        "IQ-TREE fixture model-selection log\\nBEST SCORE FOUND : -123.456\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".model").write_text(
        "Best-fit model: GTR+G\\n",
        encoding="utf-8",
    )
    raise SystemExit(0)
if "-bb" in args:
    support_tree = "((A:0.1,B:0.1)95:0.2,(C:0.1,D:0.1)88:0.2);\\n"
    prefix.with_suffix(".treefile").write_text(
        support_tree,
        encoding="utf-8",
    )
    prefix.with_suffix(".contree").write_text(support_tree, encoding="utf-8")
    prefix.with_suffix(".ufboot").write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".iqtree").write_text(
        "Best-fit model: GTR+G\\nLog-likelihood of the tree: -222.222\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".log").write_text(
        "IQ-TREE fixture bootstrap log\\nBEST SCORE FOUND : -222.222\\n",
        encoding="utf-8",
    )
    raise SystemExit(0)
if "-m" in args:
    prefix.with_suffix(".treefile").write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".iqtree").write_text(
        "Best-fit model: GTR+G\\nLog-likelihood of the tree: -200.000\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".log").write_text(
        "IQ-TREE fixture ml log\\nBEST SCORE FOUND : -200.000\\n",
        encoding="utf-8",
    )
    raise SystemExit(0)
raise SystemExit(2)
""",
    )


@pytest.mark.slow
def test_phylo_replay_cli_reuses_recorded_workflow_executables(
    tmp_path: Path, capsys
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    workflow = run_fasta_to_tree_workflow(
        fixture("alignments/example_sequences_raw.fasta"),
        out_dir=tmp_path / "fasta-to-tree",
        prefix="example",
        mafft_executable=mafft,
        trimal_executable=trimal,
        iqtree_executable=iqtree,
        bootstrap_replicates=1000,
    )

    exit_code = main(
        [
            "phylo",
            "replay",
            str(workflow.manifest_path),
            "--out-dir",
            str(tmp_path / "fasta-to-tree-replay"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["workflow"] == "fasta-to-tree"
    assert payload["metrics"]["engine_version_drift_count"] == 0
    assert payload["metrics"]["outputs_equivalent"] is True
