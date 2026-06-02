from __future__ import annotations

from pathlib import Path
import sys


def write_executable(path: Path, body: str) -> Path:
    normalized_body = body
    if normalized_body.startswith("#!/usr/bin/env python3\n"):
        normalized_body = f"#!{sys.executable}\n" + normalized_body.removeprefix(
            "#!/usr/bin/env python3\n"
        )
    path.write_text(normalized_body, encoding="utf-8")
    path.chmod(0o755)
    return path


def fake_mafft(path: Path, *, version: str = "7.999") -> Path:
    return write_executable(
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


def fake_trimal(path: Path, *, version: str = "2.0") -> Path:
    return write_executable(
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


def fake_iqtree(path: Path, *, version: str = "2.9.9") -> Path:
    return write_executable(
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


def fake_fasttree(path: Path, *, version: str = "2.2") -> Path:
    return write_executable(
        path,
        f"""#!/usr/bin/env python3
import sys

args = sys.argv[1:]
if not args or "-help" in args:
    print("FastTree Version {version} fixture")
    raise SystemExit(0)

print("((A:0.1,B:0.1)0.98:0.3,(C:0.1,D:0.1)0.62:0.3);")
print("warning: fasttree fixture approximate support only", file=sys.stderr)
""",
    )
