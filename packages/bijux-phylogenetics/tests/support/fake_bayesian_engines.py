from __future__ import annotations

from pathlib import Path
import sys


def write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)
    return path


def fake_beast(path: Path, *, version_text: str = "BEAST v2.7.7 fixture") -> Path:
    return write_executable(
        path,
        f"""#!{sys.executable}
import sys
from pathlib import Path

args = sys.argv[1:]
if "-version" in args:
    print({version_text!r})
    raise SystemExit(0)

xml_path = Path(args[-1])
seed = args[args.index("-seed") + 1]
log_path = xml_path.with_name(f"{{xml_path.stem}}.{{seed}}.log")
tree_path = xml_path.with_name(f"{{xml_path.stem}}.{{seed}}.trees")
log_path.write_text(
    "Sample\\tposterior\\tlikelihood\\tprior\\ttreeHeight\\tclockRate\\tbirthRate\\n"
    "0\\t-120.0\\t-80.0\\t-40.0\\t1.1\\t0.01\\t0.2\\n"
    "20\\t-118.0\\t-79.0\\t-39.0\\t1.0\\t0.011\\t0.21\\n",
    encoding="utf-8",
)
tree_path.write_text(
    "#NEXUS\\n"
    "Begin trees;\\n"
    "tree STATE_0 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "tree STATE_20 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "End;\\n",
    encoding="utf-8",
)
print("warning: beast fixture posterior run", file=sys.stderr)
""",
    )


def fake_beast_timeout(path: Path) -> Path:
    return write_executable(
        path,
        f"""#!{sys.executable}
import sys
import time

if "-version" in sys.argv[1:]:
    print("BEAST v2.7.7 fixture")
    raise SystemExit(0)

time.sleep(1.0)
""",
    )


def fake_beast_killed(path: Path) -> Path:
    return write_executable(
        path,
        f"""#!{sys.executable}
import os
import signal
import sys

if "-version" in sys.argv[1:]:
    print("BEAST v2.7.7 fixture")
    raise SystemExit(0)

os.kill(os.getpid(), signal.SIGTERM)
""",
    )


def fake_beast_malformed_outputs(path: Path) -> Path:
    return write_executable(
        path,
        f"""#!{sys.executable}
import sys
from pathlib import Path

args = sys.argv[1:]
if "-version" in args:
    print("BEAST v2.7.7 fixture")
    raise SystemExit(0)

xml_path = Path(args[-1])
seed = args[args.index("-seed") + 1]
log_path = xml_path.with_name(f"{{xml_path.stem}}.{{seed}}.log")
tree_path = xml_path.with_name(f"{{xml_path.stem}}.{{seed}}.trees")
log_path.write_text(
    "Sample\\tposterior\\tlikelihood\\tprior\\ttreeHeight\\tclockRate\\tbirthRate\\n"
    "0\\t-120.0\\t-80.0\\t-40.0\\t1.1\\t0.01\\t0.2\\n"
    "20\\tbad\\t-79.0\\t-39.0\\t1.0\\t0.011\\t0.21\\n",
    encoding="utf-8",
)
tree_path.write_text(
    "#NEXUS\\n"
    "Begin trees;\\n"
    "tree STATE_0 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "End;\\n",
    encoding="utf-8",
)
""",
    )


def fake_mrbayes(
    path: Path,
    *,
    version_text: str = "MrBayes v3.2.7a fixture",
) -> Path:
    return write_executable(
        path,
        f"""#!{sys.executable}
import sys
from pathlib import Path

if "--version" in sys.argv[1:] or "-v" in sys.argv[1:]:
    print({version_text!r})
    raise SystemExit(0)

nexus_path = Path(sys.argv[1])
trace_path = Path(f"{{nexus_path}}.run1.p")
tree_path = Path(f"{{nexus_path}}.run1.t")
mcmc_path = Path(f"{{nexus_path}}.mcmc")
consensus_path = Path(f"{{nexus_path}}.con.tre")
trace_path.write_text(
    "Gen\\tLnL\\tTL\\talpha\\n"
    "0\\t-110.0\\t0.40\\t0.90\\n"
    "100\\t-108.0\\t0.41\\t0.95\\n"
    "200\\t-107.0\\t0.42\\t1.00\\n"
    "300\\t-106.5\\t0.43\\t1.05\\n",
    encoding="utf-8",
)
tree_path.write_text(
    "#NEXUS\\n"
    "begin trees;\\n"
    "tree gen1 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "tree gen2 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "tree gen3 = [&R] ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);\\n"
    "tree gen4 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "end;\\n",
    encoding="utf-8",
)
mcmc_path.write_text(
    "[ID: 1]\\n"
    "[   Gen -- Generation]\\n"
    "Gen\\tMove$acc_run1\\tSwap(1<>2)$acc(1)\\tAvgStdDev(s)\\n"
    "100\\t0.5\\t0.75\\t0.20\\n"
    "200\\tNA\\t1.0\\t0.10\\n",
    encoding="utf-8",
)
consensus_path.write_text(
    "#NEXUS\\n"
    "begin trees;\\n"
    "tree con_50_majrule = [&R] ((A[&prob=1.0,prob(percent)=\\\"100\\\"]:0.1,B[&prob=1.0,prob(percent)=\\\"100\\\"]:0.1)"
    "[&prob=0.75,prob(percent)=\\\"75\\\"]:0.2,(C[&prob=1.0,prob(percent)=\\\"100\\\"]:0.1,D[&prob=1.0,prob(percent)=\\\"100\\\"]:0.1)"
    "[&prob=0.5,prob(percent)=\\\"50\\\"]:0.2);\\n"
    "end;\\n",
    encoding="utf-8",
)
print("warning: mrbayes fixture posterior run", file=sys.stderr)
""",
    )


def fake_mrbayes_timeout(path: Path) -> Path:
    return write_executable(
        path,
        f"""#!{sys.executable}
import sys
import time

if "--version" in sys.argv[1:] or "-v" in sys.argv[1:]:
    print("MrBayes v3.2.7a fixture")
    raise SystemExit(0)

time.sleep(1.0)
""",
    )


def fake_mrbayes_killed(path: Path) -> Path:
    return write_executable(
        path,
        f"""#!{sys.executable}
import os
import signal
import sys

if "--version" in sys.argv[1:] or "-v" in sys.argv[1:]:
    print("MrBayes v3.2.7a fixture")
    raise SystemExit(0)

os.kill(os.getpid(), signal.SIGTERM)
""",
    )


def fake_mrbayes_malformed_outputs(path: Path) -> Path:
    return write_executable(
        path,
        f"""#!{sys.executable}
import sys
from pathlib import Path

if "--version" in sys.argv[1:] or "-v" in sys.argv[1:]:
    print("MrBayes v3.2.7a fixture")
    raise SystemExit(0)

nexus_path = Path(sys.argv[1])
trace_path = Path(f"{{nexus_path}}.run1.p")
tree_path = Path(f"{{nexus_path}}.run1.t")
mcmc_path = Path(f"{{nexus_path}}.mcmc")
consensus_path = Path(f"{{nexus_path}}.con.tre")
trace_path.write_text(
    "Gen\\tLnL\\tTL\\talpha\\n"
    "0\\t-110.0\\t0.40\\t0.90\\n"
    "100\\tbad\\t0.41\\t0.95\\n",
    encoding="utf-8",
)
tree_path.write_text(
    "#NEXUS\\n"
    "begin trees;\\n"
    "tree gen1 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "end;\\n",
    encoding="utf-8",
)
mcmc_path.write_text(
    "Gen\\tMove$acc_run1\\tSwap(1<>2)$acc(1)\\tAvgStdDev(s)\\n"
    "100\\t0.5\\t0.75\\t0.20\\n",
    encoding="utf-8",
)
consensus_path.write_text(
    "#NEXUS\\n"
    "begin trees;\\n"
    "tree con_50_majrule = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "end;\\n",
    encoding="utf-8",
)
""",
    )
