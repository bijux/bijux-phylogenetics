from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.engines import (
    run_bootstrap_consensus_tree,
    run_bootstrap_support_estimation,
    run_sh_alrt_support_estimation,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError

pytestmark = pytest.mark.engine_contract

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def _write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _fake_iqtree(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    print("IQ-TREE multicore version 2.9.9")
    raise SystemExit(0)

prefix = Path(args[args.index("-pre") + 1]) if "-pre" in args else Path("iqtree")
prefix.parent.mkdir(parents=True, exist_ok=True)
if "-con" in args:
    prefix.with_suffix(".contree").write_text(
        "((A:0.1,B:0.1)90:0.2,(C:0.1,D:0.1)85:0.2);\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".log").write_text(
        "IQ-TREE fixture consensus log\\n",
        encoding="utf-8",
    )
    raise SystemExit(0)

raise SystemExit(2)
""",
    )


def _fake_iqtree_partial_bootstrap_support(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    print("IQ-TREE multicore version 2.9.9")
    raise SystemExit(0)

prefix = Path(args[args.index("-pre") + 1]) if "-pre" in args else Path("iqtree")
prefix.parent.mkdir(parents=True, exist_ok=True)
prefix.with_suffix(".contree").write_text(
    "((A:0.1,B:0.1)95:0.2,(C:0.1,D:0.1):0.2);\\n",
    encoding="utf-8",
)
prefix.with_suffix(".ufboot").write_text(
    "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
    encoding="utf-8",
)
prefix.with_suffix(".iqtree").write_text(
    "Best-fit model: GTR+G\\nLog-likelihood of the tree: -222.222\\nBootstrap analysis completed\\n",
    encoding="utf-8",
)
prefix.with_suffix(".log").write_text(
    "IQ-TREE fixture bootstrap log\\nBEST SCORE FOUND : -222.222\\n",
    encoding="utf-8",
)
""",
    )


def _fake_iqtree_bootstrap_taxa_mismatch(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    print("IQ-TREE multicore version 2.9.9")
    raise SystemExit(0)

prefix = Path(args[args.index("-pre") + 1]) if "-pre" in args else Path("iqtree")
prefix.parent.mkdir(parents=True, exist_ok=True)
prefix.with_suffix(".contree").write_text(
    "((A:0.1,B:0.1)95:0.2,(C:0.1,D:0.1)88:0.2);\\n",
    encoding="utf-8",
)
prefix.with_suffix(".ufboot").write_text(
    "((A:0.1,B:0.1):0.2,(C:0.1,E:0.1):0.2);\\n"
    "((A:0.1,B:0.1):0.2,(C:0.1,E:0.1):0.2);\\n",
    encoding="utf-8",
)
prefix.with_suffix(".iqtree").write_text(
    "Best-fit model: GTR+G\\nLog-likelihood of the tree: -222.222\\nBootstrap analysis completed\\n",
    encoding="utf-8",
)
prefix.with_suffix(".log").write_text(
    "IQ-TREE fixture bootstrap log\\nBEST SCORE FOUND : -222.222\\n",
    encoding="utf-8",
)
""",
    )


def _fake_iqtree_partial_sh_alrt_support(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    print("IQ-TREE multicore version 2.9.9")
    raise SystemExit(0)

prefix = Path(args[args.index("-pre") + 1]) if "-pre" in args else Path("iqtree")
prefix.parent.mkdir(parents=True, exist_ok=True)
prefix.with_suffix(".treefile").write_text(
    "((A:0.1,B:0.1)82/97:0.2,(C:0.1,D:0.1):0.2);\\n",
    encoding="utf-8",
)
prefix.with_suffix(".ufboot").write_text(
    "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
    encoding="utf-8",
)
prefix.with_suffix(".iqtree").write_text(
    "Best-fit model: GTR+G\\nLog-likelihood of the tree: -222.222\\nSH-aLRT and ultrafast bootstrap analysis completed\\n",
    encoding="utf-8",
)
prefix.with_suffix(".log").write_text(
    "IQ-TREE fixture sh-alrt log\\nBEST SCORE FOUND : -222.222\\n",
    encoding="utf-8",
)
""",
    )


def _fake_iqtree_partial_consensus_support(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    print("IQ-TREE multicore version 2.9.9")
    raise SystemExit(0)

prefix = Path(args[args.index("-pre") + 1]) if "-pre" in args else Path("iqtree")
prefix.parent.mkdir(parents=True, exist_ok=True)
prefix.with_suffix(".contree").write_text(
    "((A:0.1,B:0.1)90:0.2,(C:0.1,D:0.1):0.2);\\n",
    encoding="utf-8",
)
prefix.with_suffix(".log").write_text(
    "IQ-TREE fixture consensus log\\n",
    encoding="utf-8",
)
""",
    )


def test_run_bootstrap_support_estimation_requires_complete_support_labels(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree_partial_bootstrap_support(
        tmp_path / "iqtree-partial-bootstrap-labels"
    )

    with pytest.raises(EngineWorkflowError) as error:
        run_bootstrap_support_estimation(
            fixture("alignments/example_alignment.fasta"),
            out_dir=tmp_path / "bootstrap",
            model="GTR+G",
            executable=executable,
            prefix="example",
            replicates=1000,
        )

    assert error.value.code == "engine_support_values_incomplete"
    assert error.value.details["workflow"] == "bootstrap-support"
    assert error.value.details["support_kind"] == "bootstrap support"


def test_run_bootstrap_support_estimation_rejects_taxa_mismatch(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree_bootstrap_taxa_mismatch(
        tmp_path / "iqtree-bootstrap-taxa-mismatch"
    )

    with pytest.raises(EngineWorkflowError) as error:
        run_bootstrap_support_estimation(
            fixture("alignments/example_alignment.fasta"),
            out_dir=tmp_path / "bootstrap",
            model="GTR+G",
            executable=executable,
            prefix="example",
            replicates=1000,
        )

    assert error.value.code == "engine_output_taxa_mismatch"
    assert error.value.details["workflow"] == "bootstrap-support"


def test_run_sh_alrt_support_estimation_requires_complete_joint_support_labels(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree_partial_sh_alrt_support(
        tmp_path / "iqtree-partial-sh-alrt-labels"
    )

    with pytest.raises(EngineWorkflowError) as error:
        run_sh_alrt_support_estimation(
            fixture("alignments/example_alignment.fasta"),
            out_dir=tmp_path / "sh-alrt",
            model="GTR+G",
            executable=executable,
            prefix="example",
            sh_alrt_replicates=1000,
            bootstrap_replicates=1000,
        )

    assert error.value.code == "engine_support_values_incomplete"
    assert error.value.details["workflow"] == "sh-alrt-support"
    assert error.value.details["support_kind"] in {
        "ultrafast bootstrap support",
        "sh-alrt support",
        "joint sh-alrt and ultrafast bootstrap support",
    }


def test_run_bootstrap_consensus_tree_requires_complete_support_labels(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree_partial_consensus_support(
        tmp_path / "iqtree-partial-consensus-support"
    )
    bootstrap_trees_path = tmp_path / "bootstrap.ufboot"
    bootstrap_trees_path.write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n"
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n",
        encoding="utf-8",
    )

    with pytest.raises(EngineWorkflowError) as error:
        run_bootstrap_consensus_tree(
            bootstrap_trees_path,
            out_dir=tmp_path / "consensus",
            executable=executable,
            prefix="example",
        )

    assert error.value.code == "engine_support_values_incomplete"
    assert error.value.details["workflow"] == "bootstrap-consensus"
    assert error.value.details["support_kind"] == "bootstrap consensus support"


def test_run_bootstrap_consensus_tree_requires_matching_input_taxa(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-consensus")
    bootstrap_trees_path = tmp_path / "bootstrap.ufboot"
    bootstrap_trees_path.write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,E:0.1):0.2);\n"
        "((A:0.1,B:0.1):0.2,(C:0.1,E:0.1):0.2);\n",
        encoding="utf-8",
    )

    with pytest.raises(EngineWorkflowError) as error:
        run_bootstrap_consensus_tree(
            bootstrap_trees_path,
            out_dir=tmp_path / "consensus",
            executable=executable,
            prefix="example",
        )

    assert error.value.code == "engine_output_taxa_mismatch"
    assert error.value.details["workflow"] == "bootstrap-consensus"
