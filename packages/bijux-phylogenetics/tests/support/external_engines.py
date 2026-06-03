from __future__ import annotations

import os
from pathlib import Path
import shutil

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]


def _configured_executable(env_var: str) -> Path | None:
    configured = os.environ.get(env_var)
    if not configured:
        return None
    candidate = Path(configured)
    if candidate.exists():
        return candidate
    return None


def _existing_candidate(candidates: list[Path]) -> Path | None:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def real_mafft_executable() -> Path | None:
    configured = _configured_executable("BIJUX_PHYLOGENETICS_MAFFT_EXECUTABLE")
    if configured is not None:
        return configured
    resolved = shutil.which("mafft")
    if resolved is not None:
        return Path(resolved)
    return _existing_candidate(
        [REPOSITORY_ROOT / "artifacts" / "mafft" / "mafft-mac" / "mafft.bat"]
    )


def real_trimal_executable() -> Path | None:
    configured = _configured_executable("BIJUX_PHYLOGENETICS_TRIMAL_EXECUTABLE")
    if configured is not None:
        return configured
    resolved = shutil.which("trimal")
    if resolved is not None:
        return Path(resolved)
    return _existing_candidate(
        [REPOSITORY_ROOT / "artifacts" / "trimal" / "trimal" / "source" / "trimal"]
    )


def real_iqtree_executable() -> Path | None:
    configured = _configured_executable("BIJUX_PHYLOGENETICS_IQTREE_EXECUTABLE")
    if configured is not None:
        return configured
    for executable_name in ("iqtree2", "iqtree", "iqtree3"):
        resolved = shutil.which(executable_name)
        if resolved is not None:
            return Path(resolved)
    return None


def real_fasttree_executable() -> Path | None:
    configured = _configured_executable("BIJUX_PHYLOGENETICS_FASTTREE_EXECUTABLE")
    if configured is not None:
        return configured
    for executable_name in ("FastTree", "fasttree", "FastTreeMP"):
        resolved = shutil.which(executable_name)
        if resolved is not None:
            return Path(resolved)
    return _existing_candidate(
        [
            REPOSITORY_ROOT / "artifacts" / "fasttree" / "FastTree",
            REPOSITORY_ROOT / "artifacts" / "fasttree" / "fasttree",
        ]
    )


def real_mrbayes_executable() -> Path | None:
    configured = _configured_executable("BIJUX_PHYLOGENETICS_MRBAYES_EXECUTABLE")
    if configured is not None:
        return configured
    resolved = shutil.which("mb")
    if resolved is not None:
        return Path(resolved)
    return None


def real_beast_executable() -> Path | None:
    configured = _configured_executable("BIJUX_PHYLOGENETICS_BEAST_EXECUTABLE")
    if configured is not None:
        return configured
    resolved = shutil.which("beast")
    if resolved is not None:
        return Path(resolved)
    cask_candidates = sorted(
        Path("/opt/homebrew/Caskroom").glob("beast2/*/BEAST */bin/beast")
    )
    applications_candidates = sorted(Path("/Applications").glob("BEAST */bin/beast"))
    return _existing_candidate(
        cask_candidates
        + applications_candidates
        + [
            REPOSITORY_ROOT
            / "artifacts"
            / "beast2-runtime"
            / "BEAST 2.7.7"
            / "bin"
            / "beast"
        ]
    )


def require_alignment_engine_executables() -> dict[str, str]:
    executables = {
        "mafft": real_mafft_executable(),
        "trimal": real_trimal_executable(),
        "iqtree2": real_iqtree_executable(),
    }
    missing = [name for name, resolved in executables.items() if resolved is None]
    if missing:
        pytest.skip(
            "real alignment and inference workflow coverage requires installed executables: "
            + ", ".join(sorted(missing))
        )
    return {
        name: str(resolved)
        for name, resolved in executables.items()
        if resolved is not None
    }


def require_alignment_validation_matrix_executables() -> dict[str, str]:
    executables = {
        "mafft": real_mafft_executable(),
        "trimal": real_trimal_executable(),
        "iqtree": real_iqtree_executable(),
        "fasttree": real_fasttree_executable(),
    }
    missing = [name for name, resolved in executables.items() if resolved is None]
    if missing:
        pytest.skip(
            "real alignment validation matrix coverage requires installed executables: "
            + ", ".join(sorted(missing))
        )
    return {
        name: str(resolved)
        for name, resolved in executables.items()
        if resolved is not None
    }


def require_bayesian_validation_matrix_executables() -> dict[str, str]:
    executable = real_mrbayes_executable()
    if executable is None:
        pytest.skip(
            "real Bayesian validation matrix coverage requires an installed MrBayes executable"
        )
    return {"mrbayes": str(executable)}
