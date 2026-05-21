from __future__ import annotations

import gzip
from math import ceil
from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmTreeRetentionFileRow
from .shared import (
    _COMPRESSED_SUFFIXES,
    _COMPRESSION_RECOMMENDED_BYTES,
    _COMPRESSION_RECOMMENDED_TREE_COUNT,
    _COMPRESSION_REQUIRED_BYTES,
    _COMPRESSION_REQUIRED_TREE_COUNT,
    _THINNING_RECOMMENDED_TREE_COUNT,
    _THINNING_REQUIRED_TREE_COUNT,
    _THINNING_TARGET_TREE_COUNT,
    _TREE_FILE_SUFFIXES,
)


def iter_tree_relative_paths(bundle_root: Path) -> tuple[Path, ...]:
    return tuple(
        relative_path
        for relative_path in sorted(
            path.relative_to(bundle_root)
            for path in bundle_root.rglob("*")
            if path.is_file() and _is_tree_bearing_file(path.relative_to(bundle_root))
        )
    )


def build_tree_retention_file_row(
    *,
    bundle_root: Path,
    relative_path: Path,
) -> RabiesMethodSensitivitySlurmTreeRetentionFileRow:
    path = bundle_root / relative_path
    byte_count = path.stat().st_size
    tree_count = _count_trees(path)
    artifact_scope = _classify_artifact_scope(relative_path)
    thinning_policy, thinning_interval, retained_tree_count = _derive_thinning_policy(
        tree_count=tree_count
    )
    compression_policy, recommended_suffix = _derive_compression_policy(
        relative_path=relative_path,
        tree_count=tree_count,
        byte_count=byte_count,
    )
    issues: list[str] = []
    if thinning_policy == "thin_recommended":
        issues.append(
            "tree-set size is large enough that interval thinning is recommended"
        )
    elif thinning_policy == "thin_required":
        issues.append(
            "tree-set size is large enough that interval thinning is required"
        )
    if compression_policy == "compress_recommended":
        issues.append(
            "tree-set size is large enough that gzip compression is recommended"
        )
    elif compression_policy == "compress_required":
        issues.append("tree-set size is large enough that gzip compression is required")
    if tree_count <= 1 and artifact_scope == "tree_artifact":
        issues.append("single-tree artifact should be kept in full without thinning")
    return RabiesMethodSensitivitySlurmTreeRetentionFileRow(
        variant_id=_variant_id_from_relative_path(relative_path),
        relative_path=relative_path.as_posix(),
        artifact_scope=artifact_scope,
        tree_count=tree_count,
        byte_count=byte_count,
        thinning_policy=thinning_policy,
        thinning_interval=thinning_interval,
        retained_tree_count=retained_tree_count,
        compression_policy=compression_policy,
        recommended_suffix=recommended_suffix,
        issue_count=len(issues),
        issues=tuple(issues),
    )


def _derive_thinning_policy(*, tree_count: int) -> tuple[str, int, int]:
    if tree_count <= 1:
        return ("not_applicable", 1, tree_count)
    if tree_count >= _THINNING_REQUIRED_TREE_COUNT:
        interval = max(2, ceil(tree_count / _THINNING_TARGET_TREE_COUNT))
        return ("thin_required", interval, ceil(tree_count / interval))
    if tree_count >= _THINNING_RECOMMENDED_TREE_COUNT:
        interval = max(2, ceil(tree_count / _THINNING_TARGET_TREE_COUNT))
        return ("thin_recommended", interval, ceil(tree_count / interval))
    return ("keep_full", 1, tree_count)


def _derive_compression_policy(
    *,
    relative_path: Path,
    tree_count: int,
    byte_count: int,
) -> tuple[str, str]:
    if any(suffix in _COMPRESSED_SUFFIXES for suffix in relative_path.suffixes):
        return ("already_compressed", "")
    if (
        tree_count >= _COMPRESSION_REQUIRED_TREE_COUNT
        or byte_count >= _COMPRESSION_REQUIRED_BYTES
    ):
        return ("compress_required", ".gz")
    if (
        tree_count >= _COMPRESSION_RECOMMENDED_TREE_COUNT
        or byte_count >= _COMPRESSION_RECOMMENDED_BYTES
    ):
        return ("compress_recommended", ".gz")
    return ("keep_plain", "")


def _classify_artifact_scope(relative_path: Path) -> str:
    name = relative_path.name.lower()
    path_text = relative_path.as_posix().lower()
    if any(
        token in path_text
        for token in ("posterior", "boottrees", "bootstrap", "treeset", "samples")
    ) or name.endswith(".trees"):
        return "posterior_sample"
    return "tree_artifact"


def _variant_id_from_relative_path(relative_path: Path) -> str:
    parts = relative_path.parts
    if len(parts) >= 3 and parts[0] == "variants":
        return parts[1]
    if len(parts) >= 3 and parts[0] == "slurm-job-evidence":
        return parts[1]
    return "workflow_shared"


def _count_trees(path: Path) -> int:
    if ".gz" in path.suffixes:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return max(1, handle.read().count(";"))
    return max(1, path.read_text(encoding="utf-8").count(";"))


def _is_tree_bearing_file(relative_path: Path) -> bool:
    if relative_path.parts and relative_path.parts[0] != "variants":
        return False
    if len(relative_path.parts) < 3:
        return False
    suffixes = relative_path.suffixes
    if not suffixes:
        return False
    terminal_suffix = suffixes[-1].lower()
    if terminal_suffix in _COMPRESSED_SUFFIXES and len(suffixes) >= 2:
        terminal_suffix = suffixes[-2].lower()
    return terminal_suffix in _TREE_FILE_SUFFIXES
