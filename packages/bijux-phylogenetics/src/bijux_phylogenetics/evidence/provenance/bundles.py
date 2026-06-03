from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import shutil

from bijux_phylogenetics.core.environment import inspect_environment
from bijux_phylogenetics.runtime.errors import EvidenceContractError


@dataclass(slots=True)
class ArtifactBundleFile:
    section: str
    relative_path: Path
    sha256: str
    size_bytes: int


@dataclass(slots=True)
class ArtifactBundleReport:
    input_roots: list[Path]
    output_roots: list[Path]
    output_root: Path
    file_count: int
    input_file_count: int
    output_file_count: int
    files: list[ArtifactBundleFile]


@dataclass(slots=True)
class ArtifactBundleMismatch:
    relative_path: Path
    reason: str


@dataclass(slots=True)
class ArtifactBundleValidationReport:
    bundle_root: Path
    file_count: int
    valid: bool
    mismatches: list[ArtifactBundleMismatch]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _copy_roots(
    *, roots: list[Path], section: str, bundle_root: Path
) -> list[ArtifactBundleFile]:
    bundled_files: list[ArtifactBundleFile] = []
    for source_root in roots:
        if not source_root.exists():
            raise EvidenceContractError(f"{section} directory not found: {source_root}")
        if not source_root.is_dir():
            raise EvidenceContractError(
                f"{section} root is not a directory: {source_root}"
            )
        for source in sorted(path for path in source_root.rglob("*") if path.is_file()):
            relative_path = Path(source_root.name) / source.relative_to(source_root)
            destination = bundle_root / section / relative_path
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            bundled_files.append(
                ArtifactBundleFile(
                    section=section,
                    relative_path=relative_path,
                    sha256=_sha256(source),
                    size_bytes=source.stat().st_size,
                )
            )
    return bundled_files


def _copy_paths(
    *, paths: list[Path], section: str, bundle_root: Path
) -> list[ArtifactBundleFile]:
    bundled_files: list[ArtifactBundleFile] = []
    for source in paths:
        if not source.exists():
            raise EvidenceContractError(f"{section} file not found: {source}")
        if not source.is_file():
            raise EvidenceContractError(f"{section} path is not a file: {source}")
        relative_path = Path(*source.resolve().parts[1:])
        destination = bundle_root / section / relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        bundled_files.append(
            ArtifactBundleFile(
                section=section,
                relative_path=relative_path,
                sha256=_sha256(source),
                size_bytes=source.stat().st_size,
            )
        )
    return bundled_files


def _write_checksums_tsv(path: Path, files: list[ArtifactBundleFile]) -> Path:
    lines = ["section\trelative_path\tsha256\tsize_bytes"]
    for file in files:
        lines.append(
            f"{file.section}\t{file.relative_path.as_posix()}\t{file.sha256}\t{file.size_bytes}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_environment_json(path: Path) -> Path:
    payload = asdict(inspect_environment())
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _write_bundle_readme(
    path: Path, *, input_roots: list[Path], output_roots: list[Path], file_count: int
) -> Path:
    lines = [
        "# Bijux Artifact Bundle",
        "",
        "This bundle captures explicit phylogenetics inputs and outputs for deterministic review.",
        "",
        "## Source Roots",
        "",
        "Inputs:",
    ]
    lines.extend(f"- `{root}`" for root in input_roots)
    lines.append("")
    lines.append("Outputs:")
    lines.extend(f"- `{root}`" for root in output_roots)
    lines.extend(
        [
            "",
            "## Contents",
            "",
            "- `manifest.json` records the copied file inventory.",
            "- `checksums.tsv` records one checksum row per copied file.",
            "- `environment.json` records runtime dependency availability.",
            "- `inputs/` contains copied input files grouped by source directory name.",
            "- `outputs/` contains copied output files grouped by source directory name.",
            "",
            f"Copied files: `{file_count}`",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def bundle_artifact_directories(
    input_roots: list[Path], output_roots: list[Path], bundle_root: Path
) -> ArtifactBundleReport:
    """Copy explicit input and output directories into a checksummed artifact bundle."""
    if not input_roots:
        raise EvidenceContractError(
            "artifact bundle requires at least one input directory"
        )
    if not output_roots:
        raise EvidenceContractError(
            "artifact bundle requires at least one output directory"
        )

    normalized_inputs = [Path(root) for root in input_roots]
    normalized_outputs = [Path(root) for root in output_roots]
    if bundle_root.exists():
        shutil.rmtree(bundle_root)
    (bundle_root / "inputs").mkdir(parents=True, exist_ok=True)
    (bundle_root / "outputs").mkdir(parents=True, exist_ok=True)

    input_files = _copy_roots(
        roots=normalized_inputs, section="inputs", bundle_root=bundle_root
    )
    output_files = _copy_roots(
        roots=normalized_outputs, section="outputs", bundle_root=bundle_root
    )
    bundled_files = input_files + output_files

    manifest_payload = {
        "input_roots": [str(root) for root in normalized_inputs],
        "output_roots": [str(root) for root in normalized_outputs],
        "bundle_root": str(bundle_root),
        "file_count": len(bundled_files),
        "files": [
            {
                "section": file.section,
                "relative_path": file.relative_path.as_posix(),
                "sha256": file.sha256,
                "size_bytes": file.size_bytes,
            }
            for file in bundled_files
        ],
    }
    (bundle_root / "manifest.json").write_text(
        json.dumps(manifest_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _write_checksums_tsv(bundle_root / "checksums.tsv", bundled_files)
    _write_environment_json(bundle_root / "environment.json")
    _write_bundle_readme(
        bundle_root / "README.md",
        input_roots=normalized_inputs,
        output_roots=normalized_outputs,
        file_count=len(bundled_files),
    )

    return ArtifactBundleReport(
        input_roots=normalized_inputs,
        output_roots=normalized_outputs,
        output_root=bundle_root,
        file_count=len(bundled_files),
        input_file_count=len(input_files),
        output_file_count=len(output_files),
        files=bundled_files,
    )


def bundle_artifact_files(
    input_paths: list[Path], output_paths: list[Path], bundle_root: Path
) -> ArtifactBundleReport:
    """Copy explicit input and output files into a checksummed artifact bundle."""
    if not input_paths:
        raise EvidenceContractError("artifact bundle requires at least one input file")
    if not output_paths:
        raise EvidenceContractError("artifact bundle requires at least one output file")

    normalized_inputs = [Path(path) for path in input_paths]
    normalized_outputs = [Path(path) for path in output_paths]
    if bundle_root.exists():
        shutil.rmtree(bundle_root)
    (bundle_root / "inputs").mkdir(parents=True, exist_ok=True)
    (bundle_root / "outputs").mkdir(parents=True, exist_ok=True)

    input_files = _copy_paths(
        paths=normalized_inputs, section="inputs", bundle_root=bundle_root
    )
    output_files = _copy_paths(
        paths=normalized_outputs, section="outputs", bundle_root=bundle_root
    )
    bundled_files = input_files + output_files

    manifest_payload = {
        "input_paths": [str(path) for path in normalized_inputs],
        "output_paths": [str(path) for path in normalized_outputs],
        "bundle_root": str(bundle_root),
        "file_count": len(bundled_files),
        "files": [
            {
                "section": file.section,
                "relative_path": file.relative_path.as_posix(),
                "sha256": file.sha256,
                "size_bytes": file.size_bytes,
            }
            for file in bundled_files
        ],
    }
    (bundle_root / "manifest.json").write_text(
        json.dumps(manifest_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _write_checksums_tsv(bundle_root / "checksums.tsv", bundled_files)
    _write_environment_json(bundle_root / "environment.json")
    _write_bundle_readme(
        bundle_root / "README.md",
        input_roots=normalized_inputs,
        output_roots=normalized_outputs,
        file_count=len(bundled_files),
    )

    return ArtifactBundleReport(
        input_roots=normalized_inputs,
        output_roots=normalized_outputs,
        output_root=bundle_root,
        file_count=len(bundled_files),
        input_file_count=len(input_files),
        output_file_count=len(output_files),
        files=bundled_files,
    )


def validate_artifact_bundle(bundle_root: Path) -> ArtifactBundleValidationReport:
    """Validate a checksummed artifact bundle created by this module."""
    manifest_path = bundle_root / "manifest.json"
    if not manifest_path.is_file():
        raise EvidenceContractError(f"evidence manifest not found: {manifest_path}")
    checksums_path = bundle_root / "checksums.tsv"
    if not checksums_path.is_file():
        raise EvidenceContractError(
            f"evidence checksums file not found: {checksums_path}"
        )
    environment_path = bundle_root / "environment.json"
    if not environment_path.is_file():
        raise EvidenceContractError(
            f"evidence environment file not found: {environment_path}"
        )
    readme_path = bundle_root / "README.md"
    if not readme_path.is_file():
        raise EvidenceContractError(f"evidence readme not found: {readme_path}")
    manifest_payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    files = manifest_payload.get("files", [])
    if not isinstance(files, list):
        raise EvidenceContractError(
            f"evidence manifest file list is invalid: {manifest_path}"
        )
    mismatches: list[ArtifactBundleMismatch] = []
    for entry in files:
        if not isinstance(entry, dict):
            raise EvidenceContractError(
                f"evidence manifest file entry is invalid: {manifest_path}"
            )
        relative_path = Path(entry["relative_path"])
        target = bundle_root / entry["section"] / relative_path
        if not target.is_file():
            mismatches.append(
                ArtifactBundleMismatch(
                    relative_path=relative_path,
                    reason="missing_file",
                )
            )
            continue
        actual_sha = _sha256(target)
        actual_size = target.stat().st_size
        if actual_sha != entry["sha256"]:
            mismatches.append(
                ArtifactBundleMismatch(
                    relative_path=relative_path,
                    reason="checksum_mismatch",
                )
            )
            continue
        if actual_size != entry["size_bytes"]:
            mismatches.append(
                ArtifactBundleMismatch(
                    relative_path=relative_path,
                    reason="size_mismatch",
                )
            )
    return ArtifactBundleValidationReport(
        bundle_root=bundle_root,
        file_count=len(files),
        valid=not mismatches,
        mismatches=mismatches,
    )
