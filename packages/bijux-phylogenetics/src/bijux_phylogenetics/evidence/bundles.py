from __future__ import annotations

import hashlib
import json
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path

from bijux_phylogenetics.core.environment import inspect_environment
from bijux_phylogenetics.errors import EvidenceContractError


@dataclass(slots=True)
class EvidenceFile:
    section: str
    relative_path: Path
    sha256: str
    size_bytes: int


@dataclass(slots=True)
class EvidenceBundleReport:
    input_roots: list[Path]
    output_roots: list[Path]
    output_root: Path
    file_count: int
    input_file_count: int
    output_file_count: int
    files: list[EvidenceFile]


@dataclass(slots=True)
class EvidenceMismatch:
    relative_path: Path
    reason: str


@dataclass(slots=True)
class EvidenceValidationReport:
    bundle_root: Path
    file_count: int
    valid: bool
    mismatches: list[EvidenceMismatch]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _copy_roots(*, roots: list[Path], section: str, bundle_root: Path) -> list[EvidenceFile]:
    bundled_files: list[EvidenceFile] = []
    for source_root in roots:
        if not source_root.exists():
            raise EvidenceContractError(f"{section} directory not found: {source_root}")
        if not source_root.is_dir():
            raise EvidenceContractError(f"{section} root is not a directory: {source_root}")
        for source in sorted(path for path in source_root.rglob("*") if path.is_file()):
            relative_path = Path(source_root.name) / source.relative_to(source_root)
            destination = bundle_root / section / relative_path
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            bundled_files.append(
                EvidenceFile(
                    section=section,
                    relative_path=relative_path,
                    sha256=_sha256(source),
                    size_bytes=source.stat().st_size,
                )
            )
    return bundled_files


def _copy_paths(*, paths: list[Path], section: str, bundle_root: Path) -> list[EvidenceFile]:
    bundled_files: list[EvidenceFile] = []
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
            EvidenceFile(
                section=section,
                relative_path=relative_path,
                sha256=_sha256(source),
                size_bytes=source.stat().st_size,
            )
        )
    return bundled_files


def _write_checksums_tsv(path: Path, files: list[EvidenceFile]) -> Path:
    lines = ["section\trelative_path\tsha256\tsize_bytes"]
    for file in files:
        lines.append(f"{file.section}\t{file.relative_path.as_posix()}\t{file.sha256}\t{file.size_bytes}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_environment_json(path: Path) -> Path:
    payload = asdict(inspect_environment())
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _write_bundle_readme(path: Path, *, input_roots: list[Path], output_roots: list[Path], file_count: int) -> Path:
    lines = [
        "# Bijux Evidence Bundle",
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


def bundle_directory(input_roots: list[Path], output_roots: list[Path], bundle_root: Path) -> EvidenceBundleReport:
    """Copy explicit input and output directories into a checksummed evidence bundle."""
    if not input_roots:
        raise EvidenceContractError("evidence bundle requires at least one input directory")
    if not output_roots:
        raise EvidenceContractError("evidence bundle requires at least one output directory")

    normalized_inputs = [Path(root) for root in input_roots]
    normalized_outputs = [Path(root) for root in output_roots]
    if bundle_root.exists():
        shutil.rmtree(bundle_root)
    (bundle_root / "inputs").mkdir(parents=True, exist_ok=True)
    (bundle_root / "outputs").mkdir(parents=True, exist_ok=True)

    input_files = _copy_roots(roots=normalized_inputs, section="inputs", bundle_root=bundle_root)
    output_files = _copy_roots(roots=normalized_outputs, section="outputs", bundle_root=bundle_root)
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
    (bundle_root / "manifest.json").write_text(json.dumps(manifest_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_checksums_tsv(bundle_root / "checksums.tsv", bundled_files)
    _write_environment_json(bundle_root / "environment.json")
    _write_bundle_readme(
        bundle_root / "README.md",
        input_roots=normalized_inputs,
        output_roots=normalized_outputs,
        file_count=len(bundled_files),
    )

    return EvidenceBundleReport(
        input_roots=normalized_inputs,
        output_roots=normalized_outputs,
        output_root=bundle_root,
        file_count=len(bundled_files),
        input_file_count=len(input_files),
        output_file_count=len(output_files),
        files=bundled_files,
    )


def bundle_file_paths(input_paths: list[Path], output_paths: list[Path], bundle_root: Path) -> EvidenceBundleReport:
    """Copy explicit input and output files into a checksummed evidence bundle."""
    if not input_paths:
        raise EvidenceContractError("evidence bundle requires at least one input file")
    if not output_paths:
        raise EvidenceContractError("evidence bundle requires at least one output file")

    normalized_inputs = [Path(path) for path in input_paths]
    normalized_outputs = [Path(path) for path in output_paths]
    if bundle_root.exists():
        shutil.rmtree(bundle_root)
    (bundle_root / "inputs").mkdir(parents=True, exist_ok=True)
    (bundle_root / "outputs").mkdir(parents=True, exist_ok=True)

    input_files = _copy_paths(paths=normalized_inputs, section="inputs", bundle_root=bundle_root)
    output_files = _copy_paths(paths=normalized_outputs, section="outputs", bundle_root=bundle_root)
    bundled_files = input_files + output_files

    manifest_payload = {
        "input_roots": [str(path) for path in normalized_inputs],
        "output_roots": [str(path) for path in normalized_outputs],
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
    (bundle_root / "manifest.json").write_text(json.dumps(manifest_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_checksums_tsv(bundle_root / "checksums.tsv", bundled_files)
    _write_environment_json(bundle_root / "environment.json")
    _write_bundle_readme(
        bundle_root / "README.md",
        input_roots=normalized_inputs,
        output_roots=normalized_outputs,
        file_count=len(bundled_files),
    )

    return EvidenceBundleReport(
        input_roots=normalized_inputs,
        output_roots=normalized_outputs,
        output_root=bundle_root,
        file_count=len(bundled_files),
        input_file_count=len(input_files),
        output_file_count=len(output_files),
        files=bundled_files,
    )


def validate_bundle(bundle_root: Path) -> EvidenceValidationReport:
    """Validate an evidence bundle against its manifest checksums."""
    manifest_path = bundle_root / "manifest.json"
    checksums_path = bundle_root / "checksums.tsv"
    environment_path = bundle_root / "environment.json"
    readme_path = bundle_root / "README.md"
    inputs_root = bundle_root / "inputs"
    outputs_root = bundle_root / "outputs"
    if not manifest_path.exists():
        raise EvidenceContractError(f"evidence manifest not found: {manifest_path}")
    if not checksums_path.exists():
        raise EvidenceContractError(f"evidence checksums file not found: {checksums_path}")
    if not environment_path.exists():
        raise EvidenceContractError(f"evidence environment file not found: {environment_path}")
    if not readme_path.exists():
        raise EvidenceContractError(f"evidence readme not found: {readme_path}")
    if not inputs_root.exists():
        raise EvidenceContractError(f"evidence inputs directory not found: {inputs_root}")
    if not outputs_root.exists():
        raise EvidenceContractError(f"evidence outputs directory not found: {outputs_root}")

    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest_files = payload.get("files")
    if not isinstance(manifest_files, list):
        raise EvidenceContractError(f"evidence manifest is missing a files list: {manifest_path}")

    mismatches: list[EvidenceMismatch] = []
    for entry in manifest_files:
        section = str(entry["section"])
        relative_path = Path(str(entry["relative_path"]))
        expected_sha = str(entry["sha256"])
        expected_size = int(entry["size_bytes"])
        if section not in {"inputs", "outputs"}:
            mismatches.append(EvidenceMismatch(relative_path=relative_path, reason="invalid_section"))
            continue
        candidate = bundle_root / section / relative_path
        if not candidate.exists():
            mismatches.append(EvidenceMismatch(relative_path=relative_path, reason="missing_file"))
            continue
        if candidate.stat().st_size != expected_size:
            mismatches.append(EvidenceMismatch(relative_path=relative_path, reason="size_mismatch"))
            continue
        if _sha256(candidate) != expected_sha:
            mismatches.append(EvidenceMismatch(relative_path=relative_path, reason="checksum_mismatch"))

    return EvidenceValidationReport(
        bundle_root=bundle_root,
        file_count=len(manifest_files),
        valid=not mismatches,
        mismatches=mismatches,
    )
