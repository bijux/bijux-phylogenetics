from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path


def sha256_file(path: Path) -> str:
    """Return the hex digest for one package artifact."""
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def to_json_ready(payload: object) -> object:
    """Normalize nested dataclass payloads to JSON-safe primitives."""
    return json.loads(json.dumps(payload, default=str))


def build_machine_manifest(
    *,
    tree_set_path: Path,
    artifact_paths: list[Path],
    reproducibility_manifest_path: Path,
    reproducibility_manifest,
    layout: str,
    plot_row_limit: int,
    summary,
    processing,
    budget_report,
    consensus,
    multimodality,
    conflicts,
    conclusions,
    methods_summary,
    audit,
    methods_summary_path: Path,
) -> dict[str, object]:
    """Build the machine-readable manifest for the uncertainty figure package."""
    return {
        "report_kind": "tree_set_uncertainty_figure_package",
        "source_path": str(tree_set_path),
        "input_checksums": {str(tree_set_path): sha256_file(tree_set_path)},
        "output_paths": [str(path) for path in artifact_paths],
        "output_checksums": {str(path): sha256_file(path) for path in artifact_paths},
        "reproducibility_manifest_path": str(reproducibility_manifest_path),
        "reproducibility_manifest_checksum": sha256_file(reproducibility_manifest_path),
        "reproducibility_manifest": reproducibility_manifest,
        "layout": layout,
        "plot_row_limit": plot_row_limit,
        "processing": asdict(processing),
        "budget": asdict(budget_report),
        "consensus": to_json_ready(asdict(consensus)),
        "multimodality": to_json_ready(asdict(multimodality)),
        "clade_conflicts": to_json_ready(asdict(conflicts)),
        "conclusions": to_json_ready(asdict(conclusions)),
        "methods_summary": to_json_ready(asdict(methods_summary)),
        "audit": to_json_ready(asdict(audit)),
        "outputs": {"methods_summary_path": str(methods_summary_path)},
        "metrics": {"methods_summary_warning_count": methods_summary.warning_count},
        "linked_artifact_count": len(artifact_paths),
        "tree_count": summary.tree_count,
    }
