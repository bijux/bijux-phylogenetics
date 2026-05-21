from __future__ import annotations

from pathlib import Path


def scan_storage_inventory(
    *,
    bundle_root: Path,
    configured_variant_ids: list[str],
) -> tuple[dict[str, dict[str, int]], dict[str, int]]:
    variant_totals = {
        variant_id: _empty_category_totals() for variant_id in configured_variant_ids
    }
    shared_totals = _empty_category_totals()

    for path in sorted(bundle_root.rglob("*")):
        if not path.is_file():
            continue
        relative_path = path.relative_to(bundle_root)
        category_id, variant_id = classify_storage_path(relative_path)
        byte_count = path.stat().st_size
        if variant_id is None:
            totals = shared_totals
        else:
            totals = variant_totals.setdefault(variant_id, _empty_category_totals())
        file_key = f"{category_id}_file_count"
        byte_key = f"{category_id}_byte_count"
        totals[file_key] += 1
        totals[byte_key] += byte_count
        totals["total_file_count"] += 1
        totals["total_byte_count"] += byte_count
    return variant_totals, shared_totals


def classify_storage_path(relative_path: Path) -> tuple[str, str | None]:
    parts = relative_path.parts
    if not parts:
        raise ValueError("relative_path must not be empty")
    top_level = parts[0]
    if top_level == "parallel-logs":
        return "logs", relative_path.stem
    if top_level == "variants" and len(parts) >= 3:
        return classify_variant_file(relative_path.name), parts[1]
    if top_level == "slurm-job-evidence" and len(parts) >= 3:
        return "reports", parts[1]
    if top_level in {"report-artifacts", "slurm-arrays"}:
        return "reports", None
    return "reports", None


def classify_variant_file(filename: str) -> str:
    normalized = filename.lower()
    if normalized.endswith((".nwk", ".tree", ".treefile", ".contree", ".ufboot")):
        return "trees"
    if _looks_like_posterior_sample(normalized):
        return "posterior_samples"
    return "outputs"


def _looks_like_posterior_sample(filename: str) -> bool:
    if filename.endswith((".trees", ".state", ".trace", ".p", ".t")):
        return True
    return any(token in filename for token in ("posterior", "mcmc", "beast", "mrbayes"))


def _empty_category_totals() -> dict[str, int]:
    return {
        "outputs_file_count": 0,
        "logs_file_count": 0,
        "trees_file_count": 0,
        "posterior_samples_file_count": 0,
        "reports_file_count": 0,
        "outputs_byte_count": 0,
        "logs_byte_count": 0,
        "trees_byte_count": 0,
        "posterior_samples_byte_count": 0,
        "reports_byte_count": 0,
        "total_file_count": 0,
        "total_byte_count": 0,
    }
