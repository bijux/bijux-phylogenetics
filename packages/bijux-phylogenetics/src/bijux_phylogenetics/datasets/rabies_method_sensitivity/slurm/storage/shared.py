from __future__ import annotations

import csv
import json
from math import ceil
from pathlib import Path

_CONFIG_FILENAME = "workflow-config.resolved.json"
_MEBIBYTE = 1024 * 1024
_STORAGE_CATEGORIES = (
    "outputs",
    "logs",
    "trees",
    "posterior_samples",
    "reports",
)
_CATEGORY_LABELS = {
    "outputs": "workflow outputs",
    "logs": "canonical logs",
    "trees": "tree artifacts",
    "posterior_samples": "posterior samples",
    "reports": "review artifacts",
}
_CATEGORY_DETAILS = {
    "outputs": (
        "Retained scientific outputs such as alignments, comparison ledgers, and "
        "variant-scoped tables that are not tree sets or reviewer reports."
    ),
    "logs": (
        "Canonical orchestration logs under parallel-logs. Copied task logs inside "
        "per-job evidence packages are counted under review artifacts to avoid "
        "double-counting the same execution record as canonical storage."
    ),
    "trees": (
        "Tree products such as rooted or unrooted Newick outputs and other retained "
        "tree-topology artifacts."
    ),
    "posterior_samples": (
        "Bayesian posterior chains or sampled tree sets. This governed workflow does "
        "not currently emit them, but the category is kept explicit so zero remains "
        "an audited, reviewer-visible statement."
    ),
    "reports": (
        "Workflow-wide summaries, manifests, HTML reports, Slurm ledgers, evidence "
        "packages, and other reviewer-facing accountability artifacts."
    ),
}


def _to_mib(byte_count: int) -> int:
    return max(0, ceil(byte_count / _MEBIBYTE))


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_tsv(
    path: Path,
    *,
    fieldnames: tuple[str, ...],
    rows: list[dict[str, object]],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path
