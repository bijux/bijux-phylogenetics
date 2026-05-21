from __future__ import annotations

from dataclasses import asdict, is_dataclass
import json
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.registry import COMMAND_SPECS
from bijux_phylogenetics.evidence.closure import (
    build_claim_reaudit,
    build_closure_criteria,
    build_completion_gates,
    build_evidence_maturity_scorecard,
)
from bijux_phylogenetics.evidence.coverage import build_evidence_coverage_gap_report
from bijux_phylogenetics.evidence.freshness import build_evidence_freshness_report
from bijux_phylogenetics.evidence.integrity import build_evidence_integrity_report
from bijux_phylogenetics.runtime.results import build_command_result


def _json_ready(value: Any) -> Any:
    if is_dataclass(value):
        return {key: _json_ready(item) for key, item in asdict(value).items()}
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    return value


def _print_result(result: Any, *, json_output: bool) -> None:
    if json_output:
        print(json.dumps(_json_ready(result), indent=2, sort_keys=True))
        return
    if isinstance(result, str):
        print(result)
        return
    print(json.dumps(_json_ready(result), indent=2, sort_keys=True))


def _tsv_field(value: Any) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return format(value, ".12g")
    return str(value)


def _write_tsv(path: Path, *, header: list[str], rows: list[list[Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)]
    lines.extend("\t".join(_tsv_field(value) for value in row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_locus_occupancy_taxa_tsv(path: Path, report: Any) -> Path:
    return _write_tsv(
        path,
        header=[
            "taxon",
            "covered_locus_count",
            "total_locus_count",
            "locus_coverage_fraction",
            "observed_site_count",
            "total_site_count",
            "site_coverage_fraction",
            "low_coverage",
        ],
        rows=[
            [
                row.taxon,
                row.covered_locus_count,
                row.total_locus_count,
                row.locus_coverage_fraction,
                row.observed_site_count,
                row.total_site_count,
                row.site_coverage_fraction,
                row.low_coverage,
            ]
            for row in report.taxa
        ],
    )


def _write_locus_occupancy_loci_tsv(path: Path, report: Any) -> Path:
    return _write_tsv(
        path,
        header=[
            "locus_name",
            "covered_taxon_count",
            "total_taxa",
            "taxon_coverage_fraction",
            "observed_site_count",
            "total_site_count",
            "site_coverage_fraction",
            "low_coverage",
        ],
        rows=[
            [
                row.locus_name,
                row.covered_taxon_count,
                row.total_taxa,
                row.taxon_coverage_fraction,
                row.observed_site_count,
                row.total_site_count,
                row.site_coverage_fraction,
                row.low_coverage,
            ]
            for row in report.loci
        ],
    )


def _write_locus_occupancy_matrix_tsv(path: Path, report: Any) -> Path:
    locus_names = [partition.name for partition in report.partitions]
    return _write_tsv(
        path,
        header=[
            "taxon",
            *locus_names,
            "covered_locus_count",
            "total_locus_count",
            "locus_coverage_fraction",
            "observed_site_count",
            "total_site_count",
            "site_coverage_fraction",
            "low_coverage",
        ],
        rows=[
            [
                row.taxon,
                *[row.occupancies[name] for name in locus_names],
                row.covered_locus_count,
                row.total_locus_count,
                row.locus_coverage_fraction,
                row.observed_site_count,
                row.total_site_count,
                row.site_coverage_fraction,
                row.low_coverage,
            ]
            for row in report.taxa
        ],
    )


def _evidence_book_metrics(repo_root: Path) -> dict[str, object]:
    freshness_report = build_evidence_freshness_report(repo_root)
    integrity_report = build_evidence_integrity_report(repo_root)
    coverage_report = build_evidence_coverage_gap_report(repo_root)
    claim_reaudit = build_claim_reaudit(repo_root)
    closure_report = build_closure_criteria(repo_root)
    scorecard = build_evidence_maturity_scorecard(repo_root)
    completion_gates = build_completion_gates(repo_root)
    freshness_counts = freshness_report["freshness_status_counts"]
    integrity_counts = integrity_report["integrity_status_counts"]
    foundational_status = next(
        criterion["current_status"]
        for criterion in closure_report["criteria"]
        if criterion["criterion_id"] == "foundational-numerical-trust"
    )
    reviewer_status = next(
        criterion["current_status"]
        for criterion in closure_report["criteria"]
        if criterion["criterion_id"] == "reviewer-readiness"
    )
    completion_state_counts = completion_gates["completion_state_counts"]
    return {
        "bundle_count": int(freshness_report["bundle_count"]),
        "freshness_current_count": int(freshness_counts.get("current", 0)),
        "freshness_stale_count": int(freshness_counts.get("stale", 0)),
        "freshness_source_unresolved_count": int(
            freshness_counts.get("source_unresolved", 0)
        ),
        "integrity_tracked_count": int(integrity_counts.get("tracked", 0)),
        "coverage_gap_count": int(coverage_report["coverage_gap_count"]),
        "family_gap_count": int(coverage_report["family_gap_count"]),
        "downgraded_claim_count": int(claim_reaudit["downgraded_claim_count"]),
        "foundational_numerical_trust_status": str(foundational_status),
        "reviewer_readiness_status": str(reviewer_status),
        "maturity_tier": str(scorecard["maturity_tier"]),
        "completion_bounded_count": int(completion_state_counts.get("bounded", 0)),
        "completion_not_ready_count": int(completion_state_counts.get("not_ready", 0)),
    }


def _print_commands(*, output_format: str) -> None:
    payload = build_command_result(
        command="commands",
        inputs=[],
        outputs=[],
        metrics={"command_count": len(COMMAND_SPECS)},
        data={"commands": list(COMMAND_SPECS)},
    )
    if output_format == "json":
        print(json.dumps(_json_ready(payload), indent=2, sort_keys=True))
        return
    for command in _json_ready(payload.data)["commands"]:
        print(f"{command['name']}: {command['domain']} - {command['summary']}")
