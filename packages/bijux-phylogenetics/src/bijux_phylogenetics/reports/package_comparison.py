from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from html import escape
import json
from pathlib import Path

from bijux_phylogenetics.compare import compare_tree_structurally, compare_tree_paths
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.io.fasta._shared import load_fasta_records
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.reports.publication_package_support import (
    SUPPORTED_PUBLICATION_PACKAGE_KIND,
    artifact_kind,
    checksum,
    entry_path,
    ignored_package_prefixes,
    mapping,
    read_manifest,
    read_tsv_rows,
    text,
)


_ARTIFACT_COLUMNS = [
    "section",
    "kind",
    "relative_path",
    "status",
    "left_sha256",
    "right_sha256",
    "left_size_bytes",
    "right_size_bytes",
    "detail",
]
_CHECK_COLUMNS = [
    "section",
    "check_id",
    "status",
    "summary",
    "evidence",
    "left_artifact_path",
    "right_artifact_path",
]


@dataclass(frozen=True, slots=True)
class PublicationPackageComparisonArtifactRow:
    """One artifact-level difference row between two stored study packages."""

    section: str
    kind: str
    relative_path: str
    status: str
    left_sha256: str | None
    right_sha256: str | None
    left_size_bytes: int | None
    right_size_bytes: int | None
    detail: str


@dataclass(frozen=True, slots=True)
class PublicationPackageComparisonCheckRow:
    """One reviewer-facing comparison decision over two package versions."""

    section: str
    check_id: str
    status: str
    summary: str
    evidence: str
    left_artifact_path: str
    right_artifact_path: str


@dataclass(slots=True)
class PublicationPackageComparisonResult:
    """Written comparison artifacts for two stored publication packages."""

    output_root: Path
    left_manifest_path: Path
    right_manifest_path: Path
    left_package_root: Path
    right_package_root: Path
    report_kind: str
    dataset_id: str
    artifact_table_path: Path
    check_table_path: Path
    summary_path: Path
    report_path: Path
    artifact_rows: list[PublicationPackageComparisonArtifactRow]
    check_rows: list[PublicationPackageComparisonCheckRow]
    same_artifact_count: int
    changed_artifact_count: int
    left_only_artifact_count: int
    right_only_artifact_count: int
    config_difference_count: int
    sequence_left_only_count: int
    sequence_right_only_count: int
    accession_left_only_count: int
    accession_right_only_count: int
    alignment_difference_count: int
    figure_or_report_difference_count: int
    scientific_finding_difference_count: int
    overall_comparison_status: str


def _status(*, blocked: bool = False, risk: bool = False) -> str:
    if blocked:
        return "blocked"
    if risk:
        return "risk"
    return "pass"


def _inventory_index(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {
        row["relative_path"]: row
        for row in rows
        if row.get("relative_path", "").strip()
    }


def _artifact_row(
    relative_path: str,
    left_row: dict[str, str] | None,
    right_row: dict[str, str] | None,
) -> PublicationPackageComparisonArtifactRow:
    left_sha = None if left_row is None else left_row.get("sha256", "").strip() or None
    right_sha = (
        None if right_row is None else right_row.get("sha256", "").strip() or None
    )
    left_size = (
        None
        if left_row is None or not left_row.get("size_bytes", "").strip()
        else int(left_row["size_bytes"])
    )
    right_size = (
        None
        if right_row is None or not right_row.get("size_bytes", "").strip()
        else int(right_row["size_bytes"])
    )
    section = (
        left_row.get("section", "").strip()
        if left_row is not None
        else right_row.get("section", "").strip()
        if right_row is not None
        else "artifact"
    )
    kind = (
        left_row.get("kind", "").strip()
        if left_row is not None
        else right_row.get("kind", "").strip()
        if right_row is not None
        else artifact_kind(relative_path)
    )
    if left_row is None:
        return PublicationPackageComparisonArtifactRow(
            section=section,
            kind=kind,
            relative_path=relative_path,
            status="right_only",
            left_sha256=None,
            right_sha256=right_sha,
            left_size_bytes=None,
            right_size_bytes=right_size,
            detail="artifact appears only in the right package version",
        )
    if right_row is None:
        return PublicationPackageComparisonArtifactRow(
            section=section,
            kind=kind,
            relative_path=relative_path,
            status="left_only",
            left_sha256=left_sha,
            right_sha256=None,
            left_size_bytes=left_size,
            right_size_bytes=None,
            detail="artifact appears only in the left package version",
        )
    if left_sha == right_sha and left_size == right_size:
        return PublicationPackageComparisonArtifactRow(
            section=section,
            kind=kind,
            relative_path=relative_path,
            status="same",
            left_sha256=left_sha,
            right_sha256=right_sha,
            left_size_bytes=left_size,
            right_size_bytes=right_size,
            detail="artifact matches across both package versions",
        )
    detail_parts: list[str] = []
    if left_sha != right_sha:
        detail_parts.append("checksum changed")
    if left_size != right_size:
        detail_parts.append(
            f"size changed from {left_size if left_size is not None else 'missing'} to {right_size if right_size is not None else 'missing'} bytes"
        )
    return PublicationPackageComparisonArtifactRow(
        section=section,
        kind=kind,
        relative_path=relative_path,
        status="changed",
        left_sha256=left_sha,
        right_sha256=right_sha,
        left_size_bytes=left_size,
        right_size_bytes=right_size,
        detail="; ".join(detail_parts),
    )


def _check_row(
    *,
    section: str,
    check_id: str,
    status: str,
    summary: str,
    evidence: str,
    left_artifact_path: str,
    right_artifact_path: str,
) -> PublicationPackageComparisonCheckRow:
    return PublicationPackageComparisonCheckRow(
        section=section,
        check_id=check_id,
        status=status,
        summary=summary,
        evidence=evidence,
        left_artifact_path=left_artifact_path,
        right_artifact_path=right_artifact_path,
    )


def _load_accession_ids(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {
            row["accession"].strip()
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("accession", "").strip()
        }


def _load_sequence_ids(path: Path) -> set[str]:
    return {record.identifier for record in load_fasta_records(path)}


def _load_scientific_findings(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {
            row["finding_id"]: row
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("finding_id", "").strip()
        }


def _inventory_rows_from_manifest(
    package_root: Path,
    manifest: dict[str, object],
) -> list[dict[str, str]]:
    inventory_entry = mapping(mapping(manifest, "package_files"), "artifact_inventory")
    inventory_path = package_root / entry_path(inventory_entry)
    return read_tsv_rows(inventory_path)


def _package_artifact_rows(
    *,
    left_inventory_rows: list[dict[str, str]],
    right_inventory_rows: list[dict[str, str]],
    left_manifest_path: Path,
    right_manifest_path: Path,
) -> list[PublicationPackageComparisonArtifactRow]:
    left_index = _inventory_index(left_inventory_rows)
    right_index = _inventory_index(right_inventory_rows)
    relative_paths = sorted(set(left_index) | set(right_index))
    rows = [
        _artifact_row(relative_path, left_index.get(relative_path), right_index.get(relative_path))
        for relative_path in relative_paths
    ]
    rows.append(
        PublicationPackageComparisonArtifactRow(
            section="package",
            kind="manifest",
            relative_path=left_manifest_path.name,
            status="same" if checksum(left_manifest_path) == checksum(right_manifest_path) else "changed",
            left_sha256=checksum(left_manifest_path),
            right_sha256=checksum(right_manifest_path),
            left_size_bytes=left_manifest_path.stat().st_size,
            right_size_bytes=right_manifest_path.stat().st_size,
            detail="package manifest checksum matches"
            if checksum(left_manifest_path) == checksum(right_manifest_path)
            else "package manifest checksum changed",
        )
    )
    return rows


def _config_differences(
    left_manifest: dict[str, object],
    right_manifest: dict[str, object],
) -> dict[str, tuple[object, object]]:
    left_config = mapping(left_manifest, "config")
    right_config = mapping(right_manifest, "config")
    differences: dict[str, tuple[object, object]] = {}
    for key in sorted(set(left_config) | set(right_config)):
        if key in {"path", "checksum"}:
            continue
        if left_config.get(key) != right_config.get(key):
            differences[key] = (left_config.get(key), right_config.get(key))
    return differences


def _finding_difference_count(
    left_rows: dict[str, dict[str, str]],
    right_rows: dict[str, dict[str, str]],
) -> int:
    count = 0
    for finding_id in sorted(set(left_rows) | set(right_rows)):
        left_row = left_rows.get(finding_id)
        right_row = right_rows.get(finding_id)
        if left_row is None or right_row is None:
            count += 1
            continue
        if {
            key: value
            for key, value in left_row.items()
            if key != "finding_id"
        } != {
            key: value
            for key, value in right_row.items()
            if key != "finding_id"
        }:
            count += 1
    return count


def _write_html_report(
    path: Path,
    *,
    result: PublicationPackageComparisonResult,
) -> Path:
    check_rows = "\n".join(
        (
            "      <tr>"
            f"<td>{escape(row.section)}</td>"
            f"<td>{escape(row.check_id)}</td>"
            f"<td>{escape(row.status)}</td>"
            f"<td>{escape(row.summary)}</td>"
            f"<td>{escape(row.evidence)}</td>"
            "</tr>"
        )
        for row in result.check_rows
    )
    artifact_rows = "\n".join(
        (
            "      <tr>"
            f"<td>{escape(row.relative_path)}</td>"
            f"<td>{escape(row.kind)}</td>"
            f"<td>{escape(row.status)}</td>"
            f"<td>{escape(row.detail)}</td>"
            "</tr>"
        )
        for row in result.artifact_rows
    )
    html = "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Publication Package Comparison</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: #f5f1ea; color: #1f2a24; }",
            "    main { max-width: 1120px; margin: 0 auto; padding: 28px; }",
            "    .cards { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; }",
            "    .card { background: #fffdf8; border: 1px solid #d7d0c3; border-radius: 14px; padding: 14px; }",
            "    .label { font-size: 12px; text-transform: uppercase; letter-spacing: 0.04em; color: #676154; }",
            "    .value { display: block; font-size: 21px; margin-top: 6px; }",
            "    .panel { background: #fffdf8; border: 1px solid #d7d0c3; border-radius: 14px; padding: 18px; margin-top: 18px; }",
            "    table { width: 100%; border-collapse: collapse; margin-top: 12px; }",
            "    th, td { text-align: left; padding: 8px 10px; border-bottom: 1px solid #e6dfd1; vertical-align: top; }",
            "    code { background: #f0eadf; padding: 1px 4px; border-radius: 4px; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Publication Package Comparison</h1>",
            "  <div class=\"cards\">",
            f"    <div class=\"card\"><span class=\"label\">dataset id</span><span class=\"value\">{escape(result.dataset_id)}</span></div>",
            f"    <div class=\"card\"><span class=\"label\">overall status</span><span class=\"value\">{escape(result.overall_comparison_status)}</span></div>",
            f"    <div class=\"card\"><span class=\"label\">changed artifacts</span><span class=\"value\">{result.changed_artifact_count}</span></div>",
            f"    <div class=\"card\"><span class=\"label\">finding differences</span><span class=\"value\">{result.scientific_finding_difference_count}</span></div>",
            "  </div>",
            "  <section class=\"panel\">",
            "    <h2>Checks</h2>",
            "    <table>",
            "      <thead><tr><th>Section</th><th>Check</th><th>Status</th><th>Summary</th><th>Evidence</th></tr></thead>",
            "      <tbody>",
            check_rows,
            "      </tbody>",
            "    </table>",
            "  </section>",
            "  <section class=\"panel\">",
            "    <h2>Artifacts</h2>",
            "    <table>",
            "      <thead><tr><th>Artifact</th><th>Kind</th><th>Status</th><th>Detail</th></tr></thead>",
            "      <tbody>",
            artifact_rows,
            "      </tbody>",
            "    </table>",
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )
    path.write_text(html + "\n", encoding="utf-8")
    return path


def write_publication_package_comparison_report(
    output_root: Path,
    left_manifest_path: Path,
    right_manifest_path: Path,
) -> PublicationPackageComparisonResult:
    """Compare two stored publication package versions for the same study."""

    output_root.mkdir(parents=True, exist_ok=True)
    left_manifest_path = left_manifest_path.resolve()
    right_manifest_path = right_manifest_path.resolve()
    left_package_root = left_manifest_path.parent.resolve()
    right_package_root = right_manifest_path.parent.resolve()
    left_manifest = read_manifest(left_manifest_path)
    right_manifest = read_manifest(right_manifest_path)

    left_report_kind = text(left_manifest.get("report_kind"))
    right_report_kind = text(right_manifest.get("report_kind"))
    report_kind = left_report_kind or right_report_kind
    left_dataset_id = text(left_manifest.get("dataset_id"))
    right_dataset_id = text(right_manifest.get("dataset_id"))
    dataset_id = left_dataset_id or right_dataset_id

    left_inventory_rows = _inventory_rows_from_manifest(left_package_root, left_manifest)
    right_inventory_rows = _inventory_rows_from_manifest(
        right_package_root,
        right_manifest,
    )
    artifact_rows = _package_artifact_rows(
        left_inventory_rows=left_inventory_rows,
        right_inventory_rows=right_inventory_rows,
        left_manifest_path=left_manifest_path,
        right_manifest_path=right_manifest_path,
    )
    check_rows: list[PublicationPackageComparisonCheckRow] = []

    supported = (
        left_report_kind == SUPPORTED_PUBLICATION_PACKAGE_KIND
        and right_report_kind == SUPPORTED_PUBLICATION_PACKAGE_KIND
    )
    check_rows.append(
        _check_row(
            section="manifest",
            check_id="supported-package-kinds",
            status=_status(blocked=not supported),
            summary="both package manifests use the governed rabies study package contract",
            evidence=(
                f"left report_kind={left_report_kind or 'missing'}; "
                f"right report_kind={right_report_kind or 'missing'}"
            ),
            left_artifact_path=left_manifest_path.name,
            right_artifact_path=right_manifest_path.name,
        )
    )
    dataset_match = left_dataset_id == right_dataset_id
    check_rows.append(
        _check_row(
            section="manifest",
            check_id="same-study-dataset-id",
            status=_status(blocked=not dataset_match),
            summary="both package versions describe the same governed study dataset",
            evidence=f"left dataset_id={left_dataset_id}; right dataset_id={right_dataset_id}",
            left_artifact_path=left_manifest_path.name,
            right_artifact_path=right_manifest_path.name,
        )
    )

    left_accession_path = left_package_root / entry_path(
        mapping(mapping(left_manifest, "dataset_files"), "source_accessions")
    )
    right_accession_path = right_package_root / entry_path(
        mapping(mapping(right_manifest, "dataset_files"), "source_accessions")
    )
    left_accessions = _load_accession_ids(left_accession_path)
    right_accessions = _load_accession_ids(right_accession_path)
    left_only_accessions = sorted(left_accessions - right_accessions)
    right_only_accessions = sorted(right_accessions - left_accessions)

    left_sequences_path = left_package_root / entry_path(
        mapping(mapping(left_manifest, "dataset_files"), "sequences")
    )
    right_sequences_path = right_package_root / entry_path(
        mapping(mapping(right_manifest, "dataset_files"), "sequences")
    )
    left_sequences = _load_sequence_ids(left_sequences_path)
    right_sequences = _load_sequence_ids(right_sequences_path)
    left_only_sequences = sorted(left_sequences - right_sequences)
    right_only_sequences = sorted(right_sequences - left_sequences)
    check_rows.append(
        _check_row(
            section="inputs",
            check_id="taxa-and-accessions",
            status=_status(
                risk=bool(
                    left_only_accessions
                    or right_only_accessions
                    or left_only_sequences
                    or right_only_sequences
                )
            ),
            summary="input accessions and sequence identifiers remain aligned across study versions",
            evidence=(
                f"left-only accessions={len(left_only_accessions)}; right-only accessions={len(right_only_accessions)}; "
                f"left-only sequences={len(left_only_sequences)}; right-only sequences={len(right_only_sequences)}"
            ),
            left_artifact_path=left_accession_path.relative_to(left_package_root).as_posix(),
            right_artifact_path=right_accession_path.relative_to(right_package_root).as_posix(),
        )
    )

    config_differences = _config_differences(left_manifest, right_manifest)
    left_config_path = text(mapping(left_manifest, "config").get("path"))
    right_config_path = text(mapping(right_manifest, "config").get("path"))
    check_rows.append(
        _check_row(
            section="config",
            check_id="workflow-config",
            status=_status(risk=bool(config_differences)),
            summary="workflow settings remain stable across study versions",
            evidence=(
                "no config differences"
                if not config_differences
                else " | ".join(
                    f"{key}: {left_value!r} -> {right_value!r}"
                    for key, (left_value, right_value) in sorted(config_differences.items())
                )
            ),
            left_artifact_path=left_config_path,
            right_artifact_path=right_config_path,
        )
    )

    alignment_rows = [
        row for row in artifact_rows if row.kind == "alignment"
    ]
    alignment_differences: list[str] = []
    for row in alignment_rows:
        if row.status != "same":
            alignment_differences.append(f"{row.relative_path}: {row.status}")
            continue
        left_summary = summarise_fasta(left_package_root / row.relative_path)
        right_summary = summarise_fasta(right_package_root / row.relative_path)
        if (
            left_summary.sequence_count != right_summary.sequence_count
            or left_summary.alignment_length != right_summary.alignment_length
        ):
            alignment_differences.append(
                f"{row.relative_path}: sequences {left_summary.sequence_count}->{right_summary.sequence_count}, "
                f"sites {left_summary.alignment_length}->{right_summary.alignment_length}"
            )
    check_rows.append(
        _check_row(
            section="alignment",
            check_id="alignment-surfaces",
            status=_status(risk=bool(alignment_differences)),
            summary="alignment inputs and intermediate alignment artifacts remain stable across study versions",
            evidence=(
                "no alignment differences"
                if not alignment_differences
                else " | ".join(alignment_differences[:10])
            ),
            left_artifact_path="dataset/sequences.fasta",
            right_artifact_path="dataset/sequences.fasta",
        )
    )

    left_tree_path = left_package_root / entry_path(
        mapping(mapping(left_manifest, "workflow_files"), "rooted_tree")
    )
    right_tree_path = right_package_root / entry_path(
        mapping(mapping(right_manifest, "workflow_files"), "rooted_tree")
    )
    topology = compare_tree_paths(left_tree_path, right_tree_path)
    structural = compare_tree_structurally(
        load_tree(left_tree_path),
        load_tree(right_tree_path),
    )
    tree_changed = (
        not topology.topology_equal
        or topology.same_taxa_different_rooting
        or topology.same_topology_different_branch_lengths
        or not structural.equivalent
    )
    check_rows.append(
        _check_row(
            section="tree",
            check_id="rooted-tree",
            status=_status(risk=tree_changed),
            summary="rooted tree topology, rooting, and branch-structure remain stable across study versions",
            evidence=(
                f"topology_equal={str(topology.topology_equal).lower()}; "
                f"same_taxa_different_rooting={str(topology.same_taxa_different_rooting).lower()}; "
                f"same_topology_different_branch_lengths={str(topology.same_topology_different_branch_lengths).lower()}; "
                f"structural_parity={str(structural.equivalent).lower()}"
            ),
            left_artifact_path=left_tree_path.relative_to(left_package_root).as_posix(),
            right_artifact_path=right_tree_path.relative_to(right_package_root).as_posix(),
        )
    )

    left_metrics = mapping(left_manifest, "metrics")
    right_metrics = mapping(right_manifest, "metrics")
    model_differences: list[str] = []
    for key in (
        "selected_model",
        "comparative_selected_model",
        "comparative_pgls_lambda",
        "comparative_pgls_r_squared",
    ):
        if left_metrics.get(key) != right_metrics.get(key):
            model_differences.append(
                f"{key}: {left_metrics.get(key)!r} -> {right_metrics.get(key)!r}"
            )
    check_rows.append(
        _check_row(
            section="models",
            check_id="inference-and-comparative-models",
            status=_status(risk=bool(model_differences)),
            summary="selected inference and comparative model surfaces remain stable across study versions",
            evidence=(
                "no model differences"
                if not model_differences
                else " | ".join(model_differences)
            ),
            left_artifact_path=text(mapping(mapping(left_manifest, "workflow_files"), "model_table").get("path")),
            right_artifact_path=text(mapping(mapping(right_manifest, "workflow_files"), "model_table").get("path")),
        )
    )

    figure_or_report_rows = [
        row
        for row in artifact_rows
        if row.kind in {"figure", "report", "markdown"}
        and row.status != "same"
        and not any(
            row.relative_path.startswith(prefix)
            for prefix in ignored_package_prefixes(report_kind)
        )
    ]
    check_rows.append(
        _check_row(
            section="review-surfaces",
            check_id="figure-and-report-artifacts",
            status=_status(risk=bool(figure_or_report_rows)),
            summary="reviewer-facing figures and report surfaces remain stable across study versions",
            evidence=(
                "no figure or report differences"
                if not figure_or_report_rows
                else " | ".join(
                    f"{row.relative_path}: {row.status}" for row in figure_or_report_rows[:10]
                )
            ),
            left_artifact_path="rabies-cross-host-geography-overview.html",
            right_artifact_path="rabies-cross-host-geography-overview.html",
        )
    )

    left_findings_path = left_package_root / entry_path(
        mapping(mapping(left_manifest, "workflow_files"), "scientific_findings")
    )
    right_findings_path = right_package_root / entry_path(
        mapping(mapping(right_manifest, "workflow_files"), "scientific_findings")
    )
    left_findings = _load_scientific_findings(left_findings_path)
    right_findings = _load_scientific_findings(right_findings_path)
    finding_difference_count = _finding_difference_count(left_findings, right_findings)
    short_answer_changed = text(left_manifest.get("short_answer")) != text(
        right_manifest.get("short_answer")
    )
    conclusion_count_differences = [
        key
        for key in (
            "conclusion_stable_count",
            "conclusion_weak_count",
            "conclusion_unstable_count",
        )
        if left_metrics.get(key) != right_metrics.get(key)
    ]
    check_rows.append(
        _check_row(
            section="conclusions",
            check_id="scientific-findings-and-summary",
            status=_status(
                risk=bool(
                    finding_difference_count
                    or short_answer_changed
                    or conclusion_count_differences
                )
            ),
            summary="biological conclusions and their reviewer-facing summaries remain stable across study versions",
            evidence=(
                f"finding differences={finding_difference_count}; "
                f"short_answer_changed={str(short_answer_changed).lower()}; "
                f"conclusion count differences={','.join(conclusion_count_differences) or 'none'}"
            ),
            left_artifact_path=left_findings_path.relative_to(left_package_root).as_posix(),
            right_artifact_path=right_findings_path.relative_to(right_package_root).as_posix(),
        )
    )

    same_artifact_count = sum(1 for row in artifact_rows if row.status == "same")
    changed_artifact_count = sum(1 for row in artifact_rows if row.status == "changed")
    left_only_artifact_count = sum(
        1 for row in artifact_rows if row.status == "left_only"
    )
    right_only_artifact_count = sum(
        1 for row in artifact_rows if row.status == "right_only"
    )
    blocked_check_count = sum(1 for row in check_rows if row.status == "blocked")
    risk_check_count = sum(1 for row in check_rows if row.status == "risk")
    overall_comparison_status = (
        "blocked"
        if blocked_check_count > 0
        else "risk"
        if risk_check_count > 0
        else "pass"
    )

    artifact_table_path = write_taxon_rows(
        output_root / "publication-package-comparison-artifacts.tsv",
        columns=_ARTIFACT_COLUMNS,
        rows=[asdict(row) for row in artifact_rows],
    )
    check_table_path = write_taxon_rows(
        output_root / "publication-package-comparison-checks.tsv",
        columns=_CHECK_COLUMNS,
        rows=[asdict(row) for row in check_rows],
    )
    summary_path = output_root / "publication-package-comparison-summary.json"
    report_path = output_root / "publication-package-comparison-report.html"

    result = PublicationPackageComparisonResult(
        output_root=output_root,
        left_manifest_path=left_manifest_path,
        right_manifest_path=right_manifest_path,
        left_package_root=left_package_root,
        right_package_root=right_package_root,
        report_kind=report_kind,
        dataset_id=dataset_id,
        artifact_table_path=artifact_table_path,
        check_table_path=check_table_path,
        summary_path=summary_path,
        report_path=report_path,
        artifact_rows=artifact_rows,
        check_rows=check_rows,
        same_artifact_count=same_artifact_count,
        changed_artifact_count=changed_artifact_count,
        left_only_artifact_count=left_only_artifact_count,
        right_only_artifact_count=right_only_artifact_count,
        config_difference_count=len(config_differences),
        sequence_left_only_count=len(left_only_sequences),
        sequence_right_only_count=len(right_only_sequences),
        accession_left_only_count=len(left_only_accessions),
        accession_right_only_count=len(right_only_accessions),
        alignment_difference_count=len(alignment_differences),
        figure_or_report_difference_count=len(figure_or_report_rows),
        scientific_finding_difference_count=finding_difference_count,
        overall_comparison_status=overall_comparison_status,
    )
    summary_payload = {
        "report_kind": report_kind,
        "dataset_id": dataset_id,
        "same_artifact_count": same_artifact_count,
        "changed_artifact_count": changed_artifact_count,
        "left_only_artifact_count": left_only_artifact_count,
        "right_only_artifact_count": right_only_artifact_count,
        "config_difference_count": len(config_differences),
        "sequence_left_only_count": len(left_only_sequences),
        "sequence_right_only_count": len(right_only_sequences),
        "accession_left_only_count": len(left_only_accessions),
        "accession_right_only_count": len(right_only_accessions),
        "alignment_difference_count": len(alignment_differences),
        "figure_or_report_difference_count": len(figure_or_report_rows),
        "scientific_finding_difference_count": finding_difference_count,
        "overall_comparison_status": overall_comparison_status,
    }
    summary_path.write_text(
        json.dumps(summary_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_html_report(report_path, result=result)
    return result
