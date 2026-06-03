from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.compare import compare_tree_paths, compare_tree_structurally
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.io.trees import load_tree

from ..support import (
    SUPPORTED_PUBLICATION_PACKAGE_KIND,
    entry_path,
    ignored_package_prefixes,
    mapping,
    read_manifest,
    text,
)
from .comparison_policy import (
    check_row,
    config_differences,
    finding_difference_count,
    status,
)
from .contracts import (
    PublicationPackageComparisonCheckRow,
    PublicationPackageComparisonResult,
)
from .inventory import (
    inventory_rows_from_manifest,
    load_accession_ids,
    load_scientific_findings,
    load_sequence_ids,
    package_artifact_rows,
)
from .presentation import write_html_report

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

    left_inventory_rows = inventory_rows_from_manifest(left_package_root, left_manifest)
    right_inventory_rows = inventory_rows_from_manifest(
        right_package_root,
        right_manifest,
    )
    artifact_rows = package_artifact_rows(
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
        check_row(
            section="manifest",
            check_id="supported-package-kinds",
            status=status(blocked=not supported),
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
        check_row(
            section="manifest",
            check_id="same-study-dataset-id",
            status=status(blocked=not dataset_match),
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
    left_accessions = load_accession_ids(left_accession_path)
    right_accessions = load_accession_ids(right_accession_path)
    left_only_accessions = sorted(left_accessions - right_accessions)
    right_only_accessions = sorted(right_accessions - left_accessions)

    left_sequences_path = left_package_root / entry_path(
        mapping(mapping(left_manifest, "dataset_files"), "sequences")
    )
    right_sequences_path = right_package_root / entry_path(
        mapping(mapping(right_manifest, "dataset_files"), "sequences")
    )
    left_sequences = load_sequence_ids(left_sequences_path)
    right_sequences = load_sequence_ids(right_sequences_path)
    left_only_sequences = sorted(left_sequences - right_sequences)
    right_only_sequences = sorted(right_sequences - left_sequences)
    check_rows.append(
        check_row(
            section="inputs",
            check_id="taxa-and-accessions",
            status=status(
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
            left_artifact_path=left_accession_path.relative_to(
                left_package_root
            ).as_posix(),
            right_artifact_path=right_accession_path.relative_to(
                right_package_root
            ).as_posix(),
        )
    )

    config_difference_map = config_differences(
        left_manifest,
        right_manifest,
        mapping=mapping,
    )
    left_config_path = text(mapping(left_manifest, "config").get("path"))
    right_config_path = text(mapping(right_manifest, "config").get("path"))
    check_rows.append(
        check_row(
            section="config",
            check_id="workflow-config",
            status=status(risk=bool(config_difference_map)),
            summary="workflow settings remain stable across study versions",
            evidence=(
                "no config differences"
                if not config_difference_map
                else " | ".join(
                    f"{key}: {left_value!r} -> {right_value!r}"
                    for key, (left_value, right_value) in sorted(
                        config_difference_map.items()
                    )
                )
            ),
            left_artifact_path=left_config_path,
            right_artifact_path=right_config_path,
        )
    )

    alignment_rows = [row for row in artifact_rows if row.kind == "alignment"]
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
        check_row(
            section="alignment",
            check_id="alignment-surfaces",
            status=status(risk=bool(alignment_differences)),
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
        check_row(
            section="tree",
            check_id="rooted-tree",
            status=status(risk=tree_changed),
            summary="rooted tree topology, rooting, and branch-structure remain stable across study versions",
            evidence=(
                f"topology_equal={str(topology.topology_equal).lower()}; "
                f"same_taxa_different_rooting={str(topology.same_taxa_different_rooting).lower()}; "
                f"same_topology_different_branch_lengths={str(topology.same_topology_different_branch_lengths).lower()}; "
                f"structural_parity={str(structural.equivalent).lower()}"
            ),
            left_artifact_path=left_tree_path.relative_to(left_package_root).as_posix(),
            right_artifact_path=right_tree_path.relative_to(
                right_package_root
            ).as_posix(),
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
        check_row(
            section="models",
            check_id="inference-and-comparative-models",
            status=status(risk=bool(model_differences)),
            summary="selected inference and comparative model surfaces remain stable across study versions",
            evidence=(
                "no model differences"
                if not model_differences
                else " | ".join(model_differences)
            ),
            left_artifact_path=text(
                mapping(mapping(left_manifest, "workflow_files"), "model_table").get(
                    "path"
                )
            ),
            right_artifact_path=text(
                mapping(mapping(right_manifest, "workflow_files"), "model_table").get(
                    "path"
                )
            ),
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
        check_row(
            section="review-surfaces",
            check_id="figure-and-report-artifacts",
            status=status(risk=bool(figure_or_report_rows)),
            summary="reviewer-facing figures and report surfaces remain stable across study versions",
            evidence=(
                "no figure or report differences"
                if not figure_or_report_rows
                else " | ".join(
                    f"{row.relative_path}: {row.status}"
                    for row in figure_or_report_rows[:10]
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
    left_findings = load_scientific_findings(left_findings_path)
    right_findings = load_scientific_findings(right_findings_path)
    scientific_finding_difference_count = finding_difference_count(
        left_findings,
        right_findings,
    )
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
        check_row(
            section="conclusions",
            check_id="scientific-findings-and-summary",
            status=status(
                risk=bool(
                    scientific_finding_difference_count
                    or short_answer_changed
                    or conclusion_count_differences
                )
            ),
            summary="biological conclusions and their reviewer-facing summaries remain stable across study versions",
            evidence=(
                f"finding differences={scientific_finding_difference_count}; "
                f"short_answer_changed={str(short_answer_changed).lower()}; "
                f"conclusion count differences={','.join(conclusion_count_differences) or 'none'}"
            ),
            left_artifact_path=left_findings_path.relative_to(
                left_package_root
            ).as_posix(),
            right_artifact_path=right_findings_path.relative_to(
                right_package_root
            ).as_posix(),
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
        config_difference_count=len(config_difference_map),
        sequence_left_only_count=len(left_only_sequences),
        sequence_right_only_count=len(right_only_sequences),
        accession_left_only_count=len(left_only_accessions),
        accession_right_only_count=len(right_only_accessions),
        alignment_difference_count=len(alignment_differences),
        figure_or_report_difference_count=len(figure_or_report_rows),
        scientific_finding_difference_count=scientific_finding_difference_count,
        overall_comparison_status=overall_comparison_status,
    )
    summary_payload = {
        "report_kind": report_kind,
        "dataset_id": dataset_id,
        "same_artifact_count": same_artifact_count,
        "changed_artifact_count": changed_artifact_count,
        "left_only_artifact_count": left_only_artifact_count,
        "right_only_artifact_count": right_only_artifact_count,
        "config_difference_count": len(config_difference_map),
        "sequence_left_only_count": len(left_only_sequences),
        "sequence_right_only_count": len(right_only_sequences),
        "accession_left_only_count": len(left_only_accessions),
        "accession_right_only_count": len(right_only_accessions),
        "alignment_difference_count": len(alignment_differences),
        "figure_or_report_difference_count": len(figure_or_report_rows),
        "scientific_finding_difference_count": scientific_finding_difference_count,
        "overall_comparison_status": overall_comparison_status,
    }
    summary_path.write_text(
        json.dumps(summary_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_html_report(report_path, result=result)
    return result
