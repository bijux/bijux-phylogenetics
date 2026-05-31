from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path

from bijux_phylogenetics.ancestral.continuous import (
    continuous_ancestral_exclusions,
    reconstruct_continuous_ancestral_states,
    summarize_continuous_ancestral_report,
    write_continuous_ancestral_exclusion_table,
    write_continuous_ancestral_summary_table,
    write_continuous_ancestral_uncertainty_table,
)
from bijux_phylogenetics.ancestral.discrete import (
    discrete_ancestral_exclusions,
    reconstruct_discrete_ancestral_states,
    summarize_discrete_ancestral_report,
    write_discrete_ancestral_exclusion_table,
    write_discrete_ancestral_probability_table,
    write_discrete_ancestral_summary_table,
)
from bijux_phylogenetics.ancestral.discrete.review import (
    summarize_ancestral_transitions,
    write_ancestral_transition_branch_table,
    write_ancestral_transition_count_table,
    write_ancestral_transition_exclusion_table,
)
from bijux_phylogenetics.ancestral.presentation.methods_text import (
    write_ancestral_methods_summary_text,
)
from bijux_phylogenetics.ancestral.presentation.report_rendering import (
    render_ancestral_state_report,
    write_ancestral_state_table,
)
from bijux_phylogenetics.ancestral.presentation.visualization import (
    render_ancestral_state_visualization,
)
from bijux_phylogenetics.reports.review import write_reviewer_audit_checklist

from .artifact_outputs import (
    write_continuous_change_branch_table,
    write_continuous_change_count_table,
)
from .continuous_changes import (
    summarize_continuous_change_branches,
    summarize_continuous_change_counts,
)
from .contracts import AncestralReportPackageResult
from .presentation import write_package_html
from .review_rows import limitations as review_limitations
from .review_rows import (
    node_table_rows,
    transition_branch_table_rows,
    transition_count_table_rows,
    uncertainty_table_rows,
)


def build_ancestral_report_package(
    *,
    tree_path: Path,
    traits_path: Path,
    trait: str,
    reconstruction_kind: str,
    out_dir: Path,
    taxon_column: str | None = None,
    model: str = "brownian",
    alpha: float = 1.0,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    compare_model: str | None = None,
    compare_tree_path: Path | None = None,
    drop_taxa: list[str] | None = None,
    coding_map: dict[str, str] | None = None,
) -> AncestralReportPackageResult:
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "ancestral-report.html"
    methods_summary_path = out_dir / "ancestral-methods-summary.md"
    reviewer_audit_checklist_path = out_dir / "reviewer-audit-checklist.tsv"
    figure_path = out_dir / "ancestral-figure.svg"
    figure_png_path = out_dir / "ancestral-figure.png"
    figure_html_path = out_dir / "ancestral-figure.html"
    summary_table_path = out_dir / "summary.tsv"
    node_table_path = out_dir / "node-table.tsv"
    uncertainty_table_path = out_dir / "uncertainty-table.tsv"
    transition_count_table_path = out_dir / "transition-counts.tsv"
    transition_branch_table_path = out_dir / "transition-branches.tsv"
    exclusion_table_path = out_dir / "exclusions.tsv"
    manifest_path = out_dir / "ancestral-report.manifest.json"

    if reconstruction_kind == "continuous":
        reconstruction = reconstruct_continuous_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
            alpha=alpha,
        )
        summary = summarize_continuous_ancestral_report(reconstruction)
        exclusions = continuous_ancestral_exclusions(reconstruction)
        transition_report = None
        transition_branch_rows = summarize_continuous_change_branches(reconstruction)
        transition_count_rows = summarize_continuous_change_counts(
            transition_branch_rows
        )
        write_continuous_ancestral_summary_table(summary_table_path, reconstruction)
        write_continuous_ancestral_uncertainty_table(
            uncertainty_table_path,
            reconstruction,
        )
        write_continuous_ancestral_exclusion_table(
            exclusion_table_path,
            reconstruction,
        )
        write_continuous_change_count_table(
            transition_count_table_path,
            transition_count_rows,
        )
        write_continuous_change_branch_table(
            transition_branch_table_path,
            transition_branch_rows,
        )
    else:
        reconstruction = reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
        )
        summary = summarize_discrete_ancestral_report(reconstruction)
        exclusions = discrete_ancestral_exclusions(reconstruction)
        transition_report = summarize_ancestral_transitions(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
        )
        transition_branch_rows = []
        transition_count_rows = []
        write_discrete_ancestral_summary_table(summary_table_path, reconstruction)
        write_discrete_ancestral_probability_table(
            uncertainty_table_path,
            reconstruction,
        )
        write_discrete_ancestral_exclusion_table(
            exclusion_table_path,
            reconstruction,
        )
        write_ancestral_transition_count_table(
            transition_count_table_path,
            transition_report,
        )
        write_ancestral_transition_branch_table(
            transition_branch_table_path,
            transition_report,
        )

    methods_summary = write_ancestral_methods_summary_text(
        methods_summary_path,
        reconstruction_kind=reconstruction_kind,
        reconstruction=reconstruction,
    )
    write_ancestral_state_table(node_table_path, reconstruction)
    report_result = render_ancestral_state_report(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        reconstruction_kind=reconstruction_kind,
        out_path=report_path,
        taxon_column=taxon_column,
        model=model,
        alpha=alpha,
        compare_model=compare_model,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
        compare_tree_path=compare_tree_path,
        drop_taxa=drop_taxa,
        coding_map=coding_map,
    )
    figure = render_ancestral_state_visualization(
        tree_path,
        reconstruction,
        out_path=figure_path,
        layout="phylogram",
        discrete_node_style="pies",
        branch_coloring="regime" if reconstruction_kind == "continuous" else "state",
    )
    figure_png = render_ancestral_state_visualization(
        tree_path,
        reconstruction,
        out_path=figure_png_path,
        layout="phylogram",
        discrete_node_style="pies",
        branch_coloring="regime" if reconstruction_kind == "continuous" else "state",
    )
    figure_html = render_ancestral_state_visualization(
        tree_path,
        reconstruction,
        out_path=figure_html_path,
        layout="phylogram",
        discrete_node_style="pies",
        branch_coloring="regime" if reconstruction_kind == "continuous" else "state",
    )
    if reconstruction_kind == "continuous":
        count_payload = transition_count_table_rows(
            reconstruction_kind,
            count_rows=transition_count_rows,
        )
        branch_payload = transition_branch_table_rows(
            reconstruction_kind,
            branch_rows=transition_branch_rows,
        )
    else:
        if transition_report is None:
            raise RuntimeError(
                "discrete ancestral report package requires transition diagnostics"
            )
        count_payload = transition_count_table_rows(
            reconstruction_kind,
            count_rows=transition_report.transition_rows,
        )
        branch_payload = transition_branch_table_rows(
            reconstruction_kind,
            branch_rows=transition_report.branch_rows,
        )
        write_ancestral_transition_exclusion_table(
            exclusion_table_path,
            transition_report,
        )

    machine_manifest = {
        "report_kind": "ancestral_report_package",
        "reconstruction_kind": reconstruction_kind,
        "input_paths": [str(tree_path), str(traits_path)],
        "input_checksums": {
            str(tree_path): _checksum(tree_path),
            str(traits_path): _checksum(traits_path),
        },
        "outputs": {
            "report_path": str(report_path),
            "methods_summary_path": str(methods_summary_path),
            "reviewer_audit_checklist_path": str(reviewer_audit_checklist_path),
            "figure_path": str(figure_path),
            "figure_png_path": str(figure_png_path),
            "figure_html_path": str(figure_html_path),
            "summary_table_path": str(summary_table_path),
            "node_table_path": str(node_table_path),
            "uncertainty_table_path": str(uncertainty_table_path),
            "transition_count_table_path": str(transition_count_table_path),
            "transition_branch_table_path": str(transition_branch_table_path),
            "exclusion_table_path": str(exclusion_table_path),
        },
        "metrics": {
            "analyzed_taxon_count": summary.analyzed_taxon_count,
            "excluded_taxon_count": summary.excluded_taxon_count,
            "warning_count": summary.warning_count,
            "methods_summary_warning_count": methods_summary.warning_count,
            "transition_count_row_count": len(count_payload),
            "transition_branch_row_count": len(branch_payload),
        },
        "machine_report_manifest": report_result.machine_manifest,
    }
    reviewer_audit_checklist = write_reviewer_audit_checklist(
        reviewer_audit_checklist_path,
        machine_manifest,
    ).checklist
    machine_manifest["reviewer_audit_checklist"] = asdict(reviewer_audit_checklist)
    manifest_path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_package_html(
        path=report_path,
        title=f"Bijux Ancestral Reconstruction Report: {trait}",
        figure_svg=figure_path.read_text(encoding="utf-8"),
        methods_summary_text=methods_summary.text,
        reconstruction_kind=reconstruction_kind,
        model=model,
        summary=summary,
        warnings=list(reconstruction.warnings),
        node_table_rows=node_table_rows(reconstruction),
        uncertainty_table_rows=uncertainty_table_rows(reconstruction),
        transition_count_rows=count_payload,
        transition_branch_rows=branch_payload,
        limitations=review_limitations(reconstruction_kind, reconstruction),
        reviewer_audit_checklist=reviewer_audit_checklist,
        manifest=machine_manifest,
    )
    return AncestralReportPackageResult(
        output_dir=out_dir,
        report_path=report_path,
        methods_summary_path=methods_summary_path,
        reviewer_audit_checklist_path=reviewer_audit_checklist_path,
        figure_path=figure_path,
        figure_png_path=figure_png_path,
        figure_html_path=figure_html_path,
        summary_table_path=summary_table_path,
        node_table_path=node_table_path,
        uncertainty_table_path=uncertainty_table_path,
        transition_count_table_path=transition_count_table_path,
        transition_branch_table_path=transition_branch_table_path,
        exclusion_table_path=exclusion_table_path,
        manifest_path=manifest_path,
        reconstruction_kind=reconstruction_kind,
        model=model,
        methods_summary=methods_summary,
        summary=summary,
        reconstruction=reconstruction,
        figure=figure,
        figure_png=figure_png,
        figure_html=figure_html,
        transition_count_rows=transition_count_rows,
        transition_branch_rows=transition_branch_rows,
        transition_report=transition_report,
        exclusions=exclusions,
        reviewer_audit_checklist=reviewer_audit_checklist,
        machine_manifest=machine_manifest,
    )


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()
