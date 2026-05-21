from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
from html import escape
import json
from pathlib import Path

from .artifact_outputs import write_continuous_change_branch_table
from .artifact_outputs import write_continuous_change_count_table
from .continuous_changes import summarize_continuous_change_branches
from .continuous_changes import summarize_continuous_change_counts
from .contracts import AncestralContinuousChangeBranchRow
from .contracts import AncestralContinuousChangeCountRow
from .contracts import AncestralReportPackageResult

from bijux_phylogenetics.ancestral.continuous import (
    ContinuousAncestralExclusion,
    ContinuousAncestralReport,
    ContinuousAncestralSummary,
    continuous_ancestral_exclusions,
    reconstruct_continuous_ancestral_states,
    summarize_continuous_ancestral_report,
    write_continuous_ancestral_exclusion_table,
    write_continuous_ancestral_summary_table,
    write_continuous_ancestral_uncertainty_table,
)
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralExclusion,
    DiscreteAncestralReport,
    DiscreteAncestralSummary,
    discrete_ancestral_exclusions,
    reconstruct_discrete_ancestral_states,
    summarize_discrete_ancestral_report,
    write_discrete_ancestral_exclusion_table,
    write_discrete_ancestral_probability_table,
    write_discrete_ancestral_summary_table,
)
from bijux_phylogenetics.ancestral.discrete.review import (
    AncestralTransitionBranchRow,
    AncestralTransitionCountRow,
    AncestralTransitionReport,
    summarize_ancestral_transitions,
    write_ancestral_transition_branch_table,
    write_ancestral_transition_count_table,
    write_ancestral_transition_exclusion_table,
)
from bijux_phylogenetics.ancestral.presentation.methods_text import (
    AncestralMethodsSummaryTextResult,
    write_ancestral_methods_summary_text,
)
from bijux_phylogenetics.ancestral.presentation.report_rendering import (
    render_ancestral_state_report,
    write_ancestral_state_table,
)
from bijux_phylogenetics.ancestral.presentation.visualization import (
    AncestralVisualizationResult,
    render_ancestral_state_visualization,
)
from bijux_phylogenetics.reports.review import (
    ReviewerAuditChecklist,
    write_reviewer_audit_checklist,
)


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_script(payload: dict[str, object]) -> str:
    serialized = json.dumps(payload, default=str, indent=2, sort_keys=True).replace(
        "</", "<\\/"
    )
    return (
        '<script id="bijux-ancestral-report-package-manifest" type="application/json">'
        f"{serialized}</script>"
    )


def _table(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{escape(header)}</th>" for header in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{escape(cell)}</td>" for cell in row) + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _list(items: list[str]) -> str:
    if not items:
        return "<p>none</p>"
    return "<ul>" + "".join(f"<li>{escape(item)}</li>" for item in items) + "</ul>"


def _reviewer_audit_table(checklist: ReviewerAuditChecklist) -> str:
    rows = []
    for item in checklist.items:
        rows.append(
            "<tr>"
            f"<td>{escape(item.section)}</td>"
            f"<td>{escape(item.status)}</td>"
            f"<td>{escape(item.summary)}</td>"
            f"<td>{escape('; '.join(item.evidence))}</td>"
            "</tr>"
        )
    return (
        "<table><thead><tr>"
        "<th>section</th><th>status</th><th>summary</th><th>evidence</th>"
        "</tr></thead><tbody>" + "".join(rows) + "</tbody></table>"
    )


def _write_package_html(
    *,
    path: Path,
    title: str,
    figure_svg: str,
    methods_summary_text: str,
    reconstruction_kind: str,
    model: str,
    summary: ContinuousAncestralSummary | DiscreteAncestralSummary,
    warnings: list[str],
    node_table_rows: list[list[str]],
    uncertainty_table_rows: list[list[str]],
    transition_count_rows: list[list[str]],
    transition_branch_rows: list[list[str]],
    limitations: list[str],
    reviewer_audit_checklist: ReviewerAuditChecklist,
    manifest: dict[str, object],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(title)}</title>
  <style>
    :root {{
      --ink: #172033;
      --muted: #52607a;
      --bg: #f7f9fc;
      --panel: #ffffff;
      --rule: #d7e0ea;
      --accent: #1d4ed8;
      --accent-soft: #dbeafe;
    }}
    body {{
      margin: 0;
      background:
        radial-gradient(circle at top right, #e4eefc 0, transparent 28rem),
        linear-gradient(180deg, #eef4fb 0%, var(--bg) 100%);
      color: var(--ink);
      font: 16px/1.5 "Iowan Old Style", "Palatino Linotype", serif;
    }}
    main {{
      max-width: 1180px;
      margin: 0 auto;
      padding: 2rem;
    }}
    .shell {{
      background: rgba(255, 255, 255, 0.94);
      border: 1px solid var(--rule);
      border-radius: 24px;
      box-shadow: 0 28px 80px rgba(15, 23, 42, 0.08);
      padding: 2rem;
    }}
    h1, h2 {{
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      letter-spacing: 0.02em;
    }}
    h1 {{
      margin: 0 0 0.4rem;
      color: var(--accent);
    }}
    .lead {{
      margin: 0 0 1.5rem;
      color: var(--muted);
    }}
    .cards {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
      gap: 0.75rem;
      margin: 1.5rem 0 2rem;
    }}
    .card {{
      background: linear-gradient(180deg, #ffffff 0%, #f8fbff 100%);
      border: 1px solid var(--rule);
      border-radius: 18px;
      padding: 1rem;
    }}
    .label {{
      color: var(--muted);
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      font-size: 0.78rem;
      text-transform: uppercase;
      letter-spacing: 0.08em;
    }}
    .value {{
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      font-size: 1.6rem;
      margin-top: 0.2rem;
    }}
    section + section {{
      margin-top: 2rem;
      padding-top: 2rem;
      border-top: 1px solid var(--rule);
    }}
    .figure-frame {{
      overflow-x: auto;
      background: #fbfdff;
      border: 1px solid var(--rule);
      border-radius: 18px;
      padding: 1rem;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      font-size: 0.92rem;
    }}
    th, td {{
      text-align: left;
      padding: 0.65rem 0.55rem;
      border-bottom: 1px solid var(--rule);
      vertical-align: top;
    }}
    th {{
      color: var(--muted);
      font-size: 0.76rem;
      text-transform: uppercase;
      letter-spacing: 0.08em;
    }}
    .note {{
      background: var(--accent-soft);
      border-left: 4px solid var(--accent);
      border-radius: 12px;
      padding: 0.9rem 1rem;
    }}
  </style>
</head>
<body>
  <main>
    <div class="shell">
      <h1>{escape(title)}</h1>
      <p class="lead">Reviewer-facing ancestral reconstruction bundle with node estimates, uncertainty, branchwise change evidence, and visualization in one package.</p>
      {_json_script(manifest)}
      <div class="cards">
        <div class="card"><div class="label">Kind</div><div class="value">{
        escape(reconstruction_kind)
    }</div></div>
        <div class="card"><div class="label">Model</div><div class="value">{
        escape(model)
    }</div></div>
        <div class="card"><div class="label">Analyzed Taxa</div><div class="value">{
        summary.analyzed_taxon_count
    }</div></div>
        <div class="card"><div class="label">Warnings</div><div class="value">{
        summary.warning_count
    }</div></div>
      </div>
      <section>
        <h2>Reviewer Summary</h2>
        {_list(warnings[:6])}
        <div class="note">
          Continuous packages preserve branch-change direction counts. Discrete packages preserve inferred state-transition counts.
        </div>
      </section>
      <section>
        <h2>Methods Summary</h2>
        <pre>{escape(methods_summary_text)}</pre>
      </section>
      <section>
        <h2>Reviewer Audit Checklist</h2>
        {_reviewer_audit_table(reviewer_audit_checklist)}
      </section>
      <section>
        <h2>Tree Visualization</h2>
        <div class="figure-frame">{figure_svg}</div>
      </section>
      <section>
        <h2>Node Table</h2>
        {
        _table(
            [
                "node",
                "node name",
                "tip",
                "descendant taxa",
                "value or state",
                "uncertainty",
            ],
            node_table_rows,
        )
    }
      </section>
      <section>
        <h2>Uncertainty Review</h2>
        {
        _table(
            ["node", "descendant taxa", "uncertainty evidence", "interpretation"],
            uncertainty_table_rows,
        )
    }
      </section>
      <section>
        <h2>Transition Review</h2>
        {
        _table(
            ["label", "count", "fraction or certainty", "detail"],
            transition_count_rows,
        )
    }
        {
        _table(
            ["parent", "child", "descendant taxa", "branch length", "change", "detail"],
            transition_branch_rows,
        )
    }
      </section>
      <section>
        <h2>Limitations</h2>
        {_list(limitations)}
      </section>
    </div>
  </main>
</body>
</html>
"""
    path.write_text(html, encoding="utf-8")
    return path


def _node_table_rows(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> list[list[str]]:
    if isinstance(report, ContinuousAncestralReport):
        return [
            [
                estimate.node,
                "" if estimate.node_name is None else estimate.node_name,
                "true" if estimate.is_tip else "false",
                ", ".join(estimate.descendant_taxa),
                format(estimate.estimate, ".6g"),
                format(estimate.standard_error, ".6g"),
            ]
            for estimate in report.estimates
            if not estimate.is_tip
        ]
    return [
        [
            estimate.node,
            "" if estimate.node_name is None else estimate.node_name,
            "true" if estimate.is_tip else "false",
            ", ".join(estimate.descendant_taxa),
            estimate.most_likely_state,
            format(estimate.confidence, ".6g"),
        ]
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def _uncertainty_rows(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> list[list[str]]:
    if isinstance(report, ContinuousAncestralReport):
        return [
            [
                estimate.node,
                ", ".join(estimate.descendant_taxa),
                (
                    f"se={format(estimate.standard_error, '.6g')}; "
                    f"95%=[{format(estimate.lower_95_interval, '.6g')}, {format(estimate.upper_95_interval, '.6g')}]"
                ),
                estimate.interpretation,
            ]
            for estimate in report.estimates
            if not estimate.is_tip
        ]
    return [
        [
            estimate.node,
            ", ".join(estimate.descendant_taxa),
            (
                f"state={estimate.most_likely_state}; "
                f"confidence={format(estimate.confidence, '.6g')}; "
                f"probabilities={json.dumps(estimate.state_probabilities, sort_keys=True)}"
            ),
            estimate.interpretation,
        ]
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def _transition_count_rows(
    reconstruction_kind: str,
    *,
    count_rows: list[AncestralContinuousChangeCountRow]
    | list[AncestralTransitionCountRow],
) -> list[list[str]]:
    if reconstruction_kind == "continuous":
        continuous_rows = count_rows
        if not all(
            isinstance(row, AncestralContinuousChangeCountRow)
            for row in continuous_rows
        ):
            raise RuntimeError(
                "continuous ancestral report package received non-continuous change counts"
            )
        return [
            [
                row.direction,
                str(row.branch_count),
                format(row.branch_fraction, ".6g"),
                (
                    f"mean_delta={format(row.mean_delta, '.6g')}; "
                    f"range=[{format(row.minimum_delta, '.6g')}, {format(row.maximum_delta, '.6g')}]"
                ),
            ]
            for row in continuous_rows
        ]
    discrete_rows = count_rows
    if not all(isinstance(row, AncestralTransitionCountRow) for row in discrete_rows):
        raise RuntimeError(
            "discrete ancestral report package received non-discrete transition counts"
        )
    return [
        [
            row.transition,
            str(row.total_change_count),
            (
                "certain"
                if row.uncertain_change_count == 0
                else "mixed"
                if row.certain_change_count > 0
                else "uncertain"
            ),
            (
                f"certain={row.certain_change_count}; "
                f"uncertain={row.uncertain_change_count}"
            ),
        ]
        for row in discrete_rows
    ]


def _transition_branch_rows(
    reconstruction_kind: str,
    *,
    branch_rows: list[AncestralContinuousChangeBranchRow]
    | list[AncestralTransitionBranchRow],
) -> list[list[str]]:
    if reconstruction_kind == "continuous":
        continuous_rows = branch_rows
        if not all(
            isinstance(row, AncestralContinuousChangeBranchRow)
            for row in continuous_rows
        ):
            raise RuntimeError(
                "continuous ancestral report package received non-continuous branch rows"
            )
        return [
            [
                row.parent_node,
                row.child_node,
                ", ".join(row.child_descendant_taxa),
                "" if row.branch_length is None else format(row.branch_length, ".6g"),
                row.direction,
                (
                    f"delta={format(row.delta, '.6g')}; "
                    f"parent={format(row.parent_estimate, '.6g')}; "
                    f"child={format(row.child_estimate, '.6g')}"
                ),
            ]
            for row in continuous_rows
        ]
    discrete_rows = branch_rows
    if not all(isinstance(row, AncestralTransitionBranchRow) for row in discrete_rows):
        raise RuntimeError(
            "discrete ancestral report package received non-discrete branch rows"
        )
    return [
        [
            row.parent_node,
            row.child_node,
            ", ".join(row.child_descendant_taxa),
            "" if row.branch_length is None else format(row.branch_length, ".6g"),
            row.transition,
            row.certainty_class,
        ]
        for row in discrete_rows
        if row.changed
    ]


def _limitations(
    reconstruction_kind: str,
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> list[str]:
    limitations = list(report.warnings)
    if reconstruction_kind == "continuous":
        limitations.append(
            "continuous branch-change counts summarize direction of reconstructed value shifts, not discrete state-transition events"
        )
    else:
        limitations.append(
            "discrete transition counts are reconstructed branchwise review evidence, not stochastic mapping event histories"
        )
    return limitations


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
        count_payload = _transition_count_rows(
            reconstruction_kind,
            count_rows=transition_count_rows,
        )
        branch_payload = _transition_branch_rows(
            reconstruction_kind,
            branch_rows=transition_branch_rows,
        )
    else:
        if transition_report is None:
            raise RuntimeError(
                "discrete ancestral report package requires transition diagnostics"
            )
        count_payload = _transition_count_rows(
            reconstruction_kind,
            count_rows=transition_report.transition_rows,
        )
        branch_payload = _transition_branch_rows(
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
    _write_package_html(
        path=report_path,
        title=f"Bijux Ancestral Reconstruction Report: {trait}",
        figure_svg=figure_path.read_text(encoding="utf-8"),
        methods_summary_text=methods_summary.text,
        reconstruction_kind=reconstruction_kind,
        model=model,
        summary=summary,
        warnings=list(reconstruction.warnings),
        node_table_rows=_node_table_rows(reconstruction),
        uncertainty_table_rows=_uncertainty_rows(reconstruction),
        transition_count_rows=count_payload,
        transition_branch_rows=branch_payload,
        limitations=_limitations(reconstruction_kind, reconstruction),
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
