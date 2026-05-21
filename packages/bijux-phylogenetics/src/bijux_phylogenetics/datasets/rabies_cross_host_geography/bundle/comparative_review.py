from __future__ import annotations

from html import escape
import json
from pathlib import Path

from bijux_phylogenetics.comparative.reporting.analysis_package import (
    ComparativeAnalysisSummaryRow,
    ComparativeCoefficientTableRow,
    ComparativeInterpretationRow,
    ComparativeResidualTableRow,
    ComparativeSignalTableRow,
)

from ..models import RabiesComparativeBranchRepair
from ..shared import _checksum, _format_number, _html_list, _table


def _write_comparative_report(
    path: Path,
    *,
    summary_row: ComparativeAnalysisSummaryRow,
    coefficient_rows: list[ComparativeCoefficientTableRow],
    residual_rows: list[ComparativeResidualTableRow],
    signal_row: ComparativeSignalTableRow,
    interpretation_rows: list[ComparativeInterpretationRow],
    branch_repairs: list[RabiesComparativeBranchRepair],
) -> Path:
    key_claim = next(
        (
            row.claim
            for row in interpretation_rows
            if row.topic == "coefficient" and "nominally supported" in row.claim
        ),
        "no coefficient reached nominal support",
    )
    html = "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Rabies Comparative Report</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: linear-gradient(180deg, #f4f1ea 0%, #f3f7ef 100%); color: #173024; }",
            "    main { max-width: 1040px; margin: 0 auto; padding: 24px; }",
            "    h1, h2 { margin: 0 0 10px; }",
            "    p { line-height: 1.55; }",
            "    .panel { background: rgba(255,255,255,0.88); border: 1px solid rgba(23,48,36,0.12); border-radius: 18px; padding: 18px; margin-top: 18px; box-shadow: 0 14px 36px rgba(23,48,36,0.08); }",
            "    .cards { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 14px; margin-top: 18px; }",
            "    .card { background: rgba(255,255,255,0.88); border: 1px solid rgba(23,48,36,0.12); border-radius: 18px; padding: 16px; }",
            "    .label { color: #5f7469; font-size: 13px; text-transform: uppercase; letter-spacing: 0.04em; }",
            "    .value { display: block; font-size: 22px; margin-top: 6px; }",
            "    table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 14px; }",
            "    th, td { border-bottom: 1px solid rgba(23,48,36,0.10); padding: 8px 10px; text-align: left; vertical-align: top; }",
            "    th { color: #365443; }",
            "    ul { margin: 8px 0 0 18px; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Rabies Comparative Report</h1>",
            "  <p>This comparative section asks whether the host-associated lineages in the rabies demonstration tree are associated with a consistent eastward geographic placement when geography is summarized as regional longitude. The result is interpretive evidence, not causal proof, and it inherits the small-panel limits of this dataset.</p>",
            '  <section class="cards">',
            f'    <div class="card"><span class="label">formula</span><span class="value">{escape(summary_row.formula)}</span></div>',
            f'    <div class="card"><span class="label">analysis taxa</span><span class="value">{summary_row.analysis_taxa}</span></div>',
            f'    <div class="card"><span class="label">selected trait model</span><span class="value">{escape(summary_row.selected_model)}</span></div>',
            f'    <div class="card"><span class="label">pgls r-squared</span><span class="value">{_format_number(summary_row.pgls_r_squared)}</span></div>',
            "  </section>",
            '  <section class="panel">',
            "    <h2>Question and Answer</h2>",
            _html_list(
                [
                    "Question: does host association coincide with a consistent longitudinal shift in this rabies panel?",
                    f"Answer: {key_claim}.",
                    (
                        f"Phylogenetic signal remains strong for the response trait "
                        f"(Blomberg's K {_format_number(signal_row.blombergs_k)}, "
                        f"Pagel's lambda {_format_number(signal_row.pagels_lambda)})."
                    ),
                    (
                        "Interpret the coefficient evidence cautiously because the "
                        "residual diagnostics retain review warnings and the sample "
                        "is intentionally compact."
                    ),
                ]
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Coefficient Summary</h2>",
            _table(
                headers=[
                    "term",
                    "estimate",
                    "standard_error",
                    "p_value",
                    "significant",
                ],
                rows=[
                    [
                        row.term,
                        _format_number(row.estimate),
                        _format_number(row.standard_error),
                        _format_number(row.p_value),
                        "true" if row.significant else "false",
                    ]
                    for row in coefficient_rows
                ],
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Residual Diagnostics</h2>",
            _table(
                headers=[
                    "analysis",
                    "residual_variance",
                    "max_abs_standardized_residual",
                    "phylogenetic_residual_lambda",
                    "warnings",
                ],
                rows=[
                    [
                        row.analysis,
                        _format_number(row.residual_variance),
                        _format_number(row.max_abs_standardized_residual),
                        _format_number(row.phylogenetic_residual_lambda),
                        "; ".join(row.warnings),
                    ]
                    for row in residual_rows
                ],
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Comparative Tree Adjustments</h2>",
            _html_list(
                [
                    "The comparative fit uses the rooted demonstration tree after flooring any nonpositive nonroot branch lengths to a tiny positive value.",
                    f"Adjusted branch count: {len(branch_repairs)}",
                ]
            ),
            _table(
                headers=[
                    "node_label",
                    "original_branch_length",
                    "repaired_branch_length",
                    "reason",
                ],
                rows=[
                    [
                        row.node_label,
                        _format_number(row.original_branch_length),
                        _format_number(row.repaired_branch_length),
                        row.reason,
                    ]
                    for row in branch_repairs
                ]
                or [["", "", "", "no branch-length repair was required"]],
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Interpretation Ledger</h2>",
            _table(
                headers=["topic", "claim", "evidence", "caution"],
                rows=[
                    [row.topic, row.claim, row.evidence, row.caution]
                    for row in interpretation_rows
                ],
            ),
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )
    path.write_text(html + "\n", encoding="utf-8")
    return path


def _write_comparative_manifest(
    path: Path,
    *,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    branch_repairs: list[RabiesComparativeBranchRepair],
    output_paths: dict[str, Path],
) -> Path:
    payload = {
        "report_kind": "rabies_cross_host_geography_comparative_bundle",
        "metrics": {
            "response": comparative_summary_row.response,
            "formula": comparative_summary_row.formula,
            "analysis_taxa": comparative_summary_row.analysis_taxa,
            "selected_model": comparative_summary_row.selected_model,
            "pgls_lambda": comparative_summary_row.pgls_lambda,
            "pgls_r_squared": comparative_summary_row.pgls_r_squared,
            "branch_repair_count": len(branch_repairs),
        },
        "output_checksums": {
            key: _checksum(value) for key, value in output_paths.items()
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path
