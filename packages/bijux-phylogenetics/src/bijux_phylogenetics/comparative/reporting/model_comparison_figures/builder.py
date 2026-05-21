from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path

from bijux_phylogenetics.render.reproducibility import (
    write_figure_reproducibility_manifest,
)

from ...continuous import (
    compare_brownian_and_ou_models,
    fit_brownian_motion_model,
    fit_ornstein_uhlenbeck_model,
)
from .contracts import ComparativeModelFigurePackageResult
from .outputs import (
    write_model_figure_caption,
    write_model_figure_legend_table,
    write_model_figure_table,
)
from .presentation import build_model_figure_review_html
from .summaries import (
    build_model_criteria_rows,
    build_model_figure_audit,
    build_model_figure_caption_draft,
    build_model_figure_legend,
    build_model_fit_rows,
    build_model_likelihood_rows,
    build_model_parameter_rows,
)
from .svg import (
    write_fit_summary_svg,
    write_information_criteria_svg,
    write_likelihood_svg,
    write_parameter_svg,
)


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_ready(payload: object) -> object:
    return json.loads(json.dumps(payload, default=str))


def build_comparative_model_figure_package(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    out_dir: Path,
    taxon_column: str | None = None,
) -> ComparativeModelFigurePackageResult:
    """Build a reviewer-facing figure package for BM versus OU continuous-trait model comparison."""
    out_dir.mkdir(parents=True, exist_ok=True)
    criteria_figure_path = out_dir / "model-comparison-criteria.svg"
    likelihood_figure_path = out_dir / "model-comparison-likelihood.svg"
    parameter_figure_path = out_dir / "model-comparison-parameters.svg"
    fit_figure_path = out_dir / "model-comparison-fit-summary.svg"
    criteria_table_path = out_dir / "model-comparison-criteria.tsv"
    likelihood_table_path = out_dir / "model-comparison-likelihood.tsv"
    parameter_table_path = out_dir / "model-comparison-parameters.tsv"
    fit_table_path = out_dir / "model-comparison-fit-summary.tsv"
    legend_path = out_dir / "figure-legend.tsv"
    caption_path = out_dir / "figure-caption.md"
    review_path = out_dir / "model-comparison-review.html"
    manifest_path = out_dir / "model-comparison-package.manifest.json"
    reproducibility_manifest_path = out_dir / "figure-reproducibility.manifest.json"

    comparison_report = compare_brownian_and_ou_models(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    brownian_report = fit_brownian_motion_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    ou_report = fit_ornstein_uhlenbeck_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    criteria_rows = build_model_criteria_rows(comparison_report)
    likelihood_rows = build_model_likelihood_rows(comparison_report)
    parameter_rows = build_model_parameter_rows(brownian_report, ou_report)
    fit_rows = build_model_fit_rows(
        brownian_report,
        ou_report,
        selected_model=comparison_report.better_model,
    )
    write_model_figure_table(criteria_table_path, criteria_rows)
    write_model_figure_table(likelihood_table_path, likelihood_rows)
    write_model_figure_table(parameter_table_path, parameter_rows)
    write_model_figure_table(fit_table_path, fit_rows)
    write_information_criteria_svg(
        criteria_figure_path,
        criteria_rows,
        selected_model=comparison_report.better_model,
    )
    write_likelihood_svg(
        likelihood_figure_path,
        likelihood_rows,
        selected_model=comparison_report.better_model,
    )
    write_parameter_svg(
        parameter_figure_path,
        parameter_rows,
        selected_model=comparison_report.better_model,
    )
    write_fit_summary_svg(
        fit_figure_path,
        fit_rows,
        selected_model=comparison_report.better_model,
    )
    legend_entries = build_model_figure_legend()
    write_model_figure_legend_table(legend_path, legend_entries)
    audit = build_model_figure_audit(
        comparison_report=comparison_report,
        brownian=brownian_report,
        ou=ou_report,
        criteria_rows=criteria_rows,
        likelihood_rows=likelihood_rows,
        parameter_rows=parameter_rows,
        fit_rows=fit_rows,
        legend_entries=legend_entries,
    )
    caption_draft = build_model_figure_caption_draft(
        comparison_report=comparison_report,
        audit=audit,
        fit_rows=fit_rows,
    )
    write_model_figure_caption(caption_path, caption_draft)
    review_path.write_text(
        build_model_figure_review_html(
            criteria_figure_path=criteria_figure_path,
            likelihood_figure_path=likelihood_figure_path,
            parameter_figure_path=parameter_figure_path,
            fit_figure_path=fit_figure_path,
            criteria_table_path=criteria_table_path,
            likelihood_table_path=likelihood_table_path,
            parameter_table_path=parameter_table_path,
            fit_table_path=fit_table_path,
            legend_path=legend_path,
            caption_path=caption_path,
            audit=audit,
        ),
        encoding="utf-8",
    )
    artifact_paths = [
        criteria_figure_path,
        likelihood_figure_path,
        parameter_figure_path,
        fit_figure_path,
        criteria_table_path,
        likelihood_table_path,
        parameter_table_path,
        fit_table_path,
        legend_path,
        caption_path,
        review_path,
    ]
    reproducibility_manifest = write_figure_reproducibility_manifest(
        reproducibility_manifest_path,
        report_kind="comparative_model_figure_package",
        input_files=[
            ("tree", tree_path),
            ("traits", traits_path),
        ],
        generated_figures=[
            ("information_criteria", criteria_figure_path),
            ("likelihood", likelihood_figure_path),
            ("parameter_intervals", parameter_figure_path),
            ("fit_diagnostics", fit_figure_path),
        ],
        generated_tables=[
            ("information_criteria", criteria_table_path),
            ("likelihood", likelihood_table_path),
            ("parameter_intervals", parameter_table_path),
            ("fit_diagnostics", fit_table_path),
        ],
        filters=None,
        model={
            "kind": "comparative_model_comparison",
            "name": comparison_report.better_model,
            "candidate_models": [row.model for row in comparison_report.rows],
            "aicc_delta": audit.aicc_delta,
        },
        settings={
            "trait": trait,
            "taxon_column": taxon_column,
            "taxon_count": comparison_report.taxon_count,
        },
        linked_artifacts=[
            ("legend", legend_path),
            ("caption", caption_path),
            ("review", review_path),
        ],
    )
    machine_manifest = {
        "report_kind": "comparative_model_figure_package",
        "input_paths": [str(tree_path), str(traits_path)],
        "input_checksums": {
            str(tree_path): _checksum(tree_path),
            str(traits_path): _checksum(traits_path),
        },
        "output_paths": [str(path) for path in artifact_paths],
        "output_checksums": {str(path): _checksum(path) for path in artifact_paths},
        "reproducibility_manifest_path": str(reproducibility_manifest_path),
        "reproducibility_manifest_checksum": _checksum(reproducibility_manifest_path),
        "reproducibility_manifest": reproducibility_manifest,
        "settings": {
            "trait": trait,
            "taxon_column": taxon_column,
        },
        "metrics": {
            "taxon_count": comparison_report.taxon_count,
            "selected_model": audit.selected_model,
            "publication_ready": audit.publication_ready,
            "support_distinct": audit.support_distinct,
            "aicc_delta": audit.aicc_delta,
            "finite_aicc_model_count": audit.finite_aicc_model_count,
            "plotted_model_count": audit.plotted_model_count,
            "rendered_parameter_count": audit.rendered_parameter_count,
            "rendered_fit_row_count": audit.rendered_fit_row_count,
            "warning_count": audit.warning_count,
        },
        "comparison_report": _json_ready(asdict(comparison_report)),
        "brownian_report": _json_ready(asdict(brownian_report)),
        "ou_report": _json_ready(asdict(ou_report)),
        "criteria_rows": _json_ready([asdict(row) for row in criteria_rows]),
        "likelihood_rows": _json_ready([asdict(row) for row in likelihood_rows]),
        "parameter_rows": _json_ready([asdict(row) for row in parameter_rows]),
        "fit_rows": _json_ready([asdict(row) for row in fit_rows]),
        "audit": _json_ready(asdict(audit)),
    }
    manifest_path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return ComparativeModelFigurePackageResult(
        output_dir=out_dir,
        criteria_figure_path=criteria_figure_path,
        likelihood_figure_path=likelihood_figure_path,
        parameter_figure_path=parameter_figure_path,
        fit_figure_path=fit_figure_path,
        criteria_table_path=criteria_table_path,
        likelihood_table_path=likelihood_table_path,
        parameter_table_path=parameter_table_path,
        fit_table_path=fit_table_path,
        legend_path=legend_path,
        caption_path=caption_path,
        review_path=review_path,
        manifest_path=manifest_path,
        reproducibility_manifest_path=reproducibility_manifest_path,
        comparison_report=comparison_report,
        brownian_report=brownian_report,
        ou_report=ou_report,
        criteria_rows=criteria_rows,
        likelihood_rows=likelihood_rows,
        parameter_rows=parameter_rows,
        fit_rows=fit_rows,
        legend_entries=legend_entries,
        caption_draft=caption_draft,
        audit=audit,
        machine_manifest=machine_manifest,
    )
