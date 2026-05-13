from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path

from bijux_phylogenetics.bayesian.burnin import DEFAULT_BURNIN_FRACTIONS
from bijux_phylogenetics.bayesian.beast import (
    assess_beast_burnin_sensitivity,
    assess_beast_chain_mixing,
    assess_calibration_dominance,
    assess_time_tree_readiness,
    summarize_beast_analysis_xml,
    validate_beast_posterior_log,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.render.svg import render_tree_svg
from bijux_phylogenetics.tree_set import (
    cluster_trees_by_topology,
    compute_clade_frequency_table,
    compute_consensus_tree,
    detect_posterior_topology_multimodality,
    detect_unstable_taxa,
    summarize_clade_credibility_conflicts,
    summarize_uncertainty_aware_conclusions,
    write_topology_cluster_table,
    write_uncertainty_conclusion_table,
)


@dataclass(slots=True)
class PosteriorUncertaintyFigurePackageResult:
    output_dir: Path
    consensus_tree_path: Path
    consensus_figure_path: Path
    clade_frequency_plot_path: Path
    unstable_taxa_table_path: Path
    topology_clusters_table_path: Path
    uncertainty_conclusions_table_path: Path
    conclusion_summary_path: Path
    manifest_path: Path


@dataclass(slots=True)
class SupplementaryBayesianDiagnosticsTableResult:
    output_path: Path
    row_count: int
    chain_count: int
    warning_count: int


@dataclass(slots=True)
class BayesianMethodsSummaryTextResult:
    output_path: Path
    title: str
    warning_count: int
    text: str


@dataclass(slots=True)
class BayesianLimitationsTextResult:
    output_path: Path
    title: str
    warning_count: int
    text: str


def build_posterior_uncertainty_figure_package(
    tree_set_path: Path,
    *,
    out_dir: Path,
    layout: str = "phylogram",
    frequency_plot_limit: int = 12,
) -> PosteriorUncertaintyFigurePackageResult:
    """Build a publication-oriented posterior uncertainty package from one tree set."""
    out_dir.mkdir(parents=True, exist_ok=True)
    consensus_tree, consensus = compute_consensus_tree(tree_set_path)
    clade_frequencies = compute_clade_frequency_table(tree_set_path)
    unstable_taxa = detect_unstable_taxa(tree_set_path)
    clusters = cluster_trees_by_topology(tree_set_path)
    multimodality = detect_posterior_topology_multimodality(tree_set_path)
    conflicts = summarize_clade_credibility_conflicts(tree_set_path)
    conclusions = summarize_uncertainty_aware_conclusions(tree_set_path)

    consensus_tree_path = out_dir / "consensus-tree.nwk"
    consensus_figure_path = out_dir / "consensus-tree.svg"
    clade_frequency_plot_path = out_dir / "clade-frequency-plot.svg"
    unstable_taxa_table_path = out_dir / "unstable-taxa.tsv"
    topology_clusters_table_path = out_dir / "topology-clusters.tsv"
    uncertainty_conclusions_table_path = out_dir / "uncertainty-conclusions.tsv"
    conclusion_summary_path = out_dir / "uncertainty-summary.md"
    manifest_path = out_dir / "uncertainty-package-manifest.json"

    write_newick(consensus_tree_path, consensus_tree)
    render_tree_svg(
        consensus_tree_path,
        out_path=consensus_figure_path,
        layout=layout,
        show_support_values=True,
    )
    _write_clade_frequency_plot(
        clade_frequency_plot_path,
        clade_frequencies=clade_frequencies,
        limit=frequency_plot_limit,
    )
    write_taxon_rows(
        unstable_taxa_table_path,
        columns=[
            "taxon",
            "unique_placements",
            "dominant_frequency",
            "instability_score",
            "placement_signatures",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "unique_placements": str(row.unique_placements),
                "dominant_frequency": format(row.dominant_frequency, ".15g"),
                "instability_score": format(row.instability_score, ".15g"),
                "placement_signatures": "; ".join(
                    f"{placement.signature} ({format(placement.frequency, '.15g')})"
                    for placement in row.placements
                ),
            }
            for row in unstable_taxa.taxa
        ],
    )
    write_topology_cluster_table(topology_clusters_table_path, clusters)
    write_uncertainty_conclusion_table(uncertainty_conclusions_table_path, conclusions)
    conclusion_summary_path.write_text(
        _uncertainty_summary_markdown(
            consensus_newick=consensus.consensus_newick,
            multimodality=multimodality,
            conflicts=conflicts,
            conclusions=conclusions,
        ),
        encoding="utf-8",
    )
    manifest_path.write_text(
        json.dumps(
            {
                "tree_set_path": str(tree_set_path),
                "layout": layout,
                "artifacts": {
                    "consensus_tree": str(consensus_tree_path),
                    "consensus_figure": str(consensus_figure_path),
                    "clade_frequency_plot": str(clade_frequency_plot_path),
                    "unstable_taxa_table": str(unstable_taxa_table_path),
                    "topology_clusters_table": str(topology_clusters_table_path),
                    "uncertainty_conclusions_table": str(
                        uncertainty_conclusions_table_path
                    ),
                    "uncertainty_summary": str(conclusion_summary_path),
                },
                "consensus": asdict(consensus),
                "multimodality": asdict(multimodality),
                "clade_conflicts": asdict(conflicts),
                "conclusions": asdict(conclusions),
            },
            indent=2,
            sort_keys=True,
            default=str,
        )
        + "\n",
        encoding="utf-8",
    )
    return PosteriorUncertaintyFigurePackageResult(
        output_dir=out_dir,
        consensus_tree_path=consensus_tree_path,
        consensus_figure_path=consensus_figure_path,
        clade_frequency_plot_path=clade_frequency_plot_path,
        unstable_taxa_table_path=unstable_taxa_table_path,
        topology_clusters_table_path=topology_clusters_table_path,
        uncertainty_conclusions_table_path=uncertainty_conclusions_table_path,
        conclusion_summary_path=conclusion_summary_path,
        manifest_path=manifest_path,
    )


def write_supplementary_bayesian_diagnostics_table(
    path: Path,
    *,
    posterior_tree_path: Path,
    primary_log_path: Path,
    additional_log_paths: list[Path] | None = None,
    burnin_fractions: tuple[float, ...] = DEFAULT_BURNIN_FRACTIONS,
    required_columns: tuple[str, ...] = ("posterior", "likelihood"),
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
    cross_chain_mean_shift_threshold: float = 0.75,
) -> SupplementaryBayesianDiagnosticsTableResult:
    """Write a reviewer-facing supplementary diagnostics table for Bayesian posterior analysis."""
    validation = validate_beast_posterior_log(
        primary_log_path, required_columns=required_columns
    )
    burnin = assess_beast_burnin_sensitivity(
        posterior_tree_path,
        log_path=primary_log_path,
        burnin_fractions=burnin_fractions,
    )
    log_paths = [primary_log_path, *(additional_log_paths or [])]
    mixing = assess_beast_chain_mixing(
        log_paths,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
        cross_chain_mean_shift_threshold=cross_chain_mean_shift_threshold,
    )
    rows: list[dict[str, str]] = []
    rows.append(
        {
            "row_kind": "log-validation",
            "chain": "primary",
            "parameter": "",
            "burnin_fraction": "",
            "effective_sample_size": "",
            "mean": "",
            "posterior_mean": "",
            "likelihood_mean": "",
            "tree_height_mean": "",
            "warning_code": ",".join(issue.code for issue in validation.issues),
            "details": f"valid={str(validation.valid).lower()}; missing_columns={','.join(validation.missing_columns)}",
        }
    )
    for burnin_slice in burnin.slices:
        rows.append(
            {
                "row_kind": "burnin-summary",
                "chain": "primary",
                "parameter": "",
                "burnin_fraction": format(burnin_slice.burnin_fraction, ".15g"),
                "effective_sample_size": "",
                "mean": "",
                "posterior_mean": ""
                if burnin_slice.posterior_mean is None
                else format(burnin_slice.posterior_mean, ".15g"),
                "likelihood_mean": ""
                if burnin_slice.likelihood_mean is None
                else format(burnin_slice.likelihood_mean, ".15g"),
                "tree_height_mean": ""
                if burnin_slice.tree_height_mean is None
                else format(burnin_slice.tree_height_mean, ".15g"),
                "warning_code": "",
                "details": (
                    f"selected_tree_index={burnin_slice.selected_tree_index}; "
                    f"rooted_topology_count={burnin_slice.rooted_topology_count}; "
                    f"kept_tree_count={burnin_slice.kept_tree_count}"
                ),
            }
        )
    for chain_index, summary in enumerate(mixing.chain_summaries, start=1):
        warning_codes = ",".join(str(warning["code"]) for warning in summary.warnings)
        for parameter_summary in summary.parameter_summaries:
            rows.append(
                {
                    "row_kind": "chain-parameter",
                    "chain": f"chain_{chain_index}",
                    "parameter": str(parameter_summary["parameter"]),
                    "burnin_fraction": "",
                    "effective_sample_size": format(
                        float(parameter_summary["effective_sample_size"]), ".15g"
                    ),
                    "mean": format(float(parameter_summary["mean"]), ".15g"),
                    "posterior_mean": "",
                    "likelihood_mean": "",
                    "tree_height_mean": "",
                    "warning_code": warning_codes,
                    "details": (
                        f"sample_count={parameter_summary['sample_count']}; "
                        f"standardized_mean_shift={format(float(parameter_summary['standardized_mean_shift']), '.15g')}"
                    ),
                }
            )
    write_taxon_rows(
        path,
        columns=[
            "row_kind",
            "chain",
            "parameter",
            "burnin_fraction",
            "effective_sample_size",
            "mean",
            "posterior_mean",
            "likelihood_mean",
            "tree_height_mean",
            "warning_code",
            "details",
        ],
        rows=rows,
    )
    return SupplementaryBayesianDiagnosticsTableResult(
        output_path=path,
        row_count=len(rows),
        chain_count=len(log_paths),
        warning_count=len(validation.issues)
        + len(burnin.warnings)
        + len(mixing.issues),
    )


def write_bayesian_methods_summary_text(
    path: Path,
    *,
    posterior_tree_path: Path,
    primary_log_path: Path,
    additional_log_paths: list[Path] | None = None,
    analysis_xml_path: Path | None = None,
    tree_prior: str = "unspecified",
    clock_model: str = "unspecified",
    burnin_fractions: tuple[float, ...] = DEFAULT_BURNIN_FRACTIONS,
    required_columns: tuple[str, ...] = ("posterior", "likelihood"),
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
    cross_chain_mean_shift_threshold: float = 0.75,
    calibration_path: Path | None = None,
    tip_dates_path: Path | None = None,
) -> BayesianMethodsSummaryTextResult:
    """Write reviewer-facing methods text for one Bayesian posterior analysis."""
    validation = validate_beast_posterior_log(
        primary_log_path, required_columns=required_columns
    )
    burnin = assess_beast_burnin_sensitivity(
        posterior_tree_path,
        log_path=primary_log_path,
        burnin_fractions=burnin_fractions,
    )
    mixing = assess_beast_chain_mixing(
        [primary_log_path, *(additional_log_paths or [])],
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
        cross_chain_mean_shift_threshold=cross_chain_mean_shift_threshold,
    )
    minimum_ess = min(
        float(parameter_summary["effective_sample_size"])
        for chain in mixing.chain_summaries
        for parameter_summary in chain.parameter_summaries
    )
    retained_range = (
        min(summary.kept_tree_count for summary in burnin.slices),
        max(summary.kept_tree_count for summary in burnin.slices),
    )
    resolved_analysis_xml_path = _resolve_beast_analysis_xml_path(
        analysis_xml_path=analysis_xml_path,
        primary_log_path=primary_log_path,
    )
    analysis_summary = (
        None
        if resolved_analysis_xml_path is None
        else summarize_beast_analysis_xml(resolved_analysis_xml_path)
    )
    resolved_clock_model = clock_model
    resolved_tree_prior = tree_prior
    chain_settings_text = "Chain settings were not available because no analysis XML was supplied."
    analysis_xml_text = "No analysis XML was supplied, so XML-derived BEAST assumptions were unavailable."
    if analysis_summary is not None:
        if analysis_summary.clock_model is not None:
            resolved_clock_model = analysis_summary.clock_model
        if analysis_summary.tree_prior is not None:
            resolved_tree_prior = analysis_summary.tree_prior
        chain_settings_text = (
            f"The prepared BEAST XML `{resolved_analysis_xml_path.name}` declares chain length "
            f"`{analysis_summary.chain_length}` with posterior log `{analysis_summary.posterior_log_path.name if analysis_summary.posterior_log_path is not None else 'unspecified'}` "
            f"and posterior tree log `{analysis_summary.posterior_tree_path.name if analysis_summary.posterior_tree_path is not None else 'unspecified'}`."
        )
        analysis_xml_text = (
            f"The prepared BEAST XML `{resolved_analysis_xml_path.name}` declares `{analysis_summary.starting_tree_source}` starting trees, "
            f"`{analysis_summary.substitution_model}` substitution, `{resolved_clock_model}` clock, `{resolved_tree_prior}` tree prior, "
            f"`{analysis_summary.calibration_count}` calibration prior(s), and `{analysis_summary.tip_date_count}` dated tip(s)."
        )
    calibration_text = "No calibration table was provided."
    if calibration_path is not None:
        calibration_text = (
            f"Calibration constraints were supplied from `{calibration_path.name}`."
        )
    tip_date_text = "No tip-date table was supplied."
    if tip_dates_path is not None:
        tip_date_text = f"Tip dates were supplied from `{tip_dates_path.name}`."
    warning_count = len(validation.issues) + len(burnin.warnings) + len(mixing.issues)
    text = (
        "# Bayesian Analysis Methods Summary\n\n"
        f"We evaluated the posterior tree sample in `{posterior_tree_path.name}` against the primary log `{primary_log_path.name}`. "
        f"The analysis is described here as using a `{resolved_clock_model}` clock model and a `{resolved_tree_prior}` tree prior. "
        f"{calibration_text} {tip_date_text}\n\n"
        f"{analysis_xml_text} {chain_settings_text}\n\n"
        f"The posterior log was checked for required columns `{', '.join(required_columns)}`, monotonic state progression, and numeric parameter values; "
        f"{len(validation.issues)} validation issue(s) were detected. "
        f"Burn-in sensitivity was reviewed across fractions `{', '.join(format(value, '.15g') for value in burnin_fractions)}`, retaining between "
        f"`{retained_range[0]}` and `{retained_range[1]}` trees after filtering. "
        f"Across `{len(mixing.chain_summaries)}` chain log(s), the minimum observed effective sample size was `{format(minimum_ess, '.15g')}` "
        f"under thresholds ESS >= `{format(ess_threshold, '.15g')}` and standardized mean shift <= `{format(mean_shift_threshold, '.15g')}`. "
        f"Chain-mixing diagnostics yielded `{len(mixing.issues)}` issue(s), and the full diagnostics workflow produced `{warning_count}` warning or issue record(s).\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return BayesianMethodsSummaryTextResult(
        output_path=path,
        title="Bayesian Analysis Methods Summary",
        warning_count=warning_count,
        text=text,
    )


def _resolve_beast_analysis_xml_path(
    *,
    analysis_xml_path: Path | None,
    primary_log_path: Path,
) -> Path | None:
    if analysis_xml_path is not None:
        return analysis_xml_path
    direct_candidate = primary_log_path.with_suffix(".xml")
    if direct_candidate.exists():
        return direct_candidate
    seeded_candidate = primary_log_path.with_name(
        _strip_beast_output_suffix(primary_log_path.name) + ".xml"
    )
    if seeded_candidate.exists():
        return seeded_candidate
    return None


def _strip_beast_output_suffix(file_name: str) -> str:
    for suffix in (".log", ".trees"):
        if not file_name.endswith(suffix):
            continue
        stem = file_name[: -len(suffix)]
        if stem.endswith(".$(seed)"):
            return stem[: -len(".$(seed)")]
        parts = stem.rsplit(".", 1)
        if len(parts) == 2 and parts[1].isdigit():
            return parts[0]
        return stem
    return file_name


def write_bayesian_limitations_text(
    path: Path,
    *,
    posterior_tree_path: Path,
    primary_log_path: Path,
    additional_log_paths: list[Path] | None = None,
    tree_path: Path | None = None,
    calibration_path: Path | None = None,
    tip_dates_path: Path | None = None,
    alignment_path: Path | None = None,
    burnin_fractions: tuple[float, ...] = DEFAULT_BURNIN_FRACTIONS,
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
    cross_chain_mean_shift_threshold: float = 0.75,
) -> BayesianLimitationsTextResult:
    """Write reviewer-facing Bayesian limitations text from current diagnostics surfaces."""
    burnin = assess_beast_burnin_sensitivity(
        posterior_tree_path,
        log_path=primary_log_path,
        burnin_fractions=burnin_fractions,
    )
    mixing = assess_beast_chain_mixing(
        [primary_log_path, *(additional_log_paths or [])],
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
        cross_chain_mean_shift_threshold=cross_chain_mean_shift_threshold,
    )
    readiness = (
        assess_time_tree_readiness(
            tree_path,
            calibration_path=calibration_path,
            tip_dates_path=tip_dates_path,
            alignment_path=alignment_path,
        )
        if tree_path is not None
        else None
    )
    dominance = (
        assess_calibration_dominance(tree_path, calibration_path)
        if tree_path is not None and calibration_path is not None
        else None
    )
    limitations: list[str] = [
        "posterior clade support reflects the supplied model, priors, and sampled trees rather than direct biological truth",
        "burn-in choice can change retained topologies and parameter means, so posterior interpretation should be checked against burn-in sensitivity",
    ]
    if burnin.warnings:
        limitations.extend(burnin.warnings)
    if mixing.issues:
        limitations.append(
            "one or more chains show ESS, drift, or cross-chain mixing concerns that limit confidence in posterior summaries"
        )
    if dominance is not None and dominance.warnings:
        limitations.extend(dominance.warnings)
    if readiness is not None:
        if readiness.blockers:
            limitations.extend(readiness.blockers)
        if readiness.warnings:
            limitations.extend(readiness.warnings)
    text = (
        "# Bayesian Analysis Limitations\n\n"
        + "\n".join(f"- {line}" for line in sorted(dict.fromkeys(limitations)))
        + "\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return BayesianLimitationsTextResult(
        output_path=path,
        title="Bayesian Analysis Limitations",
        warning_count=len(sorted(dict.fromkeys(limitations))),
        text=text,
    )


def _write_clade_frequency_plot(
    path: Path,
    *,
    clade_frequencies,
    limit: int,
) -> Path:
    top_rows = sorted(
        clade_frequencies.clade_frequencies,
        key=lambda row: (-row.frequency, row.clade),
    )[: max(limit, 1)]
    width = 900
    height = 80 + 42 * max(len(top_rows), 1)
    bar_left = 250
    bar_width = 560
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="24" y="36" font-family="Avenir Next, Segoe UI, sans-serif" font-size="22" fill="#0f766e">Posterior Clade Frequencies</text>',
    ]
    for index, row in enumerate(top_rows):
        y = 70 + index * 42
        filled_width = round(bar_width * row.frequency, 3)
        lines.extend(
            [
                f'<text x="24" y="{y + 18}" font-family="SFMono-Regular, Consolas, monospace" font-size="14" fill="#1f2937">{row.clade}</text>',
                f'<rect x="{bar_left}" y="{y}" width="{bar_width}" height="22" rx="8" fill="#e2e8f0"/>',
                f'<rect x="{bar_left}" y="{y}" width="{filled_width}" height="22" rx="8" fill="#0f766e"/>',
                f'<text x="{bar_left + bar_width + 12}" y="{y + 17}" font-family="SFMono-Regular, Consolas, monospace" font-size="13" fill="#1f2937">{format(row.frequency, ".3f")}</text>',
            ]
        )
    lines.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _uncertainty_summary_markdown(
    *,
    consensus_newick: str,
    multimodality,
    conflicts,
    conclusions,
) -> str:
    robust = ", ".join(row.clade for row in conclusions.robust_clades) or "none"
    uncertain = ", ".join(row.clade for row in conclusions.uncertain_clades) or "none"
    conflicting = (
        ", ".join(row.clade for row in conclusions.conflicting_clades) or "none"
    )
    return (
        "# Posterior Uncertainty Summary\n\n"
        f"- Consensus tree: `{consensus_newick}`\n"
        f"- Topology modes detected: `{multimodality.mode_count}`\n"
        f"- High-credibility clade conflicts: `{conflicts.conflict_count}`\n"
        f"- Robust clades: `{robust}`\n"
        f"- Uncertain clades: `{uncertain}`\n"
        f"- Conflict-prone clades: `{conflicting}`\n"
    )
