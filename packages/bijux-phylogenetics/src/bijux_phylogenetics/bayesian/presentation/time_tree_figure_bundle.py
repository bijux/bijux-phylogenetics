from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.bayesian.beast.validation import assess_time_tree_readiness
from bijux_phylogenetics.bayesian.posterior_sets.tree_sets import (
    PosteriorTreeSubsamplingReport,
    subsample_beast_posterior_tree_set,
    subsample_mrbayes_posterior_tree_set,
    subsample_posterior_tree_set,
    summarize_maximum_clade_credibility_tree,
    write_posterior_tree_subsample,
)
from bijux_phylogenetics.datasets.study_inputs import (
    TaxonTable,
    load_taxon_table,
    write_taxon_rows,
)
from bijux_phylogenetics.io.newick import load_newick_tree_set
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.render.reproducibility import (
    write_figure_reproducibility_manifest,
)
from bijux_phylogenetics.render.time_tree_svg import (
    TimeTreeNodeInterval,
    TimeTreeRenderResult,
    render_time_tree_svg,
)
from bijux_phylogenetics.runtime.errors import MetadataJoinError


@dataclass(slots=True)
class TimeTreePublicationAudit:
    """Reviewer-facing publication audit for one time-tree uncertainty package."""

    publication_ready: bool
    interval_complete: bool
    ultrametric: bool
    readiness_decision: str | None
    rendered_interval_count: int
    expected_interval_count: int
    reviewer_summary: list[str]
    limitations: list[str]


@dataclass(slots=True)
class TimeTreeFigurePackageResult:
    output_dir: Path
    source_format: str
    burnin_fraction: float
    retained_tree_set_path: Path
    retained_tree_count: int
    mcc_tree_path: Path
    figure_path: Path
    interval_table_path: Path
    legend_path: Path
    caption_path: Path
    review_path: Path
    manifest_path: Path
    reproducibility_manifest_path: Path
    render: TimeTreeRenderResult
    node_intervals: list[TimeTreeNodeInterval]
    audit: TimeTreePublicationAudit


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_ready(payload: object) -> object:
    return json.loads(json.dumps(payload, default=str))


def _require_metadata_table(
    metadata_table: TaxonTable | None, *, label_column: str
) -> TaxonTable:
    if metadata_table is None:
        raise MetadataJoinError(
            f"time-tree package requires a metadata table when --label-column is set to '{label_column}'"
        )
    return metadata_table


def _build_labels(
    *,
    tree: PhyloTree,
    metadata_table: TaxonTable | None,
    label_column: str | None,
) -> dict[str, str]:
    labels = {taxon: taxon for taxon in tree.tip_names}
    if label_column is None:
        return labels
    table = _require_metadata_table(metadata_table, label_column=label_column)
    if label_column not in table.columns:
        raise MetadataJoinError(
            f"metadata table does not contain label column '{label_column}'"
        )
    for row in table.rows:
        taxon = row[table.taxon_column]
        if row[label_column]:
            labels[taxon] = row[label_column]
    return labels


def _materialize_retained_tree_set(
    source_path: Path,
    *,
    out_path: Path,
    source_format: str,
    burnin_fraction: float,
) -> PosteriorTreeSubsamplingReport:
    if source_format == "generic":
        report = subsample_posterior_tree_set(
            source_path,
            method="evenly-spaced",
            thinning_interval=1,
            burnin_fraction=burnin_fraction,
        )
    elif source_format == "beast":
        report = subsample_beast_posterior_tree_set(
            source_path,
            method="evenly-spaced",
            thinning_interval=1,
            burnin_fraction=burnin_fraction,
        )
    elif source_format == "mrbayes":
        report = subsample_mrbayes_posterior_tree_set(
            source_path,
            method="evenly-spaced",
            thinning_interval=1,
            burnin_fraction=burnin_fraction,
        )
    else:
        raise ValueError(f"unsupported posterior tree source format: {source_format}")
    write_posterior_tree_subsample(out_path, report)
    return report


def _collect_tree_depths(
    node: TreeNode,
    *,
    current_depth: float,
) -> tuple[set[str], dict[frozenset[str], float]]:
    if node.is_leaf():
        return ({node.name} if node.name is not None else set()), {}
    taxa: set[str] = set()
    age_map: dict[frozenset[str], float] = {}
    for child in node.children:
        child_taxa, child_map = _collect_tree_depths(
            child,
            current_depth=current_depth + float(child.branch_length or 0.0),
        )
        taxa.update(child_taxa)
        age_map.update(child_map)
    age_map[frozenset(taxa)] = current_depth
    return taxa, age_map


def _tip_depths(node: TreeNode, current_depth: float) -> list[float]:
    if node.is_leaf():
        return [current_depth]
    depths: list[float] = []
    for child in node.children:
        depths.extend(
            _tip_depths(child, current_depth + float(child.branch_length or 0.0))
        )
    return depths


def _quantile(sorted_values: list[float], fraction: float) -> float:
    if len(sorted_values) == 1:
        return round(sorted_values[0], 15)
    position = max(
        0,
        min(len(sorted_values) - 1, int(round(fraction * (len(sorted_values) - 1)))),
    )
    return round(sorted_values[position], 15)


def _summarize_time_tree_intervals(tree_set_path: Path) -> list[TimeTreeNodeInterval]:
    trees = load_newick_tree_set(tree_set_path)
    taxa_sets = {frozenset(tree.tip_names) for tree in trees}
    if len(taxa_sets) != 1:
        raise ValueError(
            "time-tree publication package requires all posterior trees to share the same taxon set"
        )
    age_rows: dict[frozenset[str], list[float]] = {}
    for tree in trees:
        tip_depths = _tip_depths(tree.root, 0.0)
        root_age = max(tip_depths, default=0.0)
        _taxa, depth_map = _collect_tree_depths(tree.root, current_depth=0.0)
        for clade, current_depth in depth_map.items():
            age_rows.setdefault(clade, []).append(round(root_age - current_depth, 15))
    ordered_rows: list[TimeTreeNodeInterval] = []
    total_taxa = len(next(iter(taxa_sets)))
    for clade, ages in sorted(
        age_rows.items(), key=lambda item: (len(item[0]), sorted(item[0]))
    ):
        ordered = sorted(ages)
        ordered_rows.append(
            TimeTreeNodeInterval(
                clade="|".join(sorted(clade)),
                node_kind="root" if len(clade) == total_taxa else "internal",
                mean_age=round(sum(ordered) / len(ordered), 15),
                median_age=_quantile(ordered, 0.5),
                minimum_age=round(min(ordered), 15),
                maximum_age=round(max(ordered), 15),
                lower_95_credible_interval=_quantile(ordered, 0.025),
                upper_95_credible_interval=_quantile(ordered, 0.975),
                tree_count=len(ordered),
            )
        )
    return ordered_rows


def build_time_tree_figure_package(
    posterior_tree_set_path: Path,
    *,
    out_dir: Path,
    source_format: str = "generic",
    burnin_fraction: float = 0.25,
    metadata_path: Path | None = None,
    label_column: str | None = None,
    taxon_column: str | None = None,
    tip_dates_path: Path | None = None,
    calibration_path: Path | None = None,
    alignment_path: Path | None = None,
    title: str = "Bijux Time Tree Figure",
) -> TimeTreeFigurePackageResult:
    """Build a publication-oriented time-tree package with node-age uncertainty."""
    out_dir.mkdir(parents=True, exist_ok=True)
    retained_tree_set_path = out_dir / "retained-posterior-trees.nwk"
    mcc_tree_path = out_dir / "mcc-time-tree.nwk"
    figure_path = out_dir / "time-tree.svg"
    interval_table_path = out_dir / "node-age-intervals.tsv"
    legend_path = out_dir / "time-tree-legend.tsv"
    caption_path = out_dir / "figure-caption.md"
    review_path = out_dir / "time-tree-review.html"
    manifest_path = out_dir / "time-tree-package.manifest.json"
    reproducibility_manifest_path = out_dir / "figure-reproducibility.manifest.json"

    retained = _materialize_retained_tree_set(
        posterior_tree_set_path,
        out_path=retained_tree_set_path,
        source_format=source_format,
        burnin_fraction=burnin_fraction,
    )
    mcc_tree, _mcc_report = summarize_maximum_clade_credibility_tree(
        retained_tree_set_path,
        burnin_fraction=0.0,
    )
    mcc_tree_path.write_text(mcc_tree.to_newick() + "\n", encoding="utf-8")
    node_intervals = _summarize_time_tree_intervals(retained_tree_set_path)

    metadata_table = (
        load_taxon_table(metadata_path, taxon_column=taxon_column)
        if metadata_path is not None
        else None
    )
    labels = _build_labels(
        tree=mcc_tree,
        metadata_table=metadata_table,
        label_column=label_column,
    )
    render = render_time_tree_svg(
        mcc_tree_path,
        out_path=figure_path,
        node_intervals=node_intervals,
        labels=labels,
        title=title,
    )
    readiness = (
        None
        if tip_dates_path is None
        and calibration_path is None
        and alignment_path is None
        else assess_time_tree_readiness(
            mcc_tree_path,
            calibration_path=calibration_path,
            tip_dates_path=tip_dates_path,
            alignment_path=alignment_path,
            taxon_column=taxon_column,
        )
    )

    write_taxon_rows(
        interval_table_path,
        columns=[
            "clade",
            "node_kind",
            "mean_age",
            "median_age",
            "minimum_age",
            "maximum_age",
            "lower_95_credible_interval",
            "upper_95_credible_interval",
            "tree_count",
        ],
        rows=[asdict(row) for row in node_intervals],
    )
    write_taxon_rows(
        legend_path,
        columns=["surface", "label", "detail"],
        rows=[
            {
                "surface": "node-age-label",
                "label": "median age",
                "detail": "monospace labels next to internal nodes show the posterior median node age before present",
            },
            {
                "surface": "hpd-interval",
                "label": "95% HPD interval",
                "detail": "orange whiskers show the posterior 95% highest posterior density interval for each node age",
            },
            {
                "surface": "time-axis",
                "label": "age before present",
                "detail": "the horizontal axis increases from present at the tips toward older ages at the root",
            },
        ],
    )

    expected_interval_count = render.internal_node_count
    interval_complete = render.rendered_interval_count == expected_interval_count
    limitations = list(render.warnings)
    if not interval_complete:
        limitations.append(
            "one or more internal nodes are missing posterior age intervals in the rendered figure"
        )
    if readiness is not None:
        limitations.extend(readiness.blockers)
        limitations.extend(readiness.warnings)
    reviewer_summary = [
        f"retained {retained.retained_tree_count} posterior trees after burn-in filtering",
        f"rendered {render.rendered_interval_count} node-age intervals across {render.internal_node_count} internal nodes",
        "time-tree readiness inputs passed"
        if readiness is not None and readiness.decision != "blocked"
        else (
            "time-tree readiness inputs were not supplied"
            if readiness is None
            else "time-tree readiness review reported blockers"
        ),
    ]
    audit = TimeTreePublicationAudit(
        publication_ready=interval_complete
        and render.ultrametric
        and not render.warnings
        and (readiness is None or readiness.decision != "blocked"),
        interval_complete=interval_complete,
        ultrametric=render.ultrametric,
        readiness_decision=None if readiness is None else readiness.decision,
        rendered_interval_count=render.rendered_interval_count,
        expected_interval_count=expected_interval_count,
        reviewer_summary=reviewer_summary,
        limitations=limitations,
    )

    caption_path.write_text(
        (
            f"# {title}\n\n"
            "## Draft Caption\n\n"
            f"{title} shows the maximum clade credibility time tree inferred from `{retained.retained_tree_count}` retained posterior trees. "
            "Internal node labels report median node ages before present, and orange whiskers show the corresponding 95% HPD intervals. "
            "The bottom axis is scaled in age before present so the figure can be read directly as a dated phylogeny. "
            + (
                "No additional time-tree readiness blockers were supplied for this package."
                if readiness is None or readiness.decision != "blocked"
                else "Time-tree readiness blockers remain and are summarized in the reviewer report."
            )
            + "\n\n## Figure Specifications\n\n"
            f"- Source format: `{source_format}`\n"
            f"- Burn-in fraction: `{burnin_fraction}`\n"
            f"- Retained posterior trees: `{retained.retained_tree_count}`\n"
            f"- Root age: `{format(render.root_age, '.6g')}`\n"
            f"- Rendered node-age intervals: `{render.rendered_interval_count}`\n"
            f"- Publication ready: `{audit.publication_ready}`\n"
        ),
        encoding="utf-8",
    )
    reproducibility_manifest = write_figure_reproducibility_manifest(
        reproducibility_manifest_path,
        report_kind="time_tree_package",
        input_files=[
            ("posterior_tree_set", posterior_tree_set_path),
            *([("metadata", metadata_path)] if metadata_path is not None else []),
            *([("tip_dates", tip_dates_path)] if tip_dates_path is not None else []),
            *(
                [("calibration", calibration_path)]
                if calibration_path is not None
                else []
            ),
            *([("alignment", alignment_path)] if alignment_path is not None else []),
        ],
        generated_figures=[
            ("time_tree", figure_path),
        ],
        generated_tables=[
            ("node_age_intervals", interval_table_path),
            ("legend", legend_path),
        ],
        filters=None,
        model={
            "kind": "time_tree",
            "name": source_format,
            "readiness_decision": None if readiness is None else readiness.decision,
        },
        settings={
            "burnin_fraction": burnin_fraction,
            "title": title,
            "label_column": label_column,
            "taxon_column": taxon_column,
            "retained_tree_count": retained.retained_tree_count,
        },
        linked_artifacts=[
            ("retained_tree_set", retained_tree_set_path),
            ("mcc_tree", mcc_tree_path),
            ("caption", caption_path),
        ],
    )

    manifest = {
        "report_kind": "time_tree_package",
        "title": title,
        "posterior_tree_set_path": str(posterior_tree_set_path),
        "source_format": source_format,
        "burnin_fraction": burnin_fraction,
        "metadata_path": None if metadata_path is None else str(metadata_path),
        "label_column": label_column,
        "tip_dates_path": None if tip_dates_path is None else str(tip_dates_path),
        "calibration_path": None if calibration_path is None else str(calibration_path),
        "alignment_path": None if alignment_path is None else str(alignment_path),
        "input_checksums": {
            str(path): _sha256(path)
            for path in (
                posterior_tree_set_path,
                metadata_path,
                tip_dates_path,
                calibration_path,
                alignment_path,
            )
            if path is not None
        },
        "retained_tree_set_path": str(retained_tree_set_path),
        "retained_tree_set_checksum": _sha256(retained_tree_set_path),
        "mcc_tree_path": str(mcc_tree_path),
        "mcc_tree_checksum": _sha256(mcc_tree_path),
        "figure_path": str(figure_path),
        "figure_checksum": _sha256(figure_path),
        "interval_table_path": str(interval_table_path),
        "interval_table_checksum": _sha256(interval_table_path),
        "legend_path": str(legend_path),
        "legend_checksum": _sha256(legend_path),
        "caption_path": str(caption_path),
        "caption_checksum": _sha256(caption_path),
        "review_path": str(review_path),
        "reproducibility_manifest_path": str(reproducibility_manifest_path),
        "reproducibility_manifest_checksum": _sha256(reproducibility_manifest_path),
        "reproducibility_manifest": reproducibility_manifest,
        "retained_tree_count": retained.retained_tree_count,
        "render": asdict(render),
        "node_intervals": [asdict(row) for row in node_intervals],
        "audit": asdict(audit),
        "time_tree_readiness": None if readiness is None else asdict(readiness),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )

    sections = [
        ("reviewer summary", "\n".join(f"- {line}" for line in reviewer_summary)),
        (
            "publication limitations",
            "none"
            if not limitations
            else "\n".join(f"- {line}" for line in limitations),
        ),
        (
            "node-age intervals",
            "\n".join(
                (
                    f"{row.clade} [{row.node_kind}] median={format(row.median_age, '.6g')} "
                    f"hpd95={format(row.lower_95_credible_interval, '.6g')}..{format(row.upper_95_credible_interval, '.6g')} "
                    f"trees={row.tree_count}"
                )
                for row in node_intervals
            ),
        ),
    ]
    if readiness is not None:
        sections.append(
            (
                "time-tree readiness",
                json.dumps(
                    _json_ready(asdict(readiness)),
                    indent=2,
                    sort_keys=True,
                ),
            )
        )
    write_html_report(
        title=title,
        out_path=review_path,
        embedded_json=_json_ready(manifest),
        summary_metrics=[
            ("publication_ready", audit.publication_ready),
            ("retained_tree_count", retained.retained_tree_count),
            ("root_age", format(render.root_age, ".6g")),
            ("interval_complete", audit.interval_complete),
            ("ultrametric", audit.ultrametric),
            (
                "readiness_decision",
                (
                    "not_supplied"
                    if audit.readiness_decision is None
                    else audit.readiness_decision
                ),
            ),
        ],
        artifact_links=[
            (
                "time-tree figure",
                figure_path.name,
                "dated SVG with median ages and HPD intervals",
            ),
            ("retained posterior trees", retained_tree_set_path.name, None),
            ("mcc tree", mcc_tree_path.name, None),
            ("node-age ledger", interval_table_path.name, None),
            ("legend ledger", legend_path.name, None),
            ("caption draft", caption_path.name, None),
            ("package manifest", manifest_path.name, None),
        ],
        sections=sections,
    )

    return TimeTreeFigurePackageResult(
        output_dir=out_dir,
        source_format=source_format,
        burnin_fraction=burnin_fraction,
        retained_tree_set_path=retained_tree_set_path,
        retained_tree_count=retained.retained_tree_count,
        mcc_tree_path=mcc_tree_path,
        figure_path=figure_path,
        interval_table_path=interval_table_path,
        legend_path=legend_path,
        caption_path=caption_path,
        review_path=review_path,
        manifest_path=manifest_path,
        reproducibility_manifest_path=reproducibility_manifest_path,
        render=render,
        node_intervals=node_intervals,
        audit=audit,
    )
