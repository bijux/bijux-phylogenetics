from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative import (
    build_comparative_model_figure_package,
    build_diversification_figure_package,
)
from bijux_phylogenetics.reports.publication.alignment import (
    build_alignment_figure_package,
)
from bijux_phylogenetics.trees import (
    build_tree_set_uncertainty_method_report,
    write_tree_set_uncertainty_methods_summary_text,
)
from bijux_phylogenetics.trees.uncertainty import (
    build_tree_set_uncertainty_figure_package,
)

from .models import ReferenceValidationSuiteReport
from .shared import (
    check,
    default_fixtures_root,
    fixture,
    suite_report,
    temp_reference_dir,
)


def validate_tree_set_uncertainty_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate tree-set uncertainty figure fixtures for visible support and instability."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_set_path = fixture(root, "trees", "example_tree_set_left.nwk")
    temp_root = temp_reference_dir("bijux-tree-set-uncertainty-reference")
    temp_root.mkdir(parents=True, exist_ok=True)
    single_topology_path = temp_root / "single-topology-tree-set.nwk"
    repeated_tree = tree_set_path.read_text(encoding="utf-8").splitlines()[0]
    single_topology_path.write_text(
        "\n".join([repeated_tree, repeated_tree, repeated_tree]) + "\n",
        encoding="utf-8",
    )

    multi_topology_package = build_tree_set_uncertainty_figure_package(
        tree_set_path,
        out_dir=temp_root / "multi-topology-package",
    )
    single_topology_package = build_tree_set_uncertainty_figure_package(
        single_topology_path,
        out_dir=temp_root / "single-topology-package",
    )

    fixtures = [
        check(
            goal_id=231,
            suite="tree-set-uncertainty-publication-reference",
            name="multi_topology_tree_set_uncertainty_visible",
            fixture_paths=[tree_set_path],
            expected={
                "publication_ready": True,
                "support_labels_validated": True,
                "unstable_taxa_visible": True,
                "topology_clusters_visible": True,
                "plotted_topology_cluster_count": 2,
            },
            observed={
                "publication_ready": multi_topology_package.audit.publication_ready,
                "support_labels_validated": (
                    multi_topology_package.audit.support_labels_validated
                ),
                "unstable_taxa_visible": (
                    multi_topology_package.audit.unstable_taxa_visible
                ),
                "topology_clusters_visible": (
                    multi_topology_package.audit.topology_clusters_visible
                ),
                "plotted_topology_cluster_count": (
                    multi_topology_package.audit.plotted_topology_cluster_count
                ),
            },
            notes=[
                "the governed multi-topology fixture must keep consensus support, unstable taxa, and topology clusters visible on the figure surfaces"
            ],
        ),
        check(
            goal_id=231,
            suite="tree-set-uncertainty-publication-reference",
            name="single_topology_tree_set_keeps_empty_instability_panel",
            fixture_paths=[single_topology_path],
            expected={
                "publication_ready": True,
                "unstable_taxon_count": 0,
                "plotted_unstable_taxon_count": 0,
                "unstable_taxa_visible": True,
                "plotted_topology_cluster_count": 1,
            },
            observed={
                "publication_ready": single_topology_package.audit.publication_ready,
                "unstable_taxon_count": (
                    single_topology_package.audit.unstable_taxon_count
                ),
                "plotted_unstable_taxon_count": (
                    single_topology_package.audit.plotted_unstable_taxon_count
                ),
                "unstable_taxa_visible": (
                    single_topology_package.audit.unstable_taxa_visible
                ),
                "plotted_topology_cluster_count": (
                    single_topology_package.audit.plotted_topology_cluster_count
                ),
            },
            notes=[
                "the package must keep the instability panel explicit even when all trees share one topology and no taxon changes placement"
            ],
        ),
    ]
    return suite_report(
        goal_id=231,
        suite="tree-set-uncertainty-publication-reference",
        reviewer_goal="Tree-set uncertainty figure fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins one governed multi-topology tree set where consensus support, unstable taxa, clade support, and topology clusters all remain visible on the publication package surface",
            "proves that a single-topology tree set still emits an explicit empty instability panel instead of silently dropping that uncertainty surface",
        ],
        limitations=[
            "tree-set publication readiness is governed through rendered support counts, plot-presence checks, and explicit empty-state review rather than screenshot goldens",
        ],
    )


def validate_tree_set_uncertainty_methods_summary_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate tree-set uncertainty methods-summary fixtures for multi-topology and single-topology lanes."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_set_path = fixture(root, "trees", "example_tree_set_left.nwk")
    temp_root = temp_reference_dir("bijux-tree-set-uncertainty-methods-reference")
    temp_root.mkdir(parents=True, exist_ok=True)
    single_topology_path = temp_root / "single-topology-tree-set.nwk"
    repeated_tree = tree_set_path.read_text(encoding="utf-8").splitlines()[0]
    single_topology_path.write_text(
        "\n".join([repeated_tree, repeated_tree, repeated_tree]) + "\n",
        encoding="utf-8",
    )

    multi_topology_summary = write_tree_set_uncertainty_methods_summary_text(
        temp_root / "multi-topology-methods-summary.md",
        build_tree_set_uncertainty_method_report(tree_set_path),
    )
    single_topology_summary = write_tree_set_uncertainty_methods_summary_text(
        temp_root / "single-topology-methods-summary.md",
        build_tree_set_uncertainty_method_report(single_topology_path),
    )

    fixtures = [
        check(
            goal_id=253,
            suite="tree-set-uncertainty-methods-reference",
            name="multi_topology_methods_summary_keeps_support_and_instability",
            fixture_paths=[tree_set_path],
            expected={
                "warning_count": 4,
                "topology_cluster_count": 2,
                "unstable_taxon_count": 4,
                "multimodal": True,
                "output_present": True,
            },
            observed={
                "warning_count": multi_topology_summary.warning_count,
                "topology_cluster_count": (
                    multi_topology_summary.topology_cluster_count
                ),
                "unstable_taxon_count": multi_topology_summary.unstable_taxon_count,
                "multimodal": (multi_topology_summary.report.multimodality.multimodal),
                "output_present": multi_topology_summary.output_path.exists(),
            },
            notes=[
                "the governed multi-topology fixture must keep consensus, support dispersion, and instability caveats explicit in reviewer-facing methods text"
            ],
        ),
        check(
            goal_id=253,
            suite="tree-set-uncertainty-methods-reference",
            name="single_topology_methods_summary_keeps_clean_lane_explicit",
            fixture_paths=[single_topology_path],
            expected={
                "warning_count": 0,
                "topology_cluster_count": 1,
                "unstable_taxon_count": 0,
                "multimodal": False,
                "output_present": True,
            },
            observed={
                "warning_count": single_topology_summary.warning_count,
                "topology_cluster_count": (
                    single_topology_summary.topology_cluster_count
                ),
                "unstable_taxon_count": single_topology_summary.unstable_taxon_count,
                "multimodal": (single_topology_summary.report.multimodality.multimodal),
                "output_present": single_topology_summary.output_path.exists(),
            },
            notes=[
                "the clean single-topology lane still has to describe the consensus and zero-instability state instead of dropping the methods surface when uncertainty is low"
            ],
        ),
    ]
    return suite_report(
        goal_id=253,
        suite="tree-set-uncertainty-methods-reference",
        reviewer_goal="Tree-set uncertainty methods-summary fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins one governed multi-topology tree set where the methods summary keeps support dispersion, multimodality, and instability warnings explicit",
            "proves that a single-topology tree set still emits a real methods summary with a clean zero-warning lane instead of collapsing into an omitted artifact",
        ],
        limitations=[
            "tree-set methods-summary governance is based on output presence and exact uncertainty counts rather than text snapshot goldens",
        ],
    )


def validate_alignment_figure_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate alignment figure fixtures for visible quality surfaces and blocked suspicious cases."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    clean_alignment = fixture(root, "alignments", "example_alignment.fasta")
    missing_alignment = fixture(
        root, "alignments", "example_alignment_missingness.fasta"
    )
    temp_root = temp_reference_dir("bijux-alignment-figure-reference")
    temp_root.mkdir(parents=True, exist_ok=True)

    clean_package = build_alignment_figure_package(
        clean_alignment,
        out_dir=temp_root / "clean-alignment-package",
    )
    missing_package = build_alignment_figure_package(
        missing_alignment,
        out_dir=temp_root / "missingness-alignment-package",
    )

    fixtures = [
        check(
            goal_id=232,
            suite="alignment-figure-publication-reference",
            name="clean_alignment_quality_figures_ready",
            fixture_paths=[clean_alignment],
            expected={
                "publication_ready": True,
                "heatmap_visible": True,
                "site_summary_visible": True,
                "sequence_panel_visible": True,
                "heatmap_bin_count": 8,
                "plotted_window_count": 1,
                "plotted_sequence_count": 4,
            },
            observed={
                "publication_ready": clean_package.audit.publication_ready,
                "heatmap_visible": clean_package.audit.heatmap_visible,
                "site_summary_visible": clean_package.audit.site_summary_visible,
                "sequence_panel_visible": clean_package.audit.sequence_panel_visible,
                "heatmap_bin_count": clean_package.audit.heatmap_bin_count,
                "plotted_window_count": clean_package.audit.plotted_window_count,
                "plotted_sequence_count": clean_package.audit.plotted_sequence_count,
            },
            notes=[
                "the clean alignment fixture must keep the missingness heatmap, site-quality summary, and sequence-quality panel all visible on reviewer surfaces"
            ],
        ),
        check(
            goal_id=232,
            suite="alignment-figure-publication-reference",
            name="missingness_alignment_quality_figures_blocked",
            fixture_paths=[missing_alignment],
            expected={
                "publication_ready": False,
                "suspicious_alignment": True,
                "heatmap_visible": True,
                "site_summary_visible": True,
                "sequence_panel_visible": True,
                "heatmap_bin_count": 6,
                "plotted_window_count": 1,
                "plotted_sequence_count": 3,
            },
            observed={
                "publication_ready": missing_package.audit.publication_ready,
                "suspicious_alignment": missing_package.audit.suspicious_alignment,
                "heatmap_visible": missing_package.audit.heatmap_visible,
                "site_summary_visible": missing_package.audit.site_summary_visible,
                "sequence_panel_visible": missing_package.audit.sequence_panel_visible,
                "heatmap_bin_count": missing_package.audit.heatmap_bin_count,
                "plotted_window_count": missing_package.audit.plotted_window_count,
                "plotted_sequence_count": missing_package.audit.plotted_sequence_count,
            },
            notes=[
                "publication readiness remains blocked for the missingness-heavy fixture even though all three figure surfaces still render for review"
            ],
        ),
    ]
    return suite_report(
        goal_id=232,
        suite="alignment-figure-publication-reference",
        reviewer_goal="Alignment quality figure publication fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins one clean alignment figure package where missingness, site windows, and sequence burden all remain visible and publication-ready",
            "proves that suspicious alignments still render the full figure package while the publication audit remains blocked instead of silently claiming readiness",
        ],
        limitations=[
            "alignment figure publication readiness is governed through visible figure-surface counts, suspicious-alignment review, and machine-readable audits rather than screenshot goldens",
        ],
    )


def validate_diversification_figure_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate diversification figure fixtures for visible macroevolution surfaces and sampling-aware blocking."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_path = fixture(root, "trees", "example_tree.nwk")
    complete_sampling = fixture(root, "metadata", "example_sampling_fractions.tsv")
    incomplete_sampling = fixture(
        root, "metadata", "example_sampling_fractions_incomplete.tsv"
    )
    temp_root = temp_reference_dir("bijux-diversification-figure-reference")
    temp_root.mkdir(parents=True, exist_ok=True)

    complete_package = build_diversification_figure_package(
        tree_path,
        metadata_path=complete_sampling,
        out_dir=temp_root / "complete-diversification-package",
    )
    blocked_package = build_diversification_figure_package(
        tree_path,
        metadata_path=incomplete_sampling,
        out_dir=temp_root / "incomplete-sampling-diversification-package",
    )

    fixtures = [
        check(
            goal_id=233,
            suite="diversification-figure-publication-reference",
            name="sampling_complete_diversification_figures_ready",
            fixture_paths=[tree_path, complete_sampling],
            expected={
                "publication_ready": True,
                "sampling_metadata_complete": True,
                "plotted_ltt_point_count": 4,
                "plotted_clade_count": 3,
                "highlighted_outlier_count": 2,
                "plotted_model_count": 2,
                "better_model": "yule",
                "methods_summary_present": True,
                "methods_summary_warning_count": 2,
            },
            observed={
                "publication_ready": complete_package.audit.publication_ready,
                "sampling_metadata_complete": (
                    complete_package.audit.sampling_metadata_complete
                ),
                "plotted_ltt_point_count": (
                    complete_package.audit.plotted_ltt_point_count
                ),
                "plotted_clade_count": complete_package.audit.plotted_clade_count,
                "highlighted_outlier_count": (
                    complete_package.audit.highlighted_outlier_count
                ),
                "plotted_model_count": complete_package.audit.plotted_model_count,
                "better_model": complete_package.audit.better_model,
                "methods_summary_present": complete_package.methods_summary_path.exists(),
                "methods_summary_warning_count": (
                    complete_package.methods_summary.warning_count
                ),
            },
            notes=[
                "the clean diversification fixture must keep the lineage-through-time curve, clade-rate outlier panel, and model-comparison surface all visible while remaining sampling-complete",
                "the publication bundle must also materialize a diversification methods summary instead of leaving the reviewer-facing package to figures and ledgers alone",
            ],
        ),
        check(
            goal_id=233,
            suite="diversification-figure-publication-reference",
            name="incomplete_sampling_blocks_diversification_readiness",
            fixture_paths=[tree_path, incomplete_sampling],
            expected={
                "publication_ready": False,
                "sampling_metadata_complete": False,
                "lineage_curve_visible": True,
                "clade_outlier_surface_visible": True,
                "model_comparison_visible": True,
                "better_model": "yule",
                "methods_summary_present": True,
                "methods_summary_warning_count": 4,
            },
            observed={
                "publication_ready": blocked_package.audit.publication_ready,
                "sampling_metadata_complete": (
                    blocked_package.audit.sampling_metadata_complete
                ),
                "lineage_curve_visible": blocked_package.audit.lineage_curve_visible,
                "clade_outlier_surface_visible": (
                    blocked_package.audit.clade_outlier_surface_visible
                ),
                "model_comparison_visible": (
                    blocked_package.audit.model_comparison_visible
                ),
                "better_model": blocked_package.audit.better_model,
                "methods_summary_present": blocked_package.methods_summary_path.exists(),
                "methods_summary_warning_count": (
                    blocked_package.methods_summary.warning_count
                ),
            },
            notes=[
                "publication readiness remains blocked when sampling metadata is incomplete even though the full diversification figure package still renders for review",
                "the blocked lane still has to keep the methods summary artifact, so reviewer-facing caveats remain explicit instead of disappearing with readiness",
            ],
        ),
    ]
    return suite_report(
        goal_id=233,
        suite="diversification-figure-publication-reference",
        reviewer_goal="Diversification figure publication fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins one sampling-complete diversification package where lineage accumulation, clade outliers, and model support all remain explicit on reviewer-facing figures",
            "proves that incomplete sampling metadata still allows full figure rendering while keeping publication readiness blocked instead of silently claiming a corrected diversification review",
        ],
        limitations=[
            "diversification figure publication readiness is governed through rendered surface counts, sampling-metadata completeness, and machine-readable audits rather than screenshot goldens",
        ],
    )


def validate_comparative_model_figure_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate comparative model-comparison figure fixtures for decisive and ambiguous support lanes."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    positive_tree = fixture(
        root, "trees", "example_tree_phytools_ultrametric_twenty_four_taxa.nwk"
    )
    positive_traits = fixture(
        root, "metadata", "example_traits_phytools_signal_twenty_four_taxa.tsv"
    )
    ambiguous_tree = fixture(
        root,
        "trees",
        "example_tree_phytools_ultrametric_one_hundred_twenty_eight_taxa.nwk",
    )
    ambiguous_traits = fixture(
        root,
        "metadata",
        "example_traits_phytools_signal_one_hundred_twenty_eight_taxa.tsv",
    )
    temp_root = temp_reference_dir("bijux-comparative-model-figure-reference")
    temp_root.mkdir(parents=True, exist_ok=True)

    positive_package = build_comparative_model_figure_package(
        positive_tree,
        positive_traits,
        trait="signal_strong",
        out_dir=temp_root / "signal-strong-model-package",
    )
    ambiguous_package = build_comparative_model_figure_package(
        ambiguous_tree,
        ambiguous_traits,
        trait="signal_strong",
        out_dir=temp_root / "signal-strong-ambiguous-model-package",
    )

    fixtures = [
        check(
            goal_id=234,
            suite="comparative-model-figure-reference",
            name="signal_strong_twenty_four_taxa_model_figures_ready",
            fixture_paths=[positive_tree, positive_traits],
            expected={
                "publication_ready": True,
                "selected_model": "brownian",
                "support_distinct": True,
                "finite_aicc_model_count": 2,
                "plotted_model_count": 2,
                "rendered_parameter_count": 5,
                "rendered_fit_row_count": 2,
                "aicc_delta": 2.364966153428,
            },
            observed={
                "publication_ready": positive_package.audit.publication_ready,
                "selected_model": positive_package.audit.selected_model,
                "support_distinct": positive_package.audit.support_distinct,
                "finite_aicc_model_count": (
                    positive_package.audit.finite_aicc_model_count
                ),
                "plotted_model_count": positive_package.audit.plotted_model_count,
                "rendered_parameter_count": (
                    positive_package.audit.rendered_parameter_count
                ),
                "rendered_fit_row_count": (
                    positive_package.audit.rendered_fit_row_count
                ),
                "aicc_delta": round(positive_package.audit.aicc_delta, 12),
            },
            notes=[
                "the clean comparative fixture must keep information criteria, likelihood, parameters, and fit diagnostics all explicit while remaining publication-ready"
            ],
        ),
        check(
            goal_id=234,
            suite="comparative-model-figure-reference",
            name="signal_strong_one_hundred_twenty_eight_taxa_support_ambiguous",
            fixture_paths=[ambiguous_tree, ambiguous_traits],
            expected={
                "publication_ready": False,
                "selected_model": "ou",
                "support_distinct": False,
                "finite_aicc_model_count": 2,
                "criteria_surface_visible": True,
                "likelihood_surface_visible": True,
                "parameter_surface_visible": True,
                "fit_surface_visible": True,
                "aicc_delta": 0.152810034867,
            },
            observed={
                "publication_ready": ambiguous_package.audit.publication_ready,
                "selected_model": ambiguous_package.audit.selected_model,
                "support_distinct": ambiguous_package.audit.support_distinct,
                "finite_aicc_model_count": (
                    ambiguous_package.audit.finite_aicc_model_count
                ),
                "criteria_surface_visible": (
                    ambiguous_package.audit.criteria_surface_visible
                ),
                "likelihood_surface_visible": (
                    ambiguous_package.audit.likelihood_surface_visible
                ),
                "parameter_surface_visible": (
                    ambiguous_package.audit.parameter_surface_visible
                ),
                "fit_surface_visible": ambiguous_package.audit.fit_surface_visible,
                "aicc_delta": round(ambiguous_package.audit.aicc_delta, 12),
            },
            notes=[
                "publication readiness remains blocked when the model-support gap is too small even though the full comparative model package still renders for review"
            ],
        ),
    ]
    return suite_report(
        goal_id=234,
        suite="comparative-model-figure-reference",
        reviewer_goal="Comparative model-comparison publication fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins one clean comparative model package where information criteria, likelihood, parameters, and fit diagnostics all remain publication-ready",
            "proves that weak AICc separation still allows full figure rendering while keeping publication readiness blocked instead of silently claiming a decisive winner",
        ],
        limitations=[
            "comparative model figure publication readiness is governed through rendered surface counts, finite information criteria, and explicit AICc support thresholds rather than screenshot goldens",
        ],
    )
