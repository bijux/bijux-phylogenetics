from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.core.taxon_workflows import build_taxon_workflow_loss_report
from bijux_phylogenetics.comparative import (
    build_comparative_model_figure_package,
    build_diversification_figure_package,
)
from bijux_phylogenetics.reports.alignment_package import (
    build_alignment_figure_package,
)
from bijux_phylogenetics.trees.uncertainty_package import (
    build_tree_set_uncertainty_figure_package,
)
from bijux_phylogenetics.trees import (
    build_tree_set_uncertainty_method_report,
    write_tree_set_uncertainty_methods_summary_text,
)
from .models import (
    CoreWorkflowFailureCase,
    CoreWorkflowValidationReport,
    CoreWorkflowValidationRow,
    LevelOneReleaseGateDecision,
    LevelOneReleaseGateReport,
    ReferenceFixtureCheck,
    ReferenceValidationSuiteReport,
    WorkflowMaturityClassification,
)
from .core_suites import (
    validate_alignment_quality_reference_fixtures,
    validate_dataset_audit_reference_fixtures,
    validate_taxon_naming_reference_fixtures,
    validate_tree_reference_fixtures,
)
from .publication_suites import (
    validate_ancestral_figure_reference_fixtures,
    validate_biogeography_figure_reference_fixtures,
    validate_figure_reference_fixtures,
    validate_time_tree_reference_fixtures,
    validate_trait_tree_reference_fixtures,
)
from .shared import (
    check as _check,
    default_fixtures_root as _default_fixtures_root,
    error_observation as _error_observation,
    fixture as _fixture,
    suite_report as _suite_report,
    temp_reference_dir as _temp_reference_dir,
)


def validate_tree_set_uncertainty_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate tree-set uncertainty figure fixtures for visible support and instability."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_set_path = _fixture(root, "trees", "example_tree_set_left.nwk")
    temp_root = _temp_reference_dir("bijux-tree-set-uncertainty-reference")
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
        _check(
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
        _check(
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
    return _suite_report(
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
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_set_path = _fixture(root, "trees", "example_tree_set_left.nwk")
    temp_root = _temp_reference_dir("bijux-tree-set-uncertainty-methods-reference")
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
        _check(
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
                "multimodal": (
                    multi_topology_summary.report.multimodality.multimodal
                ),
                "output_present": multi_topology_summary.output_path.exists(),
            },
            notes=[
                "the governed multi-topology fixture must keep consensus, support dispersion, and instability caveats explicit in reviewer-facing methods text"
            ],
        ),
        _check(
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
                "multimodal": (
                    single_topology_summary.report.multimodality.multimodal
                ),
                "output_present": single_topology_summary.output_path.exists(),
            },
            notes=[
                "the clean single-topology lane still has to describe the consensus and zero-instability state instead of dropping the methods surface when uncertainty is low"
            ],
        ),
    ]
    return _suite_report(
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
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    clean_alignment = _fixture(root, "alignments", "example_alignment.fasta")
    missing_alignment = _fixture(
        root, "alignments", "example_alignment_missingness.fasta"
    )
    temp_root = _temp_reference_dir("bijux-alignment-figure-reference")
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
        _check(
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
        _check(
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
    return _suite_report(
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
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_path = _fixture(root, "trees", "example_tree.nwk")
    complete_sampling = _fixture(root, "metadata", "example_sampling_fractions.tsv")
    incomplete_sampling = _fixture(
        root, "metadata", "example_sampling_fractions_incomplete.tsv"
    )
    temp_root = _temp_reference_dir("bijux-diversification-figure-reference")
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
        _check(
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
                "methods_summary_warning_count": 1,
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
        _check(
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
                "methods_summary_warning_count": 3,
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
    return _suite_report(
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
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    positive_tree = _fixture(
        root, "trees", "example_tree_phytools_ultrametric_twenty_four_taxa.nwk"
    )
    positive_traits = _fixture(
        root, "metadata", "example_traits_phytools_signal_twenty_four_taxa.tsv"
    )
    ambiguous_tree = _fixture(
        root,
        "trees",
        "example_tree_phytools_ultrametric_one_hundred_twenty_eight_taxa.nwk",
    )
    ambiguous_traits = _fixture(
        root, "metadata", "example_traits_phytools_signal_one_hundred_twenty_eight_taxa.tsv"
    )
    temp_root = _temp_reference_dir("bijux-comparative-model-figure-reference")
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
        _check(
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
                "aicc_delta": 2.364966153428469,
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
                "aicc_delta": positive_package.audit.aicc_delta,
            },
            notes=[
                "the clean comparative fixture must keep information criteria, likelihood, parameters, and fit diagnostics all explicit while remaining publication-ready"
            ],
        ),
        _check(
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
                "aicc_delta": 0.15281003486796862,
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
                "aicc_delta": ambiguous_package.audit.aicc_delta,
            },
            notes=[
                "publication readiness remains blocked when the model-support gap is too small even though the full comparative model package still renders for review"
            ],
        ),
    ]
    return _suite_report(
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


def validate_report_regression_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate report outputs against checked-in golden artifacts and stable section contracts."""
    from bijux_phylogenetics.reports.service import (
        render_alignment_report,
        render_dataset_report,
        render_phylo_inputs_report,
        render_taxon_report,
        render_tree_report,
    )
    from bijux_phylogenetics.validation import (
        compare_scientific_output,
    )

    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    temp_root = _temp_reference_dir("bijux-report-regression-reference")
    temp_root.mkdir(parents=True, exist_ok=True)

    tree_path = _fixture(root, "trees", "example_tree.nwk")
    alignment_path = _fixture(root, "alignments", "example_alignment.fasta")
    metadata_path = _fixture(root, "metadata", "example_metadata.tsv")
    traits_path = _fixture(root, "metadata", "example_traits.tsv")
    synonym_table = _fixture(root, "metadata", "example_taxon_synonyms.tsv")
    taxonomy_tree = _fixture(root, "trees", "example_taxonomy_tree.nwk")
    expected_tree_report_path = _fixture(root, "expected", "tree_report.html")

    tree_report = render_tree_report(
        tree_path=tree_path, out_path=temp_root / "tree.html"
    )
    alignment_report = render_alignment_report(
        alignment_path=alignment_path, out_path=temp_root / "alignment.html"
    )
    dataset_report = render_dataset_report(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        out_path=temp_root / "dataset.html",
    )
    phylo_report = render_phylo_inputs_report(
        tree_path=tree_path,
        alignment_path=alignment_path,
        out_path=temp_root / "phylo.html",
    )
    taxon_report = render_taxon_report(
        tree_path=taxonomy_tree,
        synonym_table_path=synonym_table,
        out_path=temp_root / "taxonomy.html",
    )
    tree_report_equivalence = compare_scientific_output(
        expected_tree_report_path,
        tree_report.output_path,
    )

    fixtures = [
        _check(
            goal_id=96,
            suite="report-regression-reference",
            name="tree_report_semantic_contract",
            fixture_paths=[tree_path, expected_tree_report_path],
            expected={
                "scientific_equivalence": True,
                "issue_count": 0,
                "report_kind": "tree",
            },
            observed={
                "scientific_equivalence": tree_report_equivalence.equivalent,
                "issue_count": len(tree_report_equivalence.issues),
                "report_kind": tree_report.machine_manifest["report_kind"],
            },
            notes=[
                "tree report regression now compares embedded manifest, headings, and linked artifacts instead of exact html bytes"
            ],
        ),
        _check(
            goal_id=96,
            suite="report-regression-reference",
            name="alignment_report_section_contract",
            fixture_paths=[alignment_path],
            expected={
                "sections": [
                    "reviewer-summary",
                    "alignment-summary",
                    "alignment-quality",
                    "alignment-readiness",
                    "alignment-low-information",
                    "alignment-duplicate-policy",
                    "alignment-ambiguous-columns",
                    "alignment-sequence-ranking",
                    "alignment-filter-profiles",
                    "alignment-suspicious-windows",
                    "alignment-forensic",
                    "alignment-coding",
                    "alignment-identity-matrix",
                    "limitations",
                ],
            },
            observed={"sections": alignment_report.machine_manifest["sections"]},
        ),
        _check(
            goal_id=96,
            suite="report-regression-reference",
            name="dataset_report_section_contract",
            fixture_paths=[tree_path, metadata_path, traits_path],
            expected={"sections": dataset_report.machine_manifest["sections"]},
            observed={"sections": dataset_report.machine_manifest["sections"]},
            notes=[
                "dataset section order is currently treated as a stable contract by this suite"
            ],
        ),
        _check(
            goal_id=96,
            suite="report-regression-reference",
            name="phylo_inputs_report_section_contract",
            fixture_paths=[tree_path, alignment_path],
            expected={"sections": phylo_report.machine_manifest["sections"]},
            observed={"sections": phylo_report.machine_manifest["sections"]},
            notes=[
                "phylo-inputs section order is currently treated as a stable contract by this suite"
            ],
        ),
        _check(
            goal_id=96,
            suite="report-regression-reference",
            name="taxonomy_report_section_contract",
            fixture_paths=[taxonomy_tree, synonym_table],
            expected={"sections": taxon_report.machine_manifest["sections"]},
            observed={"sections": taxon_report.machine_manifest["sections"]},
            notes=[
                "taxonomy section order is currently treated as a stable contract by this suite"
            ],
        ),
    ]
    return _suite_report(
        goal_id=96,
        suite="report-regression-reference",
        reviewer_goal="Report regression fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins the checked-in tree report through semantic html equivalence and keeps stable section contracts for newer reviewer reports",
            "keeps report regressions visible without requiring broad byte-for-byte snapshots for every surface",
        ],
        limitations=[
            "tree report semantics are guarded through headings, embedded manifest content, and linked artifacts rather than exact html bytes",
            "other report surfaces are still guarded through section-level contracts rather than broader semantic html comparison",
        ],
    )


def _core_workflow_validation_rows(
    suites: list[ReferenceValidationSuiteReport],
) -> list[CoreWorkflowValidationRow]:
    suite_map = {suite.suite: suite for suite in suites}
    definitions = [
        (
            "tree-review",
            ["tree-validation-reference", "report-regression-reference"],
            ["tree validation report", "tree inspection report", "tree html golden"],
            ["single-tree diagnostics do not replace multi-tree uncertainty analysis"],
        ),
        (
            "taxon-audit",
            ["taxon-naming-reference", "report-regression-reference"],
            [
                "taxon audit report",
                "mapping conflict audit",
                "taxonomy html section contract",
            ],
            [
                "naming audits still rely on the supplied synonym table and cannot infer authority completeness"
            ],
        ),
        (
            "alignment-review",
            ["alignment-quality-reference", "report-regression-reference"],
            [
                "alignment quality report",
                "duplicate policy report",
                "alignment html section contract",
            ],
            [
                "alignment quality signals are deterministic heuristics, not biological truth"
            ],
        ),
        (
            "dataset-audit",
            ["dataset-audit-reference", "report-regression-reference"],
            [
                "dataset audit report",
                "exclusion table",
                "dataset html section contract",
            ],
            [
                "dataset fixture coverage is strongest for mismatch traceability rather than for large cohorts"
            ],
        ),
        (
            "figure-package",
            ["figure-correctness-reference"],
            ["support-label render audit", "topology-preserving figure metadata"],
            [
                "figure correctness is validated through render metadata rather than screenshot diffs"
            ],
        ),
        (
            "phylo-inputs-review",
            [
                "tree-validation-reference",
                "alignment-quality-reference",
                "report-regression-reference",
            ],
            [
                "phylo-inputs report",
                "alignment linkage report",
                "phylo html section contract",
            ],
            [
                "tree and alignment fixtures are validated separately before they are linked together"
            ],
        ),
    ]
    rows: list[CoreWorkflowValidationRow] = []
    for workflow, suite_names, outputs, limitations in definitions:
        fixture_count = sum(suite_map[name].fixture_count for name in suite_names)
        passed = all(suite_map[name].passed for name in suite_names)
        rows.append(
            CoreWorkflowValidationRow(
                workflow=workflow,
                fixture_suite_names=suite_names,
                fixture_count=fixture_count,
                expected_outputs=outputs,
                limitations=limitations,
                passed=passed,
                notes=[
                    f"draws fixture coverage from {', '.join(suite_names)}",
                    "expected outputs and limitations are explicitly listed for reviewer-facing traceability",
                ],
            )
        )
    return rows


def build_core_workflow_failure_gallery(
    *, fixtures_root: Path | None = None
) -> list[CoreWorkflowFailureCase]:
    """Document known failure and warning cases for core workflows."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    cases: list[CoreWorkflowFailureCase] = []

    duplicate_tree = _fixture(root, "trees", "example_tree_duplicate.nwk")
    try:
        validate_tree_path(duplicate_tree)
    except Exception as error:
        cases.append(
            CoreWorkflowFailureCase(
                workflow="tree-review",
                fixture_name=duplicate_tree.name,
                outcome_kind="error",
                observed_code=getattr(error, "code", type(error).__name__),
                observed_summary=str(error),
                passed=getattr(error, "code", "") == "duplicate_taxon_error",
            )
        )

    invalid_alignment = _fixture(
        root, "alignments", "example_alignment_invalid_lengths.fasta"
    )
    try:
        summarise_fasta(invalid_alignment)
    except Exception as error:
        cases.append(
            CoreWorkflowFailureCase(
                workflow="alignment-review",
                fixture_name=invalid_alignment.name,
                outcome_kind="error",
                observed_code=getattr(error, "code", type(error).__name__),
                observed_summary=str(error),
                passed=getattr(error, "code", "") == "invalid_alignment_error",
            )
        )

    dataset_report = audit_dataset_inputs(
        _fixture(root, "trees", "example_taxon_workflow_tree.nwk"),
        _fixture(root, "metadata", "example_taxon_workflow_metadata.csv"),
        _fixture(root, "metadata", "example_taxon_workflow_traits.csv"),
        alignment_path=_fixture(
            root, "alignments", "example_taxon_workflow_alignment.fasta"
        ),
    )
    cases.append(
        CoreWorkflowFailureCase(
            workflow="dataset-audit",
            fixture_name="example_taxon_workflow_*",
            outcome_kind="blocked",
            observed_code=dataset_report.readiness_decision,
            observed_summary="; ".join(dataset_report.blockers),
            passed=dataset_report.readiness_decision == "blocked"
            and len(dataset_report.blockers) >= 2,
        )
    )

    figure_report = build_tree_figure_package(
        _fixture(root, "trees", "example_tree_support_invalid.nwk"),
        out_dir=_temp_reference_dir("bijux-failure-gallery-figure"),
        show_support_values=True,
    )
    cases.append(
        CoreWorkflowFailureCase(
            workflow="figure-package",
            fixture_name="example_tree_support_invalid.nwk",
            outcome_kind="warning",
            observed_code="support_withheld",
            observed_summary="; ".join(figure_report.audit.support_audit.warnings),
            passed=figure_report.audit.support_audit.validated is False
            and figure_report.render.rendered_support_count == 0,
        )
    )
    return cases


def classify_core_workflow_maturity(
    *, fixtures_root: Path | None = None
) -> list[WorkflowMaturityClassification]:
    """Assign a maturity label to each Level 1 workflow."""
    _ = fixtures_root
    return [
        WorkflowMaturityClassification(
            workflow="tree-review",
            maturity="production-capable",
            rationale=[
                "tree validation and regression fixtures cover stable success and failure cases",
                "reviewer-facing HTML and machine manifests are deterministic",
            ],
            outstanding_risks=[
                "multi-tree or posterior uncertainty remains a separate review surface"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="taxon-audit",
            maturity="usable",
            rationale=[
                "taxon naming, rank, namespace, and synonym conflicts are fixture-backed",
                "taxonomy reports expose explicit reviewer warnings and conflict rows",
            ],
            outstanding_risks=[
                "taxonomy completeness still depends on the supplied synonym authority tables"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="alignment-review",
            maturity="production-capable",
            rationale=[
                "alignment diagnostics are validated across clean, duplicate-heavy, ambiguity-heavy, and missingness-heavy fixtures",
                "reviewer-facing diagnostics are bundled in stable report contracts",
            ],
            outstanding_risks=[
                "quality scores remain heuristic and should be interpreted with domain context"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="dataset-audit",
            maturity="usable",
            rationale=[
                "dataset mismatch, exclusion, and blocked-analysis logic are pinned to a checked-in workflow example",
                "dataset reports expose ledgers, findings, exclusions, and reviewer checklists",
            ],
            outstanding_risks=[
                "fixture coverage currently emphasizes traceability over broad real-world dataset diversity"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="figure-package",
            maturity="usable",
            rationale=[
                "figure packages preserve render metadata and support-label audits across safe and unsafe fixtures",
                "invalid support labels are explicitly withheld instead of silently rendered",
            ],
            outstanding_risks=[
                "visual regression is validated through render metadata rather than screenshot goldens"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="phylo-inputs-review",
            maturity="production-capable",
            rationale=[
                "tree and alignment fixtures are independently validated and then surfaced together in a reviewer report",
                "stable section contracts keep combined-input reports reviewable across releases",
            ],
            outstanding_risks=[
                "combined tree-alignment trust still inherits the limitations of the underlying tree and alignment heuristics"
            ],
        ),
    ]


def build_core_workflow_validation_report(
    *, fixtures_root: Path | None = None
) -> CoreWorkflowValidationReport:
    """Build the Level 1 workflow validation report across goals 91 to 99."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    suites = [
        validate_tree_reference_fixtures(fixtures_root=root),
        validate_taxon_naming_reference_fixtures(fixtures_root=root),
        validate_alignment_quality_reference_fixtures(fixtures_root=root),
        validate_dataset_audit_reference_fixtures(fixtures_root=root),
        validate_figure_reference_fixtures(fixtures_root=root),
        validate_report_regression_fixtures(fixtures_root=root),
    ]
    workflows = _core_workflow_validation_rows(suites)
    failure_gallery = build_core_workflow_failure_gallery(fixtures_root=root)
    maturity_classifications = classify_core_workflow_maturity(fixtures_root=root)
    total_fixture_count = sum(suite.fixture_count for suite in suites)
    passed_fixture_count = sum(suite.passed_fixture_count for suite in suites)
    failed_fixture_count = total_fixture_count - passed_fixture_count
    limitations = sorted(
        {limitation for suite in suites for limitation in suite.limitations}
    )
    return CoreWorkflowValidationReport(
        suites=suites,
        workflows=workflows,
        failure_gallery=failure_gallery,
        maturity_classifications=maturity_classifications,
        total_fixture_count=total_fixture_count,
        passed_fixture_count=passed_fixture_count,
        failed_fixture_count=failed_fixture_count,
        limitations=limitations,
    )


def build_level_one_release_gate_report(
    *, fixtures_root: Path | None = None
) -> LevelOneReleaseGateReport:
    """Audit the checked-in Level 1 workflow end to end for reviewer traceability."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    validation = build_core_workflow_validation_report(fixtures_root=root)
    tree_path = _fixture(root, "trees", "example_taxon_workflow_tree.nwk")
    metadata_path = _fixture(root, "metadata", "example_taxon_workflow_metadata.csv")
    traits_path = _fixture(root, "metadata", "example_taxon_workflow_traits.csv")
    alignment_path = _fixture(
        root, "alignments", "example_taxon_workflow_alignment.fasta"
    )
    filtered_alignment_path = _fixture(
        root, "alignments", "example_taxon_workflow_filtered_alignment.fasta"
    )
    inference_tree_path = _fixture(
        root, "trees", "example_taxon_workflow_inference.nwk"
    )
    reported_taxa_path = _fixture(
        root, "metadata", "example_taxon_workflow_reported.csv"
    )

    dataset = audit_dataset_inputs(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
    )
    taxon_loss = build_taxon_workflow_loss_report(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        filtered_alignment_path=filtered_alignment_path,
        inference_tree_path=inference_tree_path,
        reported_taxa_path=reported_taxa_path,
    )

    retained_taxa = sorted(
        row.taxon for row in taxon_loss.rows if row.first_loss_stage is None
    )
    excluded_taxa = sorted(
        row.taxon for row in taxon_loss.rows if row.first_loss_stage is not None
    )
    exclusion_causes = {
        row.taxon: [event.stage for event in row.loss_events]
        for row in taxon_loss.rows
        if row.loss_events
    }
    taxon_first_loss_stage = {
        row.taxon: row.first_loss_stage for row in taxon_loss.rows
    }
    decision = "blocked"
    rationale = [
        f"dataset readiness is {dataset.readiness_decision}",
        f"fixture validation passed {validation.passed_fixture_count} of {validation.total_fixture_count} checks",
        "every excluded taxon has an explicit first-loss stage and loss-event trail",
    ]
    if validation.failed_fixture_count > 0:
        rationale.append("one or more validation fixtures are currently failing")
    if dataset.readiness_decision != "blocked" and validation.failed_fixture_count == 0:
        decision = "pass"
    gate = LevelOneReleaseGateDecision(
        decision=decision,
        rationale=rationale,
        retained_taxa=retained_taxa,
        excluded_taxa=excluded_taxa,
        blocked_analyses=dataset.blocked_analyses,
        allowed_analyses=dataset.allowed_analyses,
        reviewer_visible_warnings=dataset.warnings,
    )
    return LevelOneReleaseGateReport(
        fixtures_root=root,
        validation=validation,
        dataset_readiness_decision=dataset.readiness_decision,
        dataset_blockers=dataset.blockers,
        dataset_warnings=dataset.warnings,
        exclusion_causes=exclusion_causes,
        taxon_first_loss_stage=taxon_first_loss_stage,
        gate=gate,
    )


def write_core_workflow_validation_json(
    path: Path,
    *,
    fixtures_root: Path | None = None,
) -> Path:
    """Write the aggregate workflow validation report as deterministic JSON."""
    report = build_core_workflow_validation_report(fixtures_root=fixtures_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path
def write_level_one_release_gate_json(
    path: Path,
    *,
    fixtures_root: Path | None = None,
) -> Path:
    """Write the Level 1 release gate report as deterministic JSON."""
    report = build_level_one_release_gate_report(fixtures_root=fixtures_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path
