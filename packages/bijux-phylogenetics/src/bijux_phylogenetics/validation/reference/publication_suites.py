from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral import build_ancestral_figure_package
from bijux_phylogenetics.bayesian import build_time_tree_figure_package
from bijux_phylogenetics.biogeography.presentation import (
    build_biogeography_report_package,
)
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_beast_posterior_fixture,
)
from bijux_phylogenetics.render.annotated_trait_tree_package import (
    build_annotated_trait_tree_package,
)
from bijux_phylogenetics.render.tree_figure_package import build_tree_figure_package

from .models import ReferenceValidationSuiteReport
from .shared import (
    check,
    default_fixtures_root,
    fixture,
    suite_report,
    temp_reference_dir,
)


def validate_figure_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate figure correctness fixtures for topology and support-label auditing."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    valid_support_tree = fixture(root, "trees", "example_tree_support_right.nwk")
    invalid_support_tree = fixture(root, "trees", "example_tree_support_invalid.nwk")

    valid_figure = build_tree_figure_package(
        valid_support_tree,
        out_dir=temp_reference_dir("bijux-valid-figure-reference"),
        show_support_values=True,
    )
    invalid_figure = build_tree_figure_package(
        invalid_support_tree,
        out_dir=temp_reference_dir("bijux-invalid-figure-reference"),
        show_support_values=True,
    )
    fixtures = [
        check(
            goal_id=95,
            suite="figure-correctness-reference",
            name="validated_support_figure",
            fixture_paths=[valid_support_tree],
            expected={
                "support_validated": True,
                "rendered_support_count": 2,
                "table_consistent": True,
                "scale_bar_valid": True,
                "legend_complete": True,
                "caption_ready": True,
                "legible": True,
            },
            observed={
                "support_validated": valid_figure.audit.support_audit.validated,
                "rendered_support_count": valid_figure.render.rendered_support_count,
                "table_consistent": valid_figure.audit.table_consistency.consistent,
                "scale_bar_valid": valid_figure.audit.scale_bar_valid,
                "legend_complete": valid_figure.audit.legend_audit.complete,
                "caption_ready": valid_figure.caption_draft.caption_ready,
                "legible": valid_figure.legibility_audit.legible,
            },
        ),
        check(
            goal_id=95,
            suite="figure-correctness-reference",
            name="withheld_invalid_support_figure",
            fixture_paths=[invalid_support_tree],
            expected={
                "support_validated": False,
                "rendered_support_count": 0,
                "support_warning_count": 4,
                "legend_complete": True,
            },
            observed={
                "support_validated": invalid_figure.audit.support_audit.validated,
                "rendered_support_count": invalid_figure.render.rendered_support_count,
                "support_warning_count": len(
                    invalid_figure.audit.support_audit.warnings
                ),
                "legend_complete": invalid_figure.audit.legend_audit.complete,
            },
        ),
    ]
    return suite_report(
        goal_id=95,
        suite="figure-correctness-reference",
        reviewer_goal="Figure correctness reference fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "covers both safe support-label rendering and automatic withholding when support values are uninterpretable",
            "ties figure correctness to topology-preserving render audits rather than to screenshot diffs alone",
        ],
        limitations=[
            "the current suite audits figure correctness through stable render metadata instead of pixel-level visual comparison",
        ],
    )


def validate_trait_tree_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate annotated trait tree publication fixtures for coverage and readability."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_path = fixture(root, "trees", "example_tree_support_left.nwk")
    metadata_path = fixture(root, "metadata", "example_metadata.tsv")
    traits_path = fixture(root, "metadata", "example_traits_validate.tsv")
    temp_root = temp_reference_dir("bijux-trait-tree-reference")
    temp_root.mkdir(parents=True, exist_ok=True)
    incomplete_metadata_path = temp_root / "example_metadata_missing_location.tsv"
    incomplete_metadata_path.write_text(
        "\n".join(
            [
                "taxon\tspecies\tlocation",
                "A\tAlpha species\tSweden",
                "B\tBeta species\tNorway",
                "C\tGamma species\t",
                "D\tDelta species\tFinland",
                "",
            ]
        ),
        encoding="utf-8",
    )

    valid_package = build_annotated_trait_tree_package(
        tree_path,
        out_dir=temp_root / "valid-package",
        metadata_path=metadata_path,
        traits_path=traits_path,
        label_column="species",
        categorical_column="habitat",
        continuous_column="height_cm",
        metadata_strip_columns=["location"],
        heatmap_columns=["height_cm"],
    )
    incomplete_package = build_annotated_trait_tree_package(
        tree_path,
        out_dir=temp_root / "incomplete-package",
        metadata_path=incomplete_metadata_path,
        traits_path=traits_path,
        label_column="species",
        categorical_column="habitat",
        metadata_strip_columns=["location"],
    )
    incomplete_location_row = next(
        row for row in incomplete_package.coverage_rows if row.surface == "location"
    )

    fixtures = [
        check(
            goal_id=227,
            suite="trait-tree-publication-reference",
            name="complete_annotated_trait_tree",
            fixture_paths=[tree_path, metadata_path, traits_path],
            expected={
                "publication_ready": True,
                "required_surface_count": 5,
                "complete_surface_count": 5,
                "missing_surface_count": 0,
                "caption_ready": True,
                "legible": True,
            },
            observed={
                "publication_ready": valid_package.audit.publication_ready,
                "required_surface_count": valid_package.audit.required_surface_count,
                "complete_surface_count": valid_package.audit.complete_surface_count,
                "missing_surface_count": valid_package.audit.missing_surface_count,
                "caption_ready": valid_package.audit.caption_ready,
                "legible": valid_package.audit.legible,
            },
        ),
        check(
            goal_id=227,
            suite="trait-tree-publication-reference",
            name="incomplete_metadata_strip_blocks_publication",
            fixture_paths=[tree_path, incomplete_metadata_path, traits_path],
            expected={
                "publication_ready": False,
                "missing_surface_count": 1,
                "incomplete_surface": "location",
                "location_complete": False,
                "location_missing_taxa": ["C"],
            },
            observed={
                "publication_ready": incomplete_package.audit.publication_ready,
                "missing_surface_count": incomplete_package.audit.missing_surface_count,
                "incomplete_surface": incomplete_location_row.surface,
                "location_complete": incomplete_location_row.complete,
                "location_missing_taxa": incomplete_location_row.missing_taxa,
            },
            notes=[
                "publication readiness is intentionally blocked when one requested metadata strip omits a tree taxon"
            ],
        ),
    ]
    return suite_report(
        goal_id=227,
        suite="trait-tree-publication-reference",
        reviewer_goal="Annotated trait tree publication fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins a fully covered annotated trait tree package across labels, traits, metadata strips, and reviewer ledgers",
            "proves that incomplete requested annotation surfaces block publication readiness instead of silently degrading the figure package",
        ],
        limitations=[
            "trait-tree publication readiness is governed through render metadata and annotation ledgers rather than screenshot goldens",
        ],
    )


def validate_time_tree_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate time-tree publication fixtures for visible uncertainty and readiness."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    posterior_tree_path = (
        get_shared_beast_posterior_fixture(
            "strict_yule_real_posterior"
        ).posterior_trees_path
        if fixtures_root is None
        else fixture(root, "metadata", "beast2_strict_yule_posterior.trees")
    )
    metadata_path = fixture(root, "metadata", "example_metadata.tsv")
    tip_dates_path = fixture(root, "metadata", "example_tip_dates.tsv")
    invalid_tip_dates_path = fixture(root, "metadata", "example_tip_dates_invalid.tsv")
    temp_root = temp_reference_dir("bijux-time-tree-reference")
    temp_root.mkdir(parents=True, exist_ok=True)

    valid_package = build_time_tree_figure_package(
        posterior_tree_path,
        out_dir=temp_root / "valid-package",
        source_format="beast",
        burnin_fraction=0.25,
        metadata_path=metadata_path,
        label_column="species",
        tip_dates_path=tip_dates_path,
        title="Rabies time tree",
    )
    blocked_package = build_time_tree_figure_package(
        posterior_tree_path,
        out_dir=temp_root / "blocked-package",
        source_format="beast",
        burnin_fraction=0.25,
        metadata_path=metadata_path,
        label_column="species",
        tip_dates_path=invalid_tip_dates_path,
        title="Rabies time tree",
    )

    fixtures = [
        check(
            goal_id=228,
            suite="time-tree-publication-reference",
            name="visible_uncertainty_time_tree",
            fixture_paths=[posterior_tree_path, metadata_path, tip_dates_path],
            expected={
                "publication_ready": True,
                "readiness_decision": "ready",
                "interval_complete": True,
                "ultrametric": True,
                "rendered_interval_count_matches": True,
            },
            observed={
                "publication_ready": valid_package.audit.publication_ready,
                "readiness_decision": valid_package.audit.readiness_decision,
                "interval_complete": valid_package.audit.interval_complete,
                "ultrametric": valid_package.audit.ultrametric,
                "rendered_interval_count_matches": (
                    valid_package.render.rendered_interval_count
                    == valid_package.audit.expected_interval_count
                ),
            },
            notes=[
                "the governed BEAST fixture must render visible node-age labels and HPD intervals for every internal node"
            ],
        ),
        check(
            goal_id=228,
            suite="time-tree-publication-reference",
            name="invalid_tip_dates_block_publication",
            fixture_paths=[
                posterior_tree_path,
                metadata_path,
                invalid_tip_dates_path,
            ],
            expected={
                "publication_ready": False,
                "readiness_decision": "blocked",
                "interval_complete": True,
                "ultrametric": True,
                "has_tip_date_limitation": True,
            },
            observed={
                "publication_ready": blocked_package.audit.publication_ready,
                "readiness_decision": blocked_package.audit.readiness_decision,
                "interval_complete": blocked_package.audit.interval_complete,
                "ultrametric": blocked_package.audit.ultrametric,
                "has_tip_date_limitation": any(
                    "tip-date" in limitation.lower()
                    for limitation in blocked_package.audit.limitations
                ),
            },
            notes=[
                "publication readiness remains blocked when the time calibration evidence is invalid even if uncertainty intervals still render"
            ],
        ),
    ]
    return suite_report(
        goal_id=228,
        suite="time-tree-publication-reference",
        reviewer_goal="Time-tree publication fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins one governed BEAST posterior fixture where dated uncertainty stays visible through node-age labels, HPD intervals, and reviewer-facing package outputs",
            "proves that publication readiness is blocked when time-tree readiness evidence is invalid instead of silently treating a visually rendered figure as review-complete",
        ],
        limitations=[
            "time-tree publication readiness is governed through explicit interval ledgers, ultrametric checks, and readiness audits rather than screenshot goldens",
        ],
    )


def validate_ancestral_figure_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate ancestral publication fixtures for visible and interpretable uncertainty."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_path = fixture(root, "trees", "example_tree.nwk")
    continuous_traits_path = fixture(root, "metadata", "example_traits_comparative.tsv")
    discrete_traits_path = fixture(root, "metadata", "example_traits_geography.tsv")
    temp_root = temp_reference_dir("bijux-ancestral-figure-reference")
    temp_root.mkdir(parents=True, exist_ok=True)

    continuous_package = build_ancestral_figure_package(
        tree_path=tree_path,
        traits_path=continuous_traits_path,
        trait="response",
        reconstruction_kind="continuous",
        out_dir=temp_root / "continuous-package",
        model="brownian",
    )
    discrete_package = build_ancestral_figure_package(
        tree_path=tree_path,
        traits_path=discrete_traits_path,
        trait="region",
        reconstruction_kind="discrete",
        out_dir=temp_root / "discrete-package",
        model="equal-rates",
    )

    fixtures = [
        check(
            goal_id=229,
            suite="ancestral-figure-publication-reference",
            name="continuous_intervals_visible",
            fixture_paths=[tree_path, continuous_traits_path],
            expected={
                "publication_ready": True,
                "internal_state_visible": True,
                "uncertainty_visible": True,
                "rendered_internal_annotation_count_matches": True,
            },
            observed={
                "publication_ready": continuous_package.audit.publication_ready,
                "internal_state_visible": continuous_package.audit.internal_state_visible,
                "uncertainty_visible": continuous_package.audit.uncertainty_visible,
                "rendered_internal_annotation_count_matches": (
                    continuous_package.audit.rendered_internal_annotation_count
                    == continuous_package.audit.internal_node_count
                ),
            },
        ),
        check(
            goal_id=229,
            suite="ancestral-figure-publication-reference",
            name="discrete_probabilities_interpretable",
            fixture_paths=[tree_path, discrete_traits_path],
            expected={
                "publication_ready": True,
                "internal_state_visible": True,
                "uncertainty_visible": True,
                "rendered_internal_pie_count_matches": True,
            },
            observed={
                "publication_ready": discrete_package.audit.publication_ready,
                "internal_state_visible": discrete_package.audit.internal_state_visible,
                "uncertainty_visible": discrete_package.audit.uncertainty_visible,
                "rendered_internal_pie_count_matches": (
                    discrete_package.audit.rendered_internal_pie_count
                    == discrete_package.audit.internal_node_count
                ),
            },
        ),
    ]
    return suite_report(
        goal_id=229,
        suite="ancestral-figure-publication-reference",
        reviewer_goal="Ancestral figure publication fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins one continuous ancestral package where internal node labels carry explicit uncertainty and one discrete package where probability pies remain paired with readable confidence labels",
            "keeps ancestral figure publication readiness tied to rendered uncertainty visibility rather than to file-count claims or off-figure tables alone",
        ],
        limitations=[
            "ancestral figure publication readiness is audited through render metadata and reviewer ledgers rather than screenshot goldens",
        ],
    )


def validate_biogeography_figure_reference_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate biogeography publication fixtures for visible states and transitions."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_path = fixture(root, "trees", "example_tree.nwk")
    traits_path = fixture(root, "metadata", "example_traits_geography.tsv")
    centroids_path = fixture(
        root, "metadata", "example_geographic_region_centroids.tsv"
    )
    temp_root = temp_reference_dir("bijux-biogeography-figure-reference")
    temp_root.mkdir(parents=True, exist_ok=True)
    incomplete_centroids_path = (
        temp_root / "example_geographic_region_centroids_missing_island.tsv"
    )
    incomplete_centroids_path.write_text(
        "\n".join(
            [
                "region\tlatitude\tlongitude",
                "north\t59.33\t18.07",
                "south\t-33.45\t-70.66",
                "",
            ]
        ),
        encoding="utf-8",
    )

    complete_package = build_biogeography_report_package(
        tree_path=tree_path,
        traits_path=traits_path,
        centroids_path=centroids_path,
        trait="region",
        out_dir=temp_root / "complete-package",
        model="ard",
    )
    blocked_package = build_biogeography_report_package(
        tree_path=tree_path,
        traits_path=traits_path,
        centroids_path=incomplete_centroids_path,
        trait="region",
        out_dir=temp_root / "blocked-package",
        model="ard",
    )

    fixtures = [
        check(
            goal_id=230,
            suite="biogeography-figure-publication-reference",
            name="visible_state_probabilities_and_transitions",
            fixture_paths=[tree_path, traits_path, centroids_path],
            expected={
                "publication_ready": True,
                "node_probabilities_visible": True,
                "transitions_visible": True,
                "map_state_colors_complete": True,
                "rendered_internal_pie_count_matches": True,
            },
            observed={
                "publication_ready": complete_package.audit.publication_ready,
                "node_probabilities_visible": (
                    complete_package.audit.node_probabilities_visible
                ),
                "transitions_visible": complete_package.audit.transitions_visible,
                "map_state_colors_complete": (
                    complete_package.audit.map_state_colors_complete
                ),
                "rendered_internal_pie_count_matches": (
                    complete_package.audit.rendered_internal_pie_count
                    == complete_package.audit.expected_internal_node_count
                ),
            },
            notes=[
                "the governed biogeography fixture must keep internal pies, probability labels, and visible transition lines on the publication surfaces"
            ],
        ),
        check(
            goal_id=230,
            suite="biogeography-figure-publication-reference",
            name="missing_centroid_blocks_publication",
            fixture_paths=[tree_path, traits_path, incomplete_centroids_path],
            expected={
                "publication_ready": False,
                "map_state_colors_complete": False,
                "has_exclusion_limitation": True,
            },
            observed={
                "publication_ready": blocked_package.audit.publication_ready,
                "map_state_colors_complete": (
                    blocked_package.audit.map_state_colors_complete
                ),
                "has_exclusion_limitation": any(
                    "excluded" in limitation.lower()
                    for limitation in blocked_package.audit.limitations
                ),
            },
            notes=[
                "publication readiness stays blocked when one inferred region cannot be rendered on the geographic map"
            ],
        ),
    ]
    return suite_report(
        goal_id=230,
        suite="biogeography-figure-publication-reference",
        reviewer_goal="Biogeography figure publication fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins one governed biogeography figure package where the tree carries explicit node-probability pies and the map keeps shared state colors plus visible geographic transitions",
            "proves that incomplete centroid coverage blocks publication readiness instead of silently degrading the map figure",
        ],
        limitations=[
            "biogeography publication readiness is governed through shared color ledgers, render counts, and exclusion audits rather than screenshot goldens",
        ],
    )
