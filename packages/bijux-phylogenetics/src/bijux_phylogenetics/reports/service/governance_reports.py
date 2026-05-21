from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.validation import (
    build_core_workflow_validation_report,
    build_level_one_release_gate_report,
    build_production_scale_readiness_report,
    build_release_truth_report,
)

from .artifacts import report_sidecar_path, section, write_machine_manifest
from .ledger import sha256
from .models import (
    ProductionScaleReadinessReportBuildResult,
    ReleaseGateReportBuildResult,
    ReleaseTruthReportBuildResult,
    WorkflowValidationReportBuildResult,
)


def render_workflow_validation_report(
    *,
    out_path: Path,
    fixtures_root: Path | None = None,
) -> WorkflowValidationReportBuildResult:
    """Render the Level 1 workflow validation fixture report."""
    validation = build_core_workflow_validation_report(fixtures_root=fixtures_root)
    title = "Bijux Core Workflow Validation Report"
    reviewer_summary = [
        f"fixture checks passed: {validation.passed_fixture_count}/{validation.total_fixture_count}",
        f"validated workflow surfaces: {len(validation.workflows)}",
        f"known failure-gallery cases: {len(validation.failure_gallery)}",
    ]
    sections = [
        section("reviewer-summary", reviewer_summary),
        section(
            "validation-overview",
            {
                "total_fixture_count": validation.total_fixture_count,
                "passed_fixture_count": validation.passed_fixture_count,
                "failed_fixture_count": validation.failed_fixture_count,
            },
        ),
        section("validation-suites", [asdict(suite) for suite in validation.suites]),
        section("workflow-coverage", [asdict(row) for row in validation.workflows]),
        section("failure-gallery", [asdict(row) for row in validation.failure_gallery]),
        section(
            "maturity-classification",
            [asdict(row) for row in validation.maturity_classifications],
        ),
        section("limitations", validation.limitations),
    ]
    fixture_paths = sorted(
        {
            path
            for suite in validation.suites
            for fixture in suite.fixtures
            for path in fixture.fixture_paths
        }
    )
    machine_manifest = {
        "report_kind": "workflow-validation",
        "title": title,
        "input_paths": [str(path) for path in fixture_paths],
        "input_checksums": {
            str(path): sha256(path) for path in fixture_paths if path.exists()
        },
        "sections": [name for name, _ in sections],
        "metrics": {
            "total_fixture_count": validation.total_fixture_count,
            "passed_fixture_count": validation.passed_fixture_count,
            "workflow_count": len(validation.workflows),
        },
        "reviewer_summary": reviewer_summary,
        "limitations": validation.limitations,
    }
    machine_manifest_path = write_machine_manifest(
        report_sidecar_path(out_path), machine_manifest
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return WorkflowValidationReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="workflow-validation",
        title=title,
        validation=validation,
        machine_manifest=machine_manifest,
    )


def render_level_one_release_gate_report(
    *,
    out_path: Path,
    fixtures_root: Path | None = None,
) -> ReleaseGateReportBuildResult:
    """Render the Level 1 release gate for the checked-in workflow fixtures."""
    release_gate = build_level_one_release_gate_report(fixtures_root=fixtures_root)
    title = "Bijux Level 1 Release Gate Report"
    reviewer_summary = [
        f"gate decision: {release_gate.gate.decision}",
        f"dataset readiness: {release_gate.dataset_readiness_decision}",
        f"retained taxa: {len(release_gate.gate.retained_taxa)}, excluded taxa: {len(release_gate.gate.excluded_taxa)}",
    ]
    sections = [
        section("reviewer-summary", reviewer_summary),
        section("gate-decision", asdict(release_gate.gate)),
        section(
            "dataset-readiness",
            {
                "decision": release_gate.dataset_readiness_decision,
                "blockers": release_gate.dataset_blockers,
                "warnings": release_gate.dataset_warnings,
            },
        ),
        section(
            "taxon-loss-traceability",
            {
                "first_loss_stage": release_gate.taxon_first_loss_stage,
                "exclusion_causes": release_gate.exclusion_causes,
            },
        ),
        section("workflow-validation", asdict(release_gate.validation)),
        section("limitations", release_gate.validation.limitations),
    ]
    fixture_paths = sorted(
        {
            path
            for suite in release_gate.validation.suites
            for fixture in suite.fixtures
            for path in fixture.fixture_paths
        }
    )
    machine_manifest = {
        "report_kind": "release-gate",
        "title": title,
        "input_paths": [str(path) for path in fixture_paths],
        "input_checksums": {
            str(path): sha256(path) for path in fixture_paths if path.exists()
        },
        "sections": [name for name, _ in sections],
        "metrics": {
            "retained_taxa": len(release_gate.gate.retained_taxa),
            "excluded_taxa": len(release_gate.gate.excluded_taxa),
            "blocked_analysis_count": len(release_gate.gate.blocked_analyses),
        },
        "reviewer_summary": reviewer_summary,
        "limitations": release_gate.validation.limitations,
    }
    machine_manifest_path = write_machine_manifest(
        report_sidecar_path(out_path), machine_manifest
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return ReleaseGateReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="release-gate",
        title=title,
        release_gate=release_gate,
        machine_manifest=machine_manifest,
    )


def render_release_truth_report(
    *,
    out_path: Path,
    test_report_paths: list[Path],
    real_engine_test_report_paths: list[Path],
    fixtures_root: Path | None = None,
    include_extended_parity: bool = False,
    stress_tier: str = "small",
) -> ReleaseTruthReportBuildResult:
    """Render one machine-produced report of the current release truth surface."""
    release_truth = build_release_truth_report(
        test_report_paths=test_report_paths,
        real_engine_test_report_paths=real_engine_test_report_paths,
        fixtures_root=fixtures_root,
        include_extended_parity=include_extended_parity,
        stress_tier=stress_tier,
    )
    title = "Bijux Release Truth Report"
    reviewer_summary = [
        f"total tests: {release_truth.total_tests.passed_tests} passed, {release_truth.total_tests.failed_tests} failed, {release_truth.total_tests.skipped_tests} skipped",
        f"real-engine tests: {release_truth.real_engine_tests.passed_tests} passed, {release_truth.real_engine_tests.failed_tests} failed, {release_truth.real_engine_tests.skipped_tests} skipped",
        f"supported workflows: {len(release_truth.supported_workflows)}, experimental workflows: {len(release_truth.experimental_workflows)}",
        f"flagship datasets: {len(release_truth.flagship_datasets)}, reference parity cases: {release_truth.reference_parity.case_count}, stress workloads: {len(release_truth.stress_suite.observations)}",
        f"release gate decision: {release_truth.release_gate.gate.decision}",
    ]
    sections = [
        section("reviewer-summary", reviewer_summary),
        section("total-tests", asdict(release_truth.total_tests)),
        section("real-engine-tests", asdict(release_truth.real_engine_tests)),
        section(
            "supported-workflows",
            [asdict(item) for item in release_truth.supported_workflows],
        ),
        section(
            "experimental-workflows",
            [asdict(item) for item in release_truth.experimental_workflows],
        ),
        section(
            "advisory-workflows",
            [asdict(item) for item in release_truth.advisory_workflows],
        ),
        section(
            "parser-only-workflows",
            [asdict(item) for item in release_truth.parser_only_workflows],
        ),
        section(
            "flagship-datasets",
            [asdict(item) for item in release_truth.flagship_datasets],
        ),
        section("workflow-validation", asdict(release_truth.workflow_validation)),
        section("release-gate", asdict(release_truth.release_gate)),
        section("reference-parity", asdict(release_truth.reference_parity)),
        section("stress-suite", asdict(release_truth.stress_suite)),
        section("known-limitations", release_truth.known_limitations),
    ]
    fixture_paths = sorted(
        {
            path
            for suite in release_truth.workflow_validation.suites
            for fixture in suite.fixtures
            for path in fixture.fixture_paths
        }
    )
    input_paths = [*test_report_paths, *real_engine_test_report_paths, *fixture_paths]
    machine_manifest = {
        "report_kind": "release-truth",
        "title": title,
        "input_paths": [str(path) for path in input_paths],
        "input_checksums": {
            str(path): sha256(path) for path in input_paths if path.exists()
        },
        "sections": [name for name, _ in sections],
        "metrics": {
            "total_tests": release_truth.total_tests.total_tests,
            "total_tests_passed": release_truth.total_tests.passed_tests,
            "total_tests_failed": release_truth.total_tests.failed_tests,
            "total_tests_skipped": release_truth.total_tests.skipped_tests,
            "real_engine_tests": release_truth.real_engine_tests.total_tests,
            "real_engine_tests_passed": release_truth.real_engine_tests.passed_tests,
            "real_engine_tests_failed": release_truth.real_engine_tests.failed_tests,
            "real_engine_tests_skipped": release_truth.real_engine_tests.skipped_tests,
            "supported_workflow_count": len(release_truth.supported_workflows),
            "experimental_workflow_count": len(release_truth.experimental_workflows),
            "flagship_dataset_count": len(release_truth.flagship_datasets),
            "reference_parity_case_count": release_truth.reference_parity.case_count,
            "stress_workload_count": len(release_truth.stress_suite.observations),
        },
        "reviewer_summary": reviewer_summary,
        "limitations": release_truth.known_limitations,
    }
    machine_manifest_path = write_machine_manifest(
        report_sidecar_path(out_path), machine_manifest
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return ReleaseTruthReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="release-truth",
        title=title,
        release_truth=release_truth,
        machine_manifest=machine_manifest,
    )


def render_production_scale_readiness_report(
    *,
    out_path: Path,
    replicates: int = 1,
    tree_tip_counts: list[int] | None = None,
    alignment_size_classes: list[tuple[str, int, int]] | None = None,
    tree_set_size_classes: list[tuple[str, int, int]] | None = None,
    stress_tiers: list[str] | None = None,
) -> ProductionScaleReadinessReportBuildResult:
    """Render one reviewer-facing production-scale readiness report from governed benchmark evidence."""
    production_scale_readiness = build_production_scale_readiness_report(
        replicates=replicates,
        tree_tip_counts=tree_tip_counts,
        alignment_size_classes=alignment_size_classes,
        tree_set_size_classes=tree_set_size_classes,
        stress_tiers=stress_tiers,
    )
    title = "Bijux Production-Scale Readiness Report"
    highest_ready_scale_counts = {
        scale: sum(
            1
            for entry in production_scale_readiness.entries
            if entry.highest_ready_scale == scale
        )
        for scale in sorted(
            {
                "below-small",
                *(
                    threshold.scale
                    for threshold in production_scale_readiness.scale_definitions
                ),
            }
        )
    }
    scale_coverage = [
        {
            "scale": threshold.scale,
            "description": threshold.description,
            "minimum_taxa": threshold.minimum_taxa,
            "minimum_sites": threshold.minimum_sites,
            "minimum_tree_count": threshold.minimum_tree_count,
            "minimum_posterior_size": threshold.minimum_posterior_size,
            "ready_workflow_count": sum(
                1
                for entry in production_scale_readiness.entries
                for decision in entry.scale_decisions
                if decision.scale == threshold.scale and decision.ready
            ),
            "ready_workflows": sorted(
                entry.workflow
                for entry in production_scale_readiness.entries
                for decision in entry.scale_decisions
                if decision.scale == threshold.scale and decision.ready
            ),
        }
        for threshold in production_scale_readiness.scale_definitions
    ]
    reviewer_summary = [
        f"workflow count: {len(production_scale_readiness.entries)}",
        "highest ready scale distribution: "
        + ", ".join(
            f"{scale}={count}"
            for scale, count in highest_ready_scale_counts.items()
            if count > 0
        ),
        f"stress tiers: {', '.join(production_scale_readiness.stress_tiers)}",
    ]
    sections = [
        section("reviewer-summary", reviewer_summary),
        section(
            "scale-definitions",
            [asdict(item) for item in production_scale_readiness.scale_definitions],
        ),
        section("scale-coverage", scale_coverage),
        section(
            "production-scale-readiness",
            [asdict(item) for item in production_scale_readiness.entries],
        ),
        section("known-limitations", production_scale_readiness.limitations),
    ]
    machine_manifest = {
        "report_kind": "production-scale-readiness",
        "title": title,
        "input_paths": [],
        "input_checksums": {},
        "sections": [name for name, _ in sections],
        "metrics": {
            "goal_id": production_scale_readiness.goal_id,
            "workflow_count": len(production_scale_readiness.entries),
            "replicates": production_scale_readiness.replicates,
            "stress_tier_count": len(production_scale_readiness.stress_tiers),
            "scale_definition_count": len(production_scale_readiness.scale_definitions),
            "below_small_workflow_count": highest_ready_scale_counts.get(
                "below-small", 0
            ),
            **{
                f"{threshold.scale}_ready_workflow_count": sum(
                    1
                    for entry in production_scale_readiness.entries
                    for decision in entry.scale_decisions
                    if decision.scale == threshold.scale and decision.ready
                )
                for threshold in production_scale_readiness.scale_definitions
            },
        },
        "reviewer_summary": reviewer_summary,
        "limitations": production_scale_readiness.limitations,
    }
    machine_manifest_path = write_machine_manifest(
        report_sidecar_path(out_path),
        machine_manifest,
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return ProductionScaleReadinessReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="production-scale-readiness",
        title=title,
        production_scale_readiness=production_scale_readiness,
        machine_manifest=machine_manifest,
    )
