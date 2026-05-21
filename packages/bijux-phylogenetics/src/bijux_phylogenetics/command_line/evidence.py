from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import (
    _evidence_book_metrics,
    _print_result,
)
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.evidence.book import validate_evidence_book
from bijux_phylogenetics.evidence.bundles import bundle_directory, validate_bundle
from bijux_phylogenetics.evidence.workbench import (
    DOCS_EVIDENCE_OVERVIEW,
    build_evidence_book_selection,
    build_evidence_book_study,
    list_registered_evidence_studies,
    refresh_evidence_book,
    rerun_evidence_book_selection,
)
from bijux_phylogenetics.runtime.errors import EvidenceContractError
from bijux_phylogenetics.runtime.results import build_command_result


def add_evidence_command(subparsers: Any) -> None:
    evidence = subparsers.add_parser(
        get_command_spec("evidence").name, help=get_command_spec("evidence").summary
    )
    evidence_subparsers = evidence.add_subparsers(
        dest="evidence_command", required=True
    )
    evidence_bundle = evidence_subparsers.add_parser(
        "bundle",
        help="Bundle explicit phylogenetics inputs and outputs as evidence.",
    )
    evidence_bundle.add_argument("--inputs", nargs="+", required=True, type=Path)
    evidence_bundle.add_argument("--outputs", nargs="+", required=True, type=Path)
    evidence_bundle.add_argument("--out", required=True, type=Path)
    evidence_bundle.add_argument(
        "--json", action="store_true", help="Emit the bundle report as JSON."
    )
    _add_manifest_argument(evidence_bundle)

    evidence_validate = evidence_subparsers.add_parser(
        "validate", help="Validate an existing evidence bundle."
    )
    evidence_validate.add_argument("bundle_root", type=Path)
    evidence_validate.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(evidence_validate)

    evidence_book = evidence_subparsers.add_parser(
        "book",
        help="Govern evidence-book generation, validation, and partial reruns.",
    )
    evidence_book_subparsers = evidence_book.add_subparsers(
        dest="evidence_book_command",
        required=True,
    )
    evidence_book_studies = evidence_book_subparsers.add_parser(
        "studies",
        help="List governed evidence-book studies and partial rerun capabilities.",
    )
    evidence_book_studies.add_argument(
        "--json", action="store_true", help="Emit the study registry as JSON."
    )
    _add_manifest_argument(evidence_book_studies)

    evidence_book_build = evidence_book_subparsers.add_parser(
        "build",
        help="Refresh governed evidence-book outputs or rebuild one registered study.",
    )
    evidence_book_build.add_argument(
        "study_id",
        nargs="?",
        help="Optional registered study identifier to rebuild before refreshing the evidence-book.",
    )
    evidence_book_build.add_argument(
        "--evidence-id",
        dest="evidence_ids",
        action="append",
        default=[],
        help="Optional Evidence ID to rebuild within the selected study. May be repeated.",
    )
    evidence_book_build.add_argument(
        "--json", action="store_true", help="Emit the build report as JSON."
    )
    _add_manifest_argument(evidence_book_build)

    evidence_book_validate = evidence_book_subparsers.add_parser(
        "validate",
        help="Validate the governed evidence-book surface and summarize coverage gaps.",
    )
    evidence_book_validate.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(evidence_book_validate)

    evidence_book_rerun = evidence_book_subparsers.add_parser(
        "rerun",
        help="Regenerate selected Evidence IDs for a study and refresh governed outputs.",
    )
    evidence_book_rerun.add_argument("study_id")
    evidence_book_rerun.add_argument("evidence_ids", nargs="+")
    evidence_book_rerun.add_argument(
        "--json", action="store_true", help="Emit the rerun report as JSON."
    )
    _add_manifest_argument(evidence_book_rerun)


def run_evidence_command(args: Any) -> int:
    if args.evidence_command == "bundle":
        report = bundle_directory(args.inputs, args.outputs, args.out)
        inputs = [*args.inputs, *args.outputs]
        outputs = _finalize_outputs(
            args, command="evidence", inputs=inputs, outputs=[args.out]
        )
        _print_result(
            build_command_result(
                command="evidence",
                inputs=inputs,
                outputs=outputs,
                metrics={
                    "file_count": report.file_count,
                    "input_file_count": report.input_file_count,
                    "output_file_count": report.output_file_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.evidence_command == "validate":
        report = validate_bundle(args.bundle_root)
        if not report.valid:
            raise EvidenceContractError(
                f"evidence bundle validation failed with {len(report.mismatches)} mismatch(es)"
            )
        outputs = _finalize_outputs(args, command="evidence", inputs=[args.bundle_root])
        _print_result(
            build_command_result(
                command="evidence",
                inputs=[args.bundle_root],
                outputs=outputs,
                metrics={"file_count": report.file_count},
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    repo_root = Path.cwd()
    if args.evidence_book_command == "studies":
        studies = list_registered_evidence_studies(repo_root)
        outputs = _finalize_outputs(args, command="evidence", inputs=[])
        _print_result(
            build_command_result(
                command="evidence",
                inputs=[],
                outputs=outputs,
                metrics={
                    "study_count": len(studies),
                    "partial_rerun_capable_count": sum(
                        1 for study in studies if study.supports_partial_rerun
                    ),
                },
                data={"studies": studies},
            ),
            json_output=args.json,
        )
        return 0

    if args.evidence_book_command == "build":
        if args.study_id is None and args.evidence_ids:
            raise EvidenceContractError(
                "--evidence-id requires a study_id for evidence book build"
            )
        if args.study_id is None:
            refresh_report = refresh_evidence_book(repo_root)
            outputs = _finalize_outputs(
                args,
                command="evidence",
                inputs=[],
                outputs=refresh_report.updated_paths,
            )
            metrics: dict[str, object] = {
                "reviewer_summary_count": refresh_report.reviewer_summary_count,
                "updated_path_count": len(refresh_report.updated_paths),
                **_evidence_book_metrics(repo_root),
            }
            _print_result(
                build_command_result(
                    command="evidence",
                    inputs=[],
                    outputs=outputs,
                    metrics=metrics,
                    data=refresh_report,
                ),
                json_output=args.json,
            )
            return 0
        if args.evidence_ids:
            report = build_evidence_book_selection(
                repo_root,
                args.study_id,
                args.evidence_ids,
            )
            outputs = _finalize_outputs(
                args,
                command="evidence",
                inputs=[],
                outputs=report.refresh_report.updated_paths,
            )
            metrics = {
                "selected_study_count": 1,
                "selected_evidence_count": len(report.selected_evidence_ids),
                "updated_path_count": len(report.refresh_report.updated_paths),
                "reviewer_summary_count": report.refresh_report.reviewer_summary_count,
                **_evidence_book_metrics(repo_root),
            }
            _print_result(
                build_command_result(
                    command="evidence",
                    inputs=[],
                    outputs=outputs,
                    metrics=metrics,
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        report = build_evidence_book_study(repo_root, args.study_id)
        build_inputs = (
            []
            if report.study_report.build_script_path is None
            else [Path(report.study_report.build_script_path)]
        )
        outputs = _finalize_outputs(
            args,
            command="evidence",
            inputs=build_inputs,
            outputs=report.refresh_report.updated_paths,
        )
        metrics = {
            "selected_study_count": 1,
            "updated_path_count": len(report.refresh_report.updated_paths),
            "reviewer_summary_count": report.refresh_report.reviewer_summary_count,
            **_evidence_book_metrics(repo_root),
        }
        _print_result(
            build_command_result(
                command="evidence",
                inputs=build_inputs,
                outputs=outputs,
                metrics=metrics,
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.evidence_book_command == "validate":
        report = validate_evidence_book(repo_root)
        if not report.valid:
            raise EvidenceContractError(
                f"evidence-book validation failed with {len(report.issues)} issue(s)"
            )
        outputs = _finalize_outputs(
            args,
            command="evidence",
            inputs=[repo_root / "evidence-book"],
            outputs=[
                repo_root / "evidence-book" / "index" / "coverage-gaps.json",
                repo_root / "evidence-book" / "index" / "freshness-report.json",
                repo_root / "evidence-book" / "index" / "integrity-report.json",
                repo_root / DOCS_EVIDENCE_OVERVIEW,
            ],
        )
        _print_result(
            build_command_result(
                command="evidence",
                inputs=[repo_root / "evidence-book"],
                outputs=outputs,
                metrics={
                    "issue_count": len(report.issues),
                    **_evidence_book_metrics(repo_root),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    report = rerun_evidence_book_selection(repo_root, args.study_id, args.evidence_ids)
    outputs = _finalize_outputs(
        args,
        command="evidence",
        inputs=[],
        outputs=report.refresh_report.updated_paths,
    )
    _print_result(
        build_command_result(
            command="evidence",
            inputs=[],
            outputs=outputs,
            metrics={
                "selected_evidence_count": len(
                    report.rerun_report.selected_evidence_ids
                ),
                "updated_path_count": len(report.refresh_report.updated_paths),
                **_evidence_book_metrics(repo_root),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
