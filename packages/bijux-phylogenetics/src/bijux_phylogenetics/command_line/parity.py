from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.parity import (
    run_ape_parity_cases,
    run_geiger_parity_cases,
    run_phytools_parity_cases,
    validate_reference_parity_examples,
    write_ape_parity_observation_table,
    write_ape_parity_summary_table,
    write_geiger_boundary_warning_table,
    write_geiger_likelihood_policy_table,
    write_geiger_model_confidence_table,
    write_geiger_optimizer_triage_table,
    write_geiger_parameterization_registry_table,
    write_geiger_parity_observation_table,
    write_geiger_parity_summary_table,
    write_generated_geiger_parity_report_json,
    write_generated_geiger_parity_report_markdown,
    write_phytools_parity_observation_table,
    write_phytools_parity_summary_table,
    write_reference_parity_observation_table,
    write_reference_parity_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_parity_command(subparsers: Any) -> None:
    parity = subparsers.add_parser(
        get_command_spec("parity").name, help=get_command_spec("parity").summary
    )
    parity.add_argument(
        "--reference-source",
        choices=("checked-fixture", "ape-live", "phytools-live", "geiger-live"),
        default="checked-fixture",
        help="Choose the checked fixture parity suite or one of the live external execution harnesses.",
    )
    parity.add_argument(
        "--extended",
        action="store_true",
        help="Include the optional larger posterior tree-set parity fixtures.",
    )
    parity.add_argument(
        "--ape-case",
        action="append",
        dest="ape_cases",
        help="Restrict the live ape parity harness to one or more governed case ids.",
    )
    parity.add_argument(
        "--ape-rscript-executable",
        default="Rscript",
        help="Executable used to launch the live ape parity runner.",
    )
    parity.add_argument(
        "--ape-failure-root",
        type=Path,
        help="Directory for reproducible live ape mismatch and skip artifacts.",
    )
    parity.add_argument(
        "--phytools-case",
        action="append",
        dest="phytools_cases",
        help="Restrict the live phytools parity harness to one or more governed case ids.",
    )
    parity.add_argument(
        "--phytools-rscript-executable",
        default="Rscript",
        help="Executable used to launch the live phytools parity runner.",
    )
    parity.add_argument(
        "--phytools-failure-root",
        type=Path,
        help="Directory for reproducible live phytools mismatch and skip artifacts.",
    )
    parity.add_argument(
        "--geiger-case",
        action="append",
        dest="geiger_cases",
        help="Restrict the live geiger parity harness to one or more governed case ids.",
    )
    parity.add_argument(
        "--geiger-rscript-executable",
        default="Rscript",
        help="Executable used to launch the live geiger parity runner.",
    )
    parity.add_argument(
        "--geiger-failure-root",
        type=Path,
        help="Directory for reproducible live geiger mismatch and skip artifacts.",
    )
    parity.add_argument("--summary-out", type=Path)
    parity.add_argument("--observations-out", type=Path)
    parity.add_argument("--optimizer-triage-out", type=Path)
    parity.add_argument("--boundary-warning-out", type=Path)
    parity.add_argument("--likelihood-policy-out", type=Path)
    parity.add_argument("--model-confidence-out", type=Path)
    parity.add_argument("--parameterization-registry-out", type=Path)
    parity.add_argument("--generated-report-out", type=Path)
    parity.add_argument("--generated-report-json-out", type=Path)
    parity.add_argument(
        "--json", action="store_true", help="Emit the parity report as JSON."
    )
    _add_manifest_argument(parity)


def _command_line_api() -> Any:
    import bijux_phylogenetics.command_line as command_line_api

    return command_line_api


def run_parity_command(args: Any) -> int:
    command_line_api = _command_line_api()
    if (
        args.reference_source in {"ape-live", "phytools-live", "geiger-live"}
        and args.extended
    ):
        raise ValueError(
            "--extended is only supported for the checked fixture parity suite"
        )
    if args.reference_source == "ape-live":
        report = run_ape_parity_cases(
            case_ids=args.ape_cases,
            rscript_executable=args.ape_rscript_executable,
            failure_root=args.ape_failure_root,
        )
        output_paths: list[Path | str] = []
        summary_path = None
        observation_path = None
        if args.summary_out is not None:
            summary_path = write_ape_parity_summary_table(
                args.summary_out,
                report,
            )
            output_paths.append(summary_path)
        if args.observations_out is not None:
            observation_path = write_ape_parity_observation_table(
                args.observations_out,
                report,
            )
            output_paths.append(observation_path)
        outputs = _finalize_outputs(
            args,
            command="parity",
            inputs=[],
            outputs=output_paths,
        )
        _print_result(
            build_command_result(
                command="parity",
                inputs=[],
                outputs=outputs,
                metrics={
                    "all_passed": report.all_passed,
                    "case_count": report.case_count,
                    "function_count": len(report.summary_rows),
                    "failed_case_count": report.failed_case_count,
                    "skipped_case_count": report.skipped_case_count,
                    "reference_source": args.reference_source,
                },
                data={
                    "report": report,
                    "summary_table": summary_path,
                    "observation_table": observation_path,
                },
            ),
            json_output=args.json,
        )
        return 0
    if args.reference_source == "phytools-live":
        report = run_phytools_parity_cases(
            case_ids=args.phytools_cases,
            rscript_executable=args.phytools_rscript_executable,
            failure_root=args.phytools_failure_root,
        )
        output_paths: list[Path | str] = []
        summary_path = None
        observation_path = None
        if args.summary_out is not None:
            summary_path = write_phytools_parity_summary_table(
                args.summary_out,
                report,
            )
            output_paths.append(summary_path)
        if args.observations_out is not None:
            observation_path = write_phytools_parity_observation_table(
                args.observations_out,
                report,
            )
            output_paths.append(observation_path)
        outputs = _finalize_outputs(
            args,
            command="parity",
            inputs=[],
            outputs=output_paths,
        )
        _print_result(
            build_command_result(
                command="parity",
                inputs=[],
                outputs=outputs,
                metrics={
                    "all_passed": report.all_passed,
                    "case_count": report.case_count,
                    "function_count": len(report.summary_rows),
                    "failed_case_count": report.failed_case_count,
                    "skipped_case_count": report.skipped_case_count,
                    "reference_source": args.reference_source,
                },
                data={
                    "report": report,
                    "summary_table": summary_path,
                    "observation_table": observation_path,
                },
            ),
            json_output=args.json,
        )
        return 0
    if args.reference_source == "geiger-live":
        report = run_geiger_parity_cases(
            case_ids=args.geiger_cases,
            rscript_executable=args.geiger_rscript_executable,
            failure_root=args.geiger_failure_root,
        )
        output_paths: list[Path | str] = []
        summary_path = None
        observation_path = None
        optimizer_triage_path = None
        boundary_warning_path = None
        likelihood_policy_path = None
        model_confidence_path = None
        parameterization_registry_path = None
        generated_report_path = None
        generated_report_json_path = None
        if args.summary_out is not None:
            summary_path = write_geiger_parity_summary_table(
                args.summary_out,
                report,
            )
            output_paths.append(summary_path)
        if args.observations_out is not None:
            observation_path = write_geiger_parity_observation_table(
                args.observations_out,
                report,
            )
            output_paths.append(observation_path)
        if args.optimizer_triage_out is not None:
            optimizer_triage_path = write_geiger_optimizer_triage_table(
                args.optimizer_triage_out,
                report,
            )
            output_paths.append(optimizer_triage_path)
        if args.boundary_warning_out is not None:
            boundary_warning_path = write_geiger_boundary_warning_table(
                args.boundary_warning_out,
                report,
            )
            output_paths.append(boundary_warning_path)
        if args.likelihood_policy_out is not None:
            likelihood_policy_path = write_geiger_likelihood_policy_table(
                args.likelihood_policy_out,
                report,
            )
            output_paths.append(likelihood_policy_path)
        if args.model_confidence_out is not None:
            model_confidence_path = write_geiger_model_confidence_table(
                args.model_confidence_out,
                report,
            )
            output_paths.append(model_confidence_path)
        if args.parameterization_registry_out is not None:
            parameterization_registry_path = (
                write_geiger_parameterization_registry_table(
                    args.parameterization_registry_out,
                    report,
                )
            )
            output_paths.append(parameterization_registry_path)
        generated_report = None
        if (
            args.generated_report_out is not None
            or args.generated_report_json_out is not None
        ):
            generated_report = command_line_api.build_generated_geiger_parity_report(
                parity_report=report,
            )
        if args.generated_report_out is not None:
            generated_report_path = write_generated_geiger_parity_report_markdown(
                args.generated_report_out,
                generated_report,
            )
            output_paths.append(generated_report_path)
        if args.generated_report_json_out is not None:
            generated_report_json_path = write_generated_geiger_parity_report_json(
                args.generated_report_json_out,
                generated_report,
            )
            output_paths.append(generated_report_json_path)
        outputs = _finalize_outputs(
            args,
            command="parity",
            inputs=[],
            outputs=output_paths,
        )
        _print_result(
            build_command_result(
                command="parity",
                inputs=[],
                outputs=outputs,
                metrics={
                    "all_passed": report.all_passed,
                    "case_count": report.case_count,
                    "function_count": len(report.summary_rows),
                    "failed_case_count": report.failed_case_count,
                    "skipped_case_count": report.skipped_case_count,
                    "boundary_warning_row_count": len(report.boundary_warning_rows),
                    "model_confidence_row_count": len(report.model_confidence_rows),
                    "generated_report_written": generated_report_path is not None
                    or generated_report_json_path is not None,
                    "reference_source": args.reference_source,
                },
                data={
                    "report": report,
                    "summary_table": summary_path,
                    "observation_table": observation_path,
                    "optimizer_triage_table": optimizer_triage_path,
                    "boundary_warning_table": boundary_warning_path,
                    "likelihood_policy_table": likelihood_policy_path,
                    "model_confidence_table": model_confidence_path,
                    "parameterization_registry_table": parameterization_registry_path,
                    "generated_report_markdown": generated_report_path,
                    "generated_report_json": generated_report_json_path,
                    "generated_report": generated_report,
                },
            ),
            json_output=args.json,
        )
        return 0

    report = validate_reference_parity_examples(include_extended=args.extended)
    output_paths: list[Path | str] = []
    summary_path = None
    observation_path = None
    if args.summary_out is not None:
        summary_path = write_reference_parity_summary_table(
            args.summary_out,
            report,
        )
        output_paths.append(summary_path)
    if args.observations_out is not None:
        observation_path = write_reference_parity_observation_table(
            args.observations_out,
            report,
        )
        output_paths.append(observation_path)
    outputs = _finalize_outputs(
        args,
        command="parity",
        inputs=[],
        outputs=output_paths,
    )
    _print_result(
        build_command_result(
            command="parity",
            inputs=[],
            outputs=outputs,
            metrics={
                "all_passed": report.all_passed,
                "case_count": report.case_count,
                "method_count": len(report.covered_methods),
                "failed_case_count": report.failed_case_count,
                "reference_source": args.reference_source,
                "extended": args.extended,
            },
            data={
                "report": report,
                "summary_table": summary_path,
                "observation_table": observation_path,
            },
        ),
        json_output=args.json,
    )
    return 0
