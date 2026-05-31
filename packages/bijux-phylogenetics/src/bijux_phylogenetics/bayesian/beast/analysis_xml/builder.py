from __future__ import annotations

from .._shared import (
    ElementTree,
    EngineWorkflowError,
    Path,
    _ensure_inference_ready_alignment,
    deepcopy,
    infer_alignment_alphabet,
    load_fasta_alignment,
)
from ..models import BeastCalibration, BeastPreparationReport
from ..validation import validate_fossil_calibration_table, validate_tip_dating_metadata
from .document_components import (
    _append_sequence_alignment,
    _append_starting_tree,
    _append_substitution_and_site_model,
    _append_tree_likelihood,
)
from .inference_models import (
    _append_clock_model,
    _append_tree_prior,
    _translate_calibration_distribution,
)
from .summary import validate_beast_analysis_xml
from .xml_primitives import (
    _beast_data_type,
    _default_beast_substitution_model,
    _validate_tree_taxa_against_alignment,
    _xml_element,
    _xml_identifier,
)


def prepare_beast_time_tree_analysis(
    alignment_path: Path,
    output_path: Path,
    *,
    tree_path: Path | None = None,
    calibration_path: Path | None = None,
    tip_dates_path: Path | None = None,
    clock_model: str = "strict",
    tree_prior: str = "yule",
    chain_length: int = 1000000,
    log_every: int = 1000,
    taxon_column: str | None = None,
    date_column: str = "date",
) -> BeastPreparationReport:
    """Prepare a deterministic BEAST2 XML configuration from alignment and dating inputs."""
    _ensure_inference_ready_alignment(alignment_path)
    records = load_fasta_alignment(alignment_path)
    inferred_alphabet = infer_alignment_alphabet(records)
    if calibration_path is not None and tree_path is None:
        raise ValueError(
            "BEAST preparation requires tree_path when calibration_path is provided"
        )
    if tip_dates_path is not None and tree_path is None:
        raise ValueError(
            "BEAST preparation requires tree_path when tip_dates_path is provided"
        )
    beast_data_type = _beast_data_type(inferred_alphabet)
    substitution_model = _default_beast_substitution_model(beast_data_type)

    alignment_taxa = {record.identifier for record in records}
    if tree_path is not None:
        _validate_tree_taxa_against_alignment(
            tree_path=tree_path,
            alignment_taxa=alignment_taxa,
        )

    calibration_report = (
        validate_fossil_calibration_table(tree_path, calibration_path)
        if calibration_path is not None and tree_path is not None
        else None
    )
    if calibration_report is not None and calibration_report.invalid_calibration_count:
        raise ValueError(
            "BEAST preparation requires all fossil calibrations to validate successfully"
        )
    tip_date_report = (
        validate_tip_dating_metadata(
            tree_path,
            tip_dates_path,
            alignment_path=alignment_path,
            taxon_column=taxon_column,
            date_column=date_column,
        )
        if tip_dates_path is not None and tree_path is not None
        else None
    )
    if tip_date_report is not None and tip_date_report.invalid_tip_count:
        raise ValueError(
            "BEAST preparation requires all tip dates to validate successfully"
        )

    root = ElementTree.Element(
        "beast",
        {
            "version": "2.7",
            "namespace": (
                "beast.pkgmgmt:beast.base.core:beast.base.inference:"
                "beast.base.evolution.alignment:beast.base.evolution.tree:"
                "beast.base.evolution.tree.coalescent:beast.base.evolution.speciation:"
                "beast.base.evolution.branchratemodel:beast.base.evolution.operator:"
                "beast.base.inference.operator:beast.base.evolution.sitemodel:"
                "beast.base.evolution.substitutionmodel:beast.base.evolution.likelihood:"
                "beast.base.inference.parameter:beast.base.inference.distribution:"
                "beast.base.math:beast.base.math.distributions:beast.base.util:"
                "beast.evolution:beast.evolution.nuc"
            ),
        },
    )

    _append_sequence_alignment(root, records=records, beast_data_type=beast_data_type)
    (
        site_state_ids,
        site_prior_elements,
        site_operator_elements,
        site_logger_elements,
    ) = _append_substitution_and_site_model(root, beast_data_type=beast_data_type)
    starting_tree_source = _append_starting_tree(
        root,
        tree_path=tree_path,
        tip_date_report=tip_date_report,
    )
    (
        clock_state_ids,
        clock_prior_elements,
        clock_operator_elements,
        clock_logger_elements,
    ) = _append_clock_model(root, clock_model=clock_model, taxon_count=len(records))
    (
        tree_state_ids,
        tree_prior_elements,
        tree_operator_elements,
        tree_logger_elements,
    ) = _append_tree_prior(root, tree_prior=tree_prior)
    _append_tree_likelihood(root)

    run = ElementTree.SubElement(
        root,
        "run",
        {
            "spec": "MCMC",
            "id": "mcmc",
            "chainLength": str(chain_length),
        },
    )
    state = ElementTree.SubElement(run, "state")
    for node_id in dict.fromkeys(
        [*site_state_ids, *clock_state_ids, *tree_state_ids, "tree"]
    ):
        ElementTree.SubElement(state, "stateNode", {"idref": node_id})

    posterior = ElementTree.SubElement(
        run, "distribution", {"spec": "CompoundDistribution", "id": "posterior"}
    )
    prior = ElementTree.SubElement(
        posterior, "distribution", {"spec": "CompoundDistribution", "id": "prior"}
    )
    for element in [*site_prior_elements, *clock_prior_elements, *tree_prior_elements]:
        prior.append(element)

    warnings: list[str] = []
    if tip_date_report is not None and tree_prior.strip().lower() == "birth-death":
        warnings.append(
            "tip-dated analyses with the standard birth-death tree prior are exploratory in this template; "
            "BEAST warns that this prior is not serial-sampling aware"
        )
    beast_calibrations: list[BeastCalibration] = []
    if calibration_report is not None:
        for calibration in calibration_report.calibrations:
            prior_element, beast_calibration, warning = (
                _translate_calibration_distribution(calibration)
            )
            prior.append(prior_element)
            beast_calibrations.append(beast_calibration)
            if warning is not None:
                warnings.append(warning)

    ElementTree.SubElement(posterior, "distribution", {"idref": "likelihood"})

    generic_tree_operators = [
        _xml_element(
            "operator",
            {
                "id": "treeScaler",
                "spec": "ScaleOperator",
                "scaleFactor": "0.75",
                "weight": "3",
                "tree": "@tree",
            },
        ),
        _xml_element("operator", {"spec": "Uniform", "weight": "30", "tree": "@tree"}),
        _xml_element(
            "operator",
            {
                "spec": "SubtreeSlide",
                "weight": "15",
                "gaussian": "true",
                "size": "0.5",
                "tree": "@tree",
            },
        ),
        _xml_element(
            "operator",
            {
                "id": "narrowExchange",
                "spec": "Exchange",
                "isNarrow": "true",
                "weight": "15",
                "tree": "@tree",
            },
        ),
        _xml_element(
            "operator",
            {
                "id": "wideExchange",
                "spec": "Exchange",
                "isNarrow": "false",
                "weight": "3",
                "tree": "@tree",
            },
        ),
        _xml_element(
            "operator",
            {
                "id": "wilsonBalding",
                "spec": "WilsonBalding",
                "weight": "3",
                "tree": "@tree",
            },
        ),
    ]
    for operator in [
        *site_operator_elements,
        *clock_operator_elements,
        *tree_operator_elements,
        *generic_tree_operators,
    ]:
        run.append(operator)

    log_path = output_path.with_name(f"{output_path.stem}.$(seed).log")
    tree_log_path = output_path.with_name(f"{output_path.stem}.$(seed).trees")
    file_logger = ElementTree.SubElement(
        run,
        "logger",
        {"logEvery": str(log_every), "fileName": log_path.name},
    )
    ElementTree.SubElement(file_logger, "model", {"idref": "posterior"})
    for log_id in ["posterior", "prior", "likelihood"]:
        ElementTree.SubElement(file_logger, "log", {"idref": log_id})
    ElementTree.SubElement(
        file_logger,
        "log",
        {"spec": "beast.base.evolution.tree.TreeHeightLogger", "tree": "@tree"},
    )
    for element in [
        *tree_logger_elements,
        *site_logger_elements,
        *clock_logger_elements,
    ]:
        file_logger.append(element)
    for calibration in beast_calibrations:
        ElementTree.SubElement(
            file_logger, "log", {"idref": _xml_identifier(calibration.calibration_id)}
        )

    tree_logger = ElementTree.SubElement(
        run,
        "logger",
        {"logEvery": str(log_every), "fileName": tree_log_path.name},
    )
    ElementTree.SubElement(tree_logger, "log", {"idref": "tree"})

    screen_logger = ElementTree.SubElement(
        run,
        "logger",
        {"logEvery": str(max(log_every, 10000))},
    )
    ElementTree.SubElement(screen_logger, "model", {"idref": "posterior"})
    for log_id in ["posterior", "prior", "likelihood"]:
        ElementTree.SubElement(screen_logger, "log", {"idref": log_id})
    ElementTree.SubElement(
        screen_logger,
        "log",
        {"spec": "beast.base.evolution.tree.TreeHeightLogger", "tree": "@tree"},
    )
    for element in [
        *tree_logger_elements,
        *site_logger_elements,
        *clock_logger_elements,
    ]:
        cloned = deepcopy(element)
        if cloned.get("id") == "rateStatistic":
            cloned.set("id", "rateStatistic.screen")
        screen_logger.append(cloned)

    xml_tree = ElementTree.ElementTree(root)
    ElementTree.indent(xml_tree, space="    ")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    xml_tree.write(output_path, encoding="utf-8", xml_declaration=True)
    output_path.write_text(
        output_path.read_text(encoding="utf-8") + "\n", encoding="utf-8"
    )
    xml_report = validate_beast_analysis_xml(output_path)
    if not xml_report.valid:
        messages = "; ".join(issue.message for issue in xml_report.issues)
        raise EngineWorkflowError(
            f"generated BEAST analysis XML failed structural validation: {messages}"
        )
    return BeastPreparationReport(
        alignment_path=alignment_path,
        output_path=output_path,
        tree_path=tree_path,
        calibration_path=calibration_path,
        tip_dates_path=tip_dates_path,
        taxon_count=len(records),
        character_count=len(records[0].sequence),
        inferred_alphabet=inferred_alphabet,
        beast_data_type=beast_data_type,
        substitution_model=substitution_model,
        clock_model=clock_model,
        tree_prior=tree_prior,
        starting_tree_source=starting_tree_source,
        chain_length=chain_length,
        log_every=log_every,
        calibration_count=0
        if calibration_report is None
        else calibration_report.calibration_count,
        tip_date_count=0
        if tip_date_report is None
        else tip_date_report.valid_tip_count,
        warning_count=len(warnings),
        warnings=warnings,
        log_path=log_path,
        tree_log_path=tree_log_path,
        calibrations=beast_calibrations,
    )
