from __future__ import annotations

from .._shared import (
    DefusedXmlException,
    Path,
    SafeXmlET,
    XmlElement,
    _beast_artifact_error,
)
from ..models import (
    BeastAnalysisXmlIssue,
    BeastAnalysisXmlLogger,
    BeastAnalysisXmlReport,
)


def summarize_beast_analysis_xml(path: Path) -> BeastAnalysisXmlReport:
    """Summarize one prepared BEAST analysis XML into reviewer-facing assumptions."""
    if not path.exists():
        raise _beast_artifact_error(
            f"BEAST analysis XML file was not found: {path}",
            code="beast_xml_missing_file",
            path=path,
            artifact_kind="beast-analysis-xml",
        )
    issues: list[BeastAnalysisXmlIssue] = []
    try:
        root = SafeXmlET.parse(path).getroot()
    except (SafeXmlET.ParseError, DefusedXmlException) as error:
        issues.append(
            BeastAnalysisXmlIssue(
                code="invalid-xml",
                message=f"BEAST analysis XML is not well formed: {error}",
            )
        )
        return BeastAnalysisXmlReport(
            path=path,
            beast_version=None,
            beast_namespace=None,
            taxon_count=0,
            character_count=0,
            beast_data_type=None,
            substitution_model=None,
            clock_model=None,
            tree_prior=None,
            starting_tree_source=None,
            chain_length=None,
            state_node_count=0,
            logger_count=0,
            posterior_log_path=None,
            posterior_tree_path=None,
            calibration_count=0,
            calibration_ids=[],
            tip_date_count=0,
            tip_date_units=None,
            tip_date_direction=None,
            issues=issues,
            valid=False,
        )

    if root.tag != "beast":
        issues.append(
            BeastAnalysisXmlIssue(
                code="unexpected-root",
                message="BEAST analysis XML must use a top-level <beast> element",
            )
        )

    alignment = root.find("./data[@id='alignment']")
    if alignment is None:
        issues.append(
            BeastAnalysisXmlIssue(
                code="missing-alignment",
                message="BEAST analysis XML must define one alignment data block",
            )
        )
        taxon_count = 0
        character_count = 0
        beast_data_type = None
    else:
        sequences = alignment.findall("./sequence")
        taxon_count = len(sequences)
        sequence_lengths = {len(sequence.text or "") for sequence in sequences}
        character_count = 0 if not sequence_lengths else max(sequence_lengths)
        if len(sequence_lengths) > 1:
            issues.append(
                BeastAnalysisXmlIssue(
                    code="ragged-alignment",
                    message="BEAST analysis XML alignment sequences must share one common length",
                )
            )
        beast_data_type = alignment.get("dataType")
        if beast_data_type is None:
            issues.append(
                BeastAnalysisXmlIssue(
                    code="missing-data-type",
                    message="BEAST analysis XML alignment must declare its dataType",
                )
            )

    run = root.find("./run[@id='mcmc']")
    if run is None:
        issues.append(
            BeastAnalysisXmlIssue(
                code="missing-run",
                message="BEAST analysis XML must define one MCMC run with id 'mcmc'",
            )
        )
        chain_length = None
        state_node_count = 0
        loggers: list[BeastAnalysisXmlLogger] = []
    else:
        chain_length = _safe_int_attribute(
            run,
            "chainLength",
            issues=issues,
            issue_code="missing-chain-length",
            issue_message="BEAST analysis XML run must declare a numeric chainLength",
        )
        state_node_count = len(run.findall("./state/stateNode"))
        if state_node_count == 0:
            issues.append(
                BeastAnalysisXmlIssue(
                    code="missing-state-nodes",
                    message="BEAST analysis XML run must declare at least one stateNode",
                )
            )
        loggers = _collect_beast_analysis_xml_loggers(run)
        if not loggers:
            issues.append(
                BeastAnalysisXmlIssue(
                    code="missing-loggers",
                    message="BEAST analysis XML run must declare posterior logging outputs",
                )
            )

    substitution_model = _summarize_beast_xml_substitution_model(root)
    if substitution_model is None:
        issues.append(
            BeastAnalysisXmlIssue(
                code="missing-substitution-model",
                message="BEAST analysis XML must declare a substitution model in the site model",
            )
        )
    clock_model = _summarize_beast_xml_clock_model(root)
    if clock_model is None:
        issues.append(
            BeastAnalysisXmlIssue(
                code="missing-clock-model",
                message="BEAST analysis XML must declare one branch-rate model",
            )
        )
    tree_prior = _summarize_beast_xml_tree_prior(root)
    if tree_prior is None:
        issues.append(
            BeastAnalysisXmlIssue(
                code="missing-tree-prior",
                message="BEAST analysis XML must declare one tree prior",
            )
        )
    starting_tree_source = _summarize_beast_xml_starting_tree_source(root)
    if starting_tree_source is None:
        issues.append(
            BeastAnalysisXmlIssue(
                code="missing-starting-tree",
                message="BEAST analysis XML must declare a starting tree source",
            )
        )

    posterior_log_path = _beast_xml_logged_output_path(
        loggers, logger_kind="posterior-log"
    )
    posterior_tree_path = _beast_xml_logged_output_path(
        loggers, logger_kind="posterior-trees"
    )
    if posterior_log_path is None:
        issues.append(
            BeastAnalysisXmlIssue(
                code="missing-posterior-log",
                message="BEAST analysis XML must define one posterior parameter log output",
            )
        )
    if posterior_tree_path is None:
        issues.append(
            BeastAnalysisXmlIssue(
                code="missing-posterior-trees",
                message="BEAST analysis XML must define one posterior tree log output",
            )
        )

    calibration_distributions = root.findall(
        ".//distribution[@spec='beast.base.evolution.tree.MRCAPrior']"
    )
    calibration_ids = [
        distribution.get("id", "")
        for distribution in calibration_distributions
        if distribution.get("id")
    ]
    tip_trait = root.find(".//*[@traitname='date-forward']")
    tip_date_units = None if tip_trait is None else tip_trait.get("units")
    tip_date_direction = None if tip_trait is None else tip_trait.get("traitname")
    tip_date_count = 0
    if tip_trait is not None:
        tip_value = tip_trait.get("value", "")
        tip_date_count = len([part for part in tip_value.split(",") if part.strip()])

    return BeastAnalysisXmlReport(
        path=path,
        beast_version=root.get("version"),
        beast_namespace=root.get("namespace"),
        taxon_count=taxon_count,
        character_count=character_count,
        beast_data_type=beast_data_type,
        substitution_model=substitution_model,
        clock_model=clock_model,
        tree_prior=tree_prior,
        starting_tree_source=starting_tree_source,
        chain_length=chain_length,
        state_node_count=state_node_count,
        logger_count=len(loggers),
        posterior_log_path=posterior_log_path,
        posterior_tree_path=posterior_tree_path,
        calibration_count=len(calibration_ids),
        calibration_ids=calibration_ids,
        tip_date_count=tip_date_count,
        tip_date_units=tip_date_units,
        tip_date_direction=tip_date_direction,
        issues=issues,
        valid=not issues,
    )


def validate_beast_analysis_xml(path: Path) -> BeastAnalysisXmlReport:
    """Validate that a prepared BEAST analysis XML is structurally complete."""
    return summarize_beast_analysis_xml(path)


def _collect_beast_analysis_xml_loggers(
    run: XmlElement,
) -> list[BeastAnalysisXmlLogger]:
    loggers: list[BeastAnalysisXmlLogger] = []
    for logger in run.findall("./logger"):
        file_name = logger.get("fileName")
        log_every = None
        raw_log_every = logger.get("logEvery")
        if raw_log_every is not None:
            try:
                log_every = int(raw_log_every)
            except ValueError:
                log_every = None
        loggers.append(
            BeastAnalysisXmlLogger(
                logger_kind=_classify_beast_analysis_xml_logger(logger),
                file_name=file_name,
                log_every=log_every,
            )
        )
    return loggers


def _classify_beast_analysis_xml_logger(logger: XmlElement) -> str:
    file_name = logger.get("fileName")
    has_tree_log = logger.find("./log[@idref='tree']") is not None
    has_posterior_log = logger.find("./log[@idref='posterior']") is not None
    if file_name is None:
        return "screen"
    if has_tree_log and file_name.endswith(".trees"):
        return "posterior-trees"
    if has_posterior_log:
        return "posterior-log"
    return "other-file"


def _beast_xml_logged_output_path(
    loggers: list[BeastAnalysisXmlLogger], *, logger_kind: str
) -> Path | None:
    for logger in loggers:
        if logger.logger_kind == logger_kind and logger.file_name is not None:
            return Path(logger.file_name)
    return None


def _summarize_beast_xml_substitution_model(root: XmlElement) -> str | None:
    substitution_model = root.find("./input[@id='siteModel']/substModel")
    if substitution_model is None:
        return None
    spec = substitution_model.get("spec")
    if spec:
        return spec.split(".")[-1]
    id_ref = substitution_model.get("idref")
    if id_ref == "hky":
        return "HKY"
    if id_ref is not None:
        return id_ref
    return None


def _summarize_beast_xml_clock_model(root: XmlElement) -> str | None:
    branch_rates = root.find("./input[@id='branchRates']")
    if branch_rates is None:
        return None
    spec = branch_rates.get("spec")
    if spec is None:
        return None
    if spec.endswith("StrictClockModel"):
        return "strict"
    if spec.endswith("UCRelaxedClockModel"):
        return "relaxed-lognormal"
    return spec.split(".")[-1]


def _summarize_beast_xml_tree_prior(root: XmlElement) -> str | None:
    tree_prior = root.find("./input[@id='treePrior']")
    if tree_prior is None:
        return None
    spec = tree_prior.get("spec")
    if spec is None:
        return None
    if spec.endswith("YuleModel"):
        return "yule"
    if spec.endswith("BirthDeathGernhard08Model"):
        return "birth-death"
    return spec.split(".")[-1]


def _summarize_beast_xml_starting_tree_source(root: XmlElement) -> str | None:
    tree = root.find("./tree[@id='tree']")
    if tree is not None:
        return "provided-tree"
    cluster_tree = root.find("./input[@id='tree']")
    if cluster_tree is None:
        return None
    if cluster_tree.get("spec", "").endswith("ClusterTree"):
        return cluster_tree.get("clusterType", "cluster")
    return cluster_tree.get("spec")


def _safe_int_attribute(
    element: XmlElement,
    attribute: str,
    *,
    issues: list[BeastAnalysisXmlIssue],
    issue_code: str,
    issue_message: str,
) -> int | None:
    raw = element.get(attribute)
    if raw is None:
        issues.append(BeastAnalysisXmlIssue(code=issue_code, message=issue_message))
        return None
    try:
        return int(raw)
    except ValueError:
        issues.append(BeastAnalysisXmlIssue(code=issue_code, message=issue_message))
        return None
