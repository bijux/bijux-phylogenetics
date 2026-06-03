from __future__ import annotations

from .._shared import ElementTree, XmlElement, _format_decimal
from ..models import BeastCalibration, ValidatedCalibration
from .xml_primitives import _xml_element, _xml_identifier


def _append_clock_model(
    root: XmlElement,
    *,
    clock_model: str,
    taxon_count: int,
) -> tuple[list[str], list[XmlElement], list[XmlElement], list[XmlElement]]:
    normalized = clock_model.strip().lower()
    state_node_ids: list[str] = []
    prior_elements: list[XmlElement] = []
    operator_elements: list[XmlElement] = []
    logger_elements: list[XmlElement] = []
    if normalized == "strict":
        strict = ElementTree.SubElement(
            root,
            "input",
            {
                "spec": "beast.base.evolution.branchratemodel.StrictClockModel",
                "id": "branchRates",
            },
        )
        ElementTree.SubElement(strict, "clock.rate", {"idref": "clockRate"})
        ElementTree.SubElement(
            root,
            "parameter",
            {"id": "clockRate", "value": "0.001", "lower": "0.0", "upper": "100.0"},
        )
        state_node_ids.append("clockRate")
        prior_elements.append(
            _xml_element(
                "distribution",
                {
                    "id": "clockRate.prior",
                    "spec": "beast.base.inference.distribution.Prior",
                    "x": "@clockRate",
                },
                children=(
                    _xml_element(
                        "distr",
                        {
                            "id": "clockRate.uniform",
                            "lower": "0.0",
                            "upper": "100.0",
                            "spec": "beast.base.inference.distribution.Uniform",
                        },
                    ),
                ),
            )
        )
        operator_elements.append(
            _xml_element(
                "operator",
                {
                    "id": "clockRateScaler",
                    "spec": "ScaleOperator",
                    "scaleFactor": "0.75",
                    "weight": "3",
                    "parameter": "@clockRate",
                },
            )
        )
        logger_elements.append(_xml_element("log", {"idref": "clockRate"}))
        return state_node_ids, prior_elements, operator_elements, logger_elements
    if normalized != "relaxed-lognormal":
        raise ValueError(
            "BEAST preparation supports clock_model values 'strict' and 'relaxed-lognormal'"
        )
    relaxed = ElementTree.SubElement(
        root,
        "input",
        {
            "spec": "beast.base.evolution.branchratemodel.UCRelaxedClockModel",
            "id": "branchRates",
        },
    )
    distr = ElementTree.SubElement(
        relaxed,
        "distr",
        {
            "id": "ucld.lognormal",
            "spec": "beast.base.inference.distribution.LogNormalDistributionModel",
            "meanInRealSpace": "true",
        },
    )
    ElementTree.SubElement(
        distr,
        "parameter",
        {"name": "M", "id": "ucld.mean", "value": "1.0", "lower": "0.0"},
    )
    ElementTree.SubElement(
        distr,
        "parameter",
        {
            "name": "S",
            "id": "ucld.stdev",
            "value": "0.3333333333333333",
            "lower": "0.0",
        },
    )
    ElementTree.SubElement(
        relaxed,
        "parameter",
        {
            "spec": "IntegerParameter",
            "name": "rateCategories",
            "id": "rateCategories",
            "dimension": str(max(1, (2 * taxon_count) - 1)),
            "value": "1",
        },
    )
    ElementTree.SubElement(relaxed, "tree", {"idref": "tree"})
    state_node_ids.extend(["ucld.mean", "ucld.stdev", "rateCategories"])
    prior_elements.append(
        _xml_element(
            "distribution",
            {
                "id": "ucld.stdev.prior",
                "spec": "beast.base.inference.distribution.Prior",
                "x": "@ucld.stdev",
            },
            children=(
                _xml_element(
                    "distr",
                    {
                        "id": "ucld.stdev.exponential",
                        "mean": "0.3333333333333333",
                        "spec": "beast.base.inference.distribution.Exponential",
                    },
                ),
            ),
        )
    )
    operator_elements.extend(
        [
            _xml_element(
                "operator",
                {
                    "id": "ucldMeanScaler",
                    "spec": "ScaleOperator",
                    "scaleFactor": "0.75",
                    "weight": "1",
                    "parameter": "@ucld.mean",
                },
            ),
            _xml_element(
                "operator",
                {
                    "id": "ucldStdevScaler",
                    "spec": "ScaleOperator",
                    "scaleFactor": "0.75",
                    "weight": "3",
                    "parameter": "@ucld.stdev",
                },
            ),
            _xml_element(
                "operator",
                {
                    "id": "rateCategoriesRandomWalk",
                    "spec": "IntRandomWalkOperator",
                    "windowSize": "1",
                    "weight": "10",
                    "parameter": "@rateCategories",
                },
            ),
            _xml_element(
                "operator",
                {
                    "id": "rateCategoriesSwap",
                    "spec": "SwapOperator",
                    "howMany": "1",
                    "weight": "10",
                    "intparameter": "@rateCategories",
                },
            ),
            _xml_element(
                "operator",
                {
                    "id": "rateCategoriesUniform",
                    "spec": "UniformOperator",
                    "weight": "10",
                    "parameter": "@rateCategories",
                },
            ),
        ]
    )
    logger_elements.extend(
        [
            _xml_element("log", {"idref": "ucld.mean"}),
            _xml_element("log", {"idref": "ucld.stdev"}),
            _xml_element(
                "log",
                {
                    "id": "rateStatistic",
                    "spec": "beast.base.evolution.RateStatistic",
                    "tree": "@tree",
                    "branchratemodel": "@branchRates",
                },
            ),
            _xml_element("log", {"idref": "rateCategories"}),
        ]
    )
    return state_node_ids, prior_elements, operator_elements, logger_elements


def _append_tree_prior(
    root: XmlElement,
    *,
    tree_prior: str,
) -> tuple[list[str], list[XmlElement], list[XmlElement], list[XmlElement]]:
    normalized = tree_prior.strip().lower()
    state_node_ids: list[str] = ["birthRate"]
    prior_elements: list[XmlElement] = []
    operator_elements: list[XmlElement] = [
        _xml_element(
            "operator",
            {
                "id": "birthRateScaler",
                "spec": "ScaleOperator",
                "scaleFactor": "0.75",
                "weight": "1",
                "parameter": "@birthRate",
            },
        )
    ]
    logger_elements: list[XmlElement] = [_xml_element("log", {"idref": "birthRate"})]
    if normalized == "yule":
        yule = ElementTree.SubElement(
            root,
            "input",
            {"spec": "beast.base.evolution.speciation.YuleModel", "id": "treePrior"},
        )
        ElementTree.SubElement(yule, "birthDiffRate", {"idref": "birthRate"})
        ElementTree.SubElement(yule, "tree", {"idref": "tree"})
        ElementTree.SubElement(
            root,
            "parameter",
            {"id": "birthRate", "value": "1.0", "lower": "0.0", "upper": "100.0"},
        )
        prior_elements.extend(
            [
                _xml_element(
                    "distribution",
                    {"id": "treePrior.distribution", "idref": "treePrior"},
                ),
                _xml_element(
                    "distribution",
                    {
                        "id": "birthRate.prior",
                        "spec": "beast.base.inference.distribution.Prior",
                        "x": "@birthRate",
                    },
                    children=(
                        _xml_element(
                            "distr",
                            {
                                "id": "birthRate.oneOnX",
                                "offset": "0.0",
                                "spec": "beast.base.inference.distribution.OneOnX",
                            },
                        ),
                    ),
                ),
            ]
        )
        logger_elements.insert(0, _xml_element("log", {"idref": "treePrior"}))
        return state_node_ids, prior_elements, operator_elements, logger_elements
    if normalized != "birth-death":
        raise ValueError(
            "BEAST preparation supports tree_prior values 'yule' and 'birth-death'"
        )
    birth_death = ElementTree.SubElement(
        root,
        "input",
        {
            "spec": "beast.base.evolution.speciation.BirthDeathGernhard08Model",
            "id": "treePrior",
        },
    )
    ElementTree.SubElement(birth_death, "birthDiffRate", {"idref": "birthRate"})
    ElementTree.SubElement(
        birth_death, "relativeDeathRate", {"idref": "relativeDeathRate"}
    )
    ElementTree.SubElement(
        birth_death, "sampleProbability", {"idref": "sampleProbability"}
    )
    ElementTree.SubElement(birth_death, "tree", {"idref": "tree"})
    ElementTree.SubElement(
        root,
        "parameter",
        {"id": "birthRate", "value": "1.0", "lower": "0.0", "upper": "1000000.0"},
    )
    ElementTree.SubElement(
        root,
        "parameter",
        {"id": "relativeDeathRate", "value": "0.5", "lower": "0.0", "upper": "1.0"},
    )
    ElementTree.SubElement(
        root, "parameter", {"id": "sampleProbability", "value": "1.0"}
    )
    state_node_ids.append("relativeDeathRate")
    prior_elements.extend(
        [
            _xml_element(
                "distribution",
                {"id": "treePrior.distribution", "idref": "treePrior"},
            ),
            _xml_element(
                "distribution",
                {
                    "id": "birthRate.prior",
                    "spec": "beast.base.inference.distribution.Prior",
                    "x": "@birthRate",
                },
                children=(
                    _xml_element(
                        "distr",
                        {
                            "id": "birthRate.oneOnX",
                            "offset": "0.0",
                            "spec": "beast.base.inference.distribution.OneOnX",
                        },
                    ),
                ),
            ),
        ]
    )
    logger_elements.extend(
        [
            _xml_element("log", {"idref": "treePrior"}),
            _xml_element("log", {"idref": "relativeDeathRate"}),
            _xml_element("log", {"idref": "sampleProbability"}),
        ]
    )
    operator_elements.append(
        _xml_element(
            "operator",
            {
                "id": "relativeDeathRateScaler",
                "spec": "ScaleOperator",
                "scaleFactor": "0.75",
                "weight": "1",
                "parameter": "@relativeDeathRate",
            },
        )
    )
    return state_node_ids, prior_elements, operator_elements, logger_elements


def _translate_calibration_distribution(
    calibration: ValidatedCalibration,
) -> tuple[XmlElement, BeastCalibration, str | None]:
    calibration_id = _xml_identifier(calibration.calibration_id)
    mrca = ElementTree.Element(
        "distribution",
        {
            "spec": "beast.base.evolution.tree.MRCAPrior",
            "tree": "@tree",
            "id": calibration_id,
            "monophyletic": "true",
        },
    )
    taxonset = ElementTree.SubElement(
        mrca,
        "taxonset",
        {"spec": "TaxonSet", "id": f"{calibration_id}.taxa"},
    )
    for taxon in calibration.taxa:
        ElementTree.SubElement(taxonset, "taxon", {"spec": "Taxon", "id": taxon})

    requested = calibration.distribution.strip().lower() or "uniform"
    lower = calibration.minimum_age
    upper = calibration.maximum_age
    warning: str | None = None
    translated = False
    beast_distribution = "Uniform"
    if upper is not None and lower is None:
        raise ValueError(
            f"BEAST preparation does not support upper-bound-only calibration '{calibration.calibration_id}'"
        )
    if lower is not None and upper is not None:
        ElementTree.SubElement(
            mrca,
            "distr",
            {
                "id": f"{calibration_id}.bounds",
                "spec": "beast.base.inference.distribution.Uniform",
                "lower": _format_decimal(lower),
                "upper": _format_decimal(upper),
            },
        )
        if requested != "uniform":
            translated = True
            warning = (
                f"{calibration.calibration_id}: preserved the supplied hard bounds as a BEAST uniform prior "
                f"because the template generator does not infer parametric {requested} shape parameters automatically"
            )
    else:
        if lower is None:
            raise RuntimeError(
                f"BEAST lower-bound calibration unexpectedly missing for '{calibration.calibration_id}'"
            )
        if requested == "lognormal":
            beast_distribution = "LogNormalDistributionModel"
            ElementTree.SubElement(
                mrca,
                "distr",
                {
                    "id": f"{calibration_id}.lognormal",
                    "spec": "beast.base.inference.distribution.LogNormalDistributionModel",
                    "offset": _format_decimal(lower),
                    "M": "1.0",
                    "S": "1.25",
                    "meanInRealSpace": "false",
                },
            )
            translated = True
            warning = (
                f"{calibration.calibration_id}: translated a lower-bound-only calibration into an offset lognormal prior "
                "with default broad shape parameters"
            )
        else:
            beast_distribution = "Exponential"
            ElementTree.SubElement(
                mrca,
                "distr",
                {
                    "id": f"{calibration_id}.exponential",
                    "spec": "beast.base.inference.distribution.Exponential",
                    "offset": _format_decimal(lower),
                    "mean": _format_decimal(max(1.0, lower * 0.25)),
                },
            )
            translated = True
            warning = (
                f"{calibration.calibration_id}: translated a lower-bound-only {requested} calibration into an offset exponential prior "
                "because BEAST2 requires an explicit density above the hard minimum bound"
            )

    return (
        mrca,
        BeastCalibration(
            calibration_id=calibration.calibration_id,
            beast_distribution=beast_distribution,
            target_label=calibration.target_label,
            lower_bound=lower,
            upper_bound=upper,
            translated=translated,
            translation_note=warning,
        ),
        warning,
    )
