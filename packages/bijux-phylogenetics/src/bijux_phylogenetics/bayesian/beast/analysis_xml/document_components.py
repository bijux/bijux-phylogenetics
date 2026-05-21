from __future__ import annotations

from .._shared import ElementTree, Path, XmlElement, _tip_date_trait_value
from ..models import TipDatingValidationReport
from .xml_primitives import _read_newick_text, _xml_element


def _append_sequence_alignment(
    root: XmlElement,
    *,
    records,
    beast_data_type: str,
) -> None:
    data = ElementTree.SubElement(
        root,
        "data",
        {
            "id": "alignment",
            "dataType": beast_data_type,
        },
    )
    for record in records:
        sequence = ElementTree.SubElement(
            data, "sequence", {"taxon": record.identifier}
        )
        sequence.text = record.sequence


def _append_substitution_and_site_model(
    root: XmlElement,
    *,
    beast_data_type: str,
) -> tuple[list[str], list[XmlElement], list[XmlElement], list[XmlElement]]:
    state_node_ids: list[str] = []
    prior_elements: list[XmlElement] = []
    operator_elements: list[XmlElement] = []
    logger_elements: list[XmlElement] = []

    if beast_data_type == "nucleotide":
        hky = ElementTree.SubElement(root, "input", {"spec": "HKY", "id": "hky"})
        ElementTree.SubElement(
            hky, "parameter", {"name": "kappa", "idref": "hky.kappa"}
        )
        ElementTree.SubElement(
            hky,
            "input",
            {
                "id": "freqs",
                "name": "frequencies",
                "spec": "Frequencies",
                "frequencies": "@hky.frequencies",
            },
        )
        site_model = ElementTree.SubElement(
            root, "input", {"spec": "SiteModel", "id": "siteModel"}
        )
        ElementTree.SubElement(site_model, "substModel", {"idref": "hky"})
        ElementTree.SubElement(
            root,
            "parameter",
            {"id": "hky.kappa", "value": "2.0", "lower": "0.0"},
        )
        ElementTree.SubElement(
            root,
            "parameter",
            {"id": "hky.frequencies", "value": "0.25", "dimension": "4"},
        )
        state_node_ids.extend(["hky.kappa", "hky.frequencies"])

        prior_elements.extend(
            [
                _xml_element(
                    "distribution",
                    {
                        "id": "hky.kappa.prior",
                        "spec": "beast.base.inference.distribution.Prior",
                        "x": "@hky.kappa",
                    },
                    children=(
                        _xml_element(
                            "distr",
                            {
                                "id": "hky.kappa.lognormal",
                                "M": "1.0",
                                "S": "1.25",
                                "meanInRealSpace": "false",
                                "spec": "beast.base.inference.distribution.LogNormalDistributionModel",
                            },
                        ),
                    ),
                ),
                _xml_element(
                    "distribution",
                    {
                        "id": "hky.frequencies.prior",
                        "spec": "beast.base.inference.distribution.Prior",
                        "x": "@hky.frequencies",
                    },
                    children=(
                        _xml_element(
                            "distr",
                            {
                                "id": "hky.frequencies.uniform",
                                "lower": "0.0",
                                "upper": "1.0",
                                "spec": "beast.base.inference.distribution.Uniform",
                            },
                        ),
                    ),
                ),
            ]
        )
        operator_elements.extend(
            [
                _xml_element(
                    "operator",
                    {
                        "id": "kappaScaler",
                        "spec": "ScaleOperator",
                        "scaleFactor": "0.75",
                        "weight": "0.1",
                        "parameter": "@hky.kappa",
                    },
                ),
                _xml_element(
                    "operator",
                    {
                        "id": "frequenciesDelta",
                        "spec": "DeltaExchangeOperator",
                        "delta": "0.01",
                        "weight": "0.1",
                        "parameter": "@hky.frequencies",
                    },
                ),
            ]
        )
        logger_elements.extend(
            [
                _xml_element("log", {"idref": "hky.kappa"}),
                _xml_element("log", {"idref": "hky.frequencies"}),
            ]
        )
        return state_node_ids, prior_elements, operator_elements, logger_elements

    site_model = ElementTree.SubElement(
        root, "input", {"spec": "SiteModel", "id": "siteModel"}
    )
    ElementTree.SubElement(site_model, "substModel", {"spec": "JTT", "id": "jtt"})
    return state_node_ids, prior_elements, operator_elements, logger_elements


def _append_starting_tree(
    root: XmlElement,
    *,
    tree_path: Path | None,
    tip_date_report: TipDatingValidationReport | None,
) -> str:
    if tree_path is not None:
        tree = ElementTree.SubElement(
            root,
            "tree",
            {
                "spec": "beast.base.evolution.tree.TreeParser",
                "id": "tree",
                "IsLabelledNewick": "true",
                "newick": _read_newick_text(tree_path),
                "taxa": "@alignment",
            },
        )
        if tip_date_report is not None:
            trait = ElementTree.SubElement(
                tree,
                "trait",
                {
                    "spec": "beast.base.evolution.tree.TraitSet",
                    "traitname": "date-forward",
                    "units": "year",
                    "value": _tip_date_trait_value(tip_date_report),
                },
            )
            ElementTree.SubElement(
                trait,
                "taxa",
                {"spec": "TaxonSet", "alignment": "@alignment"},
            )
        return "provided-tree"

    cluster_tree = ElementTree.SubElement(
        root,
        "input",
        {
            "spec": "beast.base.evolution.tree.ClusterTree",
            "id": "tree",
            "clusterType": "upgma",
        },
    )
    if tip_date_report is not None:
        trait = ElementTree.SubElement(
            cluster_tree,
            "trait",
            {
                "spec": "beast.base.evolution.tree.TraitSet",
                "traitname": "date-forward",
                "units": "year",
                "value": _tip_date_trait_value(tip_date_report),
            },
        )
        ElementTree.SubElement(
            trait,
            "taxa",
            {"spec": "TaxonSet", "alignment": "@alignment"},
        )
    ElementTree.SubElement(cluster_tree, "taxa", {"idref": "alignment"})
    return "upgma"


def _append_tree_likelihood(root: XmlElement) -> None:
    likelihood = ElementTree.SubElement(
        root,
        "distribution",
        {
            "spec": "TreeLikelihood",
            "id": "likelihood",
            "data": "@alignment",
            "tree": "@tree",
        },
    )
    ElementTree.SubElement(likelihood, "siteModel", {"idref": "siteModel"})
    ElementTree.SubElement(likelihood, "branchRateModel", {"idref": "branchRates"})
