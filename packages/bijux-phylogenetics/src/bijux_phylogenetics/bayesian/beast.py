from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from bijux_phylogenetics.bayesian.diagnostics import (
    TraceConvergenceReport,
    summarize_trace_convergence,
)
from bijux_phylogenetics.bayesian.posterior import (
    summarize_maximum_clade_credibility_tree,
)
from bijux_phylogenetics.comparative.common import descendant_taxa
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.core.tree import PhyloTree
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.engines.workflows import _ensure_inference_ready_alignment
from bijux_phylogenetics.errors import EngineWorkflowError, InvalidAlignmentError
from bijux_phylogenetics.io.biopython import loads_biophylo
from bijux_phylogenetics.io.fasta import infer_alignment_alphabet, load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.tree_set import (
    compute_clade_frequency_table,
    compute_consensus_tree,
    load_tree_set,
)

_BEAST_TREE_PATTERN = re.compile(
    r"tree\s+([^\s=]+)\s*=\s*(.+?);", flags=re.IGNORECASE | re.DOTALL
)
_BEAST_TREE_STATE_PATTERN = re.compile(r"STATE_(\d+)$", flags=re.IGNORECASE)


@dataclass(slots=True)
class CalibrationValidationIssue:
    calibration_id: str
    code: str
    message: str


@dataclass(slots=True)
class ValidatedCalibration:
    calibration_id: str
    target_kind: str
    target_label: str
    taxa: list[str]
    minimum_age: float | None
    maximum_age: float | None
    distribution: str
    valid: bool


@dataclass(slots=True)
class FossilCalibrationValidationReport:
    tree_path: Path
    calibration_path: Path
    tree_taxa: list[str]
    calibration_count: int
    valid_calibration_count: int
    invalid_calibration_count: int
    calibrations: list[ValidatedCalibration]
    issues: list[CalibrationValidationIssue]


@dataclass(slots=True)
class ImpossibleCalibrationConstraintReport:
    tree_path: Path
    calibration_path: Path
    impossible_calibration_ids: list[str]
    issues: list[CalibrationValidationIssue]


@dataclass(slots=True)
class ValidatedTipDate:
    taxon: str
    date: float | None
    valid: bool


@dataclass(slots=True)
class TipDatingValidationIssue:
    taxon: str
    code: str
    message: str


@dataclass(slots=True)
class TipDatingValidationReport:
    tree_path: Path
    tip_dates_path: Path
    alignment_path: Path | None
    taxon_column: str
    date_column: str
    valid_tip_count: int
    invalid_tip_count: int
    missing_tree_taxa: list[str]
    extra_tip_taxa: list[str]
    extra_alignment_taxa: list[str]
    tip_dates: list[ValidatedTipDate]
    issues: list[TipDatingValidationIssue]


@dataclass(slots=True)
class CalibrationDominanceObservation:
    calibration_id: str
    target_label: str
    bounded_span_fraction: float | None
    dominates_root_age: bool
    warning: str | None


@dataclass(slots=True)
class CalibrationDominanceReport:
    tree_path: Path
    calibration_path: Path
    root_age: float
    valid_calibration_count: int
    dominant_calibration_ids: list[str]
    observations: list[CalibrationDominanceObservation]
    warnings: list[str]


@dataclass(slots=True)
class TimeTreeReadinessReport:
    tree_path: Path
    calibration_path: Path | None
    tip_dates_path: Path | None
    decision: str
    rooted: bool
    ultrametric: bool
    branch_length_status: str
    blockers: list[str]
    warnings: list[str]
    calibration_report: FossilCalibrationValidationReport | None
    tip_date_report: TipDatingValidationReport | None
    calibration_dominance: CalibrationDominanceReport | None


@dataclass(slots=True)
class BeastPreparationReport:
    alignment_path: Path
    output_path: Path
    tree_path: Path | None
    calibration_path: Path | None
    tip_dates_path: Path | None
    taxon_count: int
    character_count: int
    inferred_alphabet: str
    beast_data_type: str
    substitution_model: str
    clock_model: str
    tree_prior: str
    starting_tree_source: str
    chain_length: int
    log_every: int
    calibration_count: int
    tip_date_count: int
    warning_count: int
    warnings: list[str]
    log_path: Path
    tree_log_path: Path
    calibrations: list[BeastCalibration]


@dataclass(slots=True)
class BeastCalibration:
    calibration_id: str
    beast_distribution: str
    target_label: str
    lower_bound: float | None
    upper_bound: float | None
    translated: bool
    translation_note: str | None


@dataclass(slots=True)
class BeastLogRow:
    state: int
    values: dict[str, float]


@dataclass(slots=True)
class BeastLogReport:
    path: Path
    row_count: int
    columns: list[str]
    rows: list[BeastLogRow]


@dataclass(slots=True)
class BeastLogParameterSummary:
    parameter: str
    parameter_category: str
    sample_count: int
    effective_sample_size: float
    mean: float
    minimum: float
    maximum: float
    first_half_mean: float
    second_half_mean: float
    standardized_mean_shift: float


@dataclass(slots=True)
class BeastLogSummaryReport:
    path: Path
    burnin_fraction: float
    burnin_row_count: int
    kept_row_count: int
    first_kept_state: int
    last_kept_state: int
    posterior_parameters: list[str]
    likelihood_parameters: list[str]
    prior_parameters: list[str]
    clock_parameters: list[str]
    tree_parameters: list[str]
    other_parameters: list[str]
    parameter_summaries: list[BeastLogParameterSummary]


@dataclass(slots=True)
class BeastPosteriorTreeSample:
    tree_name: str
    state: int | None
    rooted: bool
    tip_names: list[str]
    newick: str


@dataclass(slots=True)
class BeastPosteriorClade:
    clade: str
    tree_count: int
    frequency: float


@dataclass(slots=True)
class BeastPosteriorTreeSetReport:
    path: Path
    burnin_fraction: float
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    rooted_tree_count: int
    sampled_states: list[int]
    tip_names: list[str]
    clades: list[BeastPosteriorClade]
    trees: list[BeastPosteriorTreeSample]


@dataclass(slots=True)
class BeastPosteriorConsensusReport:
    source_path: Path
    retained_tree_set_path: Path
    burnin_fraction: float
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    rooted_topology_count: int
    shared_taxa: list[str]
    consensus_newick: str
    clade_frequency_count: int
    annotated_node_count: int
    minimum_posterior_probability: float | None
    maximum_posterior_probability: float | None


@dataclass(slots=True)
class BeastLogValidationIssue:
    code: str
    message: str
    row: int | None = None
    column: str | None = None


@dataclass(slots=True)
class BeastPosteriorLogValidationReport:
    path: Path
    row_count: int
    state_count: int
    required_columns: list[str]
    observed_columns: list[str]
    missing_columns: list[str]
    issues: list[BeastLogValidationIssue]
    valid: bool


@dataclass(slots=True)
class BeastBurninSensitivitySlice:
    burnin_fraction: float
    burnin_tree_count: int
    kept_tree_count: int
    rooted_topology_count: int
    selected_tree_index: int
    clade_credibility_score: float
    posterior_mean: float | None
    likelihood_mean: float | None
    tree_height_mean: float | None


@dataclass(slots=True)
class BeastBurninSensitivityReport:
    posterior_tree_path: Path
    log_path: Path | None
    slices: list[BeastBurninSensitivitySlice]
    changed_mcc_count: int
    warnings: list[str]


@dataclass(slots=True)
class BeastChainMixingIssue:
    path: Path | None
    parameter: str
    code: str
    message: str
    observed_value: float
    threshold: float


@dataclass(slots=True)
class BeastChainMixingReport:
    log_paths: list[Path]
    chain_count: int
    converged: bool
    issues: list[BeastChainMixingIssue]
    chain_summaries: list[BeastConvergenceReport]


@dataclass(slots=True)
class BeastConvergenceReport:
    path: Path
    burnin_fraction: float
    burnin_row_count: int
    sample_count: int
    converged: bool
    ess_threshold: float
    mean_shift_threshold: float
    warnings: list[dict[str, object]]
    parameter_summaries: list[dict[str, object]]


def _read_delimited_rows(path: Path) -> list[dict[str, str]]:
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"calibration table contains no header: {path}")
        return [
            {key: value or "" for key, value in row.items() if key is not None}
            for row in reader
        ]


def _named_clades(tree: PhyloTree) -> dict[str, list[str]]:
    named: dict[str, list[str]] = {}
    for node in tree.iter_nodes():
        if node.name:
            named[node.name] = descendant_taxa(node)
    return named


def _clade_taxon_sets(tree: PhyloTree) -> set[frozenset[str]]:
    clades: set[frozenset[str]] = set()
    for node in tree.iter_nodes():
        taxa = descendant_taxa(node)
        if taxa:
            clades.add(frozenset(taxa))
    return clades


def _parse_target_taxa(raw: str) -> list[str]:
    if not raw.strip():
        return []
    normalized = raw.replace(",", "|").replace(";", "|")
    return sorted({token.strip() for token in normalized.split("|") if token.strip()})


def _parse_age(
    raw: str,
    *,
    calibration_id: str,
    field_name: str,
    issues: list[CalibrationValidationIssue],
) -> float | None:
    if not raw.strip():
        return None
    try:
        value = float(raw)
    except ValueError:
        issues.append(
            CalibrationValidationIssue(
                calibration_id=calibration_id,
                code="invalid-age",
                message=f"{field_name} must be numeric when provided",
            )
        )
        return None
    return value


_XML_IDENTIFIER_PATTERN = re.compile(r"[^0-9A-Za-z._-]+")


def _xml_identifier(raw: str) -> str:
    normalized = _XML_IDENTIFIER_PATTERN.sub("_", raw.strip())
    normalized = normalized.strip("_")
    return normalized or "identifier"


def _format_decimal(value: float) -> str:
    return format(value, ".15g")


def _beast_data_type(inferred_alphabet: str) -> str:
    if inferred_alphabet in {"dna", "rna"}:
        return "nucleotide"
    if inferred_alphabet == "protein":
        return "aminoacid"
    raise ValueError(
        "BEAST preparation requires a nucleotide, RNA, or protein alignment"
    )


def _default_beast_substitution_model(beast_data_type: str) -> str:
    if beast_data_type == "nucleotide":
        return "HKY"
    return "JTT"


def _read_newick_text(path: Path) -> str:
    raw = path.read_text(encoding="utf-8").strip()
    if not raw:
        raise ValueError(f"tree file is empty: {path}")
    first_line = raw.splitlines()[0].strip()
    if not first_line.endswith(";"):
        first_line = f"{first_line};"
    return first_line


def _tip_date_trait_value(report: TipDatingValidationReport) -> str:
    parts = [
        f"{tip.taxon}={_format_decimal(tip.date)}"
        for tip in report.tip_dates
        if tip.valid and tip.date is not None
    ]
    return ",".join(parts)


def _validate_tree_taxa_against_alignment(
    *, tree_path: Path, alignment_taxa: set[str]
) -> None:
    tree_taxa = set(load_tree(tree_path).tip_names)
    missing_from_tree = sorted(alignment_taxa - tree_taxa)
    extra_in_tree = sorted(tree_taxa - alignment_taxa)
    if missing_from_tree or extra_in_tree:
        details: list[str] = []
        if missing_from_tree:
            details.append(
                "alignment taxa missing from tree: " + ", ".join(missing_from_tree)
            )
        if extra_in_tree:
            details.append("tree taxa missing from alignment: " + ", ".join(extra_in_tree))
        raise ValueError(
            "BEAST preparation requires the starting tree and alignment to contain the same taxa: "
            + "; ".join(details)
        )


def _append_sequence_alignment(
    root: ET.Element,
    *,
    records,
    beast_data_type: str,
) -> None:
    data = ET.SubElement(
        root,
        "data",
        {
            "id": "alignment",
            "dataType": beast_data_type,
        },
    )
    for record in records:
        sequence = ET.SubElement(data, "sequence", {"taxon": record.identifier})
        sequence.text = record.sequence


def _append_substitution_and_site_model(
    root: ET.Element,
    *,
    beast_data_type: str,
) -> tuple[list[str], list[ET.Element], list[ET.Element], list[ET.Element]]:
    state_node_ids: list[str] = []
    prior_elements: list[ET.Element] = []
    operator_elements: list[ET.Element] = []
    logger_elements: list[ET.Element] = []

    if beast_data_type == "nucleotide":
        hky = ET.SubElement(root, "input", {"spec": "HKY", "id": "hky"})
        ET.SubElement(hky, "parameter", {"name": "kappa", "idref": "hky.kappa"})
        ET.SubElement(
            hky,
            "input",
            {
                "id": "freqs",
                "name": "frequencies",
                "spec": "Frequencies",
                "frequencies": "@hky.frequencies",
            },
        )
        site_model = ET.SubElement(root, "input", {"spec": "SiteModel", "id": "siteModel"})
        ET.SubElement(site_model, "substModel", {"idref": "hky"})
        ET.SubElement(
            root,
            "parameter",
            {"id": "hky.kappa", "value": "2.0", "lower": "0.0"},
        )
        ET.SubElement(
            root,
            "parameter",
            {"id": "hky.frequencies", "value": "0.25", "dimension": "4"},
        )
        state_node_ids.extend(["hky.kappa", "hky.frequencies"])

        prior_elements.extend(
            [
                ET.fromstring(
                    "<distribution id='hky.kappa.prior' spec='beast.base.inference.distribution.Prior' x='@hky.kappa'>"
                    "<distr id='hky.kappa.lognormal' M='1.0' S='1.25' meanInRealSpace='false' "
                    "spec='beast.base.inference.distribution.LogNormalDistributionModel' />"
                    "</distribution>"
                ),
                ET.fromstring(
                    "<distribution id='hky.frequencies.prior' spec='beast.base.inference.distribution.Prior' x='@hky.frequencies'>"
                    "<distr id='hky.frequencies.uniform' lower='0.0' upper='1.0' "
                    "spec='beast.base.inference.distribution.Uniform' />"
                    "</distribution>"
                ),
            ]
        )
        operator_elements.extend(
            [
                ET.fromstring(
                    "<operator id='kappaScaler' spec='ScaleOperator' scaleFactor='0.75' weight='0.1' parameter='@hky.kappa' />"
                ),
                ET.fromstring(
                    "<operator id='frequenciesDelta' spec='DeltaExchangeOperator' delta='0.01' weight='0.1' parameter='@hky.frequencies' />"
                ),
            ]
        )
        logger_elements.extend(
            [
                ET.fromstring("<log idref='hky.kappa' />"),
                ET.fromstring("<log idref='hky.frequencies' />"),
            ]
        )
        return state_node_ids, prior_elements, operator_elements, logger_elements

    site_model = ET.SubElement(root, "input", {"spec": "SiteModel", "id": "siteModel"})
    ET.SubElement(site_model, "substModel", {"spec": "JTT", "id": "jtt"})
    return state_node_ids, prior_elements, operator_elements, logger_elements


def _append_starting_tree(
    root: ET.Element,
    *,
    tree_path: Path | None,
    tip_date_report: TipDatingValidationReport | None,
) -> str:
    if tree_path is not None:
        tree = ET.SubElement(
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
            trait = ET.SubElement(
                tree,
                "trait",
                {
                    "spec": "beast.base.evolution.tree.TraitSet",
                    "traitname": "date-forward",
                    "units": "year",
                    "value": _tip_date_trait_value(tip_date_report),
                },
            )
            ET.SubElement(
                trait,
                "taxa",
                {"spec": "TaxonSet", "alignment": "@alignment"},
            )
        return "provided-tree"

    cluster_tree = ET.SubElement(
        root,
        "input",
        {
            "spec": "beast.base.evolution.tree.ClusterTree",
            "id": "tree",
            "clusterType": "upgma",
        },
    )
    if tip_date_report is not None:
        trait = ET.SubElement(
            cluster_tree,
            "trait",
            {
                "spec": "beast.base.evolution.tree.TraitSet",
                "traitname": "date-forward",
                "units": "year",
                "value": _tip_date_trait_value(tip_date_report),
            },
        )
        ET.SubElement(
            trait,
            "taxa",
            {"spec": "TaxonSet", "alignment": "@alignment"},
        )
    ET.SubElement(cluster_tree, "taxa", {"idref": "alignment"})
    return "upgma"


def _append_clock_model(
    root: ET.Element,
    *,
    clock_model: str,
    taxon_count: int,
) -> tuple[list[str], list[ET.Element], list[ET.Element], list[ET.Element]]:
    normalized = clock_model.strip().lower()
    state_node_ids: list[str] = []
    prior_elements: list[ET.Element] = []
    operator_elements: list[ET.Element] = []
    logger_elements: list[ET.Element] = []
    if normalized == "strict":
        strict = ET.SubElement(
            root,
            "input",
            {
                "spec": "beast.base.evolution.branchratemodel.StrictClockModel",
                "id": "branchRates",
            },
        )
        ET.SubElement(strict, "clock.rate", {"idref": "clockRate"})
        ET.SubElement(
            root,
            "parameter",
            {"id": "clockRate", "value": "0.001", "lower": "0.0", "upper": "100.0"},
        )
        state_node_ids.append("clockRate")
        prior_elements.append(
            ET.fromstring(
                "<distribution id='clockRate.prior' spec='beast.base.inference.distribution.Prior' x='@clockRate'>"
                "<distr id='clockRate.uniform' lower='0.0' upper='100.0' spec='beast.base.inference.distribution.Uniform' />"
                "</distribution>"
            )
        )
        operator_elements.append(
            ET.fromstring(
                "<operator id='clockRateScaler' spec='ScaleOperator' scaleFactor='0.75' weight='3' parameter='@clockRate' />"
            )
        )
        logger_elements.append(ET.fromstring("<log idref='clockRate' />"))
        return state_node_ids, prior_elements, operator_elements, logger_elements
    if normalized != "relaxed-lognormal":
        raise ValueError(
            "BEAST preparation supports clock_model values 'strict' and 'relaxed-lognormal'"
        )
    relaxed = ET.SubElement(
        root,
        "input",
        {
            "spec": "beast.base.evolution.branchratemodel.UCRelaxedClockModel",
            "id": "branchRates",
        },
    )
    distr = ET.SubElement(
        relaxed,
        "distr",
        {
            "id": "ucld.lognormal",
            "spec": "beast.base.inference.distribution.LogNormalDistributionModel",
            "meanInRealSpace": "true",
        },
    )
    ET.SubElement(
        distr,
        "parameter",
        {"name": "M", "id": "ucld.mean", "value": "1.0", "lower": "0.0"},
    )
    ET.SubElement(
        distr,
        "parameter",
        {"name": "S", "id": "ucld.stdev", "value": "0.3333333333333333", "lower": "0.0"},
    )
    ET.SubElement(
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
    ET.SubElement(relaxed, "tree", {"idref": "tree"})
    state_node_ids.extend(["ucld.mean", "ucld.stdev", "rateCategories"])
    prior_elements.append(
        ET.fromstring(
            "<distribution id='ucld.stdev.prior' spec='beast.base.inference.distribution.Prior' x='@ucld.stdev'>"
            "<distr id='ucld.stdev.exponential' mean='0.3333333333333333' spec='beast.base.inference.distribution.Exponential' />"
            "</distribution>"
        )
    )
    operator_elements.extend(
        [
            ET.fromstring(
                "<operator id='ucldMeanScaler' spec='ScaleOperator' scaleFactor='0.75' weight='1' parameter='@ucld.mean' />"
            ),
            ET.fromstring(
                "<operator id='ucldStdevScaler' spec='ScaleOperator' scaleFactor='0.75' weight='3' parameter='@ucld.stdev' />"
            ),
            ET.fromstring(
                "<operator id='rateCategoriesRandomWalk' spec='IntRandomWalkOperator' windowSize='1' weight='10' parameter='@rateCategories' />"
            ),
            ET.fromstring(
                "<operator id='rateCategoriesSwap' spec='SwapOperator' howMany='1' weight='10' intparameter='@rateCategories' />"
            ),
            ET.fromstring(
                "<operator id='rateCategoriesUniform' spec='UniformOperator' weight='10' parameter='@rateCategories' />"
            ),
        ]
    )
    logger_elements.extend(
        [
            ET.fromstring("<log idref='ucld.mean' />"),
            ET.fromstring("<log idref='ucld.stdev' />"),
            ET.fromstring(
                "<log id='rateStatistic' spec='beast.base.evolution.RateStatistic' tree='@tree' branchratemodel='@branchRates' />"
            ),
            ET.fromstring("<log idref='rateCategories' />"),
        ]
    )
    return state_node_ids, prior_elements, operator_elements, logger_elements


def _append_tree_prior(
    root: ET.Element,
    *,
    tree_prior: str,
) -> tuple[list[str], list[ET.Element], list[ET.Element], list[ET.Element]]:
    normalized = tree_prior.strip().lower()
    state_node_ids: list[str] = ["birthRate"]
    prior_elements: list[ET.Element] = []
    operator_elements: list[ET.Element] = [
        ET.fromstring(
            "<operator id='birthRateScaler' spec='ScaleOperator' scaleFactor='0.75' weight='1' parameter='@birthRate' />"
        )
    ]
    logger_elements: list[ET.Element] = [ET.fromstring("<log idref='birthRate' />")]
    if normalized == "yule":
        yule = ET.SubElement(
            root,
            "input",
            {"spec": "beast.base.evolution.speciation.YuleModel", "id": "treePrior"},
        )
        ET.SubElement(yule, "birthDiffRate", {"idref": "birthRate"})
        ET.SubElement(yule, "tree", {"idref": "tree"})
        ET.SubElement(
            root,
            "parameter",
            {"id": "birthRate", "value": "1.0", "lower": "0.0", "upper": "100.0"},
        )
        prior_elements.extend(
            [
                ET.fromstring("<distribution id='treePrior.distribution' idref='treePrior' />"),
                ET.fromstring(
                    "<distribution id='birthRate.prior' spec='beast.base.inference.distribution.Prior' x='@birthRate'>"
                    "<distr id='birthRate.oneOnX' offset='0.0' spec='beast.base.inference.distribution.OneOnX' />"
                    "</distribution>"
                ),
            ]
        )
        logger_elements.insert(0, ET.fromstring("<log idref='treePrior' />"))
        return state_node_ids, prior_elements, operator_elements, logger_elements
    if normalized != "birth-death":
        raise ValueError(
            "BEAST preparation supports tree_prior values 'yule' and 'birth-death'"
        )
    birth_death = ET.SubElement(
        root,
        "input",
        {
            "spec": "beast.base.evolution.speciation.BirthDeathGernhard08Model",
            "id": "treePrior",
        },
    )
    ET.SubElement(birth_death, "birthDiffRate", {"idref": "birthRate"})
    ET.SubElement(birth_death, "relativeDeathRate", {"idref": "relativeDeathRate"})
    ET.SubElement(birth_death, "sampleProbability", {"idref": "sampleProbability"})
    ET.SubElement(birth_death, "tree", {"idref": "tree"})
    ET.SubElement(
        root,
        "parameter",
        {"id": "birthRate", "value": "1.0", "lower": "0.0", "upper": "1000000.0"},
    )
    ET.SubElement(
        root,
        "parameter",
        {"id": "relativeDeathRate", "value": "0.5", "lower": "0.0", "upper": "1.0"},
    )
    ET.SubElement(root, "parameter", {"id": "sampleProbability", "value": "1.0"})
    state_node_ids.append("relativeDeathRate")
    prior_elements.extend(
        [
            ET.fromstring("<distribution id='treePrior.distribution' idref='treePrior' />"),
            ET.fromstring(
                "<distribution id='birthRate.prior' spec='beast.base.inference.distribution.Prior' x='@birthRate'>"
                "<distr id='birthRate.oneOnX' offset='0.0' spec='beast.base.inference.distribution.OneOnX' />"
                "</distribution>"
            ),
        ]
    )
    logger_elements.extend(
        [
            ET.fromstring("<log idref='treePrior' />"),
            ET.fromstring("<log idref='relativeDeathRate' />"),
            ET.fromstring("<log idref='sampleProbability' />"),
        ]
    )
    operator_elements.append(
        ET.fromstring(
            "<operator id='relativeDeathRateScaler' spec='ScaleOperator' scaleFactor='0.75' weight='1' parameter='@relativeDeathRate' />"
        )
    )
    return state_node_ids, prior_elements, operator_elements, logger_elements


def _append_tree_likelihood(root: ET.Element) -> None:
    likelihood = ET.SubElement(
        root,
        "distribution",
        {
            "spec": "TreeLikelihood",
            "id": "likelihood",
            "data": "@alignment",
            "tree": "@tree",
        },
    )
    ET.SubElement(likelihood, "siteModel", {"idref": "siteModel"})
    ET.SubElement(likelihood, "branchRateModel", {"idref": "branchRates"})


def _translate_calibration_distribution(
    calibration: ValidatedCalibration,
) -> tuple[ET.Element, BeastCalibration, str | None]:
    calibration_id = _xml_identifier(calibration.calibration_id)
    mrca = ET.Element(
        "distribution",
        {
            "spec": "beast.base.evolution.tree.MRCAPrior",
            "tree": "@tree",
            "id": calibration_id,
            "monophyletic": "true",
        },
    )
    taxonset = ET.SubElement(
        mrca,
        "taxonset",
        {"spec": "TaxonSet", "id": f"{calibration_id}.taxa"},
    )
    for taxon in calibration.taxa:
        ET.SubElement(taxonset, "taxon", {"spec": "Taxon", "id": taxon})

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
        ET.SubElement(
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
        assert lower is not None
        if requested == "lognormal":
            beast_distribution = "LogNormalDistributionModel"
            ET.SubElement(
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
            ET.SubElement(
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


def validate_fossil_calibration_table(
    tree_path: Path, calibration_path: Path
) -> FossilCalibrationValidationReport:
    """Validate fossil calibration targets and age constraints against one tree."""
    tree = load_tree(tree_path)
    rows = _read_delimited_rows(calibration_path)
    named_clades = _named_clades(tree)
    clade_sets = _clade_taxon_sets(tree)
    tree_taxa = sorted(tree.tip_names)
    tree_taxa_set = set(tree_taxa)
    calibrations: list[ValidatedCalibration] = []
    issues: list[CalibrationValidationIssue] = []

    for index, row in enumerate(rows, start=1):
        calibration_id = row.get("calibration_id", "").strip() or f"calibration-{index}"
        clade_name = row.get("clade_name", "").strip()
        taxa = _parse_target_taxa(row.get("taxa", ""))
        distribution = row.get("distribution", "").strip() or "uniform"
        minimum_age = _parse_age(
            row.get("minimum_age", ""),
            calibration_id=calibration_id,
            field_name="minimum_age",
            issues=issues,
        )
        maximum_age = _parse_age(
            row.get("maximum_age", ""),
            calibration_id=calibration_id,
            field_name="maximum_age",
            issues=issues,
        )

        target_kind = "taxa"
        target_label = clade_name or "|".join(taxa)
        resolved_taxa: list[str] = []
        valid = True

        if clade_name:
            target_kind = "named-clade"
            resolved_taxa = named_clades.get(clade_name, [])
            if not resolved_taxa:
                valid = False
                issues.append(
                    CalibrationValidationIssue(
                        calibration_id=calibration_id,
                        code="unknown-clade-name",
                        message=f"named clade '{clade_name}' is not present in the tree",
                    )
                )
        elif taxa:
            missing_taxa = sorted(set(taxa) - tree_taxa_set)
            if missing_taxa:
                valid = False
                issues.append(
                    CalibrationValidationIssue(
                        calibration_id=calibration_id,
                        code="unknown-taxa",
                        message=f"calibration taxa are absent from the tree: {', '.join(missing_taxa)}",
                    )
                )
            elif frozenset(taxa) not in clade_sets:
                valid = False
                issues.append(
                    CalibrationValidationIssue(
                        calibration_id=calibration_id,
                        code="non-monophyletic-target",
                        message="calibration taxa do not map to a named clade, single taxon, or monophyletic descendant set",
                    )
                )
            resolved_taxa = taxa
        else:
            valid = False
            issues.append(
                CalibrationValidationIssue(
                    calibration_id=calibration_id,
                    code="missing-target",
                    message="calibration must provide either a clade_name or taxa column",
                )
            )

        if minimum_age is not None and minimum_age < 0.0:
            valid = False
            issues.append(
                CalibrationValidationIssue(
                    calibration_id=calibration_id,
                    code="negative-minimum-age",
                    message="minimum calibration age cannot be negative",
                )
            )
        if maximum_age is not None and maximum_age < 0.0:
            valid = False
            issues.append(
                CalibrationValidationIssue(
                    calibration_id=calibration_id,
                    code="negative-maximum-age",
                    message="maximum calibration age cannot be negative",
                )
            )
        if (
            minimum_age is not None
            and maximum_age is not None
            and minimum_age > maximum_age
        ):
            valid = False
            issues.append(
                CalibrationValidationIssue(
                    calibration_id=calibration_id,
                    code="minimum-exceeds-maximum",
                    message="minimum calibration age cannot exceed maximum calibration age",
                )
            )
        if minimum_age is None and maximum_age is None:
            valid = False
            issues.append(
                CalibrationValidationIssue(
                    calibration_id=calibration_id,
                    code="missing-age-bounds",
                    message="calibration must provide at least one age bound",
                )
            )

        calibrations.append(
            ValidatedCalibration(
                calibration_id=calibration_id,
                target_kind=target_kind,
                target_label=target_label,
                taxa=resolved_taxa,
                minimum_age=minimum_age,
                maximum_age=maximum_age,
                distribution=distribution,
                valid=valid,
            )
        )

    valid_count = sum(1 for calibration in calibrations if calibration.valid)
    return FossilCalibrationValidationReport(
        tree_path=tree_path,
        calibration_path=calibration_path,
        tree_taxa=tree_taxa,
        calibration_count=len(calibrations),
        valid_calibration_count=valid_count,
        invalid_calibration_count=len(calibrations) - valid_count,
        calibrations=calibrations,
        issues=issues,
    )


def detect_impossible_calibration_constraints(
    tree_path: Path, calibration_path: Path
) -> ImpossibleCalibrationConstraintReport:
    """Return calibrations with impossible target or age constraints."""
    report = validate_fossil_calibration_table(tree_path, calibration_path)
    impossible_codes = {
        "unknown-clade-name",
        "unknown-taxa",
        "non-monophyletic-target",
        "missing-target",
        "negative-minimum-age",
        "negative-maximum-age",
        "minimum-exceeds-maximum",
        "missing-age-bounds",
    }
    impossible_ids = sorted(
        {
            issue.calibration_id
            for issue in report.issues
            if issue.code in impossible_codes
        }
    )
    return ImpossibleCalibrationConstraintReport(
        tree_path=tree_path,
        calibration_path=calibration_path,
        impossible_calibration_ids=impossible_ids,
        issues=[
            issue for issue in report.issues if issue.calibration_id in impossible_ids
        ],
    )


def _tree_root_age(tree_path: Path) -> float:
    tree = load_tree(tree_path)
    lengths = [length for length in tree.root_to_tip_lengths() if length is not None]
    if not lengths:
        return 0.0
    return round(max(float(length) for length in lengths), 15)


def assess_calibration_dominance(
    tree_path: Path, calibration_path: Path
) -> CalibrationDominanceReport:
    """Flag cases where one calibration contributes a disproportionate share of the dated-tree age range."""
    report = validate_fossil_calibration_table(tree_path, calibration_path)
    root_age = _tree_root_age(tree_path)
    observations: list[CalibrationDominanceObservation] = []
    dominant_calibration_ids: list[str] = []
    warnings: list[str] = []
    for calibration in report.calibrations:
        if not calibration.valid:
            continue
        span_fraction: float | None = None
        dominates_root_age = False
        warning: str | None = None
        if (
            calibration.minimum_age is not None
            and calibration.maximum_age is not None
            and root_age > 0.0
        ):
            span_fraction = round(
                (calibration.maximum_age - calibration.minimum_age) / root_age, 15
            )
            dominates_root_age = span_fraction >= 0.5
            if dominates_root_age:
                warning = "calibration age span covers at least half of the current root-age scale"
        elif root_age > 0.0 and (
            (
                calibration.minimum_age is not None
                and calibration.minimum_age / root_age >= 0.8
            )
            or (
                calibration.maximum_age is not None
                and calibration.maximum_age / root_age >= 0.8
            )
        ):
            dominates_root_age = True
            warning = "single-sided calibration bound lies close to the current root-age scale"
        if dominates_root_age:
            dominant_calibration_ids.append(calibration.calibration_id)
            if warning is not None:
                warnings.append(f"{calibration.calibration_id}: {warning}")
        observations.append(
            CalibrationDominanceObservation(
                calibration_id=calibration.calibration_id,
                target_label=calibration.target_label,
                bounded_span_fraction=span_fraction,
                dominates_root_age=dominates_root_age,
                warning=warning,
            )
        )
    if report.valid_calibration_count == 1 and report.calibrations:
        only_calibration = next(
            (calibration for calibration in report.calibrations if calibration.valid),
            None,
        )
        if (
            only_calibration is not None
            and only_calibration.calibration_id not in dominant_calibration_ids
        ):
            dominant_calibration_ids.append(only_calibration.calibration_id)
        warnings.append(
            "only one valid calibration remains, so time estimates are effectively driven by a single calibration target"
        )
    return CalibrationDominanceReport(
        tree_path=tree_path,
        calibration_path=calibration_path,
        root_age=root_age,
        valid_calibration_count=report.valid_calibration_count,
        dominant_calibration_ids=sorted(dominant_calibration_ids),
        observations=sorted(observations, key=lambda row: row.calibration_id),
        warnings=sorted(dict.fromkeys(warnings)),
    )


def validate_tip_dating_metadata(
    tree_path: Path,
    tip_dates_path: Path,
    *,
    alignment_path: Path | None = None,
    taxon_column: str | None = None,
    date_column: str = "date",
) -> TipDatingValidationReport:
    """Validate that dated tips resolve cleanly against the tree and optional alignment."""
    tree = load_tree(tree_path)
    tree_taxa = set(tree.tip_names)
    alignment_taxa: set[str] | None = None
    if alignment_path is not None:
        _ensure_inference_ready_alignment(alignment_path)
        alignment_taxa = {
            record.identifier for record in load_fasta_alignment(alignment_path)
        }
    table = load_taxon_table(tip_dates_path, taxon_column=taxon_column)
    if date_column not in table.columns:
        raise ValueError(f"tip-dating table does not contain column '{date_column}'")

    issues: list[TipDatingValidationIssue] = []
    tip_dates: list[ValidatedTipDate] = []
    extra_tip_taxa = sorted(set(table.taxa) - tree_taxa)
    missing_tree_taxa = sorted(tree_taxa - set(table.taxa))
    extra_alignment_taxa = (
        sorted(set(table.taxa) - alignment_taxa) if alignment_taxa is not None else []
    )
    for row in table.rows:
        taxon = row[table.taxon_column]
        raw_date = row.get(date_column, "")
        valid = True
        parsed_date: float | None = None
        if taxon not in tree_taxa:
            valid = False
            issues.append(
                TipDatingValidationIssue(
                    taxon=taxon,
                    code="taxon-missing-from-tree",
                    message="dated tip is absent from the tree",
                )
            )
        if alignment_taxa is not None and taxon not in alignment_taxa:
            valid = False
            issues.append(
                TipDatingValidationIssue(
                    taxon=taxon,
                    code="taxon-missing-from-alignment",
                    message="dated tip is absent from the alignment",
                )
            )
        if not raw_date.strip():
            valid = False
            issues.append(
                TipDatingValidationIssue(
                    taxon=taxon,
                    code="missing-date",
                    message="dated tip requires a numeric date value",
                )
            )
        else:
            try:
                parsed_date = float(raw_date)
            except ValueError:
                valid = False
                issues.append(
                    TipDatingValidationIssue(
                        taxon=taxon,
                        code="invalid-date",
                        message="dated tip value must be numeric",
                    )
                )
        tip_dates.append(ValidatedTipDate(taxon=taxon, date=parsed_date, valid=valid))

    valid_tip_count = sum(1 for tip in tip_dates if tip.valid)
    return TipDatingValidationReport(
        tree_path=tree_path,
        tip_dates_path=tip_dates_path,
        alignment_path=alignment_path,
        taxon_column=table.taxon_column,
        date_column=date_column,
        valid_tip_count=valid_tip_count,
        invalid_tip_count=len(tip_dates) - valid_tip_count,
        missing_tree_taxa=missing_tree_taxa,
        extra_tip_taxa=extra_tip_taxa,
        extra_alignment_taxa=extra_alignment_taxa,
        tip_dates=tip_dates,
        issues=issues,
    )


def assess_time_tree_readiness(
    tree_path: Path,
    *,
    calibration_path: Path | None = None,
    tip_dates_path: Path | None = None,
    alignment_path: Path | None = None,
    taxon_column: str | None = None,
    date_column: str = "date",
) -> TimeTreeReadinessReport:
    """Decide whether a dataset is suitable for dated phylogenetics."""
    validation = validate_tree_path(tree_path)
    blockers: list[str] = []
    warnings: list[str] = []
    if not validation.rooted:
        blockers.append("time-tree analysis requires a rooted tree")
    if validation.ultrametric is not True:
        blockers.append("time-tree analysis requires ultrametric branch lengths")
    if validation.branch_length_status != "complete":
        blockers.append("time-tree analysis requires complete branch lengths")
    calibration_report = None
    calibration_dominance = None
    if calibration_path is not None:
        calibration_report = validate_fossil_calibration_table(
            tree_path, calibration_path
        )
        if calibration_report.invalid_calibration_count:
            blockers.append(
                "calibration table contains invalid fossil calibration targets or ages"
            )
        calibration_dominance = assess_calibration_dominance(
            tree_path, calibration_path
        )
        warnings.extend(calibration_dominance.warnings)
    tip_date_report = None
    if tip_dates_path is not None:
        tip_date_report = validate_tip_dating_metadata(
            tree_path,
            tip_dates_path,
            alignment_path=alignment_path,
            taxon_column=taxon_column,
            date_column=date_column,
        )
        if tip_date_report.invalid_tip_count:
            blockers.append(
                "tip-date table contains missing, invalid, or mismatched dated taxa"
            )
    if calibration_path is None and tip_dates_path is None:
        warnings.append(
            "no calibrations or tip dates were supplied, so the tree cannot be dated by the current workflow"
        )
    decision = "ready"
    if blockers:
        decision = "blocked"
    elif warnings:
        decision = "risky"
    return TimeTreeReadinessReport(
        tree_path=tree_path,
        calibration_path=calibration_path,
        tip_dates_path=tip_dates_path,
        decision=decision,
        rooted=validation.rooted,
        ultrametric=validation.ultrametric is True,
        branch_length_status=validation.branch_length_status,
        blockers=sorted(dict.fromkeys(blockers)),
        warnings=sorted(dict.fromkeys([*warnings, *validation.warnings])),
        calibration_report=calibration_report,
        tip_date_report=tip_date_report,
        calibration_dominance=calibration_dominance,
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

    root = ET.Element(
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

    run = ET.SubElement(
        root,
        "run",
        {
            "spec": "MCMC",
            "id": "mcmc",
            "chainLength": str(chain_length),
        },
    )
    state = ET.SubElement(run, "state")
    for node_id in dict.fromkeys(
        [*site_state_ids, *clock_state_ids, *tree_state_ids, "tree"]
    ):
        ET.SubElement(state, "stateNode", {"idref": node_id})

    posterior = ET.SubElement(run, "distribution", {"spec": "CompoundDistribution", "id": "posterior"})
    prior = ET.SubElement(posterior, "distribution", {"spec": "CompoundDistribution", "id": "prior"})
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
            prior_element, beast_calibration, warning = _translate_calibration_distribution(
                calibration
            )
            prior.append(prior_element)
            beast_calibrations.append(beast_calibration)
            if warning is not None:
                warnings.append(warning)

    ET.SubElement(posterior, "distribution", {"idref": "likelihood"})

    generic_tree_operators = [
        ET.fromstring(
            "<operator id='treeScaler' spec='ScaleOperator' scaleFactor='0.75' weight='3' tree='@tree' />"
        ),
        ET.fromstring("<operator spec='Uniform' weight='30' tree='@tree' />"),
        ET.fromstring(
            "<operator spec='SubtreeSlide' weight='15' gaussian='true' size='0.5' tree='@tree' />"
        ),
        ET.fromstring(
            "<operator id='narrowExchange' spec='Exchange' isNarrow='true' weight='15' tree='@tree' />"
        ),
        ET.fromstring(
            "<operator id='wideExchange' spec='Exchange' isNarrow='false' weight='3' tree='@tree' />"
        ),
        ET.fromstring("<operator id='wilsonBalding' spec='WilsonBalding' weight='3' tree='@tree' />"),
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
    file_logger = ET.SubElement(
        run,
        "logger",
        {"logEvery": str(log_every), "fileName": log_path.name},
    )
    ET.SubElement(file_logger, "model", {"idref": "posterior"})
    for log_id in ["posterior", "prior", "likelihood"]:
        ET.SubElement(file_logger, "log", {"idref": log_id})
    ET.SubElement(
        file_logger,
        "log",
        {"spec": "beast.base.evolution.tree.TreeHeightLogger", "tree": "@tree"},
    )
    for element in [*tree_logger_elements, *site_logger_elements, *clock_logger_elements]:
        file_logger.append(element)
    for calibration in beast_calibrations:
        ET.SubElement(file_logger, "log", {"idref": _xml_identifier(calibration.calibration_id)})

    tree_logger = ET.SubElement(
        run,
        "logger",
        {"logEvery": str(log_every), "fileName": tree_log_path.name},
    )
    ET.SubElement(tree_logger, "log", {"idref": "tree"})

    screen_logger = ET.SubElement(
        run,
        "logger",
        {"logEvery": str(max(log_every, 10000))},
    )
    ET.SubElement(screen_logger, "model", {"idref": "posterior"})
    for log_id in ["posterior", "prior", "likelihood"]:
        ET.SubElement(screen_logger, "log", {"idref": log_id})
    ET.SubElement(
        screen_logger,
        "log",
        {"spec": "beast.base.evolution.tree.TreeHeightLogger", "tree": "@tree"},
    )
    for element in [*tree_logger_elements, *site_logger_elements, *clock_logger_elements]:
        cloned = ET.fromstring(ET.tostring(element, encoding="unicode"))
        if cloned.get("id") == "rateStatistic":
            cloned.set("id", "rateStatistic.screen")
        screen_logger.append(cloned)

    xml_tree = ET.ElementTree(root)
    ET.indent(xml_tree, space="    ")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    xml_tree.write(output_path, encoding="utf-8", xml_declaration=True)
    output_path.write_text(output_path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
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


def parse_beast_log(path: Path) -> BeastLogReport:
    """Parse a BEAST-style log table into deterministic numeric rows."""
    with path.open(encoding="utf-8", newline="") as handle:
        filtered_lines = [
            line
            for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        ]
    reader = csv.DictReader(filtered_lines, delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError(f"BEAST log contains no header: {path}")
    state_field = _beast_state_field(reader.fieldnames)
    if state_field is None:
        raise ValueError(f"BEAST log lacks a state column: {path}")
    columns = [field for field in reader.fieldnames if field and field != state_field]
    rows: list[BeastLogRow] = []
    for row in reader:
        values = {
            column: float(row[column])
            for column in columns
            if row.get(column) not in {None, ""}
        }
        rows.append(BeastLogRow(state=int(float(row[state_field])), values=values))
    if not rows:
        raise ValueError(f"BEAST log contains no sampled rows: {path}")
    return BeastLogReport(path=path, row_count=len(rows), columns=columns, rows=rows)


def summarize_beast_log(
    path: Path,
    *,
    burnin_fraction: float = 0.0,
) -> BeastLogSummaryReport:
    """Summarize a BEAST posterior log after discarding an optional burn-in fraction."""
    report = parse_beast_log(path)
    burnin_row_count, kept_rows = _split_beast_log_rows(
        report, burnin_fraction=burnin_fraction
    )
    convergence = summarize_trace_convergence(
        path=path,
        rows=[row.values for row in kept_rows],
        columns=report.columns,
        ess_threshold=0.0,
        mean_shift_threshold=float("inf"),
    )
    parameter_summaries = [
        BeastLogParameterSummary(
            parameter=summary.parameter,
            parameter_category=_classify_beast_parameter(summary.parameter),
            sample_count=summary.sample_count,
            effective_sample_size=summary.effective_sample_size,
            mean=summary.mean,
            minimum=summary.minimum,
            maximum=summary.maximum,
            first_half_mean=summary.first_half_mean,
            second_half_mean=summary.second_half_mean,
            standardized_mean_shift=summary.standardized_mean_shift,
        )
        for summary in convergence.series
    ]
    return BeastLogSummaryReport(
        path=path,
        burnin_fraction=burnin_fraction,
        burnin_row_count=burnin_row_count,
        kept_row_count=len(kept_rows),
        first_kept_state=kept_rows[0].state,
        last_kept_state=kept_rows[-1].state,
        posterior_parameters=_summary_parameters_by_category(
            parameter_summaries, category="posterior"
        ),
        likelihood_parameters=_summary_parameters_by_category(
            parameter_summaries, category="likelihood"
        ),
        prior_parameters=_summary_parameters_by_category(
            parameter_summaries, category="prior"
        ),
        clock_parameters=_summary_parameters_by_category(
            parameter_summaries, category="clock"
        ),
        tree_parameters=_summary_parameters_by_category(
            parameter_summaries, category="tree"
        ),
        other_parameters=_summary_parameters_by_category(
            parameter_summaries, category="other"
        ),
        parameter_summaries=parameter_summaries,
    )


def parse_beast_posterior_tree_samples(
    path: Path,
    *,
    burnin_fraction: float = 0.0,
) -> BeastPosteriorTreeSetReport:
    """Parse a BEAST posterior tree set into state-tagged normalized trees."""
    text = path.read_text(encoding="utf-8")
    entries = _extract_beast_tree_entries(text)
    if not entries:
        raise EngineWorkflowError(f"BEAST posterior tree file contains no trees: {path}")
    burnin_tree_count, kept_entries = _split_beast_tree_entries(
        entries, burnin_fraction=burnin_fraction, path=path
    )
    translation = _parse_nexus_translate_map(text)
    samples: list[BeastPosteriorTreeSample] = []
    trees: list[PhyloTree] = []
    for tree_name, tree_text in kept_entries:
        newick, tree, rooted = _parse_beast_tree_text(
            tree_text, translation=translation
        )
        samples.append(
            BeastPosteriorTreeSample(
                tree_name=tree_name,
                state=_parse_beast_tree_state(tree_name),
                rooted=rooted if rooted is not None else True,
                tip_names=tree.tip_names,
                newick=newick,
            )
        )
        trees.append(tree)
    clades = _summarize_beast_clades(trees)
    sampled_states = [
        state for state in (sample.state for sample in samples) if state is not None
    ]
    rooted_tree_count = sum(1 for sample in samples if sample.rooted)
    return BeastPosteriorTreeSetReport(
        path=path,
        burnin_fraction=burnin_fraction,
        total_tree_count=len(entries),
        burnin_tree_count=burnin_tree_count,
        kept_tree_count=len(samples),
        rooted_tree_count=rooted_tree_count,
        sampled_states=sampled_states,
        tip_names=sorted(samples[0].tip_names),
        clades=clades,
        trees=samples,
    )


def write_beast_posterior_tree_set(
    path: Path, report: BeastPosteriorTreeSetReport
) -> Path:
    """Write a normalized Newick tree set from parsed BEAST posterior samples."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(f"{sample.newick}\n" for sample in report.trees),
        encoding="utf-8",
    )
    return path


def _annotate_consensus_tree_with_posterior_probabilities(tree: PhyloTree) -> None:
    for node in tree.iter_nodes():
        if node is tree.root or node.is_leaf() or node.name is None:
            continue
        node.name = format(float(node.name) / 100.0, ".15g")


def summarize_beast_posterior_trees(
    tree_set_path: Path,
    *,
    burnin_fraction: float = 0.25,
) -> tuple[PhyloTree, BeastPosteriorConsensusReport]:
    """Summarize BEAST posterior trees into a majority-rule consensus tree."""
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(
            f"burnin_fraction must be between 0 and 1, got {burnin_fraction}"
        )
    tree_set_report = parse_beast_posterior_tree_samples(
        tree_set_path,
        burnin_fraction=burnin_fraction,
    )
    if not tree_set_report.trees:
        raise EngineWorkflowError(
            f"BEAST posterior tree file is empty after burn-in filtering: {tree_set_path}"
        )
    retained_tree_set_path = tree_set_path.with_suffix(".postburnin.nwk")
    write_beast_posterior_tree_set(retained_tree_set_path, tree_set_report)
    summary = load_tree_set(retained_tree_set_path)
    consensus_tree, consensus = compute_consensus_tree(retained_tree_set_path)
    _annotate_consensus_tree_with_posterior_probabilities(consensus_tree)
    consensus_newick = dumps_newick(consensus_tree)
    clade_frequencies = compute_clade_frequency_table(retained_tree_set_path)
    posterior_probabilities = sorted(
        float(node.name)
        for node in consensus_tree.iter_nodes()
        if node is not consensus_tree.root and not node.is_leaf() and node.name is not None
    )
    return consensus_tree, BeastPosteriorConsensusReport(
        source_path=tree_set_path,
        retained_tree_set_path=retained_tree_set_path,
        burnin_fraction=tree_set_report.burnin_fraction,
        total_tree_count=tree_set_report.total_tree_count,
        burnin_tree_count=tree_set_report.burnin_tree_count,
        kept_tree_count=tree_set_report.kept_tree_count,
        rooted_topology_count=summary.rooted_topology_count,
        shared_taxa=summary.shared_taxa,
        consensus_newick=consensus_newick,
        clade_frequency_count=len(clade_frequencies.clade_frequencies),
        annotated_node_count=len(posterior_probabilities),
        minimum_posterior_probability=(
            None if not posterior_probabilities else posterior_probabilities[0]
        ),
        maximum_posterior_probability=(
            None if not posterior_probabilities else posterior_probabilities[-1]
        ),
    )


def write_beast_log_summary_table(path: Path, report: BeastLogSummaryReport) -> Path:
    """Write a reviewer-facing TSV summary of one BEAST posterior log."""
    rows = [
        {
            "parameter_category": summary.parameter_category,
            "parameter": summary.parameter,
            "sample_count": str(summary.sample_count),
            "effective_sample_size": format(summary.effective_sample_size, ".15g"),
            "mean": format(summary.mean, ".15g"),
            "minimum": format(summary.minimum, ".15g"),
            "maximum": format(summary.maximum, ".15g"),
            "first_half_mean": format(summary.first_half_mean, ".15g"),
            "second_half_mean": format(summary.second_half_mean, ".15g"),
            "standardized_mean_shift": format(
                summary.standardized_mean_shift, ".15g"
            ),
            "burnin_fraction": format(report.burnin_fraction, ".15g"),
            "burnin_row_count": str(report.burnin_row_count),
            "kept_row_count": str(report.kept_row_count),
            "first_kept_state": str(report.first_kept_state),
            "last_kept_state": str(report.last_kept_state),
        }
        for summary in report.parameter_summaries
    ]
    return write_taxon_rows(
        path,
        columns=[
            "parameter_category",
            "parameter",
            "sample_count",
            "effective_sample_size",
            "mean",
            "minimum",
            "maximum",
            "first_half_mean",
            "second_half_mean",
            "standardized_mean_shift",
            "burnin_fraction",
            "burnin_row_count",
            "kept_row_count",
            "first_kept_state",
            "last_kept_state",
        ],
        rows=rows,
    )


def validate_beast_posterior_log(
    path: Path,
    *,
    required_columns: tuple[str, ...] = ("posterior", "likelihood"),
) -> BeastPosteriorLogValidationReport:
    """Validate that a BEAST posterior log contains the required fields and monotonic sampled states."""
    issues: list[BeastLogValidationIssue] = []
    with path.open(encoding="utf-8", newline="") as handle:
        filtered_lines = [
            line
            for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        ]
    reader = csv.DictReader(filtered_lines, delimiter="\t")
    if reader.fieldnames is None:
        issues.append(
            BeastLogValidationIssue(
                code="missing-header", message="BEAST log contains no header row"
            )
        )
        return BeastPosteriorLogValidationReport(
            path=path,
            row_count=0,
            state_count=0,
            required_columns=list(required_columns),
            observed_columns=[],
            missing_columns=list(required_columns),
            issues=issues,
            valid=False,
        )
    state_field = _beast_state_field(reader.fieldnames)
    observed_columns = [field for field in reader.fieldnames if field]
    if state_field is None:
        issues.append(
            BeastLogValidationIssue(
                code="missing-state-column", message="BEAST log lacks a state column"
            )
        )
    missing_columns = [
        column for column in required_columns if column not in observed_columns
    ]
    for column in missing_columns:
        issues.append(
            BeastLogValidationIssue(
                code="missing-required-column",
                message=f"missing required BEAST log column '{column}'",
                column=column,
            )
        )
    previous_state: int | None = None
    row_count = 0
    state_count = 0
    for row_index, row in enumerate(reader, start=2):
        row_count += 1
        if state_field is not None:
            raw_state = row.get(state_field, "")
            if not raw_state:
                issues.append(
                    BeastLogValidationIssue(
                        code="missing-state-value",
                        message="row is missing a sampled state",
                        row=row_index,
                        column=state_field,
                    )
                )
            else:
                try:
                    state_value = int(float(raw_state))
                    state_count += 1
                    if previous_state is not None and state_value <= previous_state:
                        issues.append(
                            BeastLogValidationIssue(
                                code="nonmonotonic-state",
                                message="sampled states must increase strictly through the log",
                                row=row_index,
                                column=state_field,
                            )
                        )
                    previous_state = state_value
                except ValueError:
                    issues.append(
                        BeastLogValidationIssue(
                            code="invalid-state-value",
                            message="sampled state must be numeric",
                            row=row_index,
                            column=state_field,
                        )
                    )
        for column in observed_columns:
            if column == state_field:
                continue
            raw_value = row.get(column, "")
            if raw_value in {None, ""}:
                issues.append(
                    BeastLogValidationIssue(
                        code="missing-parameter-value",
                        message=f"missing sampled value for '{column}'",
                        row=row_index,
                        column=column,
                    )
                )
                continue
            try:
                float(raw_value)
            except ValueError:
                issues.append(
                    BeastLogValidationIssue(
                        code="invalid-parameter-value",
                        message=f"sampled value for '{column}' must be numeric",
                        row=row_index,
                        column=column,
                    )
                )
    if row_count == 0:
        issues.append(
            BeastLogValidationIssue(
                code="missing-rows", message="BEAST log contains no sampled rows"
            )
        )
    return BeastPosteriorLogValidationReport(
        path=path,
        row_count=row_count,
        state_count=state_count,
        required_columns=list(required_columns),
        observed_columns=observed_columns,
        missing_columns=missing_columns,
        issues=issues,
        valid=not issues,
    )


def assess_beast_burnin_sensitivity(
    posterior_tree_path: Path,
    *,
    log_path: Path | None = None,
    burnin_fractions: tuple[float, ...] = (0.1, 0.25, 0.5),
) -> BeastBurninSensitivityReport:
    """Compare posterior summaries across multiple BEAST burn-in fractions."""
    if not burnin_fractions:
        raise ValueError("burnin_fractions must contain at least one value")
    log_report = parse_beast_log(log_path) if log_path is not None else None
    ordered_fractions = tuple(sorted(dict.fromkeys(burnin_fractions)))
    slices: list[BeastBurninSensitivitySlice] = []
    previous_newick: str | None = None
    changed_mcc_count = 0
    for fraction in ordered_fractions:
        _, mcc_report = summarize_maximum_clade_credibility_tree(
            posterior_tree_path,
            burnin_fraction=fraction,
        )
        posterior_mean = None
        likelihood_mean = None
        tree_height_mean = None
        if log_report is not None:
            burnin_row_count = int(log_report.row_count * fraction)
            kept_rows = log_report.rows[burnin_row_count:]
            if kept_rows:
                posterior_values = [
                    row.values["posterior"]
                    for row in kept_rows
                    if "posterior" in row.values
                ]
                likelihood_values = [
                    row.values["likelihood"]
                    for row in kept_rows
                    if "likelihood" in row.values
                ]
                tree_height_values = [
                    row.values["treeHeight"]
                    for row in kept_rows
                    if "treeHeight" in row.values
                ]
                posterior_mean = _mean_or_none(posterior_values)
                likelihood_mean = _mean_or_none(likelihood_values)
                tree_height_mean = _mean_or_none(tree_height_values)
        slices.append(
            BeastBurninSensitivitySlice(
                burnin_fraction=fraction,
                burnin_tree_count=mcc_report.burnin_tree_count,
                kept_tree_count=mcc_report.kept_tree_count,
                rooted_topology_count=mcc_report.rooted_topology_count,
                selected_tree_index=mcc_report.selected_tree_index,
                clade_credibility_score=mcc_report.clade_credibility_score,
                posterior_mean=posterior_mean,
                likelihood_mean=likelihood_mean,
                tree_height_mean=tree_height_mean,
            )
        )
        if previous_newick is not None and previous_newick != mcc_report.mcc_newick:
            changed_mcc_count += 1
        previous_newick = mcc_report.mcc_newick
    warnings: list[str] = []
    if changed_mcc_count:
        warnings.append(
            "maximum clade credibility topology changes across tested burn-in fractions"
        )
    if len({row.rooted_topology_count for row in slices}) > 1:
        warnings.append(
            "rooted topology diversity changes across tested burn-in fractions"
        )
    return BeastBurninSensitivityReport(
        posterior_tree_path=posterior_tree_path,
        log_path=log_path,
        slices=slices,
        changed_mcc_count=changed_mcc_count,
        warnings=warnings,
    )


def assess_beast_chain_mixing(
    log_paths: list[Path],
    *,
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
    cross_chain_mean_shift_threshold: float = 0.75,
    stuck_parameter_span_threshold: float = 1e-9,
) -> BeastChainMixingReport:
    """Flag low ESS, mean drift, stuck parameters, and inconsistent chain means across BEAST logs."""
    if not log_paths:
        raise ValueError("assess_beast_chain_mixing requires at least one log path")
    chain_summaries = [
        assess_beast_convergence(
            path,
            ess_threshold=ess_threshold,
            mean_shift_threshold=mean_shift_threshold,
        )
        for path in log_paths
    ]
    issues: list[BeastChainMixingIssue] = []
    for summary in chain_summaries:
        for warning in summary.warnings:
            issues.append(
                BeastChainMixingIssue(
                    path=summary.path,
                    parameter=str(warning["parameter"]),
                    code=str(warning["code"]),
                    message=str(warning["message"]),
                    observed_value=float(warning["observed_value"]),
                    threshold=float(warning["threshold"]),
                )
            )
        for parameter_summary in summary.parameter_summaries:
            span = float(parameter_summary["maximum"]) - float(
                parameter_summary["minimum"]
            )
            if span <= stuck_parameter_span_threshold:
                issues.append(
                    BeastChainMixingIssue(
                        path=summary.path,
                        parameter=str(parameter_summary["parameter"]),
                        code="stuck-parameter",
                        message="parameter shows effectively no movement across the sampled chain",
                        observed_value=span,
                        threshold=stuck_parameter_span_threshold,
                    )
                )
    parameter_to_means: dict[str, list[tuple[Path, float]]] = {}
    for summary in chain_summaries:
        for parameter_summary in summary.parameter_summaries:
            parameter_to_means.setdefault(
                str(parameter_summary["parameter"]), []
            ).append((summary.path, float(parameter_summary["mean"])))
    for parameter, chain_means in parameter_to_means.items():
        if len(chain_means) < 2:
            continue
        mean_values = [value for _, value in chain_means]
        span = max(mean_values) - min(mean_values)
        if span > cross_chain_mean_shift_threshold:
            issues.append(
                BeastChainMixingIssue(
                    path=None,
                    parameter=parameter,
                    code="inconsistent-chains",
                    message="independent chains disagree more than the allowed mean-shift threshold",
                    observed_value=span,
                    threshold=cross_chain_mean_shift_threshold,
                )
            )
    return BeastChainMixingReport(
        log_paths=log_paths,
        chain_count=len(log_paths),
        converged=not issues,
        issues=issues,
        chain_summaries=chain_summaries,
    )


def assess_beast_convergence(
    path: Path,
    *,
    burnin_fraction: float = 0.0,
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
) -> BeastConvergenceReport:
    """Flag low-ESS or unstable BEAST trace parameters."""
    report = parse_beast_log(path)
    burnin_row_count, kept_rows = _split_beast_log_rows(
        report, burnin_fraction=burnin_fraction
    )
    convergence = summarize_trace_convergence(
        path=path,
        rows=[row.values for row in kept_rows],
        columns=report.columns,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
    )
    return _build_beast_convergence_report(
        convergence,
        burnin_fraction=burnin_fraction,
        burnin_row_count=burnin_row_count,
    )


def _build_beast_convergence_report(
    convergence: TraceConvergenceReport,
    *,
    burnin_fraction: float,
    burnin_row_count: int,
) -> BeastConvergenceReport:
    return BeastConvergenceReport(
        path=convergence.path,
        burnin_fraction=burnin_fraction,
        burnin_row_count=burnin_row_count,
        sample_count=convergence.sample_count,
        converged=convergence.converged,
        ess_threshold=convergence.ess_threshold,
        mean_shift_threshold=convergence.mean_shift_threshold,
        warnings=[
            {
                "parameter": warning.parameter,
                "code": warning.code,
                "message": warning.message,
                "observed_value": warning.observed_value,
                "threshold": warning.threshold,
            }
            for warning in convergence.warnings
        ],
        parameter_summaries=[
            {
                "parameter": summary.parameter,
                "sample_count": summary.sample_count,
                "effective_sample_size": summary.effective_sample_size,
                "mean": summary.mean,
                "minimum": summary.minimum,
                "maximum": summary.maximum,
                "first_half_mean": summary.first_half_mean,
                "second_half_mean": summary.second_half_mean,
                "standardized_mean_shift": summary.standardized_mean_shift,
            }
            for summary in convergence.series
        ],
    )


def _mean_or_none(values: list[float]) -> float | None:
    if not values:
        return None
    return round(sum(values) / len(values), 6)


def _split_beast_log_rows(
    report: BeastLogReport,
    *,
    burnin_fraction: float,
) -> tuple[int, list[BeastLogRow]]:
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(
            f"burnin_fraction must be between 0 and 1, got {burnin_fraction}"
        )
    burnin_row_count = int(report.row_count * burnin_fraction)
    kept_rows = report.rows[burnin_row_count:]
    if not kept_rows:
        raise ValueError(
            f"BEAST log contains no sampled rows after burn-in filtering: {report.path}"
        )
    return burnin_row_count, kept_rows


def _classify_beast_parameter(parameter: str) -> str:
    normalized = parameter.lower()
    if normalized == "posterior" or normalized.endswith(".posterior"):
        return "posterior"
    if normalized == "likelihood" or "likelihood" in normalized:
        return "likelihood"
    if normalized == "prior" or "prior" in normalized:
        return "prior"
    if any(
        token in normalized
        for token in (
            "clock",
            "clockrate",
            "ucld",
            "branchrate",
            "mutationrate",
            "meanrate",
        )
    ):
        return "clock"
    if any(
        token in normalized
        for token in (
            "tree",
            "birthrate",
            "deathrate",
            "popsize",
            "coalescent",
            "tmrca",
            "origin",
        )
    ):
        return "tree"
    return "other"


def _summary_parameters_by_category(
    parameter_summaries: list[BeastLogParameterSummary],
    *,
    category: str,
) -> list[str]:
    return [
        summary.parameter
        for summary in parameter_summaries
        if summary.parameter_category == category
    ]


def _extract_beast_tree_entries(text: str) -> list[tuple[str, str]]:
    return [
        (match.group(1), match.group(2).strip())
        for match in _BEAST_TREE_PATTERN.finditer(text)
    ]


def _split_beast_tree_entries(
    entries: list[tuple[str, str]],
    *,
    burnin_fraction: float,
    path: Path,
) -> tuple[int, list[tuple[str, str]]]:
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(
            f"burnin_fraction must be between 0 and 1, got {burnin_fraction}"
        )
    burnin_tree_count = int(len(entries) * burnin_fraction)
    kept_entries = entries[burnin_tree_count:]
    if not kept_entries:
        raise EngineWorkflowError(
            f"BEAST posterior tree file is empty after burn-in filtering: {path}"
        )
    return burnin_tree_count, kept_entries


def _parse_beast_tree_text(
    tree_text: str, *, translation: dict[str, str]
) -> tuple[str, PhyloTree, bool | None]:
    rooted = _detect_nexus_rooted_flag(tree_text)
    stripped = _strip_square_bracket_comments(tree_text).strip()
    translated = _translate_nexus_tip_labels(stripped, translation)
    tree = loads_biophylo(f"{translated};", source_format="newick")
    return dumps_newick(tree), tree, rooted


def _detect_nexus_rooted_flag(tree_text: str) -> bool | None:
    prefix = tree_text.lstrip()
    if prefix.startswith("[&R]"):
        return True
    if prefix.startswith("[&U]"):
        return False
    return None


def _parse_beast_tree_state(tree_name: str) -> int | None:
    match = _BEAST_TREE_STATE_PATTERN.search(tree_name)
    return None if match is None else int(match.group(1))


def _summarize_beast_clades(trees: list[PhyloTree]) -> list[BeastPosteriorClade]:
    if not trees:
        return []
    taxa_sets = {frozenset(tree.tip_names) for tree in trees}
    if len(taxa_sets) != 1:
        raise InvalidAlignmentError(
            "BEAST posterior tree summaries require all trees to share the exact same taxon set"
        )
    shared_taxa = set(next(iter(taxa_sets)))
    counts: dict[frozenset[str], int] = {}
    for tree in trees:
        for clade in _informative_tree_clades(tree, shared_taxa):
            counts[clade] = counts.get(clade, 0) + 1
    return [
        BeastPosteriorClade(
            clade="|".join(sorted(clade)),
            tree_count=tree_count,
            frequency=round(tree_count / len(trees), 15),
        )
        for clade, tree_count in sorted(
            counts.items(),
            key=lambda item: (-item[1], sorted(item[0])),
        )
    ]


def _informative_tree_clades(
    tree: PhyloTree, shared_taxa: set[str]
) -> list[frozenset[str]]:
    clades: list[frozenset[str]] = []
    for node in tree.iter_nodes():
        taxa = frozenset(descendant_taxa(node))
        if 1 < len(taxa) < len(shared_taxa):
            clades.append(taxa)
    return clades


def _parse_nexus_translate_map(text: str) -> dict[str, str]:
    translate_match = re.search(
        r"translate\s+(.+?);", text, flags=re.IGNORECASE | re.DOTALL
    )
    if translate_match is None:
        return {}
    mapping: dict[str, str] = {}
    for raw_line in translate_match.group(1).splitlines():
        line = raw_line.strip().rstrip(",")
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        source, label = parts
        mapping[source] = label.rstrip(",")
    return mapping


def _translate_nexus_tip_labels(newick: str, mapping: dict[str, str]) -> str:
    if not mapping:
        return newick

    def replace(match: re.Match[str]) -> str:
        token = match.group(1)
        translated = mapping.get(token, token)
        return match.group(0).replace(token, translated)

    return re.sub(r"(?<=[(,])\s*([A-Za-z0-9_.-]+)(?=\s*[:),])", replace, newick)


def _strip_square_bracket_comments(text: str) -> str:
    result: list[str] = []
    depth = 0
    for char in text:
        if char == "[":
            depth += 1
            continue
        if char == "]" and depth:
            depth -= 1
            continue
        if depth == 0:
            result.append(char)
    return "".join(result)


def _beast_state_field(fieldnames: list[str]) -> str | None:
    for candidate in ("state", "State", "Sample", "sample"):
        if candidate in fieldnames:
            return candidate
    return None
