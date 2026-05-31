from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.models import (
    NestedLikelihoodRatioModelFit,
    NestedLikelihoodRatioTestReport,
)
from bijux_phylogenetics.phylo.likelihood.substitution_parameters import (
    optimize_nucleotide_substitution_parameters,
    validate_nucleotide_substitution_optimization_model,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_BOUNDARY_CAVEAT = "no boundary parameter in declared nesting"
_P_VALUE_METHOD = "chi-square approximation"


@dataclass(frozen=True, slots=True)
class _DeclaredNestedModelPair:
    null_model_name: str
    alternative_model_name: str
    degrees_of_freedom: int
    boundary_caveat: str


@dataclass(frozen=True, slots=True)
class _LikelihoodRatioFitContext:
    fit: NestedLikelihoodRatioModelFit
    taxa: list[str]
    site_count: int
    pattern_count: int
    tree_newick: str


_DECLARED_NESTED_MODEL_PAIRS = {
    ("JC69", "K80"): _DeclaredNestedModelPair(
        null_model_name="JC69",
        alternative_model_name="K80",
        degrees_of_freedom=1,
        boundary_caveat=_BOUNDARY_CAVEAT,
    ),
    ("JC69", "F81"): _DeclaredNestedModelPair(
        null_model_name="JC69",
        alternative_model_name="F81",
        degrees_of_freedom=3,
        boundary_caveat=_BOUNDARY_CAVEAT,
    ),
    ("JC69", "HKY85"): _DeclaredNestedModelPair(
        null_model_name="JC69",
        alternative_model_name="HKY85",
        degrees_of_freedom=4,
        boundary_caveat=_BOUNDARY_CAVEAT,
    ),
    ("JC69", "GTR"): _DeclaredNestedModelPair(
        null_model_name="JC69",
        alternative_model_name="GTR",
        degrees_of_freedom=8,
        boundary_caveat=_BOUNDARY_CAVEAT,
    ),
    ("K80", "HKY85"): _DeclaredNestedModelPair(
        null_model_name="K80",
        alternative_model_name="HKY85",
        degrees_of_freedom=3,
        boundary_caveat=_BOUNDARY_CAVEAT,
    ),
    ("K80", "GTR"): _DeclaredNestedModelPair(
        null_model_name="K80",
        alternative_model_name="GTR",
        degrees_of_freedom=7,
        boundary_caveat=_BOUNDARY_CAVEAT,
    ),
    ("F81", "HKY85"): _DeclaredNestedModelPair(
        null_model_name="F81",
        alternative_model_name="HKY85",
        degrees_of_freedom=1,
        boundary_caveat=_BOUNDARY_CAVEAT,
    ),
    ("F81", "GTR"): _DeclaredNestedModelPair(
        null_model_name="F81",
        alternative_model_name="GTR",
        degrees_of_freedom=5,
        boundary_caveat=_BOUNDARY_CAVEAT,
    ),
    ("HKY85", "GTR"): _DeclaredNestedModelPair(
        null_model_name="HKY85",
        alternative_model_name="GTR",
        degrees_of_freedom=4,
        boundary_caveat=_BOUNDARY_CAVEAT,
    ),
}


def list_declared_nucleotide_likelihood_ratio_pairs() -> list[tuple[str, str]]:
    """Return the declared null-alternative nucleotide model pairs for LRT."""
    return list(_DECLARED_NESTED_MODEL_PAIRS)


def validate_declared_nucleotide_likelihood_ratio_pair(
    null_model_name: str,
    alternative_model_name: str,
) -> tuple[str, str]:
    """Validate one declared nested null and alternative nucleotide model pair."""
    normalized_null_model_name = validate_nucleotide_substitution_optimization_model(
        null_model_name
    ).upper()
    normalized_alternative_model_name = (
        validate_nucleotide_substitution_optimization_model(
            alternative_model_name
        ).upper()
    )
    key = (normalized_null_model_name, normalized_alternative_model_name)
    if key not in _DECLARED_NESTED_MODEL_PAIRS:
        declared_pairs = ", ".join(
            f"{null}->{alternative}"
            for null, alternative in list_declared_nucleotide_likelihood_ratio_pairs()
        )
        raise ValueError(
            "nested likelihood-ratio tests are only declared for these nucleotide "
            f"model pairs: {declared_pairs}"
        )
    return key


def evaluate_nucleotide_nested_likelihood_ratio_test(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    null_model_name: str,
    alternative_model_name: str,
) -> NestedLikelihoodRatioTestReport:
    """Evaluate one declared nested fixed-topology nucleotide likelihood-ratio test."""
    pair_key = validate_declared_nucleotide_likelihood_ratio_pair(
        null_model_name,
        alternative_model_name,
    )
    pair_definition = _DECLARED_NESTED_MODEL_PAIRS[pair_key]
    null_context = _fit_nucleotide_model_for_likelihood_ratio_test(
        tree,
        records,
        model_name=pair_definition.null_model_name,
    )
    alternative_context = _fit_nucleotide_model_for_likelihood_ratio_test(
        tree,
        records,
        model_name=pair_definition.alternative_model_name,
    )
    null_fit = null_context.fit
    alternative_fit = alternative_context.fit
    raw_statistic = 2.0 * (alternative_fit.log_likelihood - null_fit.log_likelihood)
    warnings: list[str] = []
    if raw_statistic < 0.0 and not math.isclose(raw_statistic, 0.0, abs_tol=1e-12):
        warnings.append(
            "alternative log likelihood fell below the null fit; statistic clipped to zero"
        )
    statistic = max(raw_statistic, 0.0)
    if (
        null_context.taxa != alternative_context.taxa
        or null_context.site_count != alternative_context.site_count
        or null_context.pattern_count != alternative_context.pattern_count
        or null_context.tree_newick != alternative_context.tree_newick
    ):
        raise ValueError(
            "nested likelihood-ratio test requires null and alternative fits on the same tree and alignment surface"
        )
    return NestedLikelihoodRatioTestReport(
        taxa=null_context.taxa,
        site_count=null_context.site_count,
        pattern_count=null_context.pattern_count,
        tree_newick=null_context.tree_newick,
        null_fit=null_fit,
        alternative_fit=alternative_fit,
        likelihood_ratio_statistic=statistic,
        degrees_of_freedom=pair_definition.degrees_of_freedom,
        p_value=_chi_square_survival(statistic, pair_definition.degrees_of_freedom),
        p_value_method=_P_VALUE_METHOD,
        boundary_caveat=pair_definition.boundary_caveat,
        warnings=warnings,
    )


def evaluate_nucleotide_nested_likelihood_ratio_test_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    null_model_name: str,
    alternative_model_name: str,
) -> NestedLikelihoodRatioTestReport:
    """Evaluate one declared nested nucleotide likelihood-ratio test from file paths."""
    return evaluate_nucleotide_nested_likelihood_ratio_test(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        null_model_name=null_model_name,
        alternative_model_name=alternative_model_name,
    )


def _fit_nucleotide_model_for_likelihood_ratio_test(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
) -> _LikelihoodRatioFitContext:
    optimization_report = optimize_nucleotide_substitution_parameters(
        tree,
        records,
        model_name=model_name,
    )
    parameter_values: dict[str, float] = {}
    if optimization_report.base_frequency_a is not None:
        parameter_values["base_frequency_a"] = optimization_report.base_frequency_a
    if optimization_report.base_frequency_c is not None:
        parameter_values["base_frequency_c"] = optimization_report.base_frequency_c
    if optimization_report.base_frequency_g is not None:
        parameter_values["base_frequency_g"] = optimization_report.base_frequency_g
    if optimization_report.base_frequency_t is not None:
        parameter_values["base_frequency_t"] = optimization_report.base_frequency_t
    parameter_values.update(optimization_report.fixed_parameter_values)
    parameter_values.update(
        {
            row.parameter_name: row.optimized_value
            for row in optimization_report.parameter_rows
        }
    )
    return _LikelihoodRatioFitContext(
        fit=NestedLikelihoodRatioModelFit(
            model_name=optimization_report.model_name,
            parameter_count=optimization_report.parameter_count,
            log_likelihood=optimization_report.optimized_log_likelihood,
            aic=optimization_report.optimized_aic,
            parameter_values=parameter_values,
            warnings=optimization_report.warnings,
        ),
        taxa=optimization_report.taxa,
        site_count=optimization_report.site_count,
        pattern_count=optimization_report.pattern_count,
        tree_newick=optimization_report.tree_newick,
    )


def _chi_square_survival(statistic: float, degrees_of_freedom: int) -> float:
    if statistic <= 0.0 or degrees_of_freedom <= 0:
        return 1.0
    if degrees_of_freedom == 1:
        return math.erfc(math.sqrt(statistic / 2.0))
    z_score = (
        ((statistic / degrees_of_freedom) ** (1.0 / 3.0))
        - (1.0 - (2.0 / (9.0 * degrees_of_freedom)))
    ) / math.sqrt(2.0 / (9.0 * degrees_of_freedom))
    return 0.5 * math.erfc(z_score / math.sqrt(2.0))
