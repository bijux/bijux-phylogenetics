"""Comparative regression workflows and diagnostics."""

from .logistic import (
    PhylogeneticLogisticCoefficient,
    PhylogeneticLogisticFittedRow,
    PhylogeneticLogisticReport,
    PhylogeneticLogisticWarning,
    summarize_phylogenetic_logistic,
    write_phylogenetic_logistic_coefficient_table,
    write_phylogenetic_logistic_excluded_taxa_table,
    write_phylogenetic_logistic_fitted_table,
)
from .residuals import (
    PhylogeneticResidualCoefficientRow,
    PhylogeneticResidualExclusion,
    PhylogeneticResidualReport,
    PhylogeneticResidualTaxonRow,
    summarize_phylogenetic_residuals,
    write_phylogenetic_residual_coefficient_table,
    write_phylogenetic_residual_exclusion_table,
    write_phylogenetic_residual_summary_table,
    write_phylogenetic_residual_taxon_table,
)

__all__ = [
    "PhylogeneticLogisticCoefficient",
    "PhylogeneticLogisticFittedRow",
    "PhylogeneticLogisticReport",
    "PhylogeneticLogisticWarning",
    "PhylogeneticResidualCoefficientRow",
    "PhylogeneticResidualExclusion",
    "PhylogeneticResidualReport",
    "PhylogeneticResidualTaxonRow",
    "summarize_phylogenetic_logistic",
    "summarize_phylogenetic_residuals",
    "write_phylogenetic_logistic_coefficient_table",
    "write_phylogenetic_logistic_excluded_taxa_table",
    "write_phylogenetic_logistic_fitted_table",
    "write_phylogenetic_residual_coefficient_table",
    "write_phylogenetic_residual_exclusion_table",
    "write_phylogenetic_residual_summary_table",
    "write_phylogenetic_residual_taxon_table",
]
