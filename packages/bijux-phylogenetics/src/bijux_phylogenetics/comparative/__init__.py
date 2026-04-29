"""Comparative-analysis methods and helpers."""

from .common import ComparativeDataset, ComparativeReadinessReport, NumericTraitSummary, summarize_numeric_trait, summarize_numeric_trait_readiness
from .pgls import PGLSCoefficient, PGLSInputReport, PGLSPredictorClassification, PGLSResult, inspect_pgls_inputs, run_pgls
from .signal import (
    BlombergKReport,
    IndependentContrast,
    IndependentContrastReport,
    PagelLambdaReport,
    PhylogeneticSignalTestReport,
    compute_blombergs_k,
    compute_phylogenetic_independent_contrasts,
    compute_phylogenetic_signal_test,
    estimate_pagels_lambda,
)

__all__ = [
    "BlombergKReport",
    "ComparativeDataset",
    "ComparativeReadinessReport",
    "IndependentContrast",
    "IndependentContrastReport",
    "NumericTraitSummary",
    "PGLSCoefficient",
    "PGLSInputReport",
    "PGLSPredictorClassification",
    "PGLSResult",
    "PagelLambdaReport",
    "PhylogeneticSignalTestReport",
    "compute_blombergs_k",
    "compute_phylogenetic_independent_contrasts",
    "compute_phylogenetic_signal_test",
    "estimate_pagels_lambda",
    "inspect_pgls_inputs",
    "run_pgls",
    "summarize_numeric_trait",
    "summarize_numeric_trait_readiness",
]
