"""Phylogenetic generalized least-squares package surface."""

from __future__ import annotations

from .design import (
    build_pgls_model_matrix as build_pgls_model_matrix,
)
from .design import (
    inspect_pgls_inputs as inspect_pgls_inputs,
)
from .design import (
    write_pgls_model_matrix_table as write_pgls_model_matrix_table,
)
from .fitting import (
    _build_pgls_diagnostics as _build_pgls_diagnostics,
)
from .fitting import (
    _fit_gls as _fit_gls,
)
from .fitting import (
    _gls_log_likelihood as _gls_log_likelihood,
)
from .fitting import (
    _grid_values as _grid_values,
)
from .fitting import (
    _lambda_log_likelihood as _lambda_log_likelihood,
)
from .fitting import (
    _quadratic_form as _quadratic_form,
)
from .fitting import (
    _resolve_lambda_fit as _resolve_lambda_fit,
)
from .fitting import (
    run_pgls as run_pgls,
)
from .formula import (
    coerce_numeric_value as _coerce_numeric_value,
)
from .formula import (
    parse_term_descriptor as _parse_term_descriptor,
)
from .formula import (
    resolve_formula_specification as _resolve_formula_specification,
)
from .models import (
    ComparativeFormulaSpecification,
    ComparativeHypothesisTestRow,
    ComparativeMultipleTestingReport,
    PGLSCoefficient,
    PGLSDiagnosticsReport,
    PGLSFittedObservation,
    PGLSFormulaAudit,
    PGLSInputReport,
    PGLSInteractionAudit,
    PGLSLambdaFitReport,
    PGLSLambdaProfileRow,
    PGLSLeverageRow,
    PGLSModelMatrixReport,
    PGLSModelMatrixRow,
    PGLSPredictorClassification,
    PGLSResidualOutlier,
    PGLSResult,
    PGLSTaxonExclusion,
)
from .multiple_testing import (
    _benjamini_hochberg_adjustment as _benjamini_hochberg_adjustment,
)
from .multiple_testing import (
    run_pgls_multiple_testing as run_pgls_multiple_testing,
)

__all__ = [
    "ComparativeFormulaSpecification",
    "ComparativeHypothesisTestRow",
    "ComparativeMultipleTestingReport",
    "PGLSCoefficient",
    "PGLSDiagnosticsReport",
    "PGLSFormulaAudit",
    "PGLSFittedObservation",
    "PGLSInputReport",
    "PGLSInteractionAudit",
    "PGLSLambdaFitReport",
    "PGLSLambdaProfileRow",
    "PGLSLeverageRow",
    "PGLSModelMatrixReport",
    "PGLSModelMatrixRow",
    "PGLSPredictorClassification",
    "PGLSResidualOutlier",
    "PGLSResult",
    "PGLSTaxonExclusion",
    "build_pgls_model_matrix",
    "inspect_pgls_inputs",
    "run_pgls",
    "run_pgls_multiple_testing",
    "write_pgls_model_matrix_table",
]
