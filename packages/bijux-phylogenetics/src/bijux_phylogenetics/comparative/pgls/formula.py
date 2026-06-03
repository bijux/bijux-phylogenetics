from __future__ import annotations

import math

from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .models import (
    ComparativeFormulaSpecification,
    PGLSPredictorClassification,
    _FormulaTermDescriptor,
)


def parse_term_descriptor(raw_term: str) -> _FormulaTermDescriptor:
    raw_term = raw_term.strip()
    if "(" not in raw_term:
        return _FormulaTermDescriptor(
            raw_term=raw_term,
            source_column=raw_term,
            transformation=None,
        )
    if not raw_term.endswith(")") or raw_term.count("(") != 1:
        raise ComparativeMethodError(
            f"unsupported comparative term syntax '{raw_term}'"
        )
    transformation, inner = raw_term.split("(", 1)
    transformation = transformation.strip()
    source_column = inner[:-1].strip()
    if transformation not in {"log", "log10", "sqrt"}:
        raise ComparativeMethodError(
            f"unsupported comparative transformation '{transformation}' in term '{raw_term}'"
        )
    if not source_column:
        raise ComparativeMethodError(
            f"comparative transformation term '{raw_term}' is missing a source column"
        )
    return _FormulaTermDescriptor(
        raw_term=raw_term,
        source_column=source_column,
        transformation=transformation,
    )


def parse_right_hand_side_terms(right_hand_side: str) -> list[tuple[str, str]]:
    terms: list[tuple[str, str]] = []
    current: list[str] = []
    current_sign = "+"
    parenthesis_depth = 0
    for character in right_hand_side:
        if character == "(":
            parenthesis_depth += 1
            current.append(character)
            continue
        if character == ")":
            parenthesis_depth -= 1
            if parenthesis_depth < 0:
                raise ComparativeMethodError(
                    "comparative formula has unmatched closing parenthesis"
                )
            current.append(character)
            continue
        if parenthesis_depth == 0 and character in {"+", "-"}:
            term = "".join(current).strip()
            if term:
                terms.append((current_sign, term))
            current = []
            current_sign = character
            continue
        current.append(character)
    if parenthesis_depth != 0:
        raise ComparativeMethodError(
            "comparative formula has unmatched opening parenthesis"
        )
    trailing_term = "".join(current).strip()
    if trailing_term:
        terms.append((current_sign, trailing_term))
    return terms


def coerce_numeric_value(
    raw_value: str, *, descriptor: _FormulaTermDescriptor
) -> float:
    value = float(raw_value)
    if descriptor.transformation is None:
        return value
    if descriptor.transformation == "log":
        if value <= 0.0:
            raise ValueError("log transformation requires strictly positive values")
        return math.log(value)
    if descriptor.transformation == "log10":
        if value <= 0.0:
            raise ValueError("log10 transformation requires strictly positive values")
        return math.log10(value)
    if descriptor.transformation == "sqrt":
        if value < 0.0:
            raise ValueError("sqrt transformation requires non-negative values")
        return math.sqrt(value)
    raise ComparativeMethodError(
        f"unsupported transformation '{descriptor.transformation}'"
    )


def column_supports_numeric_values(
    rows: list[dict[str, str]],
    column: str,
) -> bool:
    observed_numeric = False
    for row in rows:
        raw_value = row.get(column, "")
        if not raw_value:
            continue
        try:
            float(raw_value)
        except ValueError:
            return False
        observed_numeric = True
    return observed_numeric


def resolve_formula_specification(
    *,
    response: str | None,
    predictors: list[str] | None,
    formula: str | None,
) -> ComparativeFormulaSpecification:
    if formula:
        if response is not None or predictors:
            raise ComparativeMethodError(
                "provide either a formula or explicit response/predictors, not both"
            )
        return parse_formula(formula)
    if response is None:
        raise ComparativeMethodError(
            "PGLS requires a response column when no formula is provided"
        )
    requested_predictors = list(predictors or [])
    if not requested_predictors:
        raise ComparativeMethodError("PGLS requires at least one predictor column")
    return ComparativeFormulaSpecification(
        response=response,
        formula=f"{response} ~ {' + '.join(requested_predictors)}",
        predictors=requested_predictors,
        interaction_terms=[],
        include_intercept=True,
    )


def parse_formula(formula: str) -> ComparativeFormulaSpecification:
    if "~" not in formula:
        raise ComparativeMethodError("comparative formula must contain '~'")
    response, right_hand_side = [part.strip() for part in formula.split("~", 1)]
    if not response:
        raise ComparativeMethodError(
            "comparative formula requires a response on the left-hand side"
        )
    raw_terms = parse_right_hand_side_terms(right_hand_side)
    if not raw_terms:
        raise ComparativeMethodError(
            "comparative formula requires at least one predictor term"
        )
    predictors: list[str] = []
    interaction_terms: list[str] = []
    include_intercept = True
    for sign, raw_term in raw_terms:
        if raw_term in {"0", "1"}:
            if raw_term == "0" or sign == "-":
                include_intercept = False
            elif raw_term == "1":
                include_intercept = True
            continue
        if sign == "-":
            raise ComparativeMethodError(
                f"unsupported comparative formula subtraction for term '{raw_term}'"
            )
        if "*" in raw_term:
            factors = [
                factor.strip() for factor in raw_term.split("*") if factor.strip()
            ]
            if len(factors) < 2:
                raise ComparativeMethodError(
                    f"invalid interaction expansion '{raw_term}'"
                )
            for factor in factors:
                if factor not in predictors:
                    predictors.append(factor)
            interaction = ":".join(factors)
            if interaction not in interaction_terms:
                interaction_terms.append(interaction)
            continue
        if ":" in raw_term:
            factors = [
                factor.strip() for factor in raw_term.split(":") if factor.strip()
            ]
            if len(factors) < 2:
                raise ComparativeMethodError(f"invalid interaction term '{raw_term}'")
            interaction = ":".join(factors)
            if interaction not in interaction_terms:
                interaction_terms.append(interaction)
            continue
        if raw_term not in predictors:
            predictors.append(raw_term)
    if not predictors and not interaction_terms:
        raise ComparativeMethodError(
            "comparative formula requires at least one predictor term"
        )
    return ComparativeFormulaSpecification(
        response=response,
        formula=formula.strip(),
        predictors=predictors,
        interaction_terms=interaction_terms,
        include_intercept=include_intercept,
    )


def interaction_column_names(
    interaction: str,
    predictor_reports: dict[str, PGLSPredictorClassification],
) -> list[str]:
    encoded_components: list[list[str]] = []
    for factor in interaction.split(":"):
        report = predictor_reports[factor]
        if report.kind == "numeric":
            encoded_components.append([factor])
        else:
            encoded_components.append(list(report.encoded_columns or []))
    column_names = [""]
    for component_names in encoded_components:
        column_names = [
            f"{left}:{right}".strip(":")
            for left in column_names
            for right in component_names
        ]
    return column_names


def interaction_values(
    interaction: str,
    encoded_main_effects: dict[str, list[tuple[str, float]]],
) -> list[tuple[str, float]]:
    values = [("", 1.0)]
    for factor in interaction.split(":"):
        expanded: list[tuple[str, float]] = []
        for left_name, left_value in values:
            for right_name, right_value in encoded_main_effects[factor]:
                expanded.append(
                    (f"{left_name}:{right_name}".strip(":"), left_value * right_value)
                )
        values = expanded
    return values
