from __future__ import annotations

from copy import deepcopy

from tests.support.geiger_fitcontinuous_lambda_reference import (
    GEIGER_FITCONTINUOUS_LAMBDA_REFERENCE_PAYLOADS,
)


def geiger_optimizer_triage_reference_payloads() -> dict[str, dict[str, object]]:
    """Synthetic governed reference payloads for optimizer-triage regression tests."""

    same_likelihood_payload = deepcopy(
        GEIGER_FITCONTINUOUS_LAMBDA_REFERENCE_PAYLOADS[
            "fitcontinuous-lambda-strong-signal-review"
        ]
    )
    same_likelihood_payload["summary"]["parameter_value"] = 0.75
    same_likelihood_payload["rows"][-1]["value"] = 0.75

    same_parameters_payload = deepcopy(
        GEIGER_FITCONTINUOUS_LAMBDA_REFERENCE_PAYLOADS[
            "fitcontinuous-lambda-strong-signal-review"
        ]
    )
    same_parameters_payload["summary"]["log_likelihood"] = -8.197817992659278
    same_parameters_payload["summary"]["aic"] = 22.395635985318556
    same_parameters_payload["summary"]["aicc"] = 23.595635985318555
    same_parameters_payload["summary"]["optimizer_result"][
        "best_log_likelihood"
    ] = -8.197817992659278
    same_parameters_payload["rows"][2]["value"] = -8.197817992659278
    same_parameters_payload["rows"][3]["value"] = 22.395635985318556
    same_parameters_payload["rows"][4]["value"] = 23.595635985318555

    boundary_payload = deepcopy(
        GEIGER_FITCONTINUOUS_LAMBDA_REFERENCE_PAYLOADS[
            "fitcontinuous-lambda-strong-signal-review"
        ]
    )
    boundary_payload["summary"]["parameter_value"] = 0.0
    boundary_payload["summary"]["hit_lower_parameter_boundary"] = True
    boundary_payload["summary"]["hit_upper_parameter_boundary"] = False
    boundary_payload["summary"]["log_likelihood"] = -12.697817992659278
    boundary_payload["summary"]["aic"] = 31.395635985318556
    boundary_payload["summary"]["aicc"] = 32.595635985318555
    boundary_payload["summary"]["optimizer_result"][
        "best_log_likelihood"
    ] = -12.697817992659278
    boundary_payload["rows"][2]["value"] = -12.697817992659278
    boundary_payload["rows"][3]["value"] = 31.395635985318556
    boundary_payload["rows"][4]["value"] = 32.595635985318555
    boundary_payload["rows"][5]["value"] = 0.0

    return {
        "same_likelihood_different_parameters": same_likelihood_payload,
        "different_likelihood_same_parameters": same_parameters_payload,
        "boundary_solution_review": boundary_payload,
    }
