from __future__ import annotations

GEIGER_LARGE_TREE_MODEL_FITTING_REFERENCE_PAYLOADS: dict[str, dict[str, object]] = {
    "fitcontinuous-pagel-lambda-100-taxa": {
        "fit_summary": {
            "taxon_count": 100,
            "model_name": "lambda",
            "parameter_name": "lambda",
            "parameter_value": 0.9078408096138019,
            "rate": 0.8442426979533381,
            "log_likelihood": -77.31623275813774,
            "aic": 160.6324655162755,
            "aicc": 160.8824655162755,
            "runtime_seconds": 0.34,
            "optimizer_result": {
                "best_method": "subplex",
                "attempt_count": 100,
                "attempted_methods": ["L-BFGS-B", "subplex"],
                "converged_attempt_count": 100,
                "convergence_code": 0,
                "best_log_likelihood": -77.31623275813774,
            },
        }
    },
    "fitdiscrete-er-binary-100-taxa": {
        "fit_summary": {
            "taxon_count": 100,
            "model_name": "ER",
            "parameter_name": None,
            "parameter_value": None,
            "representative_rate": 0.332095429178085,
            "log_likelihood": -26.55280534511656,
            "aic": 55.10561069023312,
            "aicc": 55.14642701676373,
            "runtime_seconds": 0.31700000000000017,
            "optimizer_result": {
                "best_method": "Brent",
                "attempt_count": 100,
                "converged_attempt_count": 100,
                "convergence_code": 0,
                "best_log_likelihood": -26.55280534511656,
            },
        }
    },
    "fitcontinuous-brownian-512-taxa": {
        "fit_summary": {
            "taxon_count": 512,
            "model_name": "BM",
            "parameter_name": None,
            "parameter_value": None,
            "rate": 0.7214415560426788,
            "log_likelihood": -294.86930788906,
            "aic": 593.73861577812,
            "aicc": 593.7621914166268,
            "runtime_seconds": 0.137,
            "optimizer_result": {
                "best_method": "Brent",
                "attempt_count": 100,
            },
        }
    },
}
