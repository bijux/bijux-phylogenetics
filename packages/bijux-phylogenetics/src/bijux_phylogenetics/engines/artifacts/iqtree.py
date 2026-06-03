from __future__ import annotations

from dataclasses import dataclass, field
import gzip
from pathlib import Path
import re

_MODEL_TOKEN_PATTERN = r"[A-Za-z0-9+._-]+"  # nosec B105
_BEST_MODEL_PATTERNS = (
    re.compile(
        rf"best-fit model\s*:\s*(?P<model>{_MODEL_TOKEN_PATTERN})\s+chosen according to\s+(?P<criterion>[A-Z0-9]+)",
        re.IGNORECASE,
    ),
    re.compile(
        rf"best-fit model(?: according to (?P<criterion>[A-Z0-9]+)(?: score)?)?\s*[:=]\s*(?P<model>{_MODEL_TOKEN_PATTERN})",
        re.IGNORECASE,
    ),
)
_LOG_LIKELIHOOD_PATTERNS = (
    re.compile(
        r"(?:log-likelihood(?: of (?:the )?tree)?|log likelihood(?: of (?:the )?tree)?)\s*[:=]\s*(?P<value>-?[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)",
        re.IGNORECASE,
    ),
    re.compile(
        r"best score found\s*[:=]\s*(?P<value>-?[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)",
        re.IGNORECASE,
    ),
)
_MODEL_CANDIDATE_PATTERN = re.compile(
    r"^\s*(?P<rank>\d+)\s+"
    rf"(?P<model>{_MODEL_TOKEN_PATTERN})\s+"
    r"(?P<negative_log_likelihood>-?[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)\s+"
    r"(?P<parameter_count>\d+)\s+"
    r"(?P<aic>-?[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)\s+"
    r"(?P<aicc>-?[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)\s+"
    r"(?P<bic>-?[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)\s*$"
)
_MODEL_CANDIDATE_SORTED_PATTERN = re.compile(
    rf"^\s*(?P<model>{_MODEL_TOKEN_PATTERN})\s+"
    r"(?P<log_likelihood>-?[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)\s+"
    r"(?P<aic>-?[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)\s+[+-]\s+\S+\s+"
    r"(?P<aicc>-?[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)\s+[+-]\s+\S+\s+"
    r"(?P<bic>-?[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)\s+[+-]\s+\S+\s*$"
)
_SIDECAR_CANDIDATE_PATTERN = re.compile(
    rf"^(?P<model>{_MODEL_TOKEN_PATTERN})\s*:\s*"
    r"(?P<log_likelihood>-?[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)\s+"
    r"(?P<parameter_count>\d+)\s+"
    r"(?P<score>-?[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)\s*$"
)
_CRITERION_REPORT_PATTERNS = {
    "AIC": re.compile(
        rf"^\s*Akaike Information Criterion(?:\s*\(AIC\))?(?:\s+score)?\s*:\s*(?P<model>{_MODEL_TOKEN_PATTERN})\s*$",
        re.IGNORECASE,
    ),
    "AICc": re.compile(
        rf"^\s*Corrected Akaike Information Criterion(?:\s*\(AICc\))?(?:\s+score)?\s*:\s*(?P<model>{_MODEL_TOKEN_PATTERN})\s*$",
        re.IGNORECASE,
    ),
    "BIC": re.compile(
        rf"^\s*Bayesian Information Criterion(?:\s*\(BIC\))?(?:\s+score)?\s*:\s*(?P<model>{_MODEL_TOKEN_PATTERN})\s*$",
        re.IGNORECASE,
    ),
}
_SIDECAR_STRING_KEYS = {
    "best_model_aic": "best_model_aic",
    "best_model_aicc": "best_model_aicc",
    "best_model_bic": "best_model_bic",
}
_SIDECAR_FLOAT_KEYS = {
    "best_score_aic": "best_score_aic",
    "best_score_aicc": "best_score_aicc",
    "best_score_bic": "best_score_bic",
}


@dataclass(frozen=True, slots=True)
class IqtreeModelCandidate:
    rank: int
    model: str
    log_likelihood: float
    parameter_count: int | None
    aic: float
    aicc: float
    bic: float


@dataclass(slots=True)
class IqtreeModelSelectionSummary:
    selected_model: str | None
    selected_criterion: str | None
    best_model_aic: str | None
    best_model_aicc: str | None
    best_model_bic: str | None
    best_score_aic: float | None
    best_score_aicc: float | None
    best_score_bic: float | None
    candidate_count: int
    candidates: list[IqtreeModelCandidate] = field(default_factory=list)
    bic_near_best_models: list[str] = field(default_factory=list)


def parse_best_model_text(text: str) -> str | None:
    decision = parse_selected_model_decision_text(text)
    return None if decision is None else decision[0]


def parse_selected_model_decision_text(text: str) -> tuple[str, str | None] | None:
    for pattern in _BEST_MODEL_PATTERNS:
        match = pattern.search(text)
        if match is None:
            continue
        model = match.group("model")
        criterion = _normalize_model_selection_criterion(match.group("criterion"))
        return model, criterion
    return None


def parse_log_likelihood_text(text: str) -> float | None:
    for pattern in _LOG_LIKELIHOOD_PATTERNS:
        match = pattern.search(text)
        if match is None:
            continue
        try:
            return float(match.group("value"))
        except ValueError:
            continue
    return None


def parse_best_model_file(path: Path) -> str | None:
    if not path.exists():
        return None
    return parse_best_model_text(path.read_text(encoding="utf-8", errors="replace"))


def parse_log_likelihood_file(path: Path) -> float | None:
    if not path.exists():
        return None
    return parse_log_likelihood_text(path.read_text(encoding="utf-8", errors="replace"))


def resolve_iqtree_model_sidecar(prefix_path: Path) -> Path | None:
    for candidate in (
        prefix_path.with_suffix(".model"),
        prefix_path.with_suffix(".model.gz"),
    ):
        if candidate.exists():
            return candidate
    return None


def read_iqtree_model_sidecar(path: Path) -> str:
    if path.suffix == ".gz":
        return gzip.decompress(path.read_bytes()).decode("utf-8", errors="replace")
    return path.read_text(encoding="utf-8", errors="replace")


def parse_iqtree_model_selection_summary(
    *,
    iqtree_report_path: Path,
    model_sidecar_path: Path | None = None,
) -> IqtreeModelSelectionSummary | None:
    report_text = (
        iqtree_report_path.read_text(encoding="utf-8", errors="replace")
        if iqtree_report_path.exists()
        else ""
    )
    sidecar_text = None
    sidecar_candidates: dict[str, tuple[float, int]] = {}
    if model_sidecar_path is not None and model_sidecar_path.exists():
        sidecar_text = read_iqtree_model_sidecar(model_sidecar_path)
        sidecar_candidates = _parse_sidecar_candidates(sidecar_text)
    selected_model = None
    selected_criterion = None
    decision = parse_selected_model_decision_text(report_text)
    if decision is not None:
        selected_model, selected_criterion = decision
    criteria = _parse_report_criteria(report_text)
    candidates = _parse_report_candidates(
        report_text,
        sidecar_candidates=sidecar_candidates,
    )
    summary = IqtreeModelSelectionSummary(
        selected_model=selected_model,
        selected_criterion=selected_criterion,
        best_model_aic=criteria.get("AIC"),
        best_model_aicc=criteria.get("AICc"),
        best_model_bic=criteria.get("BIC"),
        best_score_aic=None,
        best_score_aicc=None,
        best_score_bic=None,
        candidate_count=len(candidates),
        candidates=candidates,
    )
    if sidecar_text is not None:
        _merge_model_sidecar(summary, sidecar_text)
    if summary.selected_model is None:
        summary.selected_model = (
            summary.best_model_bic or summary.best_model_aicc or summary.best_model_aic
        )
    if not summary.candidates and sidecar_candidates:
        summary.candidates = [
            IqtreeModelCandidate(
                rank=index,
                model=model,
                log_likelihood=log_likelihood,
                parameter_count=parameter_count,
                aic=float("nan"),
                aicc=float("nan"),
                bic=float("nan"),
            )
            for index, (model, (log_likelihood, parameter_count)) in enumerate(
                sidecar_candidates.items(),
                start=1,
            )
        ]
        summary.candidate_count = len(summary.candidates)
    if summary.selected_criterion is None and summary.selected_model is not None:
        summary.selected_criterion = _infer_selected_criterion(summary)
    if not _has_model_selection_data(summary):
        return None
    return summary


def write_iqtree_model_candidates_table(
    path: Path,
    summary: IqtreeModelSelectionSummary,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "rank",
        "model",
        "log_likelihood",
        "parameter_count",
        "aic",
        "aicc",
        "bic",
        "best_aic",
        "best_aicc",
        "best_bic",
        "selected_model",
    ]
    rows = ["\t".join(header)]
    for candidate in summary.candidates:
        rows.append(
            "\t".join(
                [
                    str(candidate.rank),
                    candidate.model,
                    format(candidate.log_likelihood, ".12g"),
                    (
                        ""
                        if candidate.parameter_count is None
                        else str(candidate.parameter_count)
                    ),
                    format(candidate.aic, ".12g"),
                    format(candidate.aicc, ".12g"),
                    format(candidate.bic, ".12g"),
                    "true" if candidate.model == summary.best_model_aic else "false",
                    "true" if candidate.model == summary.best_model_aicc else "false",
                    "true" if candidate.model == summary.best_model_bic else "false",
                    "true" if candidate.model == summary.selected_model else "false",
                ]
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _normalize_model_selection_criterion(value: str | None) -> str | None:
    if value is None:
        return None
    normalized = value.strip()
    if not normalized:
        return None
    upper = normalized.upper()
    if upper == "AICC":
        return "AICc"
    if upper in {"AIC", "BIC"}:
        return upper
    return normalized


def _parse_report_candidates(
    report_text: str,
    *,
    sidecar_candidates: dict[str, tuple[float, int]],
) -> list[IqtreeModelCandidate]:
    candidates: list[IqtreeModelCandidate] = []
    for line in report_text.splitlines():
        ranked_match = _MODEL_CANDIDATE_PATTERN.match(line)
        if ranked_match is not None:
            model = ranked_match.group("model")
            sidecar_metadata = sidecar_candidates.get(model)
            candidates.append(
                IqtreeModelCandidate(
                    rank=int(ranked_match.group("rank")),
                    model=model,
                    log_likelihood=-float(
                        ranked_match.group("negative_log_likelihood")
                    ),
                    parameter_count=(
                        None if sidecar_metadata is None else sidecar_metadata[1]
                    ),
                    aic=float(ranked_match.group("aic")),
                    aicc=float(ranked_match.group("aicc")),
                    bic=float(ranked_match.group("bic")),
                )
            )
            continue
        sorted_match = _MODEL_CANDIDATE_SORTED_PATTERN.match(line)
        if sorted_match is None:
            continue
        model = sorted_match.group("model")
        sidecar_metadata = sidecar_candidates.get(model)
        candidates.append(
            IqtreeModelCandidate(
                rank=len(candidates) + 1,
                model=model,
                log_likelihood=float(sorted_match.group("log_likelihood")),
                parameter_count=(
                    None if sidecar_metadata is None else sidecar_metadata[1]
                ),
                aic=float(sorted_match.group("aic")),
                aicc=float(sorted_match.group("aicc")),
                bic=float(sorted_match.group("bic")),
            )
        )
    return candidates


def _parse_sidecar_candidates(text: str) -> dict[str, tuple[float, int]]:
    candidates: dict[str, tuple[float, int]] = {}
    for line in text.splitlines():
        match = _SIDECAR_CANDIDATE_PATTERN.match(line.strip())
        if match is None:
            continue
        candidates[match.group("model")] = (
            float(match.group("log_likelihood")),
            int(match.group("parameter_count")),
        )
    return candidates


def _parse_report_criteria(report_text: str) -> dict[str, str]:
    criteria: dict[str, str] = {}
    for line in report_text.splitlines():
        for criterion, pattern in _CRITERION_REPORT_PATTERNS.items():
            match = pattern.match(line)
            if match is not None:
                criteria[criterion] = match.group("model")
    return criteria


def _merge_model_sidecar(
    summary: IqtreeModelSelectionSummary,
    sidecar_text: str,
) -> None:
    for raw_line in sidecar_text.splitlines():
        line = raw_line.strip()
        if not line or ":" not in line:
            continue
        key, value = line.split(":", 1)
        normalized_key = key.strip().lower()
        parsed_value = value.strip()
        if normalized_key in _SIDECAR_STRING_KEYS and parsed_value:
            setattr(summary, _SIDECAR_STRING_KEYS[normalized_key], parsed_value)
            continue
        if normalized_key in _SIDECAR_FLOAT_KEYS and parsed_value:
            try:
                setattr(
                    summary, _SIDECAR_FLOAT_KEYS[normalized_key], float(parsed_value)
                )
            except ValueError:
                continue
            continue
        if normalized_key == "best_model_list_bic":
            summary.bic_near_best_models = [
                token for token in parsed_value.split() if token
            ]


def _infer_selected_criterion(summary: IqtreeModelSelectionSummary) -> str | None:
    if summary.selected_model is None:
        return None
    matches: list[str] = []
    if summary.selected_model == summary.best_model_aic:
        matches.append("AIC")
    if summary.selected_model == summary.best_model_aicc:
        matches.append("AICc")
    if summary.selected_model == summary.best_model_bic:
        matches.append("BIC")
    if len(matches) == 1:
        return matches[0]
    return None


def _has_model_selection_data(summary: IqtreeModelSelectionSummary) -> bool:
    return any(
        (
            summary.selected_model is not None,
            summary.selected_criterion is not None,
            summary.best_model_aic is not None,
            summary.best_model_aicc is not None,
            summary.best_model_bic is not None,
            summary.best_score_aic is not None,
            summary.best_score_aicc is not None,
            summary.best_score_bic is not None,
            summary.candidate_count > 0,
            bool(summary.bic_near_best_models),
        )
    )
